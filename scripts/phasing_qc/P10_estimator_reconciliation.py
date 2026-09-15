#!/usr/bin/env python3
"""
Reconcile `decay.implied_design_effect` with `readlevel.design_effect`
(docs/PRIORITIES.md item 3).

On HG00146 chr15 at 500 bp the two disagree ~2x. This walks from one to the
other ONE assumption at a time, on the same windows, so each step's share of
the gap is measured rather than argued:

  A  median readlevel.design_effect, two read halves   (the reported number)
  B  ... all reads as one group                        (split artefact?)
  C  ... mean instead of median                        (right skew)
  D  per-window 1+(cbar-1)*icc_w, median / mean        (formula vs direct)
  E  icc pooled across windows, variance-weighted      (pooling weights)
  F  E de-attenuated by (1-2eps)^2                     (call error)
  G  implied_design_effect(curve, median gap, median C) (production function)
  H  pooled curve read at each READ's own CpG pairs      (what molecules observed)
  I  G with C = positions seen by >=50% of reads         (spurious-position check)

Usage: P10_estimator_reconciliation.py <extract.tsv.gz> [max_rows=12e6] [window_bp=500]
"""
import os
import sys
from collections import defaultdict

import numpy as np

_HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, _HERE)
sys.path.insert(0, os.path.dirname(os.path.dirname(_HERE)))
from P08_decay_from_extract import to_matrices  # noqa: E402
from P09_design_effect_cohort import stream, window_metrics  # noqa: E402
from ont_asm_caller.decay import comethylation_decay, deattenuate, implied_design_effect  # noqa: E402


def curve_de_per_read(curve, pos, m):
    """Mean over reads of 1 + (C_r-1) * mean corr over THAT READ's CpG pairs,
    reading the pooled curve at the real distances. Uses what each molecule
    actually observed, not the window's union of positions."""
    obs = np.isfinite(m)
    out = []
    for r in range(m.shape[0]):
        q = pos[obs[r]]
        if len(q) < 2:
            continue
        d = np.abs(q[:, None] - q[None, :])[np.triu_indices(len(q), 1)]
        corr = np.interp(d, curve.centres, curve.corr,
                         left=float(curve.corr[0]), right=float(curve.corr[-1]))
        out.append(1.0 + (len(q) - 1) * float(np.clip(corr, 0, 1).mean()))
    return float(np.mean(out)) if out else np.nan


def support(mats):
    """Fraction of a window's reads that observe each position."""
    frac = np.concatenate([np.isfinite(m).mean(axis=0) for _, m in mats])
    real = np.array([(np.isfinite(m).mean(axis=0) >= 0.5).sum() for _, m in mats])
    return frac, real


def main():
    path = sys.argv[1]
    max_rows = int(float(sys.argv[2])) if len(sys.argv) > 2 else 12_000_000
    window_bp = int(sys.argv[3]) if len(sys.argv) > 3 else 500

    stats = defaultdict(float)
    win = defaultdict(lambda: defaultdict(dict))
    for rid, pos, call in stream(path, window_bp, 1, max_rows, stats):
        win[pos // window_bp][rid][pos] = call
    mats = to_matrices(win)
    del win
    eps = stats["err_sum"] / max(stats["kept"], 1)
    att = (1 - 2 * eps) ** 2
    rows = [window_metrics(p, m) for p, m in mats]
    g = lambda k: np.array([r[k] for r in rows], float)  # noqa: E731
    fin = lambda a: a[np.isfinite(a)]                      # noqa: E731

    # E: pool icc across windows exactly as decay._accumulate does
    sxy = sxx = syy = 0.0
    for pos, m in mats:
        obs = np.isfinite(m)
        c = np.where(obs, m - np.nanmean(m, axis=0), 0.0)
        for i in range(m.shape[1] - 1):
            both = obs[:, i:i + 1] & obs[:, i + 1:]
            sxy += float((c[:, i:i + 1] * c[:, i + 1:] * both).sum())
            sxx += float(((c[:, i:i + 1] ** 2) * both).sum())
            syy += float(((c[:, i + 1:] ** 2) * both).sum())
    icc_pool = sxy / np.sqrt(sxx * syy)
    cbar_med = float(np.nanmedian(g("cbar")))
    nc_med = int(np.median(g("n_cpgs")))
    gap_med = float(np.median([np.median(np.diff(p)) for p, _ in mats if len(p) > 1]))
    gap_mean = float(np.median([(p[-1] - p[0]) / (len(p) - 1) for p, _ in mats if len(p) > 1]))

    curve = comethylation_decay(mats)
    curve_c = deattenuate(curve, eps)

    per_read = np.array([curve_de_per_read(curve_c, p, m) for p, m in mats])
    per_read_raw = np.array([curve_de_per_read(curve, p, m) for p, m in mats])
    sup, n_real = support(mats)

    print(f"{os.path.basename(path)}  first {max_rows:,} rows  {len(mats)} windows of {window_bp} bp")
    print(f"eps={eps:.4f}  attenuation={att:.3f}  median C={nc_med}  median cbar={cbar_med:.1f}  "
          f"median gap={gap_med:.0f} bp  median span/(C-1)={gap_mean:.0f} bp\n")
    L = [
        ("A  median design_effect, read halves", np.median(g("de_split"))),
        ("B  median design_effect, all reads", np.median(g("de_one"))),
        ("C  MEAN design_effect, all reads", np.mean(g("de_one"))),
        ("   median de_het (heterogeneous-mean baseline)", np.median(fin(g("de_het")))),
        ("   mean de_het", np.mean(fin(g("de_het")))),
        ("D  median 1+(cbar-1)*icc_w per window", np.median(fin(g("de_icc")))),
        ("D' mean 1+(cbar-1)*icc_w per window", np.mean(fin(g("de_icc")))),
        ("E  pooled icc_w, median cbar  (raw)", 1 + (cbar_med - 1) * icc_pool),
        ("F  E de-attenuated", 1 + (cbar_med - 1) * icc_pool / att),
        ("G  implied_design_effect, median gap, raw curve", implied_design_effect(curve, gap_med, nc_med)),
        ("G' implied_design_effect, median gap, de-attenuated", implied_design_effect(curve_c, gap_med, nc_med)),
        ("G''implied_design_effect, span/(C-1) gap, de-att.", implied_design_effect(curve_c, gap_mean, nc_med)),
        ("H  curve at each READ's own pairs, raw, median", np.nanmedian(per_read_raw)),
        ("H' ... de-attenuated, median", np.nanmedian(per_read)),
        ("H''... de-attenuated, mean", np.nanmean(per_read)),
        ("I  implied_design_effect, C = well-supported positions, de-att.",
         implied_design_effect(curve_c, gap_mean * (nc_med - 1) / max(np.median(n_real) - 1, 1),
                               int(np.median(n_real)))),
    ]
    for label, v in L:
        print(f"  {label:<55}{v:>7.2f}")
    print(f"\n  per-position read support: median {np.median(sup):.2f}; "
          f"{(sup < 0.2).mean():.1%} of positions seen by <20% of the window's reads; "
          f"{(sup >= 0.5).mean():.1%} by >=50%")
    print(f"  positions per window seen by >=50% of reads: median {np.median(n_real):.0f} "
          f"(vs n_cpgs {nc_med}, cbar {cbar_med:.1f})")
    print(f"\n  pooled icc_w={icc_pool:.4f}   median per-window icc_w={np.nanmedian(g('icc_w')):.4f}   "
          f"mean={np.nanmean(g('icc_w')):.4f}")
    print(f"  corr(de_one, de_icc) across windows = "
          f"{np.corrcoef(*[a[np.isfinite(g('de_icc'))] for a in (g('de_one'), g('de_icc'))])[0, 1]:.3f}")
    print(f"  curve: " + "  ".join(f"{c:.0f}:{r:.2f}" for c, r in zip(curve.centres, curve.corr)))

    # Is the steep short-range decay real, or sparse positions that are the SAME
    # CpG called at a shifted coordinate in a few reads? Recompute the curve on
    # positions observed by >= min_support of a window's reads.
    for min_sup in (0.2, 0.5):
        kept = []
        for pos, m in mats:
            keep = np.isfinite(m).mean(axis=0) >= min_sup
            if keep.sum() >= 2:
                kept.append((pos[keep], m[:, keep]))
        cv = comethylation_decay(kept)
        print(f"  curve, positions seen by >={min_sup:.0%} of reads: "
              + "  ".join(f"{c:.0f}:{r:.2f}" for c, r in zip(cv.centres, cv.corr)))


if __name__ == "__main__":
    main()
