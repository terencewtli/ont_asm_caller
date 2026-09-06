"""
Benchmark v2 -- on data that looks like long-read data.

Differences from compare_methods.py, both of which invalidated its numbers:

  1. READS ARE SHARED across the CpGs of a region (simulate_reads.py). The v1
     simulator drew an independent depth and an independent count for every
     CpG, so pooling read x CpG counts was legitimate there and only there.
  2. The methylome is BIMODAL, not uniformly mu=0.5. v1 used the default
     baseline_mu=0.5 everywhere, which is the one point at which the
     Pearson-chi-square dispersion estimator happens to be unbiased.

Reports FDR and power after BH at a realistic genome-wide scale.

Run:  python benchmarks/compare_methods_v2.py
"""
from __future__ import annotations

import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

import numpy as np

from ont_asm_caller.simulate_reads import simulate_read_regions, per_cpg_counts
from ont_asm_caller.dispersion import estimate_dispersion, estimate_dispersion_trend
from ont_asm_caller.model import test_locus
from ont_asm_caller.readlevel import test_region_reads, test_region_perm, design_effect


def bh(p, alpha=0.05):
    p = np.asarray(p, dtype=float)
    n = len(p)
    if n == 0:
        return np.zeros(0, bool)
    o = np.argsort(p)
    thresh = alpha * (np.arange(1, n + 1) / n)
    passed = p[o] <= thresh
    k = np.max(np.nonzero(passed)[0]) + 1 if passed.any() else 0
    out = np.zeros(n, bool)
    out[o[:k]] = True
    return out


def score(name, pvals, truth, deltas, true_deltas, alpha=0.05):
    pvals = np.asarray(pvals, float)
    truth = np.asarray(truth, bool)
    sig = bh(pvals, alpha)
    ncall = int(sig.sum())
    fdr = float((~truth[sig]).mean()) if ncall else float("nan")
    power = float(sig[truth].mean()) if truth.any() else float("nan")
    if ncall and np.isfinite(deltas[sig]).any():
        bias = float(np.nanmean(np.abs(deltas[sig]) - np.abs(true_deltas[sig])))
    else:
        bias = float("nan")
    print(f"  {name:<44s}{ncall:>7d}{fdr:>10.3f}{power:>9.3f}{bias:>11.3f}")


def run(co_methylation, n_regions=40000, pi=0.01, seed=7):
    regs = simulate_read_regions(n_regions=n_regions, pi=pi,
                                 co_methylation=co_methylation,
                                 baseline_mu="bimodal", seed=seed)
    truth = np.array([r.is_true_asm for r in regs])
    tdelta = np.array([abs(r.true_delta) for r in regs])
    de = np.mean([design_effect(r.m1, r.m2) for r in regs[:2000]])

    print(f"\n{'='*84}")
    print(f"co_methylation = {co_methylation}   regions = {n_regions}   true ASM = {truth.sum()}"
          f"   empirical design effect = {de:.2f}")
    print(f"{'='*84}")
    print(f"  {'method':<44s}{'called':>7s}{'FDR':>10s}{'power':>9s}{'|d| bias':>11s}")

    # --- dispersion, estimated per-CpG on the same data -------------------
    x1, n1, x2, n2 = per_cpg_counts(regs)
    rho_global = estimate_dispersion(x1, n1, x2, n2)
    trend = estimate_dispersion_trend(x1, n1, x2, n2, n_bins=8)

    # --- OLD PATH: pool read x CpG counts, single global rho ---------------
    pooled = [r.counts() for r in regs]
    rho_pool = estimate_dispersion(*map(np.array, zip(*pooled)))
    res = [test_locus("chr1", r.start, *c, rho_pool) for r, c in zip(regs, pooled)]
    score(f"region, pooled readxCpG, global rho={rho_pool:.3f}",
          [x.pval for x in res], truth,
          np.array([x.delta for x in res]), tdelta)

    # --- OLD PATH but with the mean-dependent rho -------------------------
    # trend must be estimated on the SAME counts it is applied to: a per-CpG
    # rho does not describe pooled region counts (an earlier version of this
    # benchmark made that mistake and produced FDR 0.67).
    trend_pool = estimate_dispersion_trend(*map(np.array, zip(*pooled)), n_bins=8)
    res_t = []
    for r, c in zip(regs, pooled):
        mu = (c[0] + c[2]) / max(c[1] + c[3], 1)
        res_t.append(test_locus("chr1", r.start, *c, float(trend_pool(mu))))
    score("region, pooled readxCpG, rho(mu) trend",
          [x.pval for x in res_t], truth,
          np.array([x.delta for x in res_t]), tdelta)

    # --- NEW PATH: read-level Welch ---------------------------------------
    rl = [test_region_reads("chr1", r.start, r.start, r.m1, r.m2) for r in regs]
    ok = np.array([x is not None for x in rl])
    score("read-level Welch (molecule = unit)",
          [x.pval for x in np.array(rl, dtype=object)[ok]], truth[ok],
          np.array([x.delta for x in np.array(rl, dtype=object)[ok]]), tdelta[ok])

    # --- NEW PATH: exact permutation, on a subsample (it is slow) ---------
    # Permutation p is floored at 1/(n_perm+1), so it CANNOT survive BH across
    # tens of thousands of regions. It is a confirmation step for screened
    # candidates, not a genome-wide test. Scored here on the true positives
    # plus a null sample only, to show calibration rather than fake a ranking.
    idx = np.random.default_rng(0).choice(len(regs), size=min(1500, len(regs)),
                                          replace=False)
    pr = [(i, test_region_perm("chr1", regs[i].start, regs[i].start,
                               regs[i].m1, regs[i].m2, n_perm=999, seed=i))
          for i in idx]
    pr = [(i, x) for i, x in pr if x is not None]
    if pr:
        ii = np.array([i for i, _ in pr])
        pv = np.array([x.pval for _, x in pr])
        print(f"  {'read-level permutation: type I @0.05':<44s}"
              f"{'':>7s}{(pv[~truth[ii]] < 0.05).mean():>10.3f}"
              f"{'(nominal 0.05, n=' + str((~truth[ii]).sum()) + ' nulls)':>21s}")
    return de


if __name__ == "__main__":
    print("Benchmark v2: shared reads + bimodal methylome. BH-FDR alpha=0.05.")
    print("co_methylation=0.0 reproduces the (unrealistic) v1 independence assumption.")
    for cm in (0.0, 0.6, 1.0):
        run(cm)
    print("\nNOTE: the honest value of co_methylation is an EMPIRICAL question.")
    print("Measure it from real `modkit extract` output with readlevel.design_effect()")
    print("before quoting any FDR number from this benchmark.")
