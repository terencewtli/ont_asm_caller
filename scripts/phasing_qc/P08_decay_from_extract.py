#!/usr/bin/env python3
"""
Measure `decay_bp` and the design effect straight from `modkit extract`,
WITHOUT a haplotagged BAM and WITHOUT the caller's region table.

WHY THIS EXISTS SEPARATELY FROM P06/P07
---------------------------------------
`decay_bp` is the number three separate analyses are conditional on, and it is
the stated blocker on quoting anything from `benchmarks/compare_pattern_methods.py`
or from `docs/2026-09-05_calibration_critique.md`. P06/P07 measure it, but they
need HP tags from the haplotagged BAM and the production region table.

Neither is actually necessary for this particular measurement: **co-methylation
decay along a molecule is a within-read quantity and does not involve
haplotypes at all.** Dropping both dependencies turns the top-priority
measurement into something runnable on an extract file alone, which is usually
the first thing available.

The design effect reported here is likewise computed per window from all reads
pooled. That is the right quantity for "is region pooling legal" -- it asks
whether CpGs on one molecule are redundant, which has nothing to do with which
haplotype the molecule came from.

WHAT IT CORRECTS FOR
--------------------
Both defects found against real extract output (see
docs/2026-09-07_cpel_review_and_pattern_tests.md):

  * strand -- extract is per-strand; a CpG's minus-strand C sits at p+1, so
    minus-strand positions are normalised by -1 or every CpG splits into two
    columns 1 bp apart sharing zero reads.
  * call confidence -- one row per (read, position, mod_code), so each CpG
    appears as both "m" and "h" and P(canonical) = 1 - p_m - p_h. Takes the
    argmax over {C, m, h}; anything below CONF_THRESH becomes a NO-CALL rather
    than a coin flip, because a wrong binary call attenuates every correlation
    by (1-2*eps)^2 and the design effect IS a correlation.

It also estimates eps from modkit's own posteriors and reports the
de-attenuated amplitude alongside the raw one, so the reader can see how much
of the answer is call noise.

Usage: P08_decay_from_extract.py <extract.tsv.gz> [max_rows] [window_bp]
"""
import gzip
import os
import sys
from collections import defaultdict

import numpy as np

_HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(os.path.dirname(_HERE)))
sys.path.insert(0, os.path.dirname(os.path.dirname(_HERE)) + "/github/ont_asm_caller")
try:
    from ont_asm_caller.decay import comethylation_decay, implied_design_effect, deattenuate
    from ont_asm_caller.readlevel import design_effect
except ImportError:                       # running from inside the repo
    sys.path.insert(0, os.path.dirname(_HERE) + "/..")
    from ont_asm_caller.decay import comethylation_decay, implied_design_effect, deattenuate
    from ont_asm_caller.readlevel import design_effect

CONF_THRESH = 0.80
H_POLICY = "nocall"
MIN_READS = 6
MIN_CPGS = 4


def stream(path, max_rows):
    """Yield (read_id, plus_strand_cpg_pos, call, confidence). call is 1/0/None."""
    stats = defaultdict(int)
    with gzip.open(path, "rt") as fh:
        fh.readline()                                    # header
        pend = {}
        for i, line in enumerate(fh):
            if max_rows and i >= max_rows:
                break
            f = line.split("\t")
            if len(f) < 12:
                continue
            pos = int(f[2])
            if pos < 0:                                  # soft-clipped, unmapped
                stats["unmapped"] += 1
                continue
            code = f[11]
            if code not in ("m", "h"):
                continue
            if f[6] == "-":
                pos -= 1                                 # normalise to plus strand
                stats["minus"] += 1
            else:
                stats["plus"] += 1
            key = (f[0], pos)
            d = pend.setdefault(key, {})
            d[code] = float(f[10])
            if len(d) < 2:
                continue
            del pend[key]
            p_m, p_h = d.get("m", 0.0), d.get("h", 0.0)
            p_c = 1.0 - p_m - p_h
            conf, which = max((p_c, "c"), (p_m, "m"), (p_h, "h"))
            stats["calls"] += 1
            if conf < CONF_THRESH:
                stats["ambiguous"] += 1
                yield f[0], pos, None, conf
            elif which == "m":
                yield f[0], pos, 1, conf
            elif which == "c":
                yield f[0], pos, 0, conf
            else:
                stats["hmc"] += 1
                yield f[0], pos, (None if H_POLICY == "nocall"
                                  else int(H_POLICY == "methylated")), conf
    stream.stats = stats


def build_windows(path, max_rows, window_bp):
    win = defaultdict(lambda: defaultdict(dict))
    confs = []
    n = 0
    for rid, pos, call, conf in stream(path, max_rows):
        n += 1
        if call is not None:
            win[pos // window_bp][rid][pos] = call
            confs.append(conf)
    return win, np.asarray(confs), n


def _merge_adjacent(positions):
    """Collapse positions 1 bp apart onto the lower one.

    Genuine CpGs CANNOT be 1 bp apart: a CpG at p means p is a C and p+1 is a G,
    so p+1 cannot start another CpG. Any 1 bp pair is therefore the SAME physical
    CpG whose two strands did not land on the same coordinate. Measured on real
    HG00146 chr15 after the -1 minus-strand shift: 5.9% of adjacent gaps were
    still exactly 1 bp, and only 37.5% of positions were supported by both
    strands. (The shift itself is definitely right -- without it the both-strand
    overlap is 0.018 instead of 0.375.)

    Returns {raw_position: canonical_position}.
    """
    out, canon = {}, None
    for q in sorted(positions):
        if canon is None or q - canon > 1:
            canon = q
        out[q] = canon
    return out


def to_matrices(win):
    out = []
    for w, reads in win.items():
        if len(reads) < MIN_READS:
            continue
        raw = sorted({p for d in reads.values() for p in d})
        cmap = _merge_adjacent(raw)
        reads = {r: {cmap[p]: c for p, c in d.items()} for r, d in reads.items()}
        positions = sorted(set(cmap.values()))
        if len(positions) < MIN_CPGS:
            continue
        idx = {p: j for j, p in enumerate(positions)}
        m = np.full((len(reads), len(positions)), np.nan)
        for r, d in enumerate(reads.values()):
            for p, c in d.items():
                m[r, idx[p]] = c
        out.append((np.asarray(positions, float), m))
    return out


def main():
    path = sys.argv[1]
    max_rows = int(float(sys.argv[2])) if len(sys.argv) > 2 else 20_000_000
    window_bp = int(sys.argv[3]) if len(sys.argv) > 3 else 2000

    name = os.path.basename(path).split("_modkit")[0]
    print(f"{name}: streaming up to {max_rows:,} rows from {path}", flush=True)
    win, confs, n_calls = build_windows(path, max_rows, window_bp)
    st = stream.stats
    tot = max(st["calls"], 1)
    print(f"  {st['calls']:,} paired calls   "
          f"({st['unmapped']:,} rows unmapped, ref_position<0)")
    print(f"  strand: {st['plus']:,} plus / {st['minus']:,} minus"
          + ("   <-- BOTH present, normalisation required" if st["plus"] and st["minus"] else ""))
    print(f"  {st['ambiguous']:,} ({st['ambiguous']/tot:.1%}) below CONF_THRESH={CONF_THRESH} -> no-call")
    print(f"  {st['hmc']:,} ({st['hmc']/tot:.1%}) 5hmC-dominant -> {H_POLICY}")

    mats = to_matrices(win)
    if not mats:
        print("no usable windows"); return
    spans = [p[-1] - p[0] for p, _ in mats]
    ncpg = [len(p) for p, _ in mats]
    nread = [m.shape[0] for _, m in mats]
    print(f"\n  {len(mats)} windows of {window_bp} bp   "
          f"median {np.median(ncpg):.0f} CpGs, {np.median(nread):.0f} reads, "
          f"span {np.median(spans):.0f} bp")

    # eps: modkit's own posterior says how often a KEPT call is wrong
    eps = float(np.mean(1.0 - confs)) if len(confs) else 0.0
    print(f"  estimated call error rate eps = {eps:.4f} "
          f"-> correlation attenuation factor (1-2eps)^2 = {(1-2*eps)**2:.3f}")

    print(f"\n  {'stratum':<26}{'n':>7}{'decay_bp':>10}{'amp':>8}{'amp_corr':>10}"
          f"{'r2':>7}{'spacing':>9}{'DE(fit)':>9}{'DE(meas)':>10}")

    def report(label, subset):
        if len(subset) < 50:
            print(f"  {label:<26}{len(subset):>7}   (too few windows)")
            return
        c = comethylation_decay(subset)
        cc = deattenuate(c, eps)
        spacing = float(np.median([np.median(np.diff(p)) for p, _ in subset
                                   if len(p) > 1]))
        nc = int(np.median([len(p) for p, _ in subset]))
        de_fit = implied_design_effect(cc, spacing, nc)
        de_meas = float(np.median([design_effect(m[:len(m)//2], m[len(m)//2:])
                                   for _, m in subset if len(m) >= 6]))
        print(f"  {label:<26}{len(subset):>7}{c.decay_bp:>10.0f}{c.amplitude:>8.3f}"
              f"{cc.amplitude:>10.3f}{c.r2:>7.3f}{spacing:>9.0f}"
              f"{de_fit:>9.2f}{de_meas:>10.2f}")

    report("ALL", mats)
    mu = np.array([np.nanmean(m) for _, m in mats])
    for lo, hi, nm in [(0, .2, "methylation < 0.2"), (.2, .8, "methylation 0.2-0.8"),
                       (.8, 1.01, "methylation > 0.8")]:
        report(nm, [x for x, u in zip(mats, mu) if np.isfinite(u) and lo <= u < hi])
    dens = np.array([len(p) / max(p[-1] - p[0], 1) * 1000 for p, _ in mats])
    q = np.nanquantile(dens, [1/3, 2/3])
    for lo, hi, nm in [(-np.inf, q[0], "CpG-sparse tercile"),
                       (q[0], q[1], "CpG-medium tercile"),
                       (q[1], np.inf, "CpG-dense tercile")]:
        report(nm, [x for x, d in zip(mats, dens) if lo <= d < hi])

    print(f"\n  NOTE: windows are {window_bp} bp, so a decay length approaching that")
    print( "  is a LOWER BOUND -- the curve cannot see past a window.")


if __name__ == "__main__":
    main()
