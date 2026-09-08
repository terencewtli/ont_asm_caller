#!/usr/bin/env python3
"""
Measure `decay_bp` -- the co-methylation decay length along a molecule -- on
real data, stratified.

THIS IS THE FIRST THING TO RUN once P06 has produced read matrices. Three
separate analyses are conditional on this one number and none of them can be
quoted before it exists:

  * whether `region.cluster_cpgs` pooling is legal at all (the design effect;
    github/ont_asm_caller/docs/2026-09-05_calibration_critique.md, "the number
    that decides everything")
  * every result in benchmarks/compare_pattern_methods.py
  * the parent project's co-methylation-decay curve, md/20260903_qc_review.md s7.3

The three must agree. This script computes the first two together and prints
them side by side, because if `implied_design_effect` and the independently
implemented `readlevel.design_effect` disagree on real data, one of them is
wrong and that is the finding, not a nuisance.

RESOLUTION CAVEAT, read before interpreting anything below. The curve only sees
CpG pairs inside a supplied region, so it cannot resolve a decay length longer
than the typical region span. The production caller caps regions at
`max_span = 1000` bp, so this script run on the caller's own regions measures
decay only out to ~1 kb -- and if the fitted `decay_bp` comes back anywhere near
the mean region span, the true value is longer and this is a lower bound. To
measure further, re-export with wider windows rather than the caller's regions.
Measured on simulation: 400 bp regions generated with decay_bp = 500 return 88.

Stratification is not optional: decay is expected to differ between CpG-dense
and CpG-sparse regions and between high- and low-methylation regions, and a
single pooled number would hide exactly the variation that decides whether the
pooled test is safe SOMEWHERE rather than everywhere.

Usage: python3 P07_comethylation_decay.py <SAMPLE> <CHROM>
"""
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, "/u/project/cluo/terencew/claude/project_ideas/asm_lr/github/ont_asm_caller")
sys.path.insert(0, "/u/project/cluo/terencew/claude/project_ideas/asm_lr/scripts/phasing_qc")

from ont_asm_caller.decay import comethylation_decay, implied_design_effect
from ont_asm_caller.readlevel import design_effect
from P06_export_read_matrices import load_read_matrices

PROJ = "/u/project/cluo/terencew/claude/project_ideas/asm_lr"


def summarise(label, regs):
    if len(regs) < 100:
        print(f"  {label:<28} (only {len(regs)} regions -- skipped)")
        return None
    c = comethylation_decay(regs)
    spacing = float(np.mean([np.diff(r.positions).mean() for r in regs
                             if len(r.positions) > 1]))
    ncpg = float(np.mean([r.n_cpgs for r in regs]))
    implied = implied_design_effect(c, spacing, int(round(ncpg)))
    measured = float(np.median([design_effect(r.m1, r.m2) for r in regs]))
    print(f"  {label:<28}{len(regs):>7d}{c.decay_bp:>10.0f}{c.amplitude:>9.3f}"
          f"{c.r2:>7.3f}{spacing:>9.0f}{implied:>10.2f}{measured:>10.2f}")
    return {"stratum": label, "n_regions": len(regs), "decay_bp": c.decay_bp,
            "amplitude": c.amplitude, "r2": c.r2, "mean_spacing": spacing,
            "mean_n_cpgs": ncpg, "implied_DE": implied, "measured_DE_median": measured}


def main():
    sample, chrom = sys.argv[1], sys.argv[2]
    regs = list(load_read_matrices(
        f"{PROJ}/tables/read_matrices/{sample}_{chrom}.npz", chrom))
    print(f"{sample} {chrom}: {len(regs)} regions loaded\n")

    print(f"  {'stratum':<28}{'n':>7}{'decay_bp':>10}{'amp':>9}{'r2':>7}"
          f"{'spacing':>9}{'DE(fit)':>10}{'DE(meas)':>10}")
    rows = [summarise("ALL", regs)]

    mu = np.array([np.nanmean(np.concatenate([r.m1.ravel(), r.m2.ravel()]))
                   for r in regs])
    for lo, hi, name in [(0.0, 0.2, "methylation < 0.2"),
                         (0.2, 0.8, "methylation 0.2-0.8"),
                         (0.8, 1.01, "methylation > 0.8")]:
        rows.append(summarise(name, [r for r, m in zip(regs, mu)
                                     if np.isfinite(m) and lo <= m < hi]))

    dens = np.array([r.n_cpgs / max(r.positions[-1] - r.positions[0], 1) * 1000
                     for r in regs])
    q = np.nanquantile(dens, [0.33, 0.67])
    for lo, hi, name in [(-np.inf, q[0], "CpG-sparse tercile"),
                         (q[0], q[1], "CpG-medium tercile"),
                         (q[1], np.inf, "CpG-dense tercile")]:
        rows.append(summarise(name, [r for r, d in zip(regs, dens)
                                     if lo <= d < hi]))

    out = pd.DataFrame([r for r in rows if r])
    path = f"{PROJ}/tables/phasing_qc/comethylation_decay_{sample}_{chrom}.tsv"
    out.to_csv(path, sep="\t", index=False)
    print(f"\nwrote {path}")
    print("\nInterpretation:")
    print("  DE(fit) vs DE(meas) must agree -- two independent estimators of the")
    print("  same quantity. If they do not, neither number is usable yet.")
    print("  DE ~ 1  -> region.cluster_cpgs pooling is legal and more powerful.")
    print("  DE >> 1 -> it cannot control FDR; use the read-level or pattern path.")
    span = float(np.mean([r.positions[-1] - r.positions[0] for r in regs]))
    print(f"\n  mean region span = {span:.0f} bp. Any fitted decay_bp within ~2x of")
    print("  that is a LOWER BOUND -- the curve cannot see past the region.")


if __name__ == "__main__":
    main()
