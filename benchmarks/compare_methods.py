#!/usr/bin/env python3
"""
Head-to-head comparison of three approaches, under a two-layer ground-truth
generative model (existence + effect size), with CpGs spatially clustered
into regions (matching a real DMR's shape) so region-level pooling can
actually be tested:

  1. Fisher-exact + fixed delta floor + BH-FDR   (original pipeline)
  2. Beta-binomial LRT, per-CpG, BH-FDR           (this package, v1)
  3. Beta-binomial LRT, region-pooled, BH-FDR     (this package, v2)

For (3), dispersion is RE-estimated on the pooled region-level counts, not
reused from the per-CpG estimate -- pooling several independent per-CpG
beta-binomial draws reduces the effective dispersion of the combined count
(averaging out site-specific noise), so the per-CpG rho is not valid at the
region level (see ont_asm_caller/region.py and tests/test_region.py for the
direct empirical check of this).

Reports EMPIRICAL FDR (using known true-ASM labels, not just detection
counts), power, and effect-size estimation bias for each.

DSS is not included: it isn't a Python tool, and its dispersion estimation
is undefined at N=1 replicate per group (see model.py docstring and
md/20260905.progress.md) -- it hung indefinitely on real project data for
this exact reason.

Usage: python benchmarks/compare_methods.py
"""
from __future__ import annotations

import numpy as np
from scipy.stats import fisher_exact

from ont_asm_caller.dispersion import estimate_dispersion
from ont_asm_caller.model import test_locus
from ont_asm_caller.region import CpG, cluster_cpgs, test_regions
from ont_asm_caller.simulate_mixture import simulate_spatial_mixture

FISHER_DELTA_FLOOR = 0.10
FDR_THRESH = 0.05


def bh_fdr(pvals: np.ndarray) -> np.ndarray:
    n = len(pvals)
    order = np.argsort(pvals)
    ranked = pvals[order]
    q = ranked * n / (np.arange(n) + 1)
    q = np.minimum.accumulate(q[::-1])[::-1]
    out = np.empty(n)
    out[order] = np.clip(q, 0, 1)
    return out


def run_fisher(cpgs) -> tuple[np.ndarray, np.ndarray]:
    pvals = np.empty(len(cpgs))
    deltas = np.empty(len(cpgs))
    for i, c in enumerate(cpgs):
        _, p = fisher_exact([[c.x1, c.n1 - c.x1], [c.x2, c.n2 - c.x2]])
        pvals[i] = p
        deltas[i] = c.x1 / c.n1 - c.x2 / c.n2
    q = bh_fdr(pvals)
    called = (q < FDR_THRESH) & (np.abs(deltas) >= FISHER_DELTA_FLOOR)
    return called, deltas


def run_betabinom_percpg(cpgs, rho: float) -> tuple[np.ndarray, np.ndarray]:
    pvals = np.empty(len(cpgs))
    deltas = np.empty(len(cpgs))
    for i, c in enumerate(cpgs):
        r = test_locus("chr1", c.pos, c.x1, c.n1, c.x2, c.n2, rho=rho)
        pvals[i] = r.pval
        deltas[i] = r.delta
    q = bh_fdr(pvals)
    return q < FDR_THRESH, deltas


def run_betabinom_region(cpgs, max_gap: int = 400):
    """Clusters CpGs into regions and RE-ESTIMATES dispersion on the pooled
    region-level counts (see module docstring for why). Returns
    (called, deltas, region_truth, region_true_delta) -- region-level
    ground truth is tracked separately since a region isn't a 1:1 stand-in
    for a single simulated CpG."""
    region_cpgs = [CpG(pos=c.pos, x1=c.x1, n1=c.n1, x2=c.x2, n2=c.n2) for c in cpgs]
    regions = cluster_cpgs("chr1", region_cpgs, max_gap=max_gap, min_cpgs=1)

    rx1 = np.array([r.x1 for r in regions]); rn1 = np.array([r.n1 for r in regions])
    rx2 = np.array([r.x2 for r in regions]); rn2 = np.array([r.n2 for r in regions])
    rho_region_hat = estimate_dispersion(rx1, rn1, rx2, rn2)

    results = test_regions(regions, rho=rho_region_hat)

    pos_to_asm = {c.pos: c.is_true_asm for c in cpgs}
    pos_to_delta = {c.pos: c.true_delta for c in cpgs}
    region_truth, region_true_delta = [], []
    for reg in regions:
        member_pos = [p for p in range(reg.start, reg.end + 1) if p in pos_to_asm]
        is_asm = any(pos_to_asm[p] for p in member_pos)
        region_truth.append(is_asm)
        region_true_delta.append(max((pos_to_delta[p] for p in member_pos), default=0.0))
    region_truth = np.array(region_truth)
    region_true_delta = np.array(region_true_delta)

    pvals = np.array([r.pval for r in results])
    deltas = np.array([r.delta for r in results])
    q = bh_fdr(pvals)
    return q < FDR_THRESH, deltas, region_truth, region_true_delta, rho_region_hat


def summarize(name, called, deltas, truth_asm, truth_delta):
    n_called = called.sum()
    n_true_pos = (called & truth_asm).sum()
    n_false_pos = (called & ~truth_asm).sum()
    empirical_fdr = n_false_pos / n_called if n_called > 0 else float("nan")
    power = n_true_pos / truth_asm.sum() if truth_asm.sum() > 0 else float("nan")

    if n_true_pos > 0:
        bias = np.mean(np.abs(deltas[called & truth_asm]) - truth_delta[called & truth_asm])
    else:
        bias = float("nan")

    print(f"{name:22} n_called={n_called:6d}  empirical_FDR={empirical_fdr:6.3f}  "
          f"power={power:6.3f}  effect-size bias={bias:+.3f}")


def main():
    rho_true = 0.05
    pi = 0.005
    n_regions = 15000  # ~4-8 CpGs/region -> comparable total CpG count to before

    print(f"=== Simulated cohort: n_regions={n_regions}, true ASM prevalence pi={pi}, "
          f"true dispersion rho={rho_true}, CpGs spatially clustered per region ===\n")

    for depth_label, depth in [
        ("uniform high depth (18x)", 18),
        ("uniform low depth (8x)", 8),
        ("realistic heterogeneous (5-40x)", (5, 40)),
    ]:
        print(f"--- depth = {depth_label} ---")
        cpgs = simulate_spatial_mixture(
            n_regions=n_regions, pi=pi, cpgs_per_region=(4, 8),
            depth=depth, rho=rho_true, seed=2024,
        )
        truth_asm_cpg = np.array([c.is_true_asm for c in cpgs])
        truth_delta_cpg = np.array([c.true_delta for c in cpgs])
        print(f"  ({len(cpgs)} CpGs total, {truth_asm_cpg.sum()} true-ASM CpGs)")

        fisher_called, fisher_deltas = run_fisher(cpgs)
        summarize("Fisher (orig.)", fisher_called, fisher_deltas, truth_asm_cpg, truth_delta_cpg)

        x1 = np.array([c.x1 for c in cpgs]); n1 = np.array([c.n1 for c in cpgs])
        x2 = np.array([c.x2 for c in cpgs]); n2 = np.array([c.n2 for c in cpgs])
        rho_cpg_hat = estimate_dispersion(x1, n1, x2, n2)
        bb_called, bb_deltas = run_betabinom_percpg(cpgs, rho=rho_cpg_hat)
        summarize(f"BetaBinom per-CpG", bb_called, bb_deltas, truth_asm_cpg, truth_delta_cpg)

        region_called, region_deltas, region_truth, region_true_delta, rho_region_hat = \
            run_betabinom_region(cpgs)
        summarize(f"BetaBinom region (rho_hat={rho_region_hat:.3f})",
                   region_called, region_deltas, region_truth, region_true_delta)
        print()


if __name__ == "__main__":
    main()
