"""
Two-layer generative model for benchmark comparisons (Fisher vs. this
package's beta-binomial test vs., where feasible, DSS): each locus
independently gets a TRUE ASM/not-ASM label, and an effect size conditional
on that label. This is what makes an actual FDR/power comparison possible --
`simulate.simulate_loci` applies one fixed delta to every locus, which is
fine for calibration/power curves at a chosen effect size, but has no
ground-truth labels to check a method's FDR against.

Layer 1 (existence):   z_l ~ Bernoulli(pi)          pi = true ASM prevalence
Layer 2 (magnitude):    delta_l | z_l=1 ~ effect_size_dist
                        delta_l | z_l=0 = 0
Then counts are generated exactly as in simulate.simulate_loci, given
depth, baseline_mu, rho, and switch_error_rate.

pi=0.005 (deCODE's validated genome-wide ASM-QTL rate, see
md/20260905.progress.md) is used as the default "realistic" base rate.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from .simulate import _rbetabinom


@dataclass
class MixtureLocus:
    n1: int
    x1: int
    n2: int
    x2: int
    is_true_asm: bool
    true_delta: float


def simulate_mixture(
    n_loci: int,
    pi: float = 0.005,
    depth=(5, 40),
    baseline_mu: float = 0.5,
    rho: float = 0.05,
    effect_size_dist="uniform_0.2_0.6",
    switch_error_rate: float = 0.0,
    seed: int | None = None,
) -> list[MixtureLocus]:
    """effect_size_dist: currently supports "uniform_LO_HI" (true |delta| for
    ASM loci drawn uniformly in [LO, HI]) -- deliberately not concentrated
    at one value, since real ASM spans a range from meQTL-scale effects to
    near-complete imprinting-scale separation."""
    rng = np.random.default_rng(seed)

    if isinstance(effect_size_dist, str) and effect_size_dist.startswith("uniform_"):
        _, lo, hi = effect_size_dist.split("_")
        lo, hi = float(lo), float(hi)
        draw_delta = lambda: rng.uniform(lo, hi)
    else:
        raise ValueError(f"unsupported effect_size_dist: {effect_size_dist}")

    loci = []
    for _ in range(n_loci):
        is_asm = bool(rng.random() < pi)
        delta = draw_delta() if is_asm else 0.0

        if isinstance(depth, tuple):
            n1 = int(rng.integers(depth[0], depth[1] + 1))
            n2 = int(rng.integers(depth[0], depth[1] + 1))
        else:
            n1 = n2 = int(depth)

        mu1 = float(np.clip(baseline_mu + delta / 2, 0.001, 0.999))
        mu2 = float(np.clip(baseline_mu - delta / 2, 0.001, 0.999))

        if switch_error_rate > 0:
            eps = switch_error_rate
            keep1 = int(rng.binomial(n1, 1 - eps)); swap1 = n1 - keep1
            keep2 = int(rng.binomial(n2, 1 - eps)); swap2 = n2 - keep2
            x1 = (int(rng.binomial(keep1, mu1)) if keep1 > 0 else 0) + \
                 (int(rng.binomial(swap2, mu2)) if swap2 > 0 else 0)
            x2 = (int(rng.binomial(keep2, mu2)) if keep2 > 0 else 0) + \
                 (int(rng.binomial(swap1, mu1)) if swap1 > 0 else 0)
            n1_obs, n2_obs = keep1 + swap2, keep2 + swap1
        else:
            n1_obs, n2_obs = n1, n2
            x1 = _rbetabinom(rng, n1, mu1, rho)
            x2 = _rbetabinom(rng, n2, mu2, rho)

        loci.append(MixtureLocus(
            n1=n1_obs, x1=x1, n2=n2_obs, x2=x2,
            is_true_asm=is_asm, true_delta=delta,
        ))
    return loci


@dataclass
class SpatialCpG:
    pos: int
    x1: int
    n1: int
    x2: int
    n2: int
    region_id: int
    is_true_asm: bool
    true_delta: float


def simulate_spatial_mixture(
    n_regions: int,
    pi: float = 0.005,
    cpgs_per_region=(3, 8),
    intra_region_span: int = 400,
    region_spacing: int = 5000,
    depth=(5, 40),
    baseline_mu: float = 0.5,
    rho: float = 0.05,
    effect_size_dist="uniform_0.2_0.6",
    seed: int | None = None,
) -> list[SpatialCpG]:
    """Like simulate_mixture, but each "locus" is a REGION spanning several
    CpGs within `intra_region_span` bp of each other (matching a real DMR),
    with ONE true ASM status and effect size shared by the whole region --
    but each CpG still gets its OWN independent beta-binomial draw (its own
    site-specific rho-driven noise), not a copy of the same count. This is
    what makes it a fair test of region.cluster_cpgs: pooling only helps if
    there's real per-CpG noise to average out, which this simulates
    honestly rather than assuming away.

    `region_spacing` (average gap between region starts) must be well above
    `intra_region_span` so that cluster_cpgs(max_gap=intra_region_span-ish)
    correctly separates regions instead of merging neighbors together --
    default 5000 vs. 400 gives ample margin.
    """
    rng = np.random.default_rng(seed)

    if isinstance(effect_size_dist, str) and effect_size_dist.startswith("uniform_"):
        _, lo, hi = effect_size_dist.split("_")
        lo, hi = float(lo), float(hi)
        draw_delta = lambda: rng.uniform(lo, hi)
    else:
        raise ValueError(f"unsupported effect_size_dist: {effect_size_dist}")

    cpgs = []
    pos = 1000
    for region_id in range(n_regions):
        is_asm = bool(rng.random() < pi)
        delta = draw_delta() if is_asm else 0.0
        mu1 = float(np.clip(baseline_mu + delta / 2, 0.001, 0.999))
        mu2 = float(np.clip(baseline_mu - delta / 2, 0.001, 0.999))

        n_cpgs = int(rng.integers(cpgs_per_region[0], cpgs_per_region[1] + 1))
        offsets = sorted(rng.integers(0, intra_region_span, size=n_cpgs))

        for off in offsets:
            if isinstance(depth, tuple):
                n1 = int(rng.integers(depth[0], depth[1] + 1))
                n2 = int(rng.integers(depth[0], depth[1] + 1))
            else:
                n1 = n2 = int(depth)
            x1 = _rbetabinom(rng, n1, mu1, rho)
            x2 = _rbetabinom(rng, n2, mu2, rho)
            cpgs.append(SpatialCpG(
                pos=pos + int(off), x1=x1, n1=n1, x2=x2, n2=n2,
                region_id=region_id, is_true_asm=is_asm, true_delta=delta,
            ))

        pos += region_spacing

    return cpgs
