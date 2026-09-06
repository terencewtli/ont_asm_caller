"""
Simulate synthetic haplotype-partitioned methylation counts under a known
ground truth, for validating the ASM test's calibration (false-positive
control) and power before trusting it on real data -- and for directly
testing hypotheses about what could have produced this project's earlier
anomalous real-data results (depth confound, phasing/switch errors).
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np


def _rbetabinom(rng: np.random.Generator, n: int, mu: float, rho: float) -> int:
    """Sample from a beta-binomial by drawing p ~ Beta(mu, rho) then
    x ~ Binomial(n, p)."""
    rho = max(rho, 1e-10)
    a = mu * (1.0 - rho) / rho
    b = (1.0 - mu) * (1.0 - rho) / rho
    p = rng.beta(a, b)
    return int(rng.binomial(n, p))


@dataclass
class SimulatedLocus:
    n1: int
    x1: int
    n2: int
    x2: int
    true_mu1: float
    true_mu2: float
    true_delta: float


def simulate_loci(
    n_loci: int,
    depth=20,
    true_delta: float = 0.0,
    baseline_mu: float = 0.5,
    rho: float = 0.05,
    switch_error_rate: float = 0.0,
    seed: int | None = None,
) -> list[SimulatedLocus]:
    """Simulate `n_loci` independent loci.

    depth: fixed int (same N for every locus/haplotype), or (low, high) for
        per-locus-per-haplotype depth drawn uniformly at random -- use a
        wide range (e.g. (5, 40)) to mimic the real cross-sample depth
        heterogeneity this project's pipeline actually has to handle.
    true_delta: |mu1 - mu2|, the real per-locus effect size. Use 0.0 to test
        calibration (false-positive rate under the null); >0 to test power.
    baseline_mu: the average of mu1, mu2 (mu1 = baseline + delta/2,
        mu2 = baseline - delta/2).
    rho: true beta-binomial dispersion used to generate the data -- pass
        different train/test values of rho to model.py's test to check
        sensitivity to dispersion misestimation, not just correct estimation.
    switch_error_rate: probability (0-1) that an individual READ's true
        haplotype origin differs from the haplotype it gets counted under --
        simulates imperfect phasing (switch errors), which no generic DMR
        tool models but is directly relevant here (see
        md/20260905.progress.md). At switch_error_rate=eps, a read
        originally on haplotype h has probability eps of being counted on
        the other haplotype instead. This ATTENUATES the observed delta
        toward zero (expected factor ~ (1 - 2*eps) for small eps on the
        counts, before binomial sampling noise) -- it does not, by itself,
        manufacture large spurious differences at null loci; it just makes
        real differences harder to see. See tests/test_simulate.py for a
        direct check of this attenuation behavior.
    """
    rng = np.random.default_rng(seed)
    loci = []
    for _ in range(n_loci):
        if isinstance(depth, tuple):
            n1 = int(rng.integers(depth[0], depth[1] + 1))
            n2 = int(rng.integers(depth[0], depth[1] + 1))
        else:
            n1 = n2 = int(depth)

        mu1 = float(np.clip(baseline_mu + true_delta / 2, 0.001, 0.999))
        mu2 = float(np.clip(baseline_mu - true_delta / 2, 0.001, 0.999))

        if switch_error_rate > 0:
            eps = switch_error_rate
            # Reallocate READS (not just methylation calls) between haplotypes
            # first, then determine each read's methylation status by its
            # TRUE haplotype's rate -- this is what a switch error actually
            # does physically (a read is assigned to the wrong haplotype bin
            # entirely, carrying whatever true methylation state it has).
            keep1 = int(rng.binomial(n1, 1 - eps))
            swap1 = n1 - keep1
            keep2 = int(rng.binomial(n2, 1 - eps))
            swap2 = n2 - keep2

            x1 = (int(rng.binomial(keep1, mu1)) if keep1 > 0 else 0) + \
                 (int(rng.binomial(swap2, mu2)) if swap2 > 0 else 0)
            x2 = (int(rng.binomial(keep2, mu2)) if keep2 > 0 else 0) + \
                 (int(rng.binomial(swap1, mu1)) if swap1 > 0 else 0)
            n1_obs, n2_obs = keep1 + swap2, keep2 + swap1
        else:
            n1_obs, n2_obs = n1, n2
            x1 = _rbetabinom(rng, n1, mu1, rho)
            x2 = _rbetabinom(rng, n2, mu2, rho)

        loci.append(SimulatedLocus(
            n1=n1_obs, x1=x1, n2=n2_obs, x2=x2,
            true_mu1=mu1, true_mu2=mu2, true_delta=abs(mu1 - mu2),
        ))
    return loci
