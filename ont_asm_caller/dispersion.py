"""
Global (genome-wide) beta-binomial overdispersion estimation.

Why global, not per-locus / per-replicate: DSS-style tools estimate
dispersion from between-REPLICATE variance at (or smoothed near) each site.
There are no replicates here -- HP1 and HP2 are the two halves of one
individual's reads, not independent biological samples, so that quantity is
undefined (see model.py docstring and md/20260905.progress.md for the fuller
argument, and the observation that DSS's dispersion-estimation step hung
indefinitely on this exact project's data, plausibly for this exact reason).

What we DO have, in abundance, is loci: millions of CpGs genome-wide. Under
the standard (and here, defensible) assumption that the large majority of
loci show no true ASM, HP1 and HP2 read counts at a given locus are, under
the null, two independent binomial draws from a SHARED rate -- i.e. they are
exchangeable in exactly the way replicates would be, just found across the
genome instead of across biological samples. The excess variance beyond pure
binomial sampling, pooled across that large set of loci, gives a
well-defined, estimable global dispersion parameter.
"""
from __future__ import annotations

import numpy as np


def estimate_dispersion(x1, n1, x2, n2, min_n: int = 5) -> float:
    """Method-of-moments estimate of the global beta-binomial dispersion
    `rho`, from the excess of a per-locus Pearson chi-square statistic
    (assuming a shared rate per locus, 1 df) over its binomial expectation.

    For locus i, with p_hat_i = (x1_i + x2_i) / (n1_i + n2_i):

        chi2_i = (x1_i - n1_i*p_hat_i)^2 / (n1_i*p_hat_i*(1-p_hat_i))
               + (x2_i - n2_i*p_hat_i)^2 / (n2_i*p_hat_i*(1-p_hat_i))

    E[chi2_i] = 1 under pure binomial sampling (Pearson stat, 1 df).
    E[chi2_i] ~= 1 + rho*(n_eff_i - 1) under overdispersion (standard
    beta-binomial variance-inflation form), where n_eff_i = (n1_i+n2_i)/2.
    `rho` is fit by (weights-free) least squares of (chi2_i - 1) against
    (n_eff_i - 1), forced through the origin.

    This deliberately INCLUDES the (rare) true-ASM loci in the pool, same as
    standard practice in DSS/edgeR-style genome-wide dispersion estimation --
    with true ASM expected at well under 1% of tested loci (c.f. deCODE's
    0.51% validated rate, md/20260905.progress.md), the bias this introduces
    is negligible; see tests/test_dispersion.py for a direct simulation check
    of that claim.
    """
    x1 = np.asarray(x1, dtype=float)
    n1 = np.asarray(n1, dtype=float)
    x2 = np.asarray(x2, dtype=float)
    n2 = np.asarray(n2, dtype=float)
    if not (len(x1) == len(n1) == len(x2) == len(n2)):
        raise ValueError("x1, n1, x2, n2 must all be the same length")

    keep = (n1 >= min_n) & (n2 >= min_n)
    x1, n1, x2, n2 = x1[keep], n1[keep], x2[keep], n2[keep]
    if len(x1) < 100:
        raise ValueError(
            f"only {len(x1)} loci with n>={min_n} on both haplotypes -- "
            "too few to estimate a global dispersion reliably (need >=100)"
        )

    p_hat = (x1 + x2) / (n1 + n2)
    p_hat = np.clip(p_hat, 1e-6, 1 - 1e-6)

    var1 = n1 * p_hat * (1 - p_hat)
    var2 = n2 * p_hat * (1 - p_hat)
    chi2 = (x1 - n1 * p_hat) ** 2 / var1 + (x2 - n2 * p_hat) ** 2 / var2

    n_eff = (n1 + n2) / 2.0
    y = chi2 - 1.0
    xw = n_eff - 1.0

    valid = xw > 0
    if valid.sum() < 100:
        return 0.0

    rho_hat = float(np.sum(xw[valid] * y[valid]) / np.sum(xw[valid] ** 2))
    return float(np.clip(rho_hat, 0.0, 0.999))


# ---------------------------------------------------------------------------
# Mean-dependent dispersion, calibrated by simulation.
#
# `estimate_dispersion` above assumes E[chi2_i] = 1 under pure binomial
# sampling. That is an ASYMPTOTIC result and it fails badly at the boundaries:
# when x1 = n1 and x2 = n2 (both haplotypes fully methylated -- the single most
# common configuration in a real bimodal methylome) p_hat -> 1, the variance
# term collapses, chi2_i -> 0, and the locus contributes y = -1, i.e. positive
# evidence for rho = 0. Measured bias on a realistic methylome is 2-3.5x too
# LOW, which makes the downstream test anti-conservative.
#
# Fix: do not assume E[chi2] = 1. Estimate it by simulating binomial data at
# the OBSERVED (n1, n2, p_hat) triples, which reproduces the same boundary
# degeneracy, and read rho off the excess over that simulated baseline. Also
# fit rho separately in bins of p_hat, because beta-binomial dispersion is
# structurally mean-dependent -- bounded near 0 and 1, largest near 0.5. A
# single global rho is simultaneously too conservative at the boundaries
# (where imprinting-style ASM lives) and too liberal in the middle.
# ---------------------------------------------------------------------------

def _mean_chi2(x1, n1, x2, n2):
    p = np.clip((x1 + x2) / (n1 + n2), 1e-9, 1 - 1e-9)
    v1, v2 = n1 * p * (1 - p), n2 * p * (1 - p)
    c = (x1 - n1 * p) ** 2 / v1 + (x2 - n2 * p) ** 2 / v2
    return float(np.mean(c))


def _simulate_mean_chi2(rng, n1, n2, p, rho, reps=1):
    """Mean Pearson chi2 for data generated at dispersion `rho` with the given
    depths and rates. rho=0 gives the exact finite-sample binomial baseline."""
    n1 = np.asarray(n1).astype(np.int64)
    n2 = np.asarray(n2).astype(np.int64)
    # beta shape params blow up at p exactly 0 or 1; clip only for simulation
    p = np.clip(np.asarray(p, dtype=float), 1e-6, 1 - 1e-6)
    out = []
    for _ in range(reps):
        if rho <= 0:
            x1 = rng.binomial(n1, p); x2 = rng.binomial(n2, p)
        else:
            a = p * (1 - rho) / rho
            b = (1 - p) * (1 - rho) / rho
            x1 = rng.binomial(n1, rng.beta(a, b))
            x2 = rng.binomial(n2, rng.beta(a, b))
        out.append(_mean_chi2(x1, n1, x2, n2))
    return float(np.mean(out))


def estimate_dispersion_trend(x1, n1, x2, n2, min_n: int = 5, n_bins: int = 8,
                               grid=None, seed: int = 0):
    """Mean-dependent dispersion rho(mu), calibrated by simulation.

    Returns a callable rho(mu) (vectorised, clipped to the fitted range), plus
    the fitted knots as `.knots` = (mu_centres, rho_values) for inspection and
    plotting -- always look at this curve before trusting a run.

    Within each bin of p_hat, rho is found by matching the observed mean
    Pearson chi2 to simulated means over a grid of candidate rho, using the
    bin's own (n1, n2, p_hat) so finite-sample and boundary effects cancel.
    """
    x1, n1, x2, n2 = (np.asarray(v, dtype=float) for v in (x1, n1, x2, n2))
    keep = (n1 >= min_n) & (n2 >= min_n)
    x1, n1, x2, n2 = x1[keep], n1[keep], x2[keep], n2[keep]
    if len(x1) < 500:
        raise ValueError(f"only {len(x1)} usable loci -- need >=500 for a trend fit")

    if grid is None:
        grid = np.concatenate([[0.0], np.geomspace(0.002, 0.6, 24)])
    rng = np.random.default_rng(seed)

    p_hat = np.clip((x1 + x2) / (n1 + n2), 0.0, 1.0)
    edges = np.quantile(p_hat, np.linspace(0, 1, n_bins + 1))
    edges[0], edges[-1] = -1e-9, 1 + 1e-9

    centres, rhos = [], []
    for i in range(n_bins):
        m = (p_hat > edges[i]) & (p_hat <= edges[i + 1])
        if m.sum() < 200:
            continue
        obs = _mean_chi2(x1[m], n1[m], x2[m], n2[m])
        sims = np.array([_simulate_mean_chi2(rng, n1[m], n2[m], p_hat[m], r)
                         for r in grid])
        # monotone in rho; take the closest match
        rhos.append(float(grid[int(np.argmin(np.abs(sims - obs)))]))
        centres.append(float(np.median(p_hat[m])))

    if not centres:
        raise ValueError("no bin had enough loci to fit a dispersion trend")

    centres = np.array(centres); rhos = np.array(rhos)
    order = np.argsort(centres)
    centres, rhos = centres[order], rhos[order]

    def rho_of_mu(mu):
        mu = np.clip(np.asarray(mu, dtype=float), centres[0], centres[-1])
        return np.interp(mu, centres, rhos)

    rho_of_mu.knots = (centres, rhos)
    return rho_of_mu
