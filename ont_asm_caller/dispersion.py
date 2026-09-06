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
