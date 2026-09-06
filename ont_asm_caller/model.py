"""
Beta-binomial two-haplotype ASM test.

Model
-----
At each candidate locus l and haplotype h in {1, 2}:

    X_{l,h} | N_{l,h}, mu_{l,h}  ~  BetaBinomial(N_{l,h}, mu_{l,h}, rho)

`rho` (0 <= rho < 1) is the beta-binomial overdispersion / intraclass
correlation parameter. It is a SINGLE nuisance value shared across all loci
for a given sample (estimated once, genome-wide -- see dispersion.py), not
estimated per locus from between-replicate variance. There is no replicate
here: a haplotype is one half of one individual's reads, not an independent
biological sample. Pooling for the dispersion estimate is done across LOCI
instead of across replicates, which is well-defined because a sample has
millions of testable CpGs even though it has exactly one observation of
each haplotype at each one.

Test
----
For each locus, a likelihood-ratio test of

    H0: mu_1 == mu_2   (no allele-specific methylation)
    H1: mu_1 != mu_2

with rho fixed at its pre-estimated global value. LR = 2*(llH1 - llH0) is
asymptotically chi-square(1 df) under H0. This is the direct fix for the
Fisher-exact-test depth confound found in this project's original pipeline
(see md/20260905.progress.md): depth (N) enters the likelihood directly, so
the same procedure and threshold can be applied to samples with very
different coverage without conflating "low power" with "no difference".
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy.optimize import minimize_scalar
from scipy.stats import betabinom, chi2

_EPS = 1e-6


def _bb_ab(mu: float, rho: float) -> tuple[float, float]:
    """Convert (mean, dispersion) to the (a, b) shape params scipy's
    betabinom uses. rho -> 0 recovers a plain binomial (a, b -> infinity
    with a/(a+b) = mu held fixed)."""
    rho = max(rho, 1e-10)
    a = mu * (1.0 - rho) / rho
    b = (1.0 - mu) * (1.0 - rho) / rho
    return a, b


def bb_logpmf(x: np.ndarray, n: np.ndarray, mu: float, rho: float) -> np.ndarray:
    """Beta-binomial log-pmf, vectorized over loci, at a single (mu, rho)."""
    a, b = _bb_ab(mu, rho)
    return betabinom.logpmf(x, n, a, b)


def _neg_ll(mu: float, x: np.ndarray, n: np.ndarray, rho: float) -> float:
    if mu <= 0.0 or mu >= 1.0:
        return np.inf
    ll = bb_logpmf(x, n, mu, rho).sum()
    return -ll if np.isfinite(ll) else np.inf


def fit_mu(x, n, rho: float) -> tuple[float, float]:
    """MLE of a single shared mean `mu` for the given (x, n) pairs at fixed
    `rho`. Returns (mu_hat, maximized log-likelihood). Accepts scalars or
    arrays for x/n (arrays let H0's "shared mu across both haplotypes" fit
    be expressed as one call over [x1, x2], [n1, n2])."""
    x = np.atleast_1d(np.asarray(x, dtype=float))
    n = np.atleast_1d(np.asarray(n, dtype=float))
    if np.any(n <= 0):
        raise ValueError("all elements of n must be > 0")

    res = minimize_scalar(
        _neg_ll, args=(x, n, rho),
        bounds=(_EPS, 1 - _EPS), method="bounded",
        options={"xatol": 1e-8},
    )
    return float(res.x), float(-res.fun)


@dataclass
class ASMResult:
    chrom: str
    pos: int
    n1: int
    x1: int
    n2: int
    x2: int
    mu1_hat: float
    mu2_hat: float
    mu_shared_hat: float
    delta: float
    lr_stat: float
    pval: float


def test_locus(chrom: str, pos: int, x1: int, n1: int, x2: int, n2: int,
                rho: float) -> ASMResult:
    """Likelihood-ratio test for HP1 vs HP2 methylation-rate equality at one
    locus, given a pre-estimated global dispersion `rho`."""
    if n1 <= 0 or n2 <= 0:
        raise ValueError("both haplotypes need at least one read")

    mu1_hat, ll1 = fit_mu(x1, n1, rho)
    mu2_hat, ll2 = fit_mu(x2, n2, rho)
    ll_h1 = ll1 + ll2

    mu_shared_hat, ll_h0 = fit_mu(np.array([x1, x2]), np.array([n1, n2]), rho)

    lr_stat = max(0.0, 2.0 * (ll_h1 - ll_h0))
    pval = float(chi2.sf(lr_stat, df=1))

    return ASMResult(
        chrom=chrom, pos=pos, n1=int(n1), x1=int(x1), n2=int(n2), x2=int(x2),
        mu1_hat=mu1_hat, mu2_hat=mu2_hat, mu_shared_hat=mu_shared_hat,
        delta=mu1_hat - mu2_hat, lr_stat=lr_stat, pval=pval,
    )


def test_loci_vectorized(chrom: np.ndarray, pos: np.ndarray,
                          x1: np.ndarray, n1: np.ndarray,
                          x2: np.ndarray, n2: np.ndarray,
                          rho: float) -> list[ASMResult]:
    """Convenience wrapper: run test_locus over parallel arrays. Not
    numerically vectorized (each locus still does its own 1-D optimization),
    but avoids Python-level boilerplate at call sites. For genome-scale runs,
    parallelize calls to this across chunks (e.g. multiprocessing.Pool),
    same pattern as the project's existing W04_dasm_loci.py."""
    return [
        test_locus(c, p, int(a), int(b), int(cc), int(d), rho)
        for c, p, a, b, cc, d in zip(chrom, pos, x1, n1, x2, n2)
    ]
