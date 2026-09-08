"""
Read-pattern simulator: the two axes of allele-specific methylation.

WHY THIS EXISTS (and why simulate_reads.py is not enough)
---------------------------------------------------------
`simulate_reads.simulate_read_regions` was written to fix the one defect it was
aimed at -- reads shared across the CpGs of a region -- and it does. But it is
NOT a usable testbed for comparing against a pattern-based method such as CPEL
(Abante, Fang, Feinberg & Goutsias, Nat Commun 2020), for two measured reasons.

1. ITS CORRELATION IS EXCHANGEABLE, NOT DISTANCE-DECAYING.
   `_molecules` copies each CpG's state from a single molecule-level latent with
   probability q, so

       corr(x_i, x_j) = q^2   for EVERY pair (i, j), independent of distance.

   Measured, at q = 0.6, lag 1..11: 0.361 0.360 0.360 0.359 ... -- flat.
   A 1-D Ising chain (what CPEL assumes, and what processivity/methylation-LD
   decay actually look like) gives 0.600 0.364 0.218 0.131 0.076 ... -- geometric
   decay. The parent project's co-methylation-decay analysis measures the DECAY;
   this simulator has none to measure. Under exchangeability the per-read
   methylation fraction is a sufficient statistic for the read, so a
   pattern-based test has, by construction, nothing to find that the existing
   read-fraction test does not already find. Benchmarking CPEL on it is rigged.

2. IT HAS ONLY ONE ASM AXIS.
   Both haplotypes get the same `co_methylation`, and differ only in mean. So
   the entropy/disorder imbalance that CPEL's NME and PDM statistics exist to
   detect has, by construction, zero true positives. A benchmark on this data
   can only ask "does CPEL match a mean test on mean differences" -- exactly the
   axis where it claims no advantage.

WHAT THIS MODULE ADDS
---------------------
Two orthogonal ASM axes, either or both of which can be switched on:

  delta_mml : mean methylation imbalance     (what every test in this package
                                              already targets)
  delta_icc : within-read CORRELATION imbalance at EQUAL MEAN -- one haplotype
              a fixed/ordered pattern, the other bistable across molecules.
              This is polymorphic/stochastic ASM. It is invisible to any test
              that compares locations, which is every test in this package.

and three generative models, deliberately including one that neither method
assumes, because benchmarking on a simulator drawn from a method's own model
tells you nothing:

  "exchangeable" -- the simulate_reads.py copy-latent model. FAVOURS the
                    read-fraction test (fraction is sufficient).
  "markov"       -- distance-decaying nearest-neighbour copying. This is a
                    1-D Ising chain reparameterised, so it FAVOURS CPEL.
  "epiallele"    -- a mixture of a few discrete latent methylation patterns with
                    haplotype-specific weights plus per-site flip noise. Neither
                    method's model. Use this one for the headline comparison and
                    the other two as the two-sided sensitivity check.

`decay_bp` is the correlation half-life along the molecule, in base pairs, and
is the parameter the parent project's `md/20260903_qc_review.md` s7.3 decay
curve actually measures. Nothing here should be trusted at a decay_bp that
measurement has not confirmed.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np


# ---------------------------------------------------------------------------
# molecule generators
# ---------------------------------------------------------------------------

def molecules_exchangeable(rng, n_reads: int, n_cpgs: int, mu: float,
                           icc: float) -> np.ndarray:
    """Copy-latent model, identical to simulate_reads._molecules but taking the
    ICC directly rather than its square root. corr(x_i, x_j) = icc at all lags."""
    q = float(np.sqrt(np.clip(icc, 0.0, 1.0)))
    latent = rng.random(n_reads) < mu
    indep = rng.random((n_reads, n_cpgs)) < mu
    copy = rng.random((n_reads, n_cpgs)) < q
    return np.where(copy, latent[:, None], indep).astype(float)


def molecules_markov(rng, n_reads: int, positions: np.ndarray, mu: float,
                     icc: float, decay_bp: float = 1000.0) -> np.ndarray:
    """Distance-decaying nearest-neighbour copying along the molecule.

    Site n copies site n-1 with probability r_n = icc * exp(-d_n / decay_bp),
    otherwise draws Bernoulli(mu). The marginal is exactly mu at every site and

        corr(x_n, x_m) = prod_{k=n+1..m} r_k

    i.e. geometric decay in CpG index and exponential decay in base pairs. This
    is an inhomogeneous 1-D Ising chain in a different parameterisation, so a
    method that assumes an Ising chain is correctly specified here.
    """
    positions = np.asarray(positions, dtype=float)
    n_cpgs = len(positions)
    x = np.empty((n_reads, n_cpgs))
    x[:, 0] = rng.random(n_reads) < mu
    d = np.diff(positions)
    r = np.clip(icc, 0.0, 1.0) * np.exp(-d / max(decay_bp, 1e-9))
    for j in range(1, n_cpgs):
        keep = rng.random(n_reads) < r[j - 1]
        fresh = rng.random(n_reads) < mu
        x[:, j] = np.where(keep, x[:, j - 1], fresh)
    return x


def molecules_epiallele(rng, n_reads: int, n_cpgs: int, patterns: np.ndarray,
                        weights: np.ndarray, flip_rate: float = 0.05
                        ) -> np.ndarray:
    """Mixture of discrete latent methylation patterns.

    `patterns` is (n_patterns, n_cpgs) of 0/1; `weights` sums to 1. Each read
    draws a pattern, then each site flips independently with `flip_rate`
    (basecalling error + genuine site-level stochasticity). Neither the
    exchangeable nor the Ising model is correct for this, which is the point.
    """
    patterns = np.asarray(patterns, dtype=float)
    w = np.asarray(weights, dtype=float)
    w = w / w.sum()
    idx = rng.choice(len(w), size=n_reads, p=w)
    x = patterns[idx].copy()
    flip = rng.random(x.shape) < flip_rate
    return np.where(flip, 1.0 - x, x)


def _epiallele_pair(rng, n_cpgs: int, mu: float, icc: float):
    """Patterns + weights realising mean `mu` and between-read correlation `icc`
    as a two-epiallele mixture: a `mu`-fraction-methylated ORDERED pattern
    (low disorder) mixed against the all-0/all-1 pair (high disorder).

    icc -> 0 : every read carries the same ordered pattern; reads agree with
               each other only as much as the pattern forces. Fraction is
               near-constant at mu across reads.
    icc -> 1 : reads are all-methylated or all-unmethylated in proportion mu.
               Fraction is bimodal at 0/1. Same mean, far more disorder.
    """
    k = int(round(mu * n_cpgs))
    ordered = np.zeros(n_cpgs)
    if k > 0:
        # spread the methylated sites out rather than blocking them, so the
        # ordered pattern is genuinely low-correlation along the molecule
        ordered[np.linspace(0, n_cpgs - 1, k).round().astype(int)] = 1.0
    hi, lo = np.ones(n_cpgs), np.zeros(n_cpgs)
    patterns = np.stack([ordered, hi, lo])
    a = float(np.clip(icc, 0.0, 1.0))
    weights = np.array([1.0 - a, a * mu, a * (1.0 - mu)])
    return patterns, weights


# ---------------------------------------------------------------------------
# region-level simulator with two orthogonal ASM axes
# ---------------------------------------------------------------------------

@dataclass
class PatternRegion:
    chrom: str
    start: int
    positions: np.ndarray
    m1: np.ndarray
    m2: np.ndarray
    is_true_asm: bool = False
    true_delta_mml: float = 0.0
    true_delta_icc: float = 0.0
    kind: str = "epiallele"

    @property
    def n_cpgs(self) -> int:
        return len(self.positions)

    def counts(self):
        """Pooled read x CpG counts -- what region.cluster_cpgs produces."""
        return (int(np.nansum(self.m1)), int(np.isfinite(self.m1).sum()),
                int(np.nansum(self.m2)), int(np.isfinite(self.m2).sum()))


def _draw_molecules(rng, kind, n_reads, positions, mu, icc, decay_bp):
    n_cpgs = len(positions)
    if kind == "exchangeable":
        return molecules_exchangeable(rng, n_reads, n_cpgs, mu, icc)
    if kind == "markov":
        return molecules_markov(rng, n_reads, positions, mu, icc, decay_bp)
    if kind == "epiallele":
        pat, w = _epiallele_pair(rng, n_cpgs, mu, icc)
        return molecules_epiallele(rng, n_reads, n_cpgs, pat, w)
    raise ValueError(f"unknown kind: {kind}")


def simulate_pattern_regions(
    n_regions: int,
    kind: str = "epiallele",
    pi_mml: float = 0.005,
    pi_icc: float = 0.005,
    cpgs_per_region=(4, 10),
    intra_region_span: int = 400,
    region_spacing: int = 5000,
    reads_per_hap=(10, 30),
    baseline_mu="bimodal",
    baseline_icc: float = 0.5,
    decay_bp: float = 1000.0,
    mml_effect=(0.2, 0.6),
    icc_effect=(0.4, 0.9),
    no_call_rate: float = 0.02,
    seed: int | None = None,
) -> list[PatternRegion]:
    """Regions carrying mean-imbalance ASM, correlation-imbalance ASM, or both.

    `pi_mml` and `pi_icc` are the independent prevalences of the two axes, so a
    region can be a true positive on one, the other, or both. Set `pi_icc=0` to
    reproduce the single-axis world simulate_reads.py lives in.

    Correlation-imbalance regions are generated at a mean HELD EQUAL between
    haplotypes (`mu1 == mu2`), so they are true positives that carry exactly
    zero mean signal. That is deliberate: it is the only way to find out whether
    a mean-based caller is blind to them or merely underpowered.
    """
    from .simulate_reads import draw_baseline_mu
    rng = np.random.default_rng(seed)
    out, pos = [], 1000

    for _ in range(n_regions):
        has_mml = bool(rng.random() < pi_mml)
        has_icc = bool(rng.random() < pi_icc)

        base_mu = draw_baseline_mu(rng, baseline_mu)
        d_mml = 0.0
        if has_mml:
            d_mml = float(rng.uniform(*mml_effect)) * (1 if rng.random() < 0.5 else -1)
            room = min(base_mu, 1 - base_mu) * 2 * 0.98
            if abs(d_mml) > room:
                d_mml = float(np.sign(d_mml) * room)

        d_icc = 0.0
        if has_icc:
            # correlation imbalance needs headroom on both sides of baseline_icc;
            # centre the pair on baseline_icc and clip the pair, not each member,
            # so the realised delta matches what is recorded as truth
            d_icc = float(rng.uniform(*icc_effect)) * (1 if rng.random() < 0.5 else -1)
            room = min(baseline_icc, 1 - baseline_icc) * 2 * 0.98
            if abs(d_icc) > room:
                d_icc = float(np.sign(d_icc) * room)

        mu1, mu2 = base_mu + d_mml / 2, base_mu - d_mml / 2
        icc1, icc2 = baseline_icc + d_icc / 2, baseline_icc - d_icc / 2

        C = int(rng.integers(cpgs_per_region[0], cpgs_per_region[1] + 1))
        offsets = np.sort(rng.integers(0, intra_region_span, size=C))
        positions = pos + offsets
        R1 = int(rng.integers(*reads_per_hap))
        R2 = int(rng.integers(*reads_per_hap))

        m1 = _draw_molecules(rng, kind, R1, positions, mu1, icc1, decay_bp)
        m2 = _draw_molecules(rng, kind, R2, positions, mu2, icc2, decay_bp)
        if no_call_rate > 0:
            m1[rng.random(m1.shape) < no_call_rate] = np.nan
            m2[rng.random(m2.shape) < no_call_rate] = np.nan

        out.append(PatternRegion(
            chrom="chr1", start=int(pos), positions=positions, m1=m1, m2=m2,
            is_true_asm=bool(has_mml or has_icc),
            true_delta_mml=d_mml, true_delta_icc=d_icc, kind=kind,
        ))
        pos += region_spacing
    return out


def per_cpg_counts(regions) -> tuple[np.ndarray, ...]:
    """Flatten to per-CpG (x1, n1, x2, n2), for dispersion estimation and for
    running the per-CpG / pooled tests on the same data."""
    x1, n1, x2, n2 = [], [], [], []
    for r in regions:
        for j in range(r.n_cpgs):
            c1, c2 = r.m1[:, j], r.m2[:, j]
            a, b = np.isfinite(c1), np.isfinite(c2)
            if a.sum() == 0 or b.sum() == 0:
                continue
            x1.append(int(np.nansum(c1))); n1.append(int(a.sum()))
            x2.append(int(np.nansum(c2))); n2.append(int(b.sum()))
    return (np.array(x1), np.array(n1), np.array(x2), np.array(n2))
