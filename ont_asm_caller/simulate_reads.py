"""
Read-level simulator: molecules, not independent per-CpG draws.

WHY THIS EXISTS
---------------
`simulate.py` and `simulate_mixture.py` draw each CpG's counts independently,
including an independent depth per CpG. That is not what long-read data looks
like. A 35 kb ONT read covers every CpG in a 500 bp region, so within a region:

  * depth is essentially IDENTICAL across CpGs (the same molecules), and
  * methylation states are CORRELATED across CpGs (the same molecules, and
    methylation is autocorrelated along a molecule -- processivity, shared
    chromatin environment, clonal inheritance).

Pooling counts across CpGs as if they were independent inflates the effective
sample size by roughly the number of CpGs pooled, which inflates the LRT
statistic and breaks type I error control. See
docs/2026-09-05_calibration_critique.md for the measured size of that effect.

This module generates data with the read-sharing structure, parameterised by:

  co_methylation : float in [0, 1]
      Probability that a CpG's state is copied from the molecule's latent
      state rather than drawn independently. 0 reproduces the old independent
      behaviour; 1 is perfect co-methylation along the molecule. Real data
      sits in between and DECAYS WITH DISTANCE -- measure it from real per-read
      extracts (`modkit extract`) before trusting any number from here.

  baseline_mu : "bimodal" | float
      "bimodal" draws a realistic methylome (70% high / 20% low / 10%
      intermediate). The old simulators hard-coded 0.5, which is both
      unrepresentative and the single most favourable point for the
      Pearson-chi-square dispersion estimator.
"""
from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np


@dataclass
class ReadRegion:
    """One region's read-level data. `m1`/`m2` are (n_reads, n_cpgs) arrays of
    0/1 methylation calls with np.nan for no-call. Rows are MOLECULES -- this
    is the independent unit of observation, which is the whole point."""
    chrom: str
    start: int
    positions: np.ndarray
    m1: np.ndarray
    m2: np.ndarray
    is_true_asm: bool = False
    true_delta: float = 0.0

    @property
    def n_cpgs(self) -> int:
        return len(self.positions)

    def counts(self) -> tuple[int, int, int, int]:
        """Pooled read x CpG counts -- what region.cluster_cpgs produces.
        Provided so the OLD path can be benchmarked on the SAME data."""
        return (int(np.nansum(self.m1)), int(np.isfinite(self.m1).sum()),
                int(np.nansum(self.m2)), int(np.isfinite(self.m2).sum()))

    def read_counts(self, threshold: float = 0.5) -> tuple[int, int, int, int]:
        """Read-level counts: one observation per molecule, the molecule
        called methylated if its mean over the region exceeds `threshold`."""
        def f(m):
            ok = np.isfinite(m).any(axis=1)
            if ok.sum() == 0:
                return 0, 0
            frac = np.nanmean(m[ok], axis=1)
            return int((frac > threshold).sum()), int(ok.sum())
        x1, n1 = f(self.m1)
        x2, n2 = f(self.m2)
        return x1, n1, x2, n2


def draw_baseline_mu(rng, spec="bimodal") -> float:
    if spec == "bimodal":
        u = rng.random()
        if u < 0.70:
            return float(rng.beta(20, 2))
        if u < 0.90:
            return float(rng.beta(2, 20))
        return float(rng.beta(3, 3))
    return float(spec)


def _molecules(rng, n_reads, n_cpgs, mu, co_methylation):
    """(n_reads, n_cpgs) 0/1 matrix with within-molecule correlation."""
    latent = rng.random(n_reads) < mu                 # molecule-level state
    indep = rng.random((n_reads, n_cpgs)) < mu        # site-level draw
    copy = rng.random((n_reads, n_cpgs)) < co_methylation
    return np.where(copy, latent[:, None], indep).astype(float)


def simulate_read_regions(
    n_regions: int,
    pi: float = 0.005,
    cpgs_per_region=(4, 8),
    intra_region_span: int = 400,
    region_spacing: int = 5000,
    reads_per_hap=(8, 25),
    baseline_mu="bimodal",
    co_methylation: float = 0.6,
    effect_size_dist="uniform_0.2_0.6",
    no_call_rate: float = 0.02,
    seed: int | None = None,
) -> list[ReadRegion]:
    """Generate regions with realistic read-sharing and a realistic methylome.

    `reads_per_hap` is drawn ONCE PER REGION per haplotype -- the same molecules
    span every CpG in the region, which is the structural fact the old
    simulators missed.
    """
    rng = np.random.default_rng(seed)
    if not (isinstance(effect_size_dist, str)
            and effect_size_dist.startswith("uniform_")):
        raise ValueError(f"unsupported effect_size_dist: {effect_size_dist}")
    _, lo, hi = effect_size_dist.split("_")
    lo, hi = float(lo), float(hi)

    out, pos = [], 1000
    for _ in range(n_regions):
        is_asm = bool(rng.random() < pi)
        delta = float(rng.uniform(lo, hi)) if is_asm else 0.0
        if is_asm and rng.random() < 0.5:
            delta = -delta

        base = draw_baseline_mu(rng, baseline_mu)
        # keep both haplotype means in range; shrink delta rather than clip,
        # so the realised effect matches `true_delta` instead of silently
        # saturating (clipping was a quiet bias source in simulate_mixture).
        room = min(base, 1 - base) * 2 * 0.98
        if abs(delta) > room:
            delta = float(np.sign(delta) * room)
        mu1, mu2 = base + delta / 2, base - delta / 2

        C = int(rng.integers(cpgs_per_region[0], cpgs_per_region[1] + 1))
        offsets = np.sort(rng.integers(0, intra_region_span, size=C))
        R1 = int(rng.integers(*reads_per_hap))
        R2 = int(rng.integers(*reads_per_hap))

        m1 = _molecules(rng, R1, C, mu1, co_methylation)
        m2 = _molecules(rng, R2, C, mu2, co_methylation)
        if no_call_rate > 0:
            m1[rng.random(m1.shape) < no_call_rate] = np.nan
            m2[rng.random(m2.shape) < no_call_rate] = np.nan

        out.append(ReadRegion(
            chrom="chr1", start=pos, positions=pos + offsets,
            m1=m1, m2=m2, is_true_asm=is_asm, true_delta=delta,
        ))
        pos += region_spacing
    return out


def per_cpg_counts(regions: list[ReadRegion]):
    """Flatten to per-CpG (x1, n1, x2, n2) arrays, for dispersion estimation
    and for running the per-CpG test on the same data."""
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
