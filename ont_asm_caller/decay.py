"""
Co-methylation decay along a molecule: measuring `decay_bp`.

WHY THIS IS THE FIRST THING TO RUN
----------------------------------
Every result in `benchmarks/compare_pattern_methods.py` is conditional on how
fast methylation correlation decays with distance along a single molecule, and
that number has never been measured on this project's data. It is also the same
quantity as, viewed from two other directions:

  * the design effect (`readlevel.design_effect`) -- DE ~= 1 + (C_bar - 1) * ICC,
    and ICC is the area under this curve over the region's CpG spacings. The
    single number that decides whether `region.cluster_cpgs` pooling is legal at
    all (docs/2026-09-05_calibration_critique.md, "the number that decides
    everything").
  * the parent project's co-methylation-decay analysis
    (md/20260903_qc_review.md s7.3).

Three analyses, one measurement. They must agree, and if they disagree that is
itself the finding.

WHAT IT IS NOT
--------------
RESOLUTION IS CAPPED BY REGION SPAN. The curve can only see CpG pairs that
co-occur inside a supplied region, so a decay length longer than the typical
region span is not measurable -- it comes back short and confidently wrong.
Measured: feeding 400 bp regions data generated with decay_bp = 500 returns
88 bp. The production caller caps regions at `max_span = 1000`
(region.cluster_cpgs), so measuring decay beyond ~1 kb requires WIDER windows
than the caller's own regions, not the same ones. Check that the fitted
decay_bp is well inside the region span before believing it.

This is a DESCRIPTIVE pooled correlation curve, not a model fit to any one
region. Decay will differ between CpG islands, shores and open sea, and between
high- and low-methylation regions -- which is the actual reason a whole
chromosome is worth extracting rather than a few Mb. Always run it stratified
before quoting a single `decay_bp`.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np


@dataclass
class DecayCurve:
    centres: np.ndarray        # distance-bin centres, bp
    corr: np.ndarray           # pooled Pearson correlation in each bin
    n_pairs: np.ndarray        # CpG pairs contributing to each bin
    decay_bp: float            # fitted half-life of exp(-d / decay_bp)
    amplitude: float           # fitted correlation at d -> 0
    r2: float                  # fit quality of the log-linear regression

    def __repr__(self):
        return (f"DecayCurve(decay_bp={self.decay_bp:.0f}, "
                f"amplitude={self.amplitude:.3f}, r2={self.r2:.3f}, "
                f"bins={len(self.centres)})")


def _accumulate(m: np.ndarray, positions: np.ndarray, edges, sxy, sxx, syy, npair):
    """Add one haplotype's CpG-pair covariances into the distance bins.

    Each column is centred by ITS OWN mean within the region before pooling.
    Pooling raw products across regions with different methylation levels would
    manufacture correlation out of between-region mean variation, which is not
    what this is measuring."""
    obs = np.isfinite(m)
    if m.shape[0] < 4:
        return
    mu = np.nanmean(np.where(obs, m, np.nan), axis=0)
    if not np.isfinite(mu).all():
        return
    c = np.where(obs, m - mu, 0.0)
    n = m.shape[1]
    for i in range(n - 1):
        d = positions[i + 1:] - positions[i]
        b = np.searchsorted(edges, d, side="right") - 1
        ok = (b >= 0) & (b < len(edges) - 1)
        if not ok.any():
            continue
        both = obs[:, i:i + 1] & obs[:, i + 1:]
        cov = (c[:, i:i + 1] * c[:, i + 1:] * both).sum(axis=0)
        vi = ((c[:, i:i + 1] ** 2) * both).sum(axis=0)
        vj = ((c[:, i + 1:] ** 2) * both).sum(axis=0)
        np.add.at(sxy, b[ok], cov[ok])
        np.add.at(sxx, b[ok], vi[ok])
        np.add.at(syy, b[ok], vj[ok])
        np.add.at(npair, b[ok], 1)


def comethylation_decay(regions, edges=None, min_pairs: int = 50,
                        fit_max_bp: float | None = None,
                        min_corr: float = 0.03) -> DecayCurve:
    """Pooled correlation vs. CpG-CpG distance, and an exponential decay fit.

    `regions` is any iterable of objects with `.positions`, `.m1` and `.m2`
    (both `simulate_patterns.PatternRegion` and `simulate_reads.ReadRegion`
    qualify), or of plain `(positions, m)` tuples for a single haplotype.

    The two haplotypes are pooled. Under the null that is correct and under real
    ASM it is mildly conservative for the decay length (an ASM region's two
    haplotypes have different means, and each is centred separately, so no
    spurious correlation is introduced).

    `min_corr` is not cosmetic. Bins far past the decay length carry a true
    correlation of ~0 and a sampling error that is not, and log(small noisy
    positive) is a large negative outlier that flattens the fitted slope and
    biases `decay_bp` UPWARD. Measured, without a floor: a true 100 bp decay was
    recovered as 117-153 bp across replicates, and the fitted amplitude fell to
    0.46-0.79 against a true 1.0. With the floor at 0.03 the same fits land
    within 5% (99-105% of truth over 100-2000 bp) and the amplitude within 4%.

    `amplitude` is the diagnostic to read alongside `decay_bp`: it should come
    back near the region's measured ICC. An amplitude well below it means the
    fit is being dragged by near-zero tail bins and the decay length is
    overstated -- raise `min_corr` or set `fit_max_bp`.
    """
    if edges is None:
        edges = np.concatenate([[0], np.geomspace(20, 5000, 24)])
    edges = np.asarray(edges, float)
    nb = len(edges) - 1
    sxy, sxx, syy = np.zeros(nb), np.zeros(nb), np.zeros(nb)
    npair = np.zeros(nb, int)

    for r in regions:
        if isinstance(r, tuple):
            pos, mats = np.asarray(r[0], float), [np.asarray(r[1], float)]
        else:
            pos, mats = np.asarray(r.positions, float), [r.m1, r.m2]
        if len(pos) < 2:
            continue
        for m in mats:
            _accumulate(np.asarray(m, float), pos, edges, sxy, sxx, syy, npair)

    keep = (npair >= min_pairs) & (sxx > 0) & (syy > 0)
    centres = np.sqrt(edges[:-1] * np.maximum(edges[1:], 1e-9))
    centres[0] = 0.5 * (edges[0] + edges[1])
    corr = np.full(nb, np.nan)
    corr[keep] = sxy[keep] / np.sqrt(sxx[keep] * syy[keep])

    # log-linear fit of corr ~ A exp(-d / L) over bins with positive correlation
    fit_mask = keep & np.isfinite(corr) & (corr > min_corr)
    if fit_max_bp is not None:
        fit_mask &= centres <= fit_max_bp
    if fit_mask.sum() >= 3:
        d = centres[fit_mask]
        y = np.log(corr[fit_mask])
        w = np.sqrt(npair[fit_mask].astype(float))
        A = np.vstack([np.ones_like(d), -d]).T
        coef, *_ = np.linalg.lstsq(A * w[:, None], y * w, rcond=None)
        amp, inv_L = float(np.exp(coef[0])), float(coef[1])
        decay = float(1.0 / inv_L) if inv_L > 1e-9 else float("inf")
        pred = coef[0] - coef[1] * d
        ss_res = float((w * (y - pred) ** 2).sum())
        ss_tot = float((w * (y - np.average(y, weights=w)) ** 2).sum())
        r2 = 1.0 - ss_res / ss_tot if ss_tot > 0 else float("nan")
    else:
        amp, decay, r2 = float("nan"), float("nan"), float("nan")

    return DecayCurve(centres[keep], corr[keep], npair[keep], decay, amp, r2)


def implied_design_effect(curve: DecayCurve, cpg_spacing_bp: float,
                          n_cpgs: int) -> float:
    """DE = 1 + (C - 1) * ICC_bar for a region of `n_cpgs` at the given mean
    spacing, read off the fitted curve. This is the bridge between this
    measurement and `readlevel.design_effect` -- if the two disagree on real
    data, one of them is wrong and the pooled region test's legality is
    unresolved."""
    if not np.isfinite(curve.decay_bp) or n_cpgs < 2:
        return 1.0
    lags = np.arange(1, n_cpgs) * cpg_spacing_bp
    weights = (n_cpgs - np.arange(1, n_cpgs))          # pairs at each lag
    icc = float((weights * curve.amplitude * np.exp(-lags / curve.decay_bp)).sum()
                / weights.sum())
    return float(1.0 + (n_cpgs - 1) * icc)
