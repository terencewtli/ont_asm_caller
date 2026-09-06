"""
Read-level ASM tests: the MOLECULE is the unit of observation.

THE PROBLEM THIS FIXES
----------------------
`region.cluster_cpgs` pools raw counts across CpGs in a window:

    x1 = sum(c.x1 for c in group);  n1 = sum(c.n1 for c in group)

With 8 CpGs and 15 reads that is N = 120, but there are only 15 independent
molecules. Long reads cover every CpG in a 500 bp window, so the pooled count
is 15 observations counted 8 times. The likelihood then behaves as if it had
8x more information than it does, the LRT statistic inflates, and type I error
breaks -- measured at 2.9x inflation at alpha=0.05 and 11x at alpha=0.001 under
strong co-methylation. See docs/2026-09-05_calibration_critique.md.

The inflation depends on how co-methylated real molecules are, which is
exactly the quantity the parent project's co-methylation-decay analysis
measures (md/20260903_qc_review.md section 7.3). Those two analyses are the
same parameter viewed from two directions; they must not disagree.

THE FIX
-------
Summarise each molecule to ONE number -- its methylation fraction over the
region -- and compare the two haplotypes' read-level distributions. Nothing
needs to model co-methylation explicitly: aggregating within a read absorbs
it, and the between-read variance that drives the test is estimated from the
reads themselves.

`test_region_reads` (Welch t on per-read fractions) is the default: analytic,
gives arbitrarily small p-values (needed for genome-wide BH), and requires no
dispersion parameter at all. `test_region_perm` is an exact randomisation test
for validation and for small read counts, permuting WHOLE READS between
haplotypes so within-molecule structure is preserved exactly under the null.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy.stats import t as student_t


@dataclass
class ReadLevelResult:
    chrom: str
    start: int
    end: int
    n_cpgs: int
    n_reads1: int
    n_reads2: int
    mu1_hat: float
    mu2_hat: float
    delta: float
    stat: float
    pval: float
    method: str


def read_fractions(m: np.ndarray, min_calls: int = 1) -> np.ndarray:
    """Per-read methylation fraction over a region. `m` is (n_reads, n_cpgs)
    of 0/1/nan. Reads with fewer than `min_calls` non-nan calls are dropped."""
    ok = np.isfinite(m).sum(axis=1) >= min_calls
    if ok.sum() == 0:
        return np.array([])
    return np.nanmean(m[ok], axis=1)


def test_region_reads(chrom, start, end, m1, m2, min_calls: int = 1,
                       min_reads: int = 3) -> ReadLevelResult | None:
    """Welch two-sample t-test on per-read methylation fractions.

    Reads are the independent units, so no dispersion parameter is needed --
    between-read variance is estimated directly. Welch (unequal variance) is
    used deliberately: the two haplotypes routinely differ in both depth and
    in spread, and the equal-variance pooled test is anti-conservative when
    the deeper haplotype is also the less variable one.

    Returns None if either haplotype has < `min_reads` usable reads.
    """
    f1, f2 = read_fractions(m1, min_calls), read_fractions(m2, min_calls)
    if len(f1) < min_reads or len(f2) < min_reads:
        return None

    m_1, m_2 = float(f1.mean()), float(f2.mean())
    n_1, n_2 = len(f1), len(f2)

    # COMPLETE SEPARATION. Every read on HP1 is on one side of every read on
    # HP2. This is the STRONGEST possible ASM signal (e.g. an imprinted DMR:
    # one haplotype fully methylated, the other fully unmethylated), and a
    # t-test cannot express it -- both within-group variances are 0, so the
    # statistic is inf and any naive p-value floor (1/(n1+n2) in an earlier
    # version of this function) is far too large to survive genome-wide BH.
    # The exact randomisation p-value is the right answer: of the
    # C(n1+n2, n1) equally likely relabellings, only the observed split and
    # its mirror are this extreme.
    if f1.min() > f2.max() or f2.min() > f1.max():
        from math import comb
        pval = min(1.0, 2.0 / comb(n_1 + n_2, n_1))
        return ReadLevelResult(chrom, start, end, m1.shape[1], n_1, n_2,
                               m_1, m_2, m_1 - m_2, np.inf, float(pval),
                               "exact_separation")

    v1, v2 = float(f1.var(ddof=1)), float(f2.var(ddof=1))
    se2 = v1 / n_1 + v2 / n_2
    if se2 <= 0:                       # identical constants; no evidence
        return ReadLevelResult(chrom, start, end, m1.shape[1], n_1, n_2,
                               m_1, m_2, m_1 - m_2, 0.0, 1.0, "degenerate")

    se = np.sqrt(se2)
    stat = (m_1 - m_2) / se
    df = se2 ** 2 / ((v1 / n_1) ** 2 / max(n_1 - 1, 1)
                     + (v2 / n_2) ** 2 / max(n_2 - 1, 1))
    p_t = float(2 * student_t.sf(abs(stat), df))

    # Read fractions are bounded in [0,1] and pile up at the boundaries under
    # strong co-methylation, where t asymptotics are poor. Mann-Whitney is
    # robust to that and to ties; take the more conservative of the two rather
    # than trusting either alone.
    try:
        from scipy.stats import mannwhitneyu
        p_u = float(mannwhitneyu(f1, f2, alternative="two-sided",
                                 method="asymptotic").pvalue)
        if not np.isfinite(p_u):
            p_u = 1.0
    except Exception:
        p_u = 1.0
    pval, method = (p_t, "welch") if p_t >= p_u else (p_u, "mannwhitney")
    return ReadLevelResult(chrom, start, end, m1.shape[1], n_1, n_2,
                           m_1, m_2, m_1 - m_2, float(stat), pval, method)


def test_region_perm(chrom, start, end, m1, m2, n_perm: int = 2000,
                      min_calls: int = 1, min_reads: int = 3,
                      seed: int | None = None) -> ReadLevelResult | None:
    """Exact randomisation test: permute whole reads between haplotypes.

    This is the reference implementation -- it makes no distributional
    assumption and preserves within-molecule co-methylation exactly, because
    a read moves between haplotypes intact. Use it to validate the analytic
    test, and for regions with few reads where t asymptotics are doubtful.

    p is bounded below by 1/(n_perm+1), so this is NOT suitable on its own for
    genome-wide BH across millions of regions. Screen with `test_region_reads`,
    then confirm survivors here.
    """
    f1, f2 = read_fractions(m1, min_calls), read_fractions(m2, min_calls)
    if len(f1) < min_reads or len(f2) < min_reads:
        return None
    rng = np.random.default_rng(seed)
    obs = abs(f1.mean() - f2.mean())
    pool = np.concatenate([f1, f2])
    n_1 = len(f1)
    ge = 0
    for _ in range(n_perm):
        rng.shuffle(pool)
        if abs(pool[:n_1].mean() - pool[n_1:].mean()) >= obs - 1e-12:
            ge += 1
    pval = (ge + 1) / (n_perm + 1)
    return ReadLevelResult(chrom, start, end, m1.shape[1], len(f1), len(f2),
                           float(f1.mean()), float(f2.mean()),
                           float(f1.mean() - f2.mean()), float(obs), pval,
                           f"perm{n_perm}")


def design_effect(m1: np.ndarray, m2: np.ndarray) -> float:
    """Empirical design effect: how much the naive read x CpG pooling inflates
    the effective sample size, for one region.

        DE = Var(pooled count) / Var(count under independence)
           ~= 1 + (C_bar - 1) * ICC

    where ICC is the intra-read correlation of methylation calls. DE = 1 means
    pooling is safe; DE = C means every CpG on a read is redundant.

    Compute this on REAL `modkit extract` output to find out where the data
    actually sits. It is the single number that decides whether the pooled
    region test in region.py is usable at all.
    """
    des = []
    for m in (m1, m2):
        ok = np.isfinite(m).sum(axis=1) >= 2
        if ok.sum() < 3:
            continue
        sub = m[ok]
        C = np.isfinite(sub).sum(axis=1).mean()
        frac = np.nanmean(sub, axis=1)
        p = float(np.nanmean(sub))
        if p <= 0 or p >= 1:
            continue
        # variance of the per-read fraction vs binomial expectation
        obs_v = float(np.var(frac, ddof=1))
        exp_v = p * (1 - p) / C
        if exp_v > 0:
            des.append(max(obs_v / exp_v, 0.0))
    return float(np.mean(des)) if des else 1.0
