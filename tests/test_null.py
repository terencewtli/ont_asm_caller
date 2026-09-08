"""The matched null's two claims: it is calibrated, and it reaches p-values
below the permutation floor that makes the existing exact tests unusable
genome-wide."""
import numpy as np
import pytest

from ont_asm_caller.null import MatchedNull, calibration_report, split_reads
from ont_asm_caller.simulate_patterns import simulate_pattern_regions
from ont_asm_caller import pattern as P


def _cheap_stat(m1, m2):
    """KS on read fractions -- fast enough to test with, and it exercises the
    same stratification and tail-fitting code the expensive statistics use."""
    ks = P.read_fraction_ks(m1, m2)
    return None if ks is None else {"ks": ks[0]}


@pytest.fixture(scope="module")
def null_and_regions():
    regs = simulate_pattern_regions(1200, pi_mml=0.0, pi_icc=0.0, seed=21)
    mn = MatchedNull(_cheap_stat, n_perm_per_region=25).build(regs, seed=1)
    return mn, regs


def test_split_reads_moves_whole_molecules_and_preserves_sizes():
    """Permuting whole reads is what keeps within-molecule co-methylation intact
    under H0. A count-level bootstrap (scripts/phasing_qc/P04) cannot do this."""
    rng = np.random.default_rng(0)
    m1 = np.arange(20).repeat(6).reshape(20, 6).astype(float)
    m2 = (np.arange(30) + 100).repeat(6).reshape(30, 6).astype(float)
    a, b = split_reads(rng, m1, m2)
    assert a.shape == m1.shape and b.shape == m2.shape
    for row in np.vstack([a, b]):
        assert len(set(row)) == 1               # rows never got mixed
    assert sorted(np.vstack([a, b])[:, 0]) == sorted(
        np.vstack([m1, m2])[:, 0])              # nothing lost or duplicated


def test_null_is_calibrated_on_h0_data(null_and_regions):
    mn, regs = null_and_regions
    rep = calibration_report(mn, regs, _cheap_stat, ["ks"], seed=99)["ks"]
    assert rep["n"] > 800
    assert rep["type1_0.05"] == pytest.approx(0.05, abs=0.025)
    assert rep["type1_0.01"] == pytest.approx(0.01, abs=0.012)


def test_tail_extrapolation_goes_below_the_permutation_floor(null_and_regions):
    """The whole point. A within-region permutation test floors p at
    1/(n_perm+1); CPEL's own bootstrap floors it at 1/(L+1) with L=1000, which
    under genome-wide BH means no region can be called unless tens of thousands
    tie at the floor. Pooling within strata and fitting the tail must reach
    well past that."""
    mn, _ = null_and_regions
    entry = mn._fits[(("*",), "ks")]
    v, u, _, fit = entry
    assert fit is not None, "GPD tail fit failed; extrapolation is inactive"
    empirical_floor = 1.0 / (len(v) + 1)
    extreme = float(v.max() * 1.6)
    p = mn.pvalue("ks", extreme, key=("*",))
    assert p < empirical_floor / 10
    assert p > 0


def test_extrapolation_is_capped(null_and_regions):
    """Regression test for a defect found by the first benchmark run: the GPD
    tail was trusted without limit and produced p-values of 1e-300 from a
    60,000-draw null. Type I error at alpha=0.05 and 0.01 still looked correct,
    so the failure was invisible in the usual calibration summary -- it showed
    up only in the extreme tail, which is exactly where genome-wide BH operates.

    A tail fitted to the top 10% of L draws supports extending the empirical
    floor by a couple of orders of magnitude, not three hundred."""
    mn, _ = null_and_regions
    v, u, _, fit = mn._fits[(("*",), "ks")]
    L = len(v)
    cap = 1.0 / (mn.max_extrapolation * (L + 1))
    # a statistic far beyond anything observable must bottom out at the cap
    assert mn.pvalue("ks", 10.0, key=("*",)) == pytest.approx(cap, rel=1e-9)
    # ... and the cap must still be small enough to survive genome-wide BH
    # across 10^5 regions, which needs roughly alpha / n = 5e-7
    assert cap < 5e-7


def test_pathological_gpd_shapes_are_rejected(null_and_regions):
    """Every statistic here is bounded in [0, 1], so a large positive GPD shape
    -- a polynomial tail with no upper endpoint -- is structurally misspecified
    and was the source of the runaway extrapolation. Fits outside a sane shape
    range must be dropped rather than used."""
    mn, _ = null_and_regions
    for (key, name), (v, u, pe, fit) in mn._fits.items():
        if fit is not None:
            c, scale = fit
            assert -0.5 <= c <= 0.25 and scale > 0


def test_extreme_tail_is_calibrated_not_just_the_usual_alphas(null_and_regions):
    """The alpha=0.05 row is not evidence about the extreme tail; check it."""
    mn, regs = null_and_regions
    rep = calibration_report(mn, regs, _cheap_stat, ["ks"], seed=123)["ks"]
    assert rep["type1_0.001"] < 0.006
    assert rep["type1_1e-4"] < 0.002


def test_pvalues_are_monotone_in_the_statistic(null_and_regions):
    mn, _ = null_and_regions
    v = mn._fits[(("*",), "ks")][0]
    grid = np.linspace(float(np.quantile(v, 0.5)), float(v.max() * 1.5), 25)
    ps = [mn.pvalue("ks", g, key=("*",)) for g in grid]
    assert all(ps[i] >= ps[i + 1] - 1e-12 for i in range(len(ps) - 1))


def test_cpel_style_null_is_conservative_at_long_read_depth():
    """CPEL generates every null draw at Cmin = 5 reads regardless of the
    region's real coverage (Supplementary s10, step 5). At ONT depth that null
    is far wider than the truth, so the same data yields systematically larger
    p-values. This is the cost the depth-matched null removes."""
    regs = simulate_pattern_regions(700, pi_mml=0.0, pi_icc=0.0, seed=31,
                                    reads_per_hap=(20, 40))
    matched = MatchedNull(_cheap_stat, n_perm_per_region=25).build(regs, seed=2)
    cpel = MatchedNull(_cheap_stat, n_perm_per_region=25, stratify=("cpgs",),
                       cpel_cmin=5).build(regs, seed=3)
    # the Cmin=5 null must sit to the RIGHT of the depth-matched one
    v_m = matched._fits[(("*",), "ks")][0]
    v_c = cpel._fits[(("*",), "ks")][0]
    assert np.quantile(v_c, 0.95) > np.quantile(v_m, 0.95) * 1.2
