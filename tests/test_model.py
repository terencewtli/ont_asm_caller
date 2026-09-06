"""
Calibration and power checks for the beta-binomial ASM test. These are the
core rigor requirement before trusting this on real data: a test that
doesn't control its false-positive rate, or whose power collapses
inexplicably at low depth, is exactly the failure mode that made the
original Fisher-exact-test pipeline untrustworthy (see
md/20260905.progress.md) -- these checks exist to catch that class of bug
directly, not just exercise the code.

Note: `ont_asm_caller.model.test_locus` is imported as `run_asm_test` here
(NOT as `test_locus`) purely to avoid pytest auto-collecting the imported
function itself as a test item -- pytest matches any module-level callable
named `test_*`, including imported ones, and its real signature doesn't
match pytest's fixture-injection expectations.
"""
import numpy as np
import pytest

from ont_asm_caller.model import test_locus as run_asm_test
from ont_asm_caller.model import fit_mu, bb_logpmf
from ont_asm_caller.simulate import simulate_loci


def test_identical_counts_give_high_pvalue():
    r = run_asm_test("chr1", 100, x1=10, n1=20, x2=10, n2=20, rho=0.05)
    assert r.pval > 0.5
    assert abs(r.delta) < 1e-6


def test_extreme_difference_gives_low_pvalue():
    r = run_asm_test("chr1", 100, x1=0, n1=30, x2=30, n2=30, rho=0.05)
    assert r.pval < 1e-6
    assert r.delta < -0.9


def test_zero_dispersion_matches_binomial_mle():
    # At rho -> 0 the beta-binomial MLE of a shared mean should reduce to
    # the pooled proportion, same as a plain binomial model would give.
    mu_hat, _ = fit_mu(np.array([7, 13]), np.array([20, 20]), rho=1e-9)
    assert mu_hat == pytest.approx(20 / 40, abs=1e-3)


def test_bb_logpmf_reduces_to_binomial_at_low_rho():
    from scipy.stats import binom
    x, n, mu = 6, 20, 0.3
    bb = bb_logpmf(np.array([x]), np.array([n]), mu, rho=1e-9)[0]
    binom_ll = binom.logpmf(x, n, mu)
    assert bb == pytest.approx(binom_ll, abs=1e-3)


@pytest.mark.parametrize("depth", [8, 20, 60])
def test_calibration_null_false_positive_rate(depth):
    """Under the null (true_delta=0), the fraction of loci called
    significant at alpha=0.05 should be close to 0.05, AT EVERY DEPTH. This
    is precisely the property the original Fisher-exact-test pipeline
    violated -- it conflated low depth with "no difference" instead of
    properly widening its uncertainty, which is what let one low-coverage
    sample's candidate rate collapse by >50,000x relative to a high-coverage
    sample using identical thresholds. If this test passes at all three
    depths, that specific failure mode is fixed."""
    rho_true = 0.05
    loci = simulate_loci(n_loci=3000, depth=depth, true_delta=0.0,
                          rho=rho_true, seed=42)
    pvals = np.array([
        run_asm_test("chr1", i, l.x1, l.n1, l.x2, l.n2, rho=rho_true).pval
        for i, l in enumerate(loci)
    ])
    fpr = (pvals < 0.05).mean()
    # Binomial 99% CI half-width for p=0.05, n=3000 is ~0.0116
    assert 0.05 - 0.02 < fpr < 0.05 + 0.02, (
        f"depth={depth}: false-positive rate {fpr:.4f}, expected ~0.05"
    )


def test_power_increases_with_depth():
    """A real effect should be easier to detect at higher depth -- sanity
    check that the test is actually depth-aware.

    Parameter choice note: at rho=0.05, the beta-binomial's per-haplotype
    variance is p(1-p)/N * (1 + (N-1)*rho) -- as N grows, this does NOT go
    to zero, it asymptotes to p(1-p)*rho. That's a real, correct property
    of overdispersed count models (genuine site-to-site variability doesn't
    vanish just because you sequence deeper), not a bug -- but it means a
    small effect size can have a hard power CEILING that no amount of depth
    clears. true_delta=0.4 here is picked so the asymptotic ceiling is
    comfortably above 0.5 (~0.79 by the same calculation), so depth=8 vs
    depth=60 cleanly demonstrates increasing-but-still-bounded power rather
    than accidentally probing a delta too small for depth to ever rescue.
    """
    rho_true = 0.05
    true_delta = 0.4

    def detection_rate(depth):
        loci = simulate_loci(n_loci=1000, depth=depth, true_delta=true_delta,
                              rho=rho_true, seed=7)
        pvals = np.array([
            run_asm_test("chr1", i, l.x1, l.n1, l.x2, l.n2, rho=rho_true).pval
            for i, l in enumerate(loci)
        ])
        return (pvals < 0.05).mean()

    low = detection_rate(8)
    high = detection_rate(60)
    assert high > low, f"power did not increase with depth: low={low:.3f} high={high:.3f}"
    assert high > 0.5, f"power at high depth unexpectedly low: {high:.3f}"


def test_uniform_depth_ratio_does_not_produce_absurd_candidate_ratio():
    """Direct regression test for the specific failure this project hit:
    running the SAME test/threshold on a ~2x-depth-different pair of samples
    should NOT produce a >100x difference in candidate rate when the
    underlying true ASM rate is identical. (The old Fisher pipeline showed a
    >50,000x difference between an 18x and 8x median-depth sample on
    identical real data.)"""
    rho_true = 0.05
    true_delta = 0.15
    n_loci = 2000

    def candidate_rate(depth):
        loci = simulate_loci(n_loci=n_loci, depth=depth, true_delta=true_delta,
                              rho=rho_true, seed=99)
        pvals = np.array([
            run_asm_test("chr1", i, l.x1, l.n1, l.x2, l.n2, rho=rho_true).pval
            for i, l in enumerate(loci)
        ])
        return (pvals < 0.05).mean()

    rate_low = candidate_rate(8)   # matches HG00126's median tested depth
    rate_high = candidate_rate(18)  # matches NA18508's median tested depth
    ratio = rate_high / max(rate_low, 1e-9)
    assert ratio < 5, (
        f"candidate-rate ratio {ratio:.1f}x between 8x and 18x depth -- "
        "still exhibiting the depth-confound failure mode this model exists to fix"
    )
