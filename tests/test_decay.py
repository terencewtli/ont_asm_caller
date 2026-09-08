"""`decay_bp` is the parameter every pattern-method result is conditional on and
the one number that has never been measured on real data. These tests establish
that the estimator recovers it, and -- more usefully -- that it agrees with the
independently written `readlevel.design_effect`, since the two are the same
quantity approached from opposite directions."""
import numpy as np
import pytest

from ont_asm_caller.decay import comethylation_decay, implied_design_effect
from ont_asm_caller.simulate_patterns import simulate_pattern_regions
from ont_asm_caller.simulate_reads import simulate_read_regions
from ont_asm_caller.readlevel import design_effect


def _markov_regions(decay_bp, n=1500, span=1500, cpgs=(8, 14), seed=2):
    return simulate_pattern_regions(
        n, kind="markov", pi_mml=0.0, pi_icc=0.0, baseline_icc=1.0,
        decay_bp=decay_bp, baseline_mu=0.5, intra_region_span=span,
        cpgs_per_region=cpgs, no_call_rate=0.0, seed=seed)


@pytest.mark.parametrize("true_L", [100.0, 200.0, 600.0, 2000.0])
def test_recovers_known_decay_length(true_L):
    c = comethylation_decay(_markov_regions(true_L))
    assert c.r2 > 0.95
    assert c.decay_bp == pytest.approx(true_L, rel=0.10)
    assert c.amplitude == pytest.approx(1.0, abs=0.08)   # true ICC is 1.0 here


def test_exchangeable_data_has_no_decay():
    """The simulate_reads.py model is exchangeable, so the curve must come back
    FLAT at q^2 and the fitted decay length effectively infinite. If this ever
    starts returning a finite length, the simulator has silently changed."""
    regs = simulate_read_regions(1500, pi=0.0, co_methylation=0.6,
                                 baseline_mu=0.5, no_call_rate=0.0, seed=3)
    c = comethylation_decay(regs)
    assert np.allclose(c.corr, 0.36, atol=0.03)
    assert c.decay_bp > 1e5 or not np.isfinite(c.decay_bp)


def test_agrees_with_independently_implemented_design_effect():
    """Two estimators written from different definitions, on the same data. If
    they disagree on real data, the legality of region.cluster_cpgs pooling is
    unresolved and neither number should be quoted."""
    for true_L in (200.0, 2000.0):
        regs = _markov_regions(true_L, n=1200, span=400, cpgs=(8, 8), seed=4)
        c = comethylation_decay(regs)
        spacing = float(np.mean([np.diff(r.positions).mean() for r in regs]))
        implied = implied_design_effect(c, spacing, 8)
        measured = float(np.mean([design_effect(r.m1, r.m2) for r in regs]))
        assert implied == pytest.approx(measured, rel=0.25)


def test_tail_bins_would_bias_the_fit_without_the_correlation_floor():
    """Documents WHY `min_corr` exists: bins past the decay length carry true
    correlation ~0 and non-zero noise, and log(small noisy positive) flattens
    the fitted slope. Dropping the floor must visibly overstate the decay."""
    regs = _markov_regions(100.0, n=1500)
    good = comethylation_decay(regs)
    bad = comethylation_decay(regs, min_corr=1e-3)
    assert good.decay_bp == pytest.approx(100.0, rel=0.10)
    assert bad.decay_bp > good.decay_bp * 1.10
    assert bad.amplitude < good.amplitude


def test_no_call_handling_does_not_bias_the_curve():
    clean = comethylation_decay(_markov_regions(600.0, n=1200))
    holey = comethylation_decay(simulate_pattern_regions(
        1200, kind="markov", pi_mml=0.0, pi_icc=0.0, baseline_icc=1.0,
        decay_bp=600.0, baseline_mu=0.5, intra_region_span=1500,
        cpgs_per_region=(8, 14), no_call_rate=0.15, seed=2))
    assert holey.decay_bp == pytest.approx(clean.decay_bp, rel=0.25)


@pytest.mark.parametrize("eps", [0.02, 0.05, 0.10, 0.20])
def test_call_error_attenuates_by_the_known_factor_and_is_correctable(eps):
    """Binary call error attenuates correlation by exactly (1-2*eps)^2. This is
    not academic: modkit extract emits a posterior per modification, and turning
    an ambiguous posterior into a hard 0/1 call injects exactly this error. The
    design effect is a correlation measure, so an uncorrected DE from noisy
    calls is a LOWER BOUND on the real one -- and DE is what decides whether
    region pooling is legal at all."""
    import copy
    from ont_asm_caller.decay import deattenuate
    base = _markov_regions(600.0, n=1200)
    clean = comethylation_decay(base)

    rng = np.random.default_rng(1)
    noisy_regs = []
    for r in base:
        r2 = copy.copy(r)
        r2.m1 = np.where(rng.random(r.m1.shape) < eps, 1 - r.m1, r.m1)
        r2.m2 = np.where(rng.random(r.m2.shape) < eps, 1 - r.m2, r.m2)
        noisy_regs.append(r2)
    noisy = comethylation_decay(noisy_regs)

    # the decay LENGTH survives call error; only the amplitude is attenuated
    assert noisy.decay_bp == pytest.approx(clean.decay_bp, rel=0.10)
    assert noisy.amplitude == pytest.approx(clean.amplitude * (1 - 2 * eps) ** 2,
                                            rel=0.05)
    # ... and correcting for a known eps recovers it
    assert deattenuate(noisy, eps).amplitude == pytest.approx(clean.amplitude,
                                                              rel=0.05)


def test_deattenuate_rejects_impossible_error_rates():
    curve = comethylation_decay(_markov_regions(600.0, n=800))
    from ont_asm_caller.decay import deattenuate
    for bad in (0.5, 0.7, -0.1):
        with pytest.raises(ValueError):
            deattenuate(curve, bad)
