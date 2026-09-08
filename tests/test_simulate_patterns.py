"""What the pattern simulator must get right for the CPEL comparison to mean
anything -- in particular, the two things simulate_reads.py does not."""
import numpy as np
import pytest

from ont_asm_caller.simulate_patterns import (
    molecules_exchangeable, molecules_markov, molecules_epiallele,
    simulate_pattern_regions, _epiallele_pair)
from ont_asm_caller.simulate_reads import _molecules


def _corr_by_lag(m):
    return np.array([np.corrcoef(m[:, 0], m[:, j])[0, 1] for j in range(1, m.shape[1])])


def test_simulate_reads_correlation_is_exchangeable_not_decaying():
    """The finding that motivates this module. simulate_reads._molecules copies
    from a single molecule-level latent, so corr(x_i, x_j) = q^2 at EVERY lag.
    A simulator with no distance decay cannot be used to benchmark a method
    whose entire premise is distance-decaying correlation."""
    rng = np.random.default_rng(0)
    q = 0.6
    m = _molecules(rng, 80000, 12, 0.5, q)
    cors = _corr_by_lag(m)
    assert np.allclose(cors, q ** 2, atol=0.02)
    assert cors[-1] / cors[0] > 0.9          # flat, no decay whatsoever


def test_markov_generator_decays_with_distance():
    rng = np.random.default_rng(0)
    positions = np.arange(12) * 100.0
    m = molecules_markov(rng, 80000, positions, mu=0.5, icc=1.0, decay_bp=300.0)
    cors = _corr_by_lag(m)
    expected = np.exp(-np.arange(1, 12) * 100.0 / 300.0)
    assert np.allclose(cors, expected, atol=0.03)
    assert cors[-1] < 0.1 * cors[0]          # decays, unlike the above


def test_marginal_mean_is_preserved_by_every_generator():
    """Correlation must be varied WITHOUT moving the mean, or the two ASM axes
    are not orthogonal and the benchmark cannot attribute power to either."""
    rng = np.random.default_rng(1)
    for mu in (0.2, 0.5, 0.85):
        for icc in (0.1, 0.5, 0.95):
            a = molecules_exchangeable(rng, 40000, 8, mu, icc)
            b = molecules_markov(rng, 40000, np.arange(8) * 50.0, mu, icc, 500.0)
            assert a.mean() == pytest.approx(mu, abs=0.02)
            assert b.mean() == pytest.approx(mu, abs=0.02)


def test_epiallele_icc_moves_read_fraction_spread_at_fixed_mean():
    rng = np.random.default_rng(2)
    n = 10
    spreads, means = [], []
    for icc in (0.05, 0.5, 0.95):
        pat, w = _epiallele_pair(rng, n, 0.5, icc)
        m = molecules_epiallele(rng, 20000, n, pat, w, flip_rate=0.02)
        means.append(m.mean())
        spreads.append(m.mean(axis=1).std())
    assert max(means) - min(means) < 0.05       # mean held
    assert spreads[0] < spreads[1] < spreads[2]  # disorder varies


def test_correlation_axis_regions_carry_zero_mean_signal():
    """A correlation-axis true positive must have mu1 == mu2 exactly, otherwise
    a mean test could pass it and the 'blind by construction' claim is false."""
    regs = simulate_pattern_regions(4000, pi_mml=0.0, pi_icc=0.5, seed=3,
                                    no_call_rate=0.0)
    icc_regs = [r for r in regs if r.true_delta_icc != 0]
    assert len(icc_regs) > 100
    d = np.array([np.nanmean(r.m1) - np.nanmean(r.m2) for r in icc_regs])
    assert abs(d.mean()) < 0.02                 # no systematic mean difference
    assert all(r.true_delta_mml == 0.0 for r in icc_regs)


def test_reads_are_shared_across_cpgs_of_a_region():
    """Inherited requirement from simulate_reads.py: one molecule spans every
    CpG in the region. Losing this reintroduces the v1 design-effect defect."""
    regs = simulate_pattern_regions(50, pi_mml=0, pi_icc=0, seed=4,
                                    no_call_rate=0.0)
    for r in regs:
        assert r.m1.shape[1] == r.n_cpgs
        assert np.isfinite(r.m1).all()          # every read observes every CpG
