"""
Sanity checks on the simulator itself -- if this is wrong, every other test
that depends on it is meaningless.
"""
import numpy as np

from ont_asm_caller.simulate import simulate_loci


def test_null_simulation_has_no_average_delta():
    loci = simulate_loci(n_loci=5000, depth=20, true_delta=0.0, rho=0.05, seed=1)
    obs_delta = np.mean([l.x1 / l.n1 - l.x2 / l.n2 for l in loci])
    assert abs(obs_delta) < 0.02


def test_nonzero_delta_simulation_shows_expected_direction():
    loci = simulate_loci(n_loci=2000, depth=30, true_delta=0.3, rho=0.02, seed=2)
    obs_delta = np.mean([l.x1 / l.n1 - l.x2 / l.n2 for l in loci])
    assert obs_delta > 0.2


def test_switch_errors_attenuate_observed_delta():
    """A haplotype-mislabeling rate should shrink the OBSERVED delta toward
    zero relative to the true delta -- this is the mechanism by which
    imperfect phasing makes real ASM harder to detect, distinct from (and
    not a source of) the huge spurious deltas seen in the original pipeline,
    which is why this alone can't explain that anomaly (see
    md/20260905.progress.md) but is still worth modeling explicitly."""
    true_delta = 0.4
    clean = simulate_loci(n_loci=3000, depth=30, true_delta=true_delta,
                           rho=0.02, switch_error_rate=0.0, seed=3)
    noisy = simulate_loci(n_loci=3000, depth=30, true_delta=true_delta,
                           rho=0.02, switch_error_rate=0.2, seed=3)

    delta_clean = np.mean([abs(l.x1 / l.n1 - l.x2 / l.n2) for l in clean])
    delta_noisy = np.mean([abs(l.x1 / l.n1 - l.x2 / l.n2) for l in noisy])

    assert delta_noisy < delta_clean
    # Expected attenuation factor ~ (1 - 2*eps) = 0.6 at eps=0.2, allow slack
    # for binomial sampling noise on top of the mixture effect.
    ratio = delta_noisy / delta_clean
    assert 0.4 < ratio < 0.85, f"attenuation ratio {ratio:.2f}, expected ~0.6"


def test_depth_tuple_produces_variable_depth():
    loci = simulate_loci(n_loci=200, depth=(5, 50), true_delta=0.0,
                          rho=0.05, seed=4)
    depths = [l.n1 for l in loci]
    assert min(depths) >= 5 and max(depths) <= 50
    assert len(set(depths)) > 10  # actually varying, not silently constant
