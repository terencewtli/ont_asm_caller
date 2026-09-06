"""
Checks that estimate_dispersion actually recovers the true dispersion used
to generate synthetic data -- and that pooling in the (rare) true-ASM loci
doesn't meaningfully bias it, which is the assumption the estimator's
docstring relies on.
"""
import numpy as np
import pytest

from ont_asm_caller.dispersion import estimate_dispersion
from ont_asm_caller.simulate import simulate_loci


@pytest.mark.parametrize("rho_true", [0.0, 0.02, 0.1, 0.2])
def test_recovers_true_dispersion_under_null(rho_true):
    loci = simulate_loci(n_loci=20000, depth=(5, 40), true_delta=0.0,
                          rho=rho_true, seed=123)
    x1 = np.array([l.x1 for l in loci])
    n1 = np.array([l.n1 for l in loci])
    x2 = np.array([l.x2 for l in loci])
    n2 = np.array([l.n2 for l in loci])

    rho_hat = estimate_dispersion(x1, n1, x2, n2)
    assert rho_hat == pytest.approx(rho_true, abs=0.03), (
        f"rho_true={rho_true}, rho_hat={rho_hat:.4f}"
    )


def test_robust_to_a_small_fraction_of_true_asm_loci():
    """~0.5% true ASM loci (matching deCODE's validated rate order of
    magnitude) mixed into an otherwise-null pool should not meaningfully
    bias the dispersion estimate -- this is the exact assumption the global
    estimator relies on (see dispersion.py docstring)."""
    rho_true = 0.05
    null_loci = simulate_loci(n_loci=19900, depth=(5, 40), true_delta=0.0,
                               rho=rho_true, seed=1)
    asm_loci = simulate_loci(n_loci=100, depth=(5, 40), true_delta=0.4,
                              rho=rho_true, seed=2)
    loci = null_loci + asm_loci

    x1 = np.array([l.x1 for l in loci])
    n1 = np.array([l.n1 for l in loci])
    x2 = np.array([l.x2 for l in loci])
    n2 = np.array([l.n2 for l in loci])

    rho_hat = estimate_dispersion(x1, n1, x2, n2)
    assert rho_hat == pytest.approx(rho_true, abs=0.03)


def test_raises_with_too_few_loci():
    with pytest.raises(ValueError):
        estimate_dispersion([5], [10], [6], [10])
