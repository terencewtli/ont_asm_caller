"""Correctness of the CPEL-style chain machinery.

Everything here is checked against brute-force enumeration of the 2^N state
space, which is the only way to know the transfer-matrix recursions are right.
"""
import itertools

import numpy as np
import pytest

from ont_asm_caller import pattern as P


def _brute(h, J):
    N = len(h)
    states = np.array(list(itertools.product([0, 1], repeat=N)))
    s = 2.0 * states - 1.0
    lp = s @ h + ((s[:, :-1] * s[:, 1:]) @ J if N > 1 else 0.0)
    Z = np.exp(lp).sum()
    return states, np.exp(lp) / Z, np.log(Z)


@pytest.mark.parametrize("seed", range(5))
def test_transfer_matrix_matches_enumeration(seed):
    rng = np.random.default_rng(seed)
    N = int(rng.integers(2, 9))
    h, J = rng.normal(0, 1.2, N), rng.normal(0, 0.8, N - 1)
    states, p, logZ = _brute(h, J)

    assert P.chain_logZ(h, J) == pytest.approx(logZ, abs=1e-9)

    p1, _ = P.chain_marginals(h, J)
    brute_marg = np.stack([[(p * (states[:, n] == 0)).sum(),
                            (p * (states[:, n] == 1)).sum()] for n in range(N)])
    assert np.abs(p1 - brute_marg).max() < 1e-9

    # entropy via H(s_1) + sum H(s_{n+1}|s_n), not via enumeration
    assert P.chain_entropy_bits(h, J) == pytest.approx(
        -(p * np.log2(p)).sum(), abs=1e-9)
    assert P.chain_mml(h, J) == pytest.approx(
        float((p[:, None] * states).sum(axis=0).mean()), abs=1e-9)
    assert np.abs(P._chain_logpmf(states, h, J) - np.log(p)).max() < 1e-9


def test_clamped_forward_gives_configuration_probability():
    """The masked forward pass is how partially-observed reads enter the
    likelihood. Fully clamping every site must reproduce log p(x) + log Z."""
    rng = np.random.default_rng(0)
    N = 6
    h, J = rng.normal(0, 1.0, N), rng.normal(0, 0.7, N - 1)
    states, p, logZ = _brute(h, J)
    NEG = -1e9
    for idx in (0, 7, 31, 63):
        x = states[idx]
        mask = np.zeros((N, 2))
        mask[np.arange(N), 1 - x] = NEG
        _, lz = P._forward(h, J, mask=mask)
        assert lz == pytest.approx(np.log(p[idx]) + logZ, abs=1e-6)


def test_batched_and_single_forward_agree():
    rng = np.random.default_rng(1)
    N = 5
    h, J = rng.normal(0, 1.0, N), rng.normal(0, 0.5, N - 1)
    masks = np.where(rng.random((7, N, 2)) < 0.3, -1e9, 0.0)
    _, batched = P._forward(h, J, mask=masks)
    for r in range(7):
        _, single = P._forward(h, J, mask=masks[r])
        assert batched[r] == pytest.approx(single, abs=1e-9)


def test_analytic_gradient_matches_finite_differences():
    """The fit uses an analytic gradient of a partially-observed likelihood.
    A wrong gradient still converges to something plausible-looking, so this
    has to be checked directly rather than inferred from the fitted values."""
    import scipy.optimize as so
    rng = np.random.default_rng(11)
    for _ in range(3):
        N, R = int(rng.integers(3, 8)), int(rng.integers(6, 22))
        x = (rng.random((R, N)) < rng.uniform(0.2, 0.8)).astype(float)
        x[rng.random(x.shape) < 0.08] = np.nan
        pos = np.sort(rng.integers(0, 900, N)).astype(float)

        captured = {}
        orig = so.minimize

        def spy(fun, x0, **kw):
            captured["fun"], captured["x0"] = fun, x0
            return orig(fun, x0, **kw)

        P.minimize = spy
        try:
            P.fit_chain(x, pos)
        finally:
            P.minimize = orig

        fun, x0 = captured["fun"], captured["x0"]
        theta = x0 + rng.normal(0, 0.4, len(x0))
        _, g = fun(theta)
        gn = np.empty_like(g)
        for i in range(len(theta)):
            e = np.zeros_like(theta); e[i] = 1e-6
            gn[i] = (fun(theta + e)[0] - fun(theta - e)[0]) / 2e-6
        assert np.abs(g - gn).max() < 1e-5 * max(1.0, np.abs(g).max())


def test_mle_recovers_known_parameters():
    rng = np.random.default_rng(3)
    for alpha, beta in [(0.8, 0.0), (0.0, 0.9), (-1.0, 0.4)]:
        h, J = np.full(8, alpha), np.full(7, beta)
        x = P._sample_chain(rng, h, J, 3000).astype(float)
        _, _, _, theta = P.fit_chain(x, np.arange(8) * 40.0)
        assert theta[0] == pytest.approx(alpha, abs=0.12)
        assert theta[-1] == pytest.approx(beta, abs=0.12)


def test_fit_is_finite_on_fully_methylated_region():
    """The single most common configuration in a real bimodal methylome. An
    unpenalised MLE diverges here; the ridge keeps it finite and the entropy
    at its correct value of zero."""
    for x in (np.ones((20, 6)), np.zeros((20, 6))):
        h, J, ll, theta = P.fit_chain(x)
        assert np.all(np.isfinite(theta))
        assert P.chain_entropy_bits(h, J) / 6 < 0.05


def test_nme_separates_disorder_at_equal_mean():
    """The claim the whole CPEL comparison rests on: NME distinguishes an
    allele whose molecules are individually noisy from one that is bistable
    across molecules, when both have mean methylation 0.5."""
    rng = np.random.default_rng(7)
    n, R = 8, 60
    iid = (rng.random((R, n)) < 0.5).astype(float)
    # exactly half the molecules fully methylated, so the two arms are matched
    # on mean BY CONSTRUCTION -- drawing the split at random leaves a ~0.065 SD
    # wobble in the bistable arm's mean that has nothing to do with entropy
    bistable = np.repeat(np.r_[np.ones(R // 2), np.zeros(R - R // 2)][:, None],
                         n, axis=1)

    def nme(x):
        h, J, _, _ = P.fit_chain(x)
        return P.chain_entropy_bits(h, J) / n

    assert abs(np.nanmean(iid) - np.nanmean(bistable)) < 0.08    # equal means
    assert nme(iid) > 0.9
    assert nme(bistable) < 0.3


def test_jsd_matches_enumeration():
    rng = np.random.default_rng(5)
    N = 5
    h1, J1 = rng.normal(0, 1.0, N), rng.normal(0, 0.6, N - 1)
    h2, J2 = rng.normal(0, 1.0, N), rng.normal(0, 0.6, N - 1)
    _, p1, _ = _brute(h1, J1)
    _, p2, _ = _brute(h2, J2)
    mix = 0.5 * (p1 + p2)
    exact = 0.5 * (p1 * np.log2(p1 / mix)).sum() + 0.5 * (p2 * np.log2(p2 / mix)).sum()
    mc = P.jsd_chains(h1, J1, h2, J2, n_mc=20000, seed=0)
    assert mc == pytest.approx(exact, abs=0.02)
    assert P.jsd_chains(h1, J1, h1, J1, n_mc=4000) == pytest.approx(0.0, abs=0.01)


def test_ks_sees_what_the_location_tests_cannot():
    """A bistable allele against a fixed-fraction allele, matched on mean. The
    existing read-level test compares locations and should find nothing; KS
    compares distributions and should."""
    from ont_asm_caller.readlevel import test_region_reads
    rng = np.random.default_rng(2)
    n, R = 8, 40
    fixed = (rng.random((R, n)) < 0.5).astype(float)                    # frac ~ 0.5
    bistable = np.repeat(np.r_[np.ones(R // 2), np.zeros(R - R // 2)][:, None],
                         n, axis=1)

    welch = test_region_reads("chr1", 0, 0, fixed, bistable)
    ks = P.read_fraction_ks(fixed, bistable)
    assert abs(welch.delta) < 0.15          # no mean signal to find
    assert welch.pval > 0.05                # and the location test finds none
    assert ks[1] < 1e-4                     # the distributional one does
