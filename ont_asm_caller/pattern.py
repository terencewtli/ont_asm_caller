"""
Read-PATTERN statistics: CPEL's device, in the form long reads make available.

WHAT CPEL DOES, AND WHY
-----------------------
Abante, Fang, Feinberg & Goutsias (Nat Commun 2020; CpelAsm.jl) -- cited in this
repo's README as "Jiang et al.", which is wrong, the first author is Abante --
fit a 1-D Ising model to each haplotype's methylation states,

    p(x) prop-to exp( sum_k alpha_k sum_{n in R_k} s_n  +  beta sum_n s_n s_{n+1} ),
    s_n = 2 x_n - 1,

with one field alpha_k per <=500 bp subregion and ONE global coupling beta, then
tests three statistics between the two haplotypes' fitted distributions:

    T_MML = |MML_1 - MML_2|      mean methylation level imbalance
    T_NME = |NME_1 - NME_2|      normalised methylation ENTROPY imbalance
    T_PDM = JSD(p_1, p_2) / N    whole-distribution imbalance

The reason the Ising model is there at all is a WGBS constraint, not a
statistical preference: a ~100 bp bisulfite read observes only a fragment of a
haplotype, so p(x) over N<=20 CpGs is never observed directly and must be
inferred, by marginalising a parametric model over everything each read failed
to see. That inference is what costs CpelAsm 48 h on 20 CPUs per sample.

WHY LONG READS CHANGE THE ANSWER
--------------------------------
A 35 kb ONT read observes EVERY CpG in the region at once. The joint pattern is
observed, not inferred. Three consequences, all implemented here:

  1. The likelihood is complete-data, so fitting is a small concave problem
     solved in milliseconds by transfer matrix + L-BFGS -- no simulated
     annealing, no marginalisation over unseen sites (`fit_chain`).
  2. A NON-PARAMETRIC alternative becomes viable that is not viable in WGBS:
     with reads observing whole patterns you can compare the two haplotypes'
     EMPIRICAL read distributions directly. CPEL's own comparator "NPD"
     (Onuchic et al.) was restricted to 4-CpG epialleles for exactly this
     reason. `ks_read_fractions` is the cheapest version of it.
  3. Coverage is 20-40x per haplotype rather than 5x, which matters because
     CPEL's null is generated at Cmin = 5 REGARDLESS of a region's actual
     coverage (Supplementary section 10, step 5). That is conservative by
     design in WGBS and would be severely conservative at ONT depth. The
     matched null in `null.py` fixes it by conditioning on depth.

WHAT IS AND IS NOT PORTED
-------------------------
Ported: the model form, MML, NME, JSD/PDM, and the read-splitting null.
Deliberately changed: complete-data fitting; optional per-edge coupling that
decays with CpG-CpG distance (CPEL's single beta is defensible over <=500 bp of
WGBS haplotype and is not defensible over a multi-kb long-read region);
depth-matched rather than Cmin-matched null.
Not ported: SNP-cluster haplotype definition. Long reads phase whole chromosome
arms, so regions are regulatory units, not SNP neighbourhoods -- which is
precisely the scope CPEL cannot reach.
"""
from __future__ import annotations

import numpy as np
from scipy.optimize import minimize

_LOG2 = np.log(2.0)


# ---------------------------------------------------------------------------
# binary chain: transfer-matrix inference. States are s in {-1, +1}, indexed
# 0 -> -1 (unmethylated), 1 -> +1 (methylated).
# ---------------------------------------------------------------------------
_S = np.array([-1.0, 1.0])


def _edge_logpot(J):
    """(N-1, 2, 2) log pairwise potentials beta_n * s_n * s_{n+1}."""
    outer = np.outer(_S, _S)                       # (2, 2)
    return np.asarray(J, float)[:, None, None] * outer[None, :, :]


def _forward(h, J, mask=None):
    """Log forward messages. `mask` is (N, 2) or a BATCH (B, N, 2) of 0/-inf
    allowing state clamping (for reads with no-calls, or to score a specific
    configuration). Batching matters: the fit loops over reads inside the
    optimiser, and doing that in Python costs ~20x -- measured 53 ms/region
    unbatched vs 2.6 ms batched, which is the difference between a null that
    takes 30 h per 100k regions and one that takes 1.5 h.

    Returns (alpha, logZ); logZ is a scalar for a single mask and shape (B,)
    for a batch."""
    h = np.asarray(h, float)
    N = len(h)
    site = h[:, None] * _S[None, :]                # (N, 2)
    if mask is None:
        batch = None
        cur = site[0].copy()
    else:
        mask = np.asarray(mask, float)
        batch = mask.ndim == 3
        site = site[None, ...] + (mask if batch else mask[None, ...])
        cur = site[:, 0].copy()                    # (B, 2)
    ep = _edge_logpot(J)
    alphas = [cur]
    for n in range(1, N):
        prev = cur[..., :, None] + ep[n - 1]       # (..., 2, 2)
        mx = prev.max(axis=-2)
        cur = (site[n] if batch is None else site[:, n]) \
            + mx + np.log(np.exp(prev - mx[..., None, :]).sum(axis=-2))
        alphas.append(cur)
    a = np.stack(alphas, axis=-2)                  # (..., N, 2)
    mx = cur.max(axis=-1)
    logZ = mx + np.log(np.exp(cur - mx[..., None]).sum(axis=-1))
    if batch is None and mask is not None:
        return a[0], float(logZ[0])
    if mask is None:
        return a, float(logZ)
    return a, logZ


def _backward(h, J):
    h = np.asarray(h, float)
    N = len(h)
    site = h[:, None] * _S[None, :]
    ep = _edge_logpot(J)
    b = np.zeros((N, 2))
    for n in range(N - 2, -1, -1):
        nxt = (b[n + 1] + site[n + 1])[None, :] + ep[n]   # (2, 2)
        mx = nxt.max(axis=1)
        b[n] = mx + np.log(np.exp(nxt - mx[:, None]).sum(axis=1))
    return b


def _backward_b(h, J, mask=None):
    """Batched log backward messages, same mask convention as `_forward`."""
    h = np.asarray(h, float)
    N = len(h)
    site = h[:, None] * _S[None, :]
    if mask is not None:
        mask = np.asarray(mask, float)
        site = site[None, ...] + (mask if mask.ndim == 3 else mask[None, ...])
        shape = (site.shape[0], N, 2)
    else:
        shape = (N, 2)
    ep = _edge_logpot(J)
    b = np.zeros(shape)
    for n in range(N - 2, -1, -1):
        s_next = site[n + 1] if mask is None else site[:, n + 1]
        b_next = b[n + 1] if mask is None else b[:, n + 1]
        nxt = (b_next + s_next)[..., None, :] + ep[n]
        mx = nxt.max(axis=-1)
        val = mx + np.log(np.exp(nxt - mx[..., None]).sum(axis=-1))
        if mask is None:
            b[n] = val
        else:
            b[:, n] = val
    return b


def _expected_stats(h, J, k_idx, K, mask=None):
    """E[T(x)] for the sufficient statistics (per-subregion field sums, coupling
    sum), summed over the batch. Used for the analytic gradient of the fit --
    a numerical gradient costs ~3x more nll evaluations and dominated the
    runtime before this was added."""
    h = np.asarray(h, float)
    N = len(h)
    a, _ = _forward(h, J, mask=mask)
    b = _backward_b(h, J, mask=mask)
    site = h[:, None] * _S[None, :]
    if mask is not None:
        mask = np.asarray(mask, float)
        site = site[None, ...] + (mask if mask.ndim == 3 else mask[None, ...])
    lp1 = a + b
    lp1 = lp1 - lp1.max(axis=-1, keepdims=True)
    p1 = np.exp(lp1); p1 /= p1.sum(axis=-1, keepdims=True)
    s_exp = p1[..., 1] - p1[..., 0]                       # E[s_n], (..., N)
    g_alpha = np.zeros(K)
    for k in range(K):
        g_alpha[k] = s_exp[..., k_idx == k].sum()
    g_beta = 0.0
    if N > 1:
        ep = _edge_logpot(J)                                  # (N-1, 2, 2)
        outer = np.outer(_S, _S)
        s_all = site if mask is None else site                # (N,2) or (B,N,2)
        rhs = (s_all[..., 1:, :] + b[..., 1:, :])             # (..., N-1, 2)
        lp = a[..., :-1, :, None] + ep + rhs[..., None, :]    # (..., N-1, 2, 2)
        lp = lp - lp.max(axis=(-2, -1), keepdims=True)
        e = np.exp(lp)
        e /= e.sum(axis=(-2, -1), keepdims=True)
        g_beta = float((e * outer).sum())
    return g_alpha, float(g_beta)


def chain_logZ(h, J) -> float:
    return _forward(h, J)[1]


def chain_marginals(h, J):
    """Site marginals (N, 2) and pairwise marginals (N-1, 2, 2)."""
    h = np.asarray(h, float)
    N = len(h)
    a, logZ = _forward(h, J)
    b = _backward(h, J)
    site = h[:, None] * _S[None, :]
    p1 = np.exp(a + b - logZ)
    p1 /= p1.sum(axis=1, keepdims=True)
    ep = _edge_logpot(J)
    p2 = np.empty((N - 1, 2, 2))
    for n in range(N - 1):
        lp = (a[n][:, None] + ep[n]
              + (site[n + 1] + b[n + 1])[None, :] - logZ)
        e = np.exp(lp)
        p2[n] = e / e.sum()
    return p1, p2


def chain_entropy_bits(h, J) -> float:
    """Exact Shannon entropy (bits) of the fitted chain. The chain is Markov, so
    H = H(s_1) + sum_n H(s_{n+1} | s_n), each term closed-form from the
    forward-backward marginals -- no enumeration of the 2^N state space."""
    p1, p2 = chain_marginals(h, J)

    def _h(p):
        p = np.clip(p, 1e-300, 1.0)
        return -float((p * np.log2(p)).sum())

    H = _h(p1[0])
    for n in range(len(p2)):
        H += _h(p2[n]) - _h(p1[n])          # H(s_n, s_{n+1}) - H(s_n)
    return H


def chain_mml(h, J) -> float:
    """E[ mean methylation ] under the fitted chain."""
    p1, _ = chain_marginals(h, J)
    return float(p1[:, 1].mean())


# ---------------------------------------------------------------------------
# fitting, from reads that may have no-calls
# ---------------------------------------------------------------------------

def _subregion_index(positions, max_bp: int = 500) -> np.ndarray:
    """CPEL partitions each haplotype into the minimum number of equally sized
    non-overlapping subregions no larger than 500 bp, one alpha each."""
    positions = np.asarray(positions, float)
    span = positions[-1] - positions[0]
    K = max(1, int(np.ceil((span + 1e-9) / max_bp)))
    width = (span + 1e-9) / K
    idx = np.floor((positions - positions[0]) / width).astype(int)
    return np.clip(idx, 0, K - 1)


def fit_chain(m: np.ndarray, positions=None, max_subregion_bp: int = 500,
              distance_coupling: bool = False, ridge: float = 1e-3,
              max_bp_scale: float = 200.0):
    """MLE of the CPEL model from a read x CpG matrix of 0/1/nan.

    Parameters are (alpha_1..alpha_K, beta): one field per <=500 bp subregion
    and one nearest-neighbour coupling, exactly as in CPEL. With
    `distance_coupling=True` the coupling on edge n is beta * exp(-d_n / scale),
    a long-read extension -- CPEL's constant beta is reasonable over a WGBS
    haplotype and is not over a multi-kb region with 10x-varying CpG spacing.

    Returns (h, J, loglik, theta). `ridge` is a small L2 penalty that keeps the
    fit finite when a region is fully methylated on every read (which is the
    single most common configuration in a real bimodal methylome, and where an
    unpenalised MLE diverges to alpha = inf).
    """
    m = np.asarray(m, float)
    R, N = m.shape
    if positions is None:
        positions = np.arange(N) * 50.0
    positions = np.asarray(positions, float)

    if N == 1:
        k_idx = np.zeros(1, int)
        d = np.zeros(0)
    else:
        k_idx = _subregion_index(positions, max_subregion_bp)
        d = np.diff(positions)
    K = int(k_idx.max()) + 1

    # per-read clamping masks: 0 where a state is allowed, -inf where excluded
    NEG = -1e9
    masks = np.zeros((R, N, 2))
    obs = np.isfinite(m)
    masks[..., 0] = np.where(obs & (m == 1), NEG, 0.0)
    masks[..., 1] = np.where(obs & (m == 0), NEG, 0.0)

    def unpack(theta):
        h = np.asarray(theta[:K])[k_idx]
        beta = theta[K]
        if N == 1:
            return h, np.zeros(0)
        J = beta * np.exp(-d / max_bp_scale) if distance_coupling \
            else np.full(N - 1, beta)
        return h, J

    usable = obs.any(axis=1)
    masks_u = masks[usable]

    Ru = int(usable.sum())

    def nll_grad(theta):
        h, J = unpack(theta)
        _, logZ = _forward(h, J)
        if Ru == 0:
            return ridge * float(np.dot(theta, theta)), 2 * ridge * theta
        _, lz = _forward(h, J, mask=masks_u)
        f = -float((lz - logZ).sum()) + ridge * float(np.dot(theta, theta))
        ga_c, gb_c = _expected_stats(h, J, k_idx, K, mask=masks_u)
        ga_f, gb_f = _expected_stats(h, J, k_idx, K, mask=None)
        g = np.empty(K + 1)
        # d/dtheta of -(sum_r logZ_r_clamped - R logZ) = -(E_clamped[T] - R E[T])
        g[:K] = -(ga_c - Ru * ga_f)
        if distance_coupling and N > 1:
            # beta enters each edge scaled by exp(-d/scale); the expected-stat
            # gradient above is w.r.t. the per-edge couplings, so the chain rule
            # factor is folded in by recomputing numerically for this one param
            eps = 1e-5
            t2 = np.array(theta, float); t2[K] += eps
            h2, J2 = unpack(t2)
            _, lz2 = _forward(h2, J2, mask=masks_u)
            _, lZ2 = _forward(h2, J2)
            f2 = -float((lz2 - lZ2).sum()) + ridge * float(np.dot(t2, t2))
            g[K] = (f2 - f) / eps
            g[:K] += 2 * ridge * theta[:K]
            return f, g
        g[K] = -(gb_c - Ru * gb_f)
        g += 2 * ridge * theta
        return f, g

    x0 = np.zeros(K + 1)
    mu0 = float(np.nanmean(m)) if np.isfinite(m).any() else 0.5
    x0[:K] = 0.5 * np.log(np.clip(mu0, 0.02, 0.98) / (1 - np.clip(mu0, 0.02, 0.98)))
    res = minimize(nll_grad, x0, method="L-BFGS-B", jac=True,
                   bounds=[(-6, 6)] * K + [(-4, 4)],
                   options={"maxiter": 200})
    h, J = unpack(res.x)
    return h, J, float(-res.fun), res.x


# ---------------------------------------------------------------------------
# CPEL's three test statistics, computed per region from read matrices
# ---------------------------------------------------------------------------

def _sample_chain(rng, h, J, n):
    """Backward-filtered forward sampling from the binary chain."""
    N = len(h)
    b = _backward(h, J)
    site = np.asarray(h, float)[:, None] * _S[None, :]
    lp0 = site[0] + b[0]
    p = np.exp(lp0 - lp0.max()); p /= p.sum()
    out = np.empty((n, N), int)
    out[:, 0] = rng.random(n) < p[1]
    ep = _edge_logpot(J)
    for k in range(1, N):
        lp = ep[k - 1] + (site[k] + b[k])[None, :]     # (2 prev, 2 cur)
        e = np.exp(lp - lp.max(axis=1, keepdims=True))
        e /= e.sum(axis=1, keepdims=True)
        out[:, k] = rng.random(n) < e[out[:, k - 1], 1]
    return out


def _chain_logpmf(x, h, J):
    """log p(x) for integer 0/1 configurations x of shape (n, N)."""
    s = 2.0 * np.asarray(x, float) - 1.0
    lp = s @ np.asarray(h, float)
    if len(J):
        lp = lp + (s[:, :-1] * s[:, 1:]) @ np.asarray(J, float)
    return lp - chain_logZ(h, J)


def jsd_chains(h1, J1, h2, J2, n_mc: int = 400, seed: int = 0) -> float:
    """Jensen-Shannon DIVERGENCE (bits, in [0,1]) between two fitted chains.

    Estimated by Monte Carlo rather than enumerating 2^N: CPEL caps N at 20
    (10^6 states) and needs the enumeration; sampling is exact in the limit and
    is what makes this affordable inside a genome-wide benchmark.
    """
    rng = np.random.default_rng(seed)
    x1 = _sample_chain(rng, h1, J1, n_mc)
    x2 = _sample_chain(rng, h2, J2, n_mc)

    def _kl_to_mix(x, ha, Ja, hb, Jb):
        la = _chain_logpmf(x, ha, Ja)
        lb = _chain_logpmf(x, hb, Jb)
        mx = np.maximum(la, lb)
        lm = mx + np.log(0.5 * np.exp(la - mx) + 0.5 * np.exp(lb - mx))
        return float(np.mean(la - lm)) / _LOG2

    return float(np.clip(0.5 * _kl_to_mix(x1, h1, J1, h2, J2)
                         + 0.5 * _kl_to_mix(x2, h2, J2, h1, J1), 0.0, 1.0))


def region_statistics(m1, m2, positions=None, n_mc: int = 400,
                      distance_coupling: bool = False, seed: int = 0) -> dict:
    """CPEL's T_MML, T_NME and T_PDM for one region, plus the fitted models.

    NME is the entropy NORMALISED BY N, so it is comparable across regions with
    different CpG counts -- that normalisation is what makes a single null
    distribution per N meaningful.
    """
    n_cpgs = m1.shape[1]
    h1, J1, ll1, _ = fit_chain(m1, positions, distance_coupling=distance_coupling)
    h2, J2, ll2, _ = fit_chain(m2, positions, distance_coupling=distance_coupling)
    mml1, mml2 = chain_mml(h1, J1), chain_mml(h2, J2)
    nme1 = chain_entropy_bits(h1, J1) / n_cpgs
    nme2 = chain_entropy_bits(h2, J2) / n_cpgs
    pdm = jsd_chains(h1, J1, h2, J2, n_mc=n_mc, seed=seed) if n_cpgs > 1 else 0.0
    return {
        "t_mml": abs(mml1 - mml2), "t_nme": abs(nme1 - nme2), "t_pdm": pdm,
        "mml1": mml1, "mml2": mml2, "nme1": nme1, "nme2": nme2,
        "n_cpgs": n_cpgs, "n_reads1": int(np.isfinite(m1).any(axis=1).sum()),
        "n_reads2": int(np.isfinite(m2).any(axis=1).sum()),
    }


# ---------------------------------------------------------------------------
# the cheap non-parametric alternative long reads make possible
# ---------------------------------------------------------------------------

def read_fraction_ks(m1, m2, min_calls: int = 1, min_reads: int = 3):
    """Two-sample Kolmogorov-Smirnov statistic on per-read methylation fractions.

    The point of including this: the read-level test already in `readlevel.py`
    compares LOCATIONS (Welch t / Mann-Whitney), so it is blind by construction
    to a haplotype pair with equal mean and different disorder -- one allele
    fixed at 50% methylated, the other bistable between fully methylated and
    fully unmethylated molecules. KS compares the whole fraction distribution
    and catches that, at the cost of one line and no model.

    If KS recovers most of what T_PDM recovers at ONT depth, the Ising machinery
    is not earning its keep in this setting. That is the comparison the
    benchmark is built to make.
    """
    from scipy.stats import ks_2samp
    from .readlevel import read_fractions
    f1, f2 = read_fractions(m1, min_calls), read_fractions(m2, min_calls)
    if len(f1) < min_reads or len(f2) < min_reads:
        return None
    r = ks_2samp(f1, f2, alternative="two-sided", method="asymp")
    return float(r.statistic), float(r.pvalue)
