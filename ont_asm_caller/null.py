"""
Depth-matched empirical null with tail extrapolation.

THE PROBLEM THIS SOLVES
-----------------------
Two tests in this package are correctly calibrated and unusable genome-wide for
the same reason:

  * `readlevel.test_region_perm` -- p >= 1/(n_perm+1).
  * CPEL's own bootstrap -- p_hat = (1 + #{t_l >= t*}) / (L+1), L >= 1000
    (Abante et al. 2020, Supplementary section 10, Eq. 73), so p >= ~1e-3.

Under BH across ~10^5-10^6 regions a p-value floored at 1e-3 can only yield a
discovery if TENS OF THOUSANDS of regions pile up at the floor simultaneously.
CPEL got away with it on 715,155 haplotypes; it also means every one of its
discoveries sits at the floor and none of them can be ranked against each other.
A permutation test is worth nothing at genome scale unless its tail can be
extended below 1/(L+1).

TWO FIXES, BOTH NEEDED
----------------------
1. POOL THE NULL ACROSS REGIONS, WITHIN A STRATUM. Under H0 a region's null
   statistic distribution depends on the region's shape, not its identity, so
   null draws from comparable regions are exchangeable and can be pooled. This
   buys L ~ 10^5 instead of 10^3 for the same compute.

   Which variables define "comparable" is the whole argument. CPEL conditions on
   the number of CpG sites N ALONE, and generates every null draw at the global
   MINIMUM coverage C_min = 5 (their step 5), a choice that is conservative in
   WGBS and would be crushing at ONT depth -- 5 reads against a region's actual
   20-40 gives a null far wider than the truth, so real signal is discarded.
   Here the strata are (n_cpgs, depth, methylation level). Depth matters for the
   obvious reason. Methylation level matters because every one of these
   statistics has a null spread that collapses as mu approaches 0 or 1 -- the
   same mean-dependence that made `dispersion.estimate_dispersion` fail at the
   boundaries (docs/2026-09-05_calibration_critique.md, Defect 2). CPEL
   conditions on neither.

2. EXTRAPOLATE THE TAIL. Fit a generalised Pareto distribution to the
   exceedances over a high threshold and read p-values off the fitted tail below
   the empirical floor. This is the standard peaks-over-threshold treatment of
   permutation p-values (Knijnenburg et al., Bioinformatics 2009, "Fewer
   permutations, more accurate P-values"), not something invented here.

WHAT IS ASSUMED, AND HOW TO CHECK IT
------------------------------------
Pooling assumes within-stratum exchangeability; extrapolation assumes the tail
is regularly varying. Both are checkable and neither is free. `calibration_report`
runs the only test that matters: are p-values uniform on data simulated under
H0? Run it before believing any q-value that came out of here.
"""
from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
from scipy.stats import genpareto


def split_reads(rng, m1: np.ndarray, m2: np.ndarray):
    """Pool the two haplotypes' reads and re-split at the observed sizes.

    WHOLE READS move between groups, so within-molecule co-methylation is
    preserved exactly under the null -- which is the entire point, and the
    reason this is valid where a count-level bootstrap (scripts/phasing_qc/
    P04_permutation_null.py) is not."""
    pool = np.vstack([m1, m2])
    idx = rng.permutation(len(pool))
    n1 = len(m1)
    return pool[idx[:n1]], pool[idx[n1:]]


def _depth_bin(n1, n2, edges):
    return int(np.searchsorted(edges, min(n1, n2), side="right"))


def _mu_bin(m1, m2, edges):
    vals = np.concatenate([m1.ravel(), m2.ravel()])
    vals = vals[np.isfinite(vals)]
    mu = float(vals.mean()) if len(vals) else 0.5
    return int(np.searchsorted(edges, mu, side="right"))


@dataclass
class MatchedNull:
    """Stratified null distributions for one or more region statistics.

    `stat_fn(m1, m2) -> dict[str, float]` is called on permuted read splits.
    """
    stat_fn: object
    n_perm_per_region: int = 20
    depth_edges: tuple = (8, 15, 25, 40)
    mu_edges: tuple = (0.15, 0.35, 0.65, 0.85)
    cpg_edges: tuple = (3, 5, 8, 12, 20)
    tail_quantile: float = 0.90
    min_stratum: int = 500
    max_extrapolation: float = 100.0   # how far past the empirical floor the
                                       # GPD tail is trusted; see `pvalue`
    stratify: tuple = ("cpgs", "depth", "mu")
    cpel_cmin: int | None = None      # emulate CPEL: force every null draw to
    cpel_dcov: int = 2                # coverage in [Cmin, Cmin+dCov], any depth
    _draws: dict = field(default_factory=dict)
    _fits: dict = field(default_factory=dict)

    # -- stratum key -------------------------------------------------------
    def key(self, m1, m2):
        k = []
        if "cpgs" in self.stratify:
            k.append(int(np.searchsorted(self.cpg_edges, m1.shape[1], side="right")))
        if "depth" in self.stratify:
            k.append(_depth_bin(len(m1), len(m2), self.depth_edges))
        if "mu" in self.stratify:
            k.append(_mu_bin(m1, m2, self.mu_edges))
        return tuple(k)

    # -- build -------------------------------------------------------------
    def build(self, regions, seed: int = 0, max_regions: int | None = None,
              progress: bool = False):
        """Accumulate null draws by permuting reads within each supplied region.

        `regions` should be objects with `.m1` / `.m2`. In real data pass the
        SAME regions being tested: true ASM is rare enough that including it
        biases the null slightly WIDE, i.e. conservative -- the same argument
        `dispersion.estimate_dispersion` already makes for pooling across loci,
        and it fails the same way if real prevalence turns out to be high."""
        rng = np.random.default_rng(seed)
        regs = regions if max_regions is None else regions[:max_regions]
        for i, r in enumerate(regs):
            if progress and i % 2000 == 0:
                print(f"    null: {i}/{len(regs)}", flush=True)
            k = self.key(r.m1, r.m2)
            bucket = self._draws.setdefault(k, {})
            for _ in range(self.n_perm_per_region):
                a, b = split_reads(rng, r.m1, r.m2)
                if self.cpel_cmin is not None:
                    # CPEL Supplementary s10 step 5: reduce BOTH groups to a
                    # coverage in [Cmin, Cmin+dCov] regardless of what the region
                    # actually has. Reproduced here to measure what it costs at
                    # ONT depth, not because it is a good idea there.
                    c = int(rng.integers(self.cpel_cmin,
                                         self.cpel_cmin + self.cpel_dcov + 1))
                    if len(a) < c or len(b) < c:
                        continue
                    a = a[rng.permutation(len(a))[:c]]
                    b = b[rng.permutation(len(b))[:c]]
                try:
                    s = self.stat_fn(a, b)
                except Exception:
                    continue
                if s is None:
                    continue
                for name, v in s.items():
                    if np.isfinite(v):
                        bucket.setdefault(name, []).append(float(v))
        self._finalise()
        return self

    def _finalise(self):
        """Sort each stratum's draws and fit a GPD to the upper tail."""
        self._fits = {}
        # fall back to the marginal (all-strata) null for thin strata
        pooled: dict = {}
        for k, bucket in self._draws.items():
            for name, vals in bucket.items():
                pooled.setdefault(name, []).extend(vals)
        self._draws[("*",)] = {n: v for n, v in pooled.items()}

        for k, bucket in self._draws.items():
            for name, vals in bucket.items():
                v = np.sort(np.asarray(vals, float))
                if len(v) < 50:
                    continue
                u = float(np.quantile(v, self.tail_quantile))
                exc = v[v > u] - u
                fit = None
                if len(exc) >= 200:
                    try:
                        c, loc, scale = genpareto.fit(exc, floc=0.0)
                        # A pathological shape parameter is the failure mode
                        # that matters. Measured before this guard: thin strata
                        # produced fits that extrapolated to p ~ 1e-300, and the
                        # extreme-tail type I error broke (min p under H0 also
                        # 1e-300) while alpha=0.05 and 0.01 still looked fine.
                        # Every statistic here is bounded in [0, 1], so a large
                        # positive shape -- a polynomial tail with no upper
                        # endpoint -- is structurally misspecified.
                        if np.isfinite(c) and scale > 0 and -0.5 <= c <= 0.25:
                            fit = (c, scale)
                    except Exception:
                        fit = None
                self._fits[(k, name)] = (v, u, float((v > u).mean()), fit)

    # -- query -------------------------------------------------------------
    def pvalue(self, name: str, stat: float, m1=None, m2=None, key=None) -> float:
        """One-sided upper-tail p-value, extrapolated below the empirical floor."""
        if key is None:
            key = self.key(m1, m2)
        entry = self._fits.get((key, name))
        if entry is None or len(entry[0]) < self.min_stratum:
            entry = self._fits.get((("*",), name))
        if entry is None:
            return float("nan")
        v, u, p_exceed, fit = entry
        L = len(v)
        if not np.isfinite(stat):
            return float("nan")
        if stat <= u or fit is None:
            # plain empirical tail, floored at 1/(L+1) as usual
            return float((1 + int((v >= stat).sum())) / (L + 1))
        c, scale = fit
        tail = float(genpareto.sf(stat - u, c, loc=0.0, scale=scale))
        # CAP THE EXTRAPOLATION. The empirical floor is 1/(L+1); a GPD fitted to
        # the top 10% of L draws supports extending that by a couple of orders
        # of magnitude, not by three hundred. Trusting it further is how the
        # extreme tail broke while alpha=0.05 still looked calibrated -- and the
        # extreme tail is exactly where genome-wide BH operates.
        floor = 1.0 / (self.max_extrapolation * (L + 1))
        return float(min(1.0, max(p_exceed * tail, floor)))

    def sizes(self):
        return {k: {n: len(v) for n, v in b.items()} for k, b in self._draws.items()}


def calibration_report(null: MatchedNull, regions, stat_fn, names,
                       seed: int = 99, max_regions: int | None = 4000) -> dict:
    """Type I error of the matched null on H0 data: permute each region once
    more (a fresh split, independent of the ones used to build the null) and
    check the resulting p-values are uniform.

    This is the check that decides whether the pooling and the tail fit are
    doing anything dishonest. A stratified null that is silently mismatched
    shows up here as a non-uniform CDF, not as an error."""
    rng = np.random.default_rng(seed)
    regs = regions if max_regions is None else regions[:max_regions]
    out = {n: [] for n in names}
    for r in regs:
        a, b = split_reads(rng, r.m1, r.m2)
        try:
            s = stat_fn(a, b)
        except Exception:
            continue
        if s is None:
            continue
        k = null.key(a, b)
        for n in names:
            if n in s and np.isfinite(s[n]):
                out[n].append(null.pvalue(n, s[n], key=k))
    rep = {}
    for n, ps in out.items():
        p = np.asarray(ps, float)
        p = p[np.isfinite(p)]
        rep[n] = {
            "n": int(len(p)),
            "type1_0.05": float((p < 0.05).mean()) if len(p) else float("nan"),
            "type1_0.01": float((p < 0.01).mean()) if len(p) else float("nan"),
            "type1_0.001": float((p < 0.001).mean()) if len(p) else float("nan"),
            # the extreme tail is where BH across 10^5-10^6 regions actually
            # operates, and it is where tail extrapolation fails first -- an
            # alpha=0.05 row that looks fine tells you nothing about it
            "type1_1e-4": float((p < 1e-4).mean()) if len(p) else float("nan"),
            "min_p": float(p.min()) if len(p) else float("nan"),
        }
    return rep
