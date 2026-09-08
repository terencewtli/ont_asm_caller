"""
Do CPEL's pattern statistics buy anything over this package's tests, on
long-read-shaped data?  And is a cheap non-parametric alternative enough?

DESIGN (the part that matters more than the numbers)
----------------------------------------------------
A benchmark run on data simulated from a method's own model tells you nothing.
`simulate_reads.py` draws exchangeable within-read correlation, under which the
per-read methylation FRACTION is a sufficient statistic -- so a pattern method
is guaranteed to find nothing extra there, and this package would "win" by
construction. Simulating from an Ising chain rigs it the other way, because that
is exactly CPEL's model.

So all three are run:

  epiallele    -- mixture of discrete methylation patterns. NEITHER method's
                  model. This is the headline scenario.
  markov       -- distance-decaying nearest-neighbour correlation = a 1-D Ising
                  chain. CPEL is correctly specified. Upper bound for CPEL.
  exchangeable -- the simulate_reads.py model. Read-fraction tests are
                  sufficient. Upper bound for this package.

and TWO orthogonal ASM axes are simulated with independent prevalence, so power
can be broken out by which axis a true positive actually carries:

  MML axis -- haplotypes differ in mean.        Every existing test targets this.
  ICC axis -- haplotypes differ in within-read CORRELATION at EQUAL MEAN
              (one allele fixed, the other bistable). No existing test in this
              package can see it, because they all compare locations.

METHODS
-------
  pooled_bb        region.py's pooled read x CpG beta-binomial LRT (v1 path)
  readfrac_welch   readlevel.py's Welch/Mann-Whitney on read fractions (v2 path)
  ks_matched       two-sample KS on read fractions + matched null  [NEW, cheap]
  t_mml / t_nme / t_pdm    CPEL's three statistics + matched null   [NEW]
  t_pdm_cpelnull   CPEL's statistic scored against CPEL'S OWN NULL PROTOCOL --
                   stratified on CpG count only and generated at Cmin=5 reads
                   regardless of the region's real depth. Included to separate
                   "the statistic is weak" from "the null is mismatched".

Run:  python benchmarks/compare_pattern_methods.py [n_regions] [kind ...]

With no `kind` all three run in sequence (~20 min each, single-core). Naming
one runs just that scenario, so the three can be run in parallel on separate
cores -- which is how the committed numbers were produced.
"""
from __future__ import annotations

import sys, os, time
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(
    globals().get("__file__", "benchmarks/x")))))

import numpy as np

from ont_asm_caller.simulate_patterns import simulate_pattern_regions, per_cpg_counts
from ont_asm_caller.dispersion import estimate_dispersion
from ont_asm_caller.model import test_locus
from ont_asm_caller.readlevel import test_region_reads
from ont_asm_caller import pattern as P
from ont_asm_caller.null import MatchedNull, calibration_report

N_MC = 200
N_PERM = 15
NULL_REGIONS = 3000


def bh(p, alpha=0.05):
    p = np.asarray(p, float)
    n = len(p)
    if n == 0:
        return np.zeros(0, bool)
    o = np.argsort(p)
    passed = p[o] <= alpha * (np.arange(1, n + 1) / n)
    k = np.max(np.nonzero(passed)[0]) + 1 if passed.any() else 0
    out = np.zeros(n, bool)
    out[o[:k]] = True
    return out


def stat_fn(m1, m2, positions=None):
    """All test statistics for one region, as one dict, so the null only has to
    permute once for all of them."""
    out = {}
    ks = P.read_fraction_ks(m1, m2)
    if ks is not None:
        out["ks"] = ks[0]
    try:
        st = P.region_statistics(m1, m2, positions, n_mc=N_MC)
        out.update({k: st[k] for k in ("t_mml", "t_nme", "t_pdm")})
    except Exception:
        pass
    return out or None


def run(kind, n_regions, pi_mml=0.01, pi_icc=0.01, seed=7, verbose=True):
    t0 = time.time()
    regs = simulate_pattern_regions(n_regions, kind=kind, pi_mml=pi_mml,
                                    pi_icc=pi_icc, seed=seed)
    has_mml = np.array([r.true_delta_mml != 0 for r in regs])
    has_icc = np.array([r.true_delta_icc != 0 for r in regs])
    truth = has_mml | has_icc

    print(f"\n{'='*96}")
    print(f"kind={kind}   regions={n_regions}   true ASM: "
          f"{truth.sum()} ({has_mml.sum()} mean-axis, {has_icc.sum()} correlation-axis, "
          f"{(has_mml & has_icc).sum()} both)")
    print(f"{'='*96}")

    # ---- observed statistics -------------------------------------------
    obs, keys = [], []
    for r in regs:
        obs.append(stat_fn(r.m1, r.m2, r.positions))
        keys.append(r)
    if verbose:
        print(f"  observed statistics: {time.time()-t0:.0f}s", flush=True)

    # ---- matched null ---------------------------------------------------
    mn = MatchedNull(lambda a, b: stat_fn(a, b), n_perm_per_region=N_PERM)
    mn.build(regs, seed=11, max_regions=NULL_REGIONS)
    cp = MatchedNull(lambda a, b: stat_fn(a, b), n_perm_per_region=N_PERM,
                     stratify=("cpgs",), cpel_cmin=5)
    cp.build(regs, seed=12, max_regions=NULL_REGIONS)
    if verbose:
        print(f"  nulls built: {time.time()-t0:.0f}s", flush=True)

    # ---- baselines that carry their own analytic p-values ---------------
    pvals = {}
    pooled = [r.counts() for r in regs]
    rho = estimate_dispersion(*map(np.array, zip(*pooled)))
    pvals["pooled_bb"] = np.array(
        [test_locus("chr1", r.start, *c, rho).pval for r, c in zip(regs, pooled)])
    wl = [test_region_reads("chr1", r.start, r.start, r.m1, r.m2) for r in regs]
    pvals["readfrac_welch"] = np.array([1.0 if x is None else x.pval for x in wl])

    # ---- null-calibrated statistics -------------------------------------
    for name, label, null in [("ks", "ks_matched", mn),
                              ("t_mml", "t_mml", mn),
                              ("t_nme", "t_nme", mn),
                              ("t_pdm", "t_pdm", mn),
                              ("t_pdm", "t_pdm_cpelnull", cp)]:
        p = np.ones(len(regs))
        for i, (s, r) in enumerate(zip(obs, keys)):
            if s is not None and name in s:
                v = null.pvalue(name, s[name], r.m1, r.m2)
                p[i] = 1.0 if not np.isfinite(v) else v
        pvals[label] = p

    # ---- score ----------------------------------------------------------
    print(f"\n  {'method':<20}{'called':>8}{'FDR':>8}{'power':>8}"
          f"{'pow(mean)':>11}{'pow(corr)':>11}{'min p':>11}")
    only_mml = has_mml & ~has_icc
    only_icc = has_icc & ~has_mml
    for label, p in pvals.items():
        sig = bh(p)
        n = int(sig.sum())
        fdr = float((~truth[sig]).mean()) if n else float("nan")
        pw = float(sig[truth].mean()) if truth.any() else float("nan")
        pm = float(sig[only_mml].mean()) if only_mml.any() else float("nan")
        pc = float(sig[only_icc].mean()) if only_icc.any() else float("nan")
        print(f"  {label:<20}{n:>8d}{fdr:>8.3f}{pw:>8.3f}{pm:>11.3f}{pc:>11.3f}"
              f"{p.min():>11.2e}")

    # ---- null calibration ------------------------------------------------
    rep = calibration_report(mn, regs, lambda a, b: stat_fn(a, b),
                             ["ks", "t_mml", "t_nme", "t_pdm"],
                             max_regions=min(2500, n_regions))
    print(f"\n  matched-null type I error on H0 splits (nominal in parentheses)")
    print(f"  {'stat':<12}{'n':>7}{'@0.05':>9}{'@0.01':>9}{'@0.001':>10}{'min p':>11}")
    for n_, d in rep.items():
        print(f"  {n_:<12}{d['n']:>7d}{d['type1_0.05']:>9.3f}{d['type1_0.01']:>9.3f}"
              f"{d['type1_0.001']:>10.4f}{d['min_p']:>11.2e}")
    print(f"\n  total {time.time()-t0:.0f}s")
    return pvals


if __name__ == "__main__":
    n = int(sys.argv[1]) if len(sys.argv) > 1 else 8000
    kinds = sys.argv[2:] or ["epiallele", "markov", "exchangeable"]
    print("Pattern-method benchmark. BH-FDR alpha=0.05.")
    print("Two ASM axes simulated independently: mean imbalance and")
    print("within-read CORRELATION imbalance at equal mean.")
    for kind in kinds:
        run(kind, n)
