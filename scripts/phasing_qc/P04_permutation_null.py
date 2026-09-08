#!/usr/bin/env python3
"""
Donor-specific parametric-bootstrap null for the beta-binomial ASM caller.

Question: given a donor's OWN already-fitted global rho_hat and the REAL
per-region depths (n1, n2) on chr1, what fraction of regions would the
caller call "significant" purely by chance, if there were truly zero ASM
and zero haplotype-labeling problems? This is a donor-specific, learned
null -- not an external population benchmark (deCODE 0.51%) or the bare
chi-square(1) asymptotic -- so it directly separates two hypotheses for a
donor with an elevated real pct_significant:

  (a) the caller's own dispersion/calibration is the problem for this
      donor's depth distribution -> the bootstrap null itself comes back
      elevated, matching the real pct_significant.
  (b) the input haplotype labels (x1/x2) are the problem (e.g. switch-error
      driven, non-random with respect to something else) -> the bootstrap
      null comes back sane (~FDR target), while the real data stays
      elevated -- the caller is fine, the inputs aren't.

Method: for each region, under H0 (shared true rate p_hat = (x1+x2)/(n1+n2),
the SAME rho_hat already estimated for that donor/chrom by the production
pipeline), simulate x1*, x2* ~ BetaBinomial(n, p_hat, rho_hat) independently
per haplotype, using the real n1/n2 depths. Re-run the exact same LRT +
BH-FDR + delta-floor significance rule from G02_betabinom_region_chr1.py.
Repeat for N_BOOT replicates and report the mean/sd null pct_significant.

Usage: python3 P04_permutation_null.py <SAMPLE> <CHROM> [n_boot]
Input:  tables/dasm_betabinom/{SAMPLE}/regions_{CHROM}.tsv + summary_{CHROM}.json
Output: tables/phasing_qc/permutation_null_{SAMPLE}_{CHROM}.json
"""
import json
import sys

import numpy as np
import pandas as pd
from scipy.stats import betabinom

sys.path.insert(0, "/u/project/cluo/terencew/claude/project_ideas/asm_lr/github/ont_asm_caller")
from ont_asm_caller.model import _bb_ab

PROJDIR = "/u/project/cluo/terencew/claude/project_ideas/asm_lr"
FDR_THRESH = 0.05
DELTA_FLOOR = 0.10
N_BOOT_DEFAULT = 10
RNG_SEED = 0


def bh_fdr(pvals):
    n = len(pvals)
    order = np.argsort(pvals)
    ranked = pvals[order]
    q = ranked * n / (np.arange(n) + 1)
    q = np.minimum.accumulate(q[::-1])[::-1]
    out = np.empty(n)
    out[order] = np.clip(q, 0, 1)
    return out


_EPS = 1e-6


def _fit_mu_vec(rho, x_a, n_a, x_b=None, n_b=None, n_iter=60):
    """Vectorized golden-section search for the per-locus mu that maximizes
    logpmf(x_a, n_a, mu, rho) [+ logpmf(x_b, n_b, mu, rho) if given], over
    arrays of loci at once. Replaces model.fit_mu's per-locus
    scipy.optimize.minimize_scalar call -- at chromosome scale (~2e5
    regions x up to 10 bootstrap replicates), the per-locus scalar
    optimizer is too slow (~1e6+ individual optimizer calls); beta-binomial
    log-likelihood in mu is unimodal here, so golden-section search
    vectorizes cleanly."""
    a_shape = mu = np.full_like(x_a, 0.5, dtype=float)  # placeholder, overwritten below
    rho_c = max(rho, 1e-10)

    def negll(m):
        a = m * (1.0 - rho_c) / rho_c
        b = (1.0 - m) * (1.0 - rho_c) / rho_c
        ll = betabinom.logpmf(x_a, n_a, a, b)
        if x_b is not None:
            ll = ll + betabinom.logpmf(x_b, n_b, a, b)
        return -ll

    lo = np.full(x_a.shape, _EPS, dtype=float)
    hi = np.full(x_a.shape, 1 - _EPS, dtype=float)
    gr = (np.sqrt(5) - 1) / 2
    for _ in range(n_iter):
        c = hi - gr * (hi - lo)
        d = lo + gr * (hi - lo)
        fc, fd = negll(c), negll(d)
        greater = fc > fd
        lo = np.where(greater, c, lo)
        hi = np.where(greater, hi, d)
    mu_hat = (lo + hi) / 2.0
    ll_hat = -negll(mu_hat)
    return mu_hat, ll_hat


def lrt_pvals_vectorized(x1, n1, x2, n2, rho):
    """Vectorized version of model.test_locus's LRT over arrays of loci."""
    from scipy.stats import chi2 as chi2dist

    _, ll1 = _fit_mu_vec(rho, x1, n1)
    _, ll2 = _fit_mu_vec(rho, x2, n2)
    _, ll0 = _fit_mu_vec(rho, x1, n1, x2, n2)
    lr = np.maximum(0.0, 2.0 * (ll1 + ll2 - ll0))
    return chi2dist.sf(lr, df=1)


def simulate_bb(n, p_hat, rho, rng):
    a, b = _bb_ab(p_hat, rho)
    # betabinom.rvs needs per-element a,b if p_hat is an array
    return betabinom.rvs(n, a, b, random_state=rng)


def main():
    sample, chrom = sys.argv[1], sys.argv[2]
    n_boot = int(sys.argv[3]) if len(sys.argv) > 3 else N_BOOT_DEFAULT

    regdir = f"{PROJDIR}/tables/dasm_betabinom/{sample}"
    regions = pd.read_csv(f"{regdir}/regions_{chrom}.tsv", sep="\t")
    summary = json.load(open(f"{regdir}/summary_{chrom}.json"))
    rho_hat = summary["rho_hat"]
    real_pct_sig = summary["pct_significant"]

    n1 = regions["n1"].to_numpy(dtype=float)
    n2 = regions["n2"].to_numpy(dtype=float)
    x1 = regions["x1"].to_numpy(dtype=float)
    x2 = regions["x2"].to_numpy(dtype=float)
    p_hat = np.clip((x1 + x2) / (n1 + n2), 1e-6, 1 - 1e-6)

    n1i = n1.astype(int)
    n2i = n2.astype(int)

    print(f"{sample} {chrom}: {len(regions):,} regions, rho_hat={rho_hat:.4f}, "
          f"real pct_sig={real_pct_sig:.3f}%  -- running {n_boot} bootstrap replicates",
          flush=True)

    rng = np.random.default_rng(RNG_SEED)
    boot_pct = []
    for b in range(n_boot):
        x1_sim = simulate_bb(n1i, p_hat, rho_hat, rng)
        x2_sim = simulate_bb(n2i, p_hat, rho_hat, rng)
        pvals = lrt_pvals_vectorized(x1_sim, n1i, x2_sim, n2i, rho_hat)
        q = bh_fdr(pvals)
        mu1_sim = x1_sim / n1
        mu2_sim = x2_sim / n2
        delta_sim = np.abs(mu1_sim - mu2_sim)
        sig = (q < FDR_THRESH) & (delta_sim >= DELTA_FLOOR)
        pct = 100.0 * sig.mean()
        boot_pct.append(pct)
        print(f"  replicate {b+1}/{n_boot}: null pct_sig={pct:.4f}%", flush=True)

    boot_pct = np.array(boot_pct)
    out = {
        "sample": sample, "chrom": chrom, "n_regions": int(len(regions)),
        "rho_hat": rho_hat,
        "real_pct_significant": real_pct_sig,
        "null_pct_significant_mean": float(boot_pct.mean()),
        "null_pct_significant_sd": float(boot_pct.std()),
        "null_pct_significant_all": boot_pct.tolist(),
        "ratio_real_over_null": float(real_pct_sig / boot_pct.mean()) if boot_pct.mean() > 0 else None,
        "n_boot": n_boot,
    }
    outdir = f"{PROJDIR}/tables/phasing_qc"
    import os
    os.makedirs(outdir, exist_ok=True)
    with open(f"{outdir}/permutation_null_{sample}_{chrom}.json", "w") as f:
        json.dump(out, f, indent=2)
    print(json.dumps(out, indent=2))


if __name__ == "__main__":
    main()
