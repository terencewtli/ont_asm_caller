# ont_asm_caller

A purpose-built statistical model for calling allele-specific methylation
(ASM) from single-individual long-read (ONT) data — built after the
project's original per-CpG Fisher-exact-test pipeline turned out to be
badly miscalibrated, and after repurposing DSS (a multi-replicate DMR tool)
for this single-individual problem turned out to be conceptually mismatched
and, on real data, ran indefinitely with no output.

## The problem

Long-read sequencing lets you assign individual reads to a haplotype (via
phased heterozygous SNPs) and ask, at each CpG: does haplotype 1's
methylation rate differ from haplotype 2's? The obvious per-locus test is a
2×2 Fisher's exact test on (methylated, unmethylated) × (HP1, HP2) read
counts. Running this on real project data (an 18-donor HPRC cohort)
surfaced a serious problem:

- On the same chromosome, with identical thresholds, one sample called
  **113,119** loci "significant" (9.8% of tested CpGs) while another called
  **2** — a >50,000× difference. The two samples' *raw* effect sizes
  (before any test) were actually similar; what diverged was purely a
  function of sequencing depth (median 18x vs. 8x per haplotype).
- Fisher's exact test doesn't model overdispersion (biological/technical
  variance beyond simple binomial sampling), so its p-values get more
  confident — not more honest — as depth increases. A test that
  conflates "low depth" with "no difference" cannot be compared across
  samples of different coverage, which real long-read cohorts always have.

Full narrative and real-data diagnosis: see the parent project's
[`md/20260905.progress.md`](https://github.com/terencewtli/asm_lr/blob/master/md/20260905.progress.md)
— this repo is the methods follow-up to that investigation.

## Why not just use DSS?

DSS is a mature, peer-reviewed beta-binomial tool for calling differentially
methylated regions — but it's built for comparing **groups of biological
replicates** (e.g. 3 tumor samples vs. 3 normal samples), and estimates its
overdispersion parameter from **between-replicate variance**. A single
individual's two haplotypes are not biological replicates — HP1 and HP2 are
the two halves of *one* sample, not independent samples of a population.
With exactly one "replicate" per group everywhere in the genome, that
variance is undefined, not just noisy. On this project's real data, DSS's
dispersion-estimation step hung indefinitely (no progress after 25+
minutes on a single chromosome) — plausibly hitting a degenerate case the
package was never built or tested for.

## The model

At each CpG *l*, haplotype *h* ∈ {1, 2}:

```
X_{l,h} | N_{l,h}, p_{l,h}  ~  Binomial(N_{l,h}, p_{l,h})
p_{l,h}                     ~  Beta(mean = μ_{l,h}, dispersion = ρ)
```

i.e. a beta-binomial with `Var(X) = Nμ(1-μ)[1+(N-1)ρ]` — variance grows with
depth beyond pure binomial sampling, capturing real per-site noise.

**The key design choice**: ρ is a single value, estimated **once, genome-wide,
by pooling across loci** — not per-locus from replicate variance (there is
none). Millions of CpGs substitute for the "replicates" a normal DMR tool
would use (`ont_asm_caller/dispersion.py`), estimated via the excess of a
per-locus Pearson χ² statistic over its binomial expectation.

Per locus, the test is a likelihood-ratio test of `H0: μ_1 = μ_2` vs.
`H1: μ_1 ≠ μ_2`, with ρ fixed at its estimate — depth (N) enters the
likelihood directly, so the same procedure and threshold apply sensibly
regardless of a sample's coverage. (`ont_asm_caller/model.py`)

### Region-level pooling (`ont_asm_caller/region.py`)

The per-CpG test above, while correctly calibrated in isolation, turned out
to have **near-zero power once BH-FDR is applied at a realistic genome-wide
scale** (~0.5% true ASM prevalence, millions of tests — see Results below).
The fix: cluster CpGs within 500bp (matching this project's existing
`GAP_MERGE_BP` convention and DSS's own `smoothing.span`) and pool their
counts before testing — fewer tests (eases the BH penalty) and genuinely
lower per-test noise (averages out site-specific dispersion across the
pooled CpGs).

**Non-obvious pitfall, caught during development**: pooling several
independent per-CpG beta-binomial draws *lowers* the effective dispersion of
the combined count. Reusing the per-CpG ρ estimate for pooled region counts
is misspecified — it makes the region test needlessly conservative, silently
defeating the point of pooling. Dispersion must be **re-estimated on the
pooled region-level counts themselves** (confirmed directly by
`test_pooled_dispersion_is_lower_than_per_cpg_dispersion` in the test suite).

## Results (simulated)

Two-layer generative model (`ont_asm_caller/simulate_mixture.py`): each
locus/region independently gets a true ASM label (Bernoulli, π=0.5% —
matching the one published long-read validated rate, deCODE's Icelandic
cohort) and, if true, an effect size (uniform 0.2–0.6). CpGs are spatially
clustered (4-8 per region, within 400bp) so region pooling can be fairly
tested. Full benchmark: `benchmarks/compare_methods.py`.

| Depth | Method | n called | Empirical FDR | Power | Effect-size bias |
|---|---|---|---|---|---|
| 18x | Fisher (original pipeline) | 365 | **0.822** | 0.168 | +0.249 |
| 18x | Beta-binomial, per-CpG | 1 | 0.000 | 0.003 | +0.479 |
| 18x | **Beta-binomial, region** | **45** | **0.111** | **0.615** | **+0.016** |
| 8x | Fisher | 0 | — | 0.000 | — |
| 8x | Beta-binomial, per-CpG | 0 | — | 0.000 | — |
| 8x | **Beta-binomial, region** | **26** | **0.077** | **0.369** | **+0.044** |
| 5-40x (heterogeneous) | Fisher | 1256 | **0.898** | 0.274 | +0.190 |
| 5-40x | Beta-binomial, per-CpG | 7 | 0.143 | 0.013 | +0.346 |
| 5-40x | **Beta-binomial, region** | **56** | **0.107** | **0.649** | **+0.012** |

**Reading these numbers honestly:**

- **Fisher is not just underpowered, it's actively unsafe**: 77–90% of
  what it calls "significant" are false positives against ground truth,
  even in a clean simulation with a known generative model. This
  quantifies (doesn't just infer) the real-data finding that motivated this
  whole package.
- **The per-CpG beta-binomial test is correctly calibrated but practically
  useless on its own** — genome-wide multiple testing at a realistic (rare)
  true-ASM prevalence crushes its power to near zero, even though it behaves
  exactly as designed in isolation (see `tests/test_model.py`'s calibration
  and power checks, which use a flat p<0.05 threshold rather than a
  realistic BH-FDR burden).
- **Region-level pooling is a real fix, not yet a finished one.** Power goes
  from ~0–1% to 37–65% across depths, and FDR (8–11%) is categorically
  better than Fisher's — but it does not yet hit the nominal 5% target
  exactly. That gap is a measured fact, not a rounding error, and likely
  reflects the LRT's chi-square asymptotics not being exact at these sample
  sizes and/or imperfect dispersion estimation. Not yet something to quote a
  specific FDR number from in a paper without more calibration work.
- **Effect-size estimates are the cleanest win**: region-level bias
  (+0.01–0.04) is far better than Fisher's systematic overestimation
  (+0.19–0.25), which matters independently of the detection-power question
  if effect sizes themselves are reported anywhere downstream.

## Literature context

The core device — beta-binomial modeling with dispersion pooled across
*loci* rather than *replicates*, because a single individual has no
replicates — is well-established in **allele-specific expression** (ASE):
MBASED (Mayba et al. 2014), QuASAR (Harvey et al. 2015), and RASQUAL
(Kumasaka et al. 2016) all solve exactly this problem for RNA-seq allelic
counts. This package ports that machinery to methylation counts.

For long-read methylation ASM specifically, the closest published
precedent — the deCODE Icelandic ONT cohort (Nat Genet 2024) — does **not**
use a statistical test for its initial haplotype-methylation calling at
all: it's a bare rule ("≥3 CpG units each with rate <0.15, within 500bp").
Their rigor comes entirely from a downstream **population-level** genotype
regression across their 7,179-genome cohort. That substitution (cohort scale
for per-individual statistical power) isn't available at the 18-donor scale
this package was built for, which is the actual reason a per-individual
statistical test matters here in a way it didn't for them.

A related but distinct problem — detecting methylation-pattern
("epiallele") heterogeneity **without** requiring a phased het SNP nearby,
via clustering per-read methylation vectors (e.g. HDBSCAN, as in
`modbamtools`) — is established as a complementary approach for
SNP-sparse regions, but is a different modeling problem (unsupervised
mixture-model estimation, with the attendant non-regular LRT asymptotics for
testing the number of components) and is not implemented here.

## Status / what's not done yet

- FDR is not yet calibrated to the nominal target at the region level (see
  Results above) — needs either a better dispersion estimator, an exact
  small-sample test, or empirical calibration before quoting a specific FDR.
- Not yet validated on real data — everything above is simulation-only.
  Next step: run on the parent project's "good donors" (confirmed-complete
  downloads, no anomalous chr1 signal) and cross-check against known
  positive controls (Zink et al. 2018 imprinted DMRs; Rosenski et al. 2025
  allele-specific methylation atlas).
- No CLI yet — use the library directly (see `ont_asm_caller/__init__.py`
  for the public API: `test_locus`, `estimate_dispersion`, `cluster_cpgs`,
  `test_regions`, plus the `simulate*` functions for further validation
  work).
- Switch-error / haplotype-mislabeling modeling exists in the simulator
  (`simulate_loci(..., switch_error_rate=...)`) to characterize how
  imperfect phasing attenuates detectable effect sizes, but isn't yet
  incorporated into the test itself.

## Development

```
pip install -e .
python -m pytest tests/          # 33 tests: calibration, power, dispersion
                                  # recovery, clustering mechanics, region
                                  # pooling's power/dispersion claims
python benchmarks/compare_methods.py   # reproduces the Results table above
```
