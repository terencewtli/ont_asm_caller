# CPEL in the long-read setting: review, applicability, and a test design

**2026-09-07.** Response to handoff item 4 ("Evaluate CPEL as a possible
off-the-shelf alternative rather than continuing to patch a custom pooled-count
model"). Covers what the paper actually does, what transfers to ONT data and
what does not, whether the existing simulator can adjudicate it (it cannot, and
why), and what was built instead.

**Citation correction first.** The README and the handoff note both cite this as
"Jiang et al." That is wrong and should be fixed before it propagates into a
manuscript. The paper is:

> Abante J, Fang Y, Feinberg AP, Goutsias J. *Detection of haplotype-dependent
> allele-specific DNA methylation in WGBS data.* Nat Commun 11:5238 (2020).
> Software: `CpelAsm.jl` (Julia), github.com/jordiabante/CpelAsm.jl

---

## 1. What CPEL actually does

Per haplotype, over N ≤ 20 CpG sites, it fits a 1-D Ising model

```
p(x) ∝ exp( Σ_k α_k Σ_{n∈R_k} s_n  +  β Σ_n s_n s_{n+1} ),      s_n = 2x_n − 1
```

with one field `α_k` per ≤500 bp subregion and **one** nearest-neighbour
coupling `β` for the whole region — so K+1 free parameters, K usually 1–2. From
the two fitted allele distributions it forms three statistics:

| statistic | definition | what it detects |
|---|---|---|
| `T_MML` | \|MML₁ − MML₂\| | mean methylation imbalance |
| `T_NME` | \|NME₁ − NME₂\| | **entropy / disorder imbalance** |
| `T_PDM` | JSD(p₁, p₂) / N | any distributional difference |

`T_NME` is the one that has no counterpart anywhere in this package.

**Significance** comes from an empirical bootstrap, not from the model
(Supplementary §10). Real reads from *non-haplotype* (homozygous) regions with
≥ N CpGs are randomly split into two equal-sized pseudo-alleles, the model is
refitted, and the statistic recorded; L ≥ 1000 such draws per N give
`p̂ = (1 + #{t_l ≥ t*}) / (L+1)`, then BH at q ≤ 0.05.

Two details of that procedure matter more than anything else in the paper:

- **The null is stratified on N alone.** Not on coverage, not on methylation
  level.
- **Every null draw is generated at the global minimum coverage `Cmin = 5`**
  (their step 5 downsamples both pseudo-alleles to 5–7 reads), *regardless of
  the region's real depth*. They are explicit that this is deliberately
  conservative.

### What is right about it

The core insight is correct and it is exactly the insight this package arrived
at independently from the other direction. CPEL's critique of marginal per-CpG
testing — that it ignores correlation between neighbouring CpGs and thereby
loses both sensitivity and specificity — is the same defect as this project's
original Fisher pipeline and the same defect as
`docs/2026-09-05_calibration_critique.md` Defect 1. Two independent literatures
converging on one diagnosis is good evidence the diagnosis is right.

Where CPEL goes further than this package: it models the correlation instead of
absorbing it. `readlevel.py` reduces each molecule to a fraction, which handles
the design effect correctly but throws away the within-read pattern. CPEL keeps
the pattern and turns it into a second, orthogonal ASM axis. **That axis is a
real gap in this package, and it is the one thing worth taking.**

### What is wrong with it, on its own terms

1. **Its p-values are floored at ~1/1001.** Across 715,155 haplotypes under BH
   at 0.05, nothing can be called unless roughly 14,000 haplotypes tie at the
   floor simultaneously. This is the *same* defect this repo already documented
   in `readlevel.test_region_perm` ("p is bounded below by 1/(n_perm+1), so this
   is NOT suitable on its own for genome-wide BH"). It means CPEL's discovery
   set is unrankable internally: the strongest and the weakest survivor carry
   identical p-values.
2. **The `Cmin = 5` null is not conservative in the way they claim once depth
   varies.** It is conservative *on average*, but it is not calibrated — it
   makes the null width independent of the region's coverage, so the effective
   threshold drifts with depth. That is the depth-confounding this whole project
   exists because of, arriving through the back door of the null instead of
   through the test statistic.
3. **No conditioning on methylation level.** Every one of these statistics has a
   null spread that collapses as μ → 0 or 1. Not conditioning on μ is the same
   error as `dispersion.estimate_dispersion`'s boundary failure (Defect 2), and
   it bites hardest at exactly the fully-methylated/fully-unmethylated regions
   where imprinting-style ASM lives.
4. **A single β for the whole region.** Defensible over a ≤500 bp WGBS
   haplotype. Not defensible over a long-read region where CpG spacing varies by
   an order of magnitude.
5. **48 h on 20 CPUs per sample**, by simulated annealing.

---

## 2. Applicability to long reads

The short version: **CPEL's statistics transfer; CPEL's machinery mostly does
not, because the constraint it was built around does not exist in ONT data.**

The Ising model is in that paper for one reason: a ~100 bp bisulfite read
observes a *fragment* of a haplotype, so `p(x)` over N ≤ 20 CpGs is never
observed and must be inferred by marginalising a parametric model over
everything each read failed to see. That is what costs 48 CPU-hours, and it is
also why their non-parametric comparator (NPD, Onuchic et al.) had to be
restricted to 4-CpG epialleles.

A 35 kb ONT read observes **every CpG in the region on a single molecule**. Five
consequences:

| | WGBS (CPEL's setting) | ONT long reads |
|---|---|---|
| observation | fragmentary; joint pattern inferred | joint pattern observed directly |
| fitting | marginalise unobserved sites; simulated annealing | complete-data concave MLE, milliseconds |
| non-parametric option | infeasible beyond 4 CpGs | **feasible** — compare empirical read distributions |
| region definition | SNP clusters ± read length | whole chromosome arms are phased; regions are regulatory units |
| depth per allele | ~5× | 20–40× |

Two of these are not incremental:

- **Region scope.** CPEL can only test where a SNP cluster is. Long reads phase
  chromosome arms, so every region is testable — including SNP-poor regions,
  which is where a large fraction of imprinted DMRs sit. This is a capability
  CPEL structurally cannot reach, and it is the strongest argument that porting
  the statistics rather than adopting the tool is the right call.
- **The null must be re-derived.** `Cmin = 5` against a real 20–40× region is
  not "conservative", it is a several-fold inflation of the null width, and it
  discards most of the power long reads were bought for.

**Recommendation: do not adopt `CpelAsm.jl`.** It takes SNPsplit-partitioned
WGBS BAMs and an N-masked reference, is stratified around a partial-observation
problem this data does not have, and would need a shim, a Julia dependency, and
a null replacement anyway. Port the two statistics that matter (`NME`, `PDM`),
fix the null, and cite the paper as the precedent it is. That is what
`ont_asm_caller/pattern.py` and `ont_asm_caller/null.py` now do.

---

## 3. Can the existing simulator adjudicate this? No — and why not is the point

`simulate_reads.simulate_read_regions` was the right fix for the defect it was
aimed at. It is **not** a valid testbed for a pattern-based method, for two
measured reasons.

### 3.1 Its correlation is exchangeable, not distance-decaying

`_molecules` copies each CpG's state from a single molecule-level latent with
probability q, so `corr(x_i, x_j) = q²` for **every** pair — no distance
dependence at all. Measured (n = 200,000 molecules, 12 CpGs, q = 0.6):

```
simulate_reads._molecules   lag 1..11:  0.361 0.360 0.360 0.359 0.359 0.357 ...
1-D Ising / Markov chain    lag 1..11:  0.600 0.364 0.218 0.131 0.076 0.043 ...
```

Two things follow. First, the parent project's co-methylation *decay* curve
(`md/20260903_qc_review.md` §7.3) has nothing to measure in this simulator —
there is no decay. Second, and worse for this comparison: **under
exchangeability the per-read methylation fraction is a sufficient statistic for
the read.** A pattern-based method is guaranteed to find nothing that
`readlevel.py` does not already find. Benchmarking CPEL on it would produce a
clean win for this package that means nothing.

The mirror-image error is just as easy: simulating from an Ising chain is
simulating from CPEL's own model, and would hand CPEL an equally meaningless
win. **Any benchmark run on one generative model is uninterpretable.** Three are
needed, and the headline one must be a model neither method assumes.

### 3.2 It has only one ASM axis

Both haplotypes get the same `co_methylation`; they differ only in mean. So the
entropy imbalance `T_NME` exists to detect has, by construction, **zero true
positives**. A benchmark on this data can only ask whether CPEL matches a mean
test on mean differences — precisely the axis where it claims no advantage.

### 3.3 What was built instead

`ont_asm_caller/simulate_patterns.py`. Two orthogonal ASM axes with independent
prevalence:

- **MML axis** — haplotypes differ in mean. What every existing test targets.
- **ICC axis** — haplotypes differ in within-read correlation **at exactly equal
  mean**: one allele carries a fixed pattern, the other is bistable between
  fully-methylated and fully-unmethylated molecules. Biologically this is
  polymorphic / stochastic ASM. `tests/test_simulate_patterns.py` asserts these
  regions carry no mean signal, so "blind by construction" is a checked claim
  and not an assertion.

and three generative models, run as a two-sided sensitivity check:

| model | correlation structure | favours |
|---|---|---|
| `epiallele` | mixture of discrete methylation patterns | **neither** — headline scenario |
| `markov` | distance-decaying, = 1-D Ising reparameterised | CPEL (correctly specified) |
| `exchangeable` | the `simulate_reads.py` model | this package (fraction sufficient) |

**Is pure simulation enough?** For the *method-selection* question — does the
pattern axis add anything, and does the cheap non-parametric version capture
it — yes, provided all three models are run and the conclusion is stable across
them. For any *prevalence* or *FDR* number that would appear in a paper, no.
`decay_bp` (the correlation half-life along a molecule) is the parameter that
sets every result here and it is currently a guess. See §6.

---

## 4. What was implemented

### `ont_asm_caller/pattern.py`
CPEL's model and statistics, in the form long reads make available.

- `fit_chain` — MLE of (α₁..α_K, β) by transfer matrix + L-BFGS with an
  **analytic gradient** (verified against finite differences to 1e-9 in
  `tests/test_pattern.py`). Complete-data where reads are complete, with
  per-read state clamping for no-calls, batched across reads. ~5 ms per
  haplotype-region, against CPEL's simulated annealing.
- `chain_entropy_bits` — exact NME via H(s₁) + Σ H(s_{n+1}|s_n) from
  forward-backward marginals, no enumeration of the 2^N state space.
- `jsd_chains` — `T_PDM` by Monte Carlo sampling from both fitted chains.
- `region_statistics` — `T_MML`, `T_NME`, `T_PDM` for one region.
- `read_fraction_ks` — **the cheap alternative**: two-sample KS on per-read
  methylation fractions. One line, no model. `readlevel.test_region_reads`
  compares *locations* (Welch / Mann-Whitney) and is therefore blind to
  equal-mean disorder imbalance; KS compares the whole fraction distribution and
  is not. If KS recovers most of what `T_PDM` recovers at ONT depth, the Ising
  machinery is not earning its keep here — which is the comparison the benchmark
  is built to make.

Every recursion is checked against brute-force enumeration of 2^N to machine
precision (`tests/test_pattern.py`).

A ridge penalty is required, not optional: on a fully-methylated region — the
single most common configuration in a real bimodal methylome — the unpenalised
MLE diverges to α = ∞.

### `ont_asm_caller/null.py`
The depth-matched null with tail extrapolation, which is the part that is not a
port.

- `split_reads` — permutes **whole molecules** between haplotypes, so
  within-read co-methylation is preserved exactly under H₀. This is why it is
  valid where `scripts/phasing_qc/P04_permutation_null.py`'s count-level
  bootstrap is not.
- `MatchedNull` — pools null draws **across regions within a stratum**, buying
  L ~ 10⁵ instead of 10³ for the same compute. Strata are
  **(n_cpgs, depth, methylation level)** — CPEL uses n_cpgs alone.
- Tail extrapolation — a generalised Pareto fit to exceedances over a high
  threshold, so p-values continue below the empirical floor. This is the
  standard peaks-over-threshold treatment of permutation p-values (Knijnenburg
  et al., *Bioinformatics* 2009), not an invention here.
- `cpel_cmin=5` — reproduces CPEL's own null protocol, so the benchmark can
  separate "the statistic is weak" from "the null is mismatched".
- `calibration_report` — the check that decides whether any of the above is
  honest: are p-values uniform on H₀ splits?

**A defect this found in its own first run, worth recording because of how it
hid.** The initial version trusted the GPD tail without limit and returned
p-values of `1e-300` from a 60,000-draw null — including on H₀ data, i.e. the
extreme tail was broken. Type I error at α = 0.05 (0.049) and α = 0.01 (0.010)
looked *correct*, so the standard calibration summary showed nothing wrong; only
the α = 10⁻³ row (0.0016, 1.6× inflated) and `min_p` gave it away. **The usual
calibration table is not evidence about the regime BH actually operates in.**

Cause: every statistic here is bounded in [0, 1], so a fitted GPD with a large
positive shape — a polynomial tail with no upper endpoint — is structurally
misspecified, and thin strata produced exactly that. Fixed by rejecting
pathological shapes, requiring 200 exceedances rather than 50, and capping
extrapolation at 100× the empirical floor. Post-fix type I error is
0.047 / 0.0088 / 0.00075 / 0.000 at α = 0.05 / 0.01 / 10⁻³ / 10⁻⁴. The cap still
leaves p ≈ 10⁻⁷ reachable, which is what BH across 10⁵ regions needs.
`calibration_report` now reports the 10⁻⁴ row by default so this cannot hide
again.

**This fix is worth more than the CPEL port.** It applies to
`readlevel.test_region_perm` unchanged, and that test is currently the only
correctly-calibrated test in the package and is unusable genome-wide for exactly
the floor reason. Fixing the floor promotes it from "confirmation step for
screened candidates" to a genome-wide test.

---

## 5. Running it on real data

Ordered, and the order matters — step 2 gates everything after it.

**Pick the donors deliberately — and note a contradiction in the repo's own
docs.** `VALIDATION_PLAN.md` §"Which chromosome and which donors" recommends
`HG00344` as one of the two first-pass donors ("complete downloads, not on the
problem-donor list"), while `scripts/phasing_qc/README.md` flags `HG00344` as one
of the **two anomalous zero-hit donors** (ρ̂ = 0.097, 5–10× the cohort, ~zero
significant regions). Both cannot be right; resolve before extracting.

The resolution is convenient rather than costly: extract **`HG00344` (anomalous)
plus one confirmed-clean donor**. That gives the decay/DE baseline *and* the
diagnostic in §8.6 in one pass.

```bash
# 1. extract once; ~tens of MB instead of the 4-5 GB modkit extract stream
python3 scripts/phasing_qc/P06_export_read_matrices.py NA19682 chr15   # clean baseline
python3 scripts/phasing_qc/P06_export_read_matrices.py HG00344 chr15   # anomalous donor

# 2. MEASURE decay_bp. Nothing below should be quoted before this exists.
#    Prints the fitted decay beside readlevel.design_effect on the same
#    regions -- two independent estimators that must agree.
python3 scripts/phasing_qc/P07_comethylation_decay.py NA19682 chr15
python3 scripts/phasing_qc/P07_comethylation_decay.py HG00344 chr15
```

```python
# 3. then, in the library, on the same loaded regions:
from ont_asm_caller import MatchedNull, region_statistics, read_fraction_ks
from scripts.phasing_qc.P06_export_read_matrices import load_read_matrices

regions = list(load_read_matrices("tables/read_matrices/NA19682_chr15.npz", "chr15"))

def stats(m1, m2):
    out = {}
    ks = read_fraction_ks(m1, m2)
    if ks: out["ks"] = ks[0]
    out.update({k: v for k, v in region_statistics(m1, m2).items()
                if k.startswith("t_")})
    return out or None

null = MatchedNull(stats, n_perm_per_region=20).build(regions, seed=0)
# ALWAYS run this before trusting a q-value from the above:
from ont_asm_caller import calibration_report
print(calibration_report(null, regions, stats, ["ks", "t_mml", "t_nme", "t_pdm"]))
```

`load_read_matrices` yields objects duck-typed to
`simulate_patterns.PatternRegion`, so there is no separate real-data code path to
keep in sync with the simulation one.

Cost, measured single-core: `region_statistics` ~10 ms per region, so ~17 min per
10⁵ regions for point estimates and ~6 h for a 20-permutation null. Both
parallelise trivially by region.

---

## 6. Benchmark

12,000 regions per scenario, π = 0.01 on each axis independently, BH α = 0.05,
10–30 reads per haplotype. Raw logs: `benchmarks/results/2026-09-07_pattern_methods.txt`.

| scenario | method | called | FDR | power | pow(mean) | pow(corr) |
|---|---|---|---|---|---|---|
| **epiallele** (DE 2.3) | pooled_bb | 21 | 0.143 | **0.073** | 0.151 | 0.000 |
| *neither model* | readfrac_welch | 4 | 0.000 | 0.016 | 0.034 | 0.000 |
| | ks_matched | 12 | 0.000 | 0.049 | 0.101 | 0.000 |
| | t_mml | 13 | 0.154 | 0.045 | 0.092 | 0.000 |
| | t_nme | 15 | 0.133 | 0.053 | 0.109 | 0.000 |
| | t_pdm | 14 | 0.357 | 0.037 | 0.050 | 0.024 |
| | **t_pdm_cpelnull** | **0** | — | **0.000** | 0.000 | 0.000 |
| **markov** (DE 2.1) | pooled_bb | 73 | 0.041 | **0.295** | 0.590 | 0.000 |
| *CPEL correctly* | readfrac_welch | 22 | 0.000 | 0.093 | 0.188 | 0.000 |
| *specified* | ks_matched | 21 | 0.048 | 0.084 | 0.162 | 0.008 |
| | t_mml | 9 | 0.000 | 0.038 | 0.077 | 0.000 |
| | t_nme | 26 | 0.077 | 0.101 | 0.205 | 0.000 |
| | t_pdm | 32 | 0.031 | 0.131 | 0.231 | 0.034 |
| | **t_pdm_cpelnull** | **0** | — | **0.000** | 0.000 | 0.000 |
| **exchangeable** (DE 3.5) | pooled_bb | 36 | 0.083 | **0.144** | 0.266 | 0.000 |
| *this package* | readfrac_welch | 2 | 0.000 | 0.009 | 0.016 | 0.000 |
| *favoured* | ks_matched | 9 | 0.111 | 0.035 | 0.065 | 0.000 |
| | t_mml | 4 | 0.000 | 0.017 | 0.032 | 0.000 |
| | t_nme | 3 | 0.000 | 0.013 | 0.024 | 0.000 |
| | t_pdm | 5 | 0.000 | 0.022 | 0.024 | 0.019 |
| | **t_pdm_cpelnull** | **0** | — | **0.000** | 0.000 | 0.000 |

Matched-null type I error, all three scenarios pooled: 0.044–0.059 at α = 0.05,
0.007–0.015 at 0.01, 0.0008–0.0024 at 0.001. The null is calibrated.

### Reading these honestly

**1. CPEL's statistics do not win anywhere — including where CPEL is correctly
specified.** On `markov`, generated from an Ising chain, the best CPEL statistic
(`t_pdm`, power 0.131) is beaten by the plain pooled beta-binomial already in
this package (0.295), and that is with CPEL's statistics given a *better* null
than CPEL's own. The added machinery does not pay for itself here.

**2. CPEL's null is the binding constraint, confirmed three times.**
`t_pdm_cpelnull` — the same statistic scored against CPEL's own protocol
(`Cmin = 5`, stratified on CpG count alone) — called **zero** regions in every
scenario, floored at p ≈ 1e-3, while the identical statistic against the
depth-matched null called 14 / 32 / 5. This is the single cleanest result here.

**3. The correlation axis is not reachable at ONT depth.** Maximum power on
equal-mean disorder imbalance across every method and every scenario: **0.034**.
`t_nme` — the statistic built for it — scored **0.000 on that axis in all three**
while retaining power on the *mean* axis (0.109 / 0.205 / 0.024), i.e. at 10–30
reads NME is largely tracking the mean rather than the entropy. That follows from
the identifiability measurement in §8.3. The second ASM axis is real and is
genuinely invisible to every existing test, and on this evidence it is **also not
detectable**, which is an argument against investing further in it rather than
for.

**4. The benchmark does not probe the regime that motivated the read-level path.**
Measured design effect in all three scenarios: **2.1 – 3.5**. `compare_methods_v2.py`
measured the pooled test's FDR against DE as 0.040 (DE 1.01) → 0.070 (DE 2.40) →
**0.695** (DE 5.49). These scenarios sit at the benign end. So "pooled_bb wins on
power everywhere" is a statement about DE ≈ 2–3.5 and says nothing about DE ≈ 5,
where the same test is already known to collapse. **The apparent winner is the
method whose failure mode this benchmark does not test**, and which regime is real
is decided by `decay_bp` (§7a). A DE sweep is the obvious missing run.

**5. The precision is poor and the method ranking is not resolved.** With 229–246
true positives and 2–73 calls per method, several FDR figures rest on fewer than
ten events; "FDR 0.000" on 4 calls means nothing. Differences among methods calling
under ~20 regions are not statistically resolvable at this scale. Treat the
qualitative findings (1)–(4) as the results and the individual numbers as
indicative.

**6. A limitation of the `epiallele` scenario itself.** It was built to be a model
neither method assumes, and parametrically it is — but its measured decay length
is infinite (amplitude 0.32), because its latent patterns span the whole region.
So it is non-parametric but still *exchangeable*, and only `markov` actually
carries distance decay. A genuinely neutral scenario with decaying correlation is
missing and should be added before this comparison is called settled.


---

## 7. What real data this needs, and how much

Three quantities need real data, with very different appetites. Conflating them
is how this ends up as a bigger data ask than it is.

### (a) `decay_bp` — the co-methylation decay length. **A few Mb is enough.**

This is the parameter every result in §6 is conditional on, and it is the same
quantity as the design effect (`readlevel.design_effect`) and as the parent
project's co-methylation-decay curve (`md/20260903_qc_review.md` §7.3). It is a
correlation-versus-distance curve pooled over millions of CpG pairs, so its
sampling error is negligible at almost any scale — 2,000 regions at 20×/haplotype
already gives third-decimal standard errors.

The reason to ask for a whole chromosome anyway is **stratification, not
volume**: decay almost certainly differs between CpG islands, shores and open
sea, and between high- and low-methylation regions, and a 5 Mb window will not
contain enough of each class. One chromosome will.

### (b) Null strata. **One chromosome is ample.**

`MatchedNull` needs regions covering the real (n_cpgs, depth, μ) space. One
chromosome gives ~10⁵ regions, i.e. thousands per stratum — far more than the
500-per-stratum minimum.

### (c) Prevalence of entropy-axis ASM — the actual scientific question.

This is the one that genuinely wants scale, because it is a rare-event rate that
has to clear genome-wide BH. **One chromosome is a legitimate pilot but not a
final answer**, and there is a specific problem with validating it:

> **chr15 validates the port but cannot validate the new axis.** The
> `VALIDATION_PLAN.md` argument for chr15 (98 of 229 Zink imprinted DMRs, vs. 16
> for the next-best chromosome) is right for `T_MML` — but imprinted DMRs are
> the canonical *mean* difference, one allele on and one off, and **both alleles
> at an imprinted DMR are LOW entropy**. They are the wrong positive control for
> `T_NME` and would make the entropy axis look useless.

There is no established truth set for allele-specific *entropy* imbalance. The
three usable substitutes, in order of strength:

1. **Cross-donor replication stratified by carrier status** — the same argument
   as `VALIDATION_PLAN.md` B2, applied to the entropy axis. Needs several donors
   on the same chromosome, not a whole genome on one.
2. **The parent project's own epiallele metrics.** `tables/W12_epiallele_metrics.csv`
   already holds `within_hp_entropy`, `hp_sep`, `within_hp_var` and
   `bimodality_coef` for **304,514 loci** on real data. That is the same axis,
   already measured descriptively — without a test, a null, or a haplotype
   contrast. **Reconcile with W12 before building anything further here.** Two
   arms of this project measuring the same quantity with different definitions
   is exactly the failure the parent `CLAUDE.md` guardrail was added for.
3. Known metastable epialleles / VMRs from the literature.

### Format matters more than volume

`modkit extract` on chr1 is 4–5 GB gzipped per sample and is **not
position-indexed**, which is why `P05_design_effect_real_data.py` streams the
whole file per window and why the handoff note flags that it "won't scale to
genome-wide". Asking for more of that format makes the problem worse.

The right ask is a **per-region read × CpG matrix dump**, produced once on the
cluster: for each region the caller already defines, one row per
(region, read_id, HP) with the packed binary call vector. At 10⁵ regions × ~40
reads × ~8 CpGs that is ~3 × 10⁷ calls — **tens of MB, not GB** — and it loads
straight into `pattern.py` and `null.py` with no reformatting, no HP-tag join
against the BAM, and no streaming.

### Bottom line

**One chromosome of read-level data for 2–3 donors, as a per-region read × CpG
matrix.** chr15 for the imprinted positive controls and the `T_MML` port check;
chr1 in addition if the depth-contrast and anomalous-donor questions
(`scripts/phasing_qc/README.md`) are to be touched at the same time. Anything
beyond that is not the binding constraint — the entropy-axis truth set is.

---

## 8. Open, in priority order

1. **`decay_bp` is a guess.** Every number in §5 is conditional on it. It is
   cheap to measure (§7a) and nothing here should be quoted before it is.
2. **Reconcile with `W12_epiallele_metrics.csv`** before extending the entropy
   axis further (§7c.2).
3. **Parameter identifiability at ONT depth.** Measured: at R = 400 reads the
   MLE recovers (α, β) = (1.2, 0.6) as (1.26, 0.58); at R = 20 it returns
   (1.81, 0.10) — the field and the coupling trade off badly and the split is
   poorly identified at realistic depth. CPEL states the same caveat. It is
   survivable because the null is built by refitting under the same conditions,
   so the bias is common to statistic and null — but it caps the power of
   `T_NME` and is a reason the non-parametric route may win.
4. **A fixed heterogeneous pattern is not representable** by K fields + one
   coupling. Measured: 25 reads all carrying the same fixed half-methylated
   pattern give NME = 0.64, not ~0. CPEL's NME therefore *overstates* disorder
   for ordered non-constant patterns. Another argument for the non-parametric
   route in long reads.
5. **Switch errors are still not modelled** in any test, here or in the
   pre-existing ones. Haplotype mislabelling moves reads between alleles, which
   attenuates `T_MML` and *inflates* `T_NME` — a mislabelled allele looks
   bistable. Given that the four flagged donors are suspected switch-error cases
   (`scripts/phasing_qc/README.md`), the entropy axis may be **more**
   switch-sensitive than the mean axis, not less. This should be simulated
   before the entropy axis is run on the flagged donors.

6. **The two unexplained-overdispersion donors are now a decidable question.**
   `HG00344` and `NA21144` carry ρ̂ = 0.097 / 0.116, 5–10× every other donor, and
   call essentially zero significant regions (`scripts/phasing_qc/README.md`;
   handoff item 3). Defect 2 predicts ρ too *low*, so neither documented defect
   explains it.

   **This is almost certainly a property of the inputs rather than a bug the
   method must fix — but the work in this commit turns it from a mystery into a
   test, and that is the reason to care.** A single global ρ is the only place a
   pooled-count model can put within-molecule correlation, so a genuinely higher
   design effect and haplotype labels that scramble reads across alleles produce
   the *same* symptom and the model cannot separate them. The read-level and
   pattern paths have no ρ parameter at all. Three checks, in order:

   1. Run `comethylation_decay` / `design_effect` on those two donors. If their
      DE really is 5–10× the cohort, the elevated ρ̂ is **correct** and the caller
      is behaving properly on genuinely noisier data.
   2. Run `test_region_reads` + `MatchedNull` on them. Normal results ⇒ it was a
      dispersion-estimation artifact. Still anomalous ⇒ the haplotype labels are
      the problem and no test will repair them.
   3. Sharpest: switch errors move reads between alleles, which makes a
      mislabelled allele look **bistable**. That should inflate `T_NME` while
      *attenuating* `T_MML`. Globally elevated `T_NME` with normal `T_MML` on
      exactly these two donors would be close to a positive identification of
      switch-error contamination — and it is the one use for the entropy axis
      that survives its poor power in §6, because it is a genome-wide aggregate
      signal, not a per-region call.
