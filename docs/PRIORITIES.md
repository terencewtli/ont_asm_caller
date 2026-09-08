# Priorities

**Living document. Updated 2026-09-07** after the first real read-level
measurement on this project's data (HG00146 chr15). Supersedes the ordering in
`2026-09-07_next_steps_for_discovery_handoff.md`, which was written before that
measurement existed.

---

## The headline: the design effect has now been measured, and it is small

`readlevel.design_effect` on **HG00146 chr15**, 12M `modkit extract` rows
(~2.5 Mb, 9,049 windows, median 68 reads and 13 CpGs per window), with
strand normalisation and confident calling applied:

| window | median DE | IQR | 90th pct | % DE>2 | % DE>5 |
|---|---|---|---|---|---|
| **500 bp** (the caller's own region scale) | **1.29** | 1.05–1.72 | 2.20 | 15% | **1%** |
| 1000 bp (`max_span`) | 1.48 | 1.17–2.06 | 2.78 | 27% | 1% |

Against `benchmarks/compare_methods_v2.py`, which measured the pooled
beta-binomial's FDR as a function of DE — 0.040 at DE 1.01, 0.070 at DE 2.40,
**0.695** at DE 5.49 — this says:

> **Region pooling is broadly defensible at the scale this caller actually uses.
> The regime where it collapses is ~1% of regions, not the typical case.**

That is close to what `P05` reported on chr1 (HG00146 median 1.39), so that
earlier number survives its call-logic bug better than expected.

**This demotes a large part of the work built in response to the calibration
critique.** The critique's headline figure — type I error 2.9× at α=0.05 and 11×
at α=0.001 — was measured at co-methylation 1.0, i.e. DE 5.5, which is the 99th
percentile of real regions rather than the typical one. The read-level path, the
pattern path and the matched null are all still correct and still worth having,
but the case for making any of them the *default* is much weaker than it looked
yesterday. The honest architecture is **DE-aware**: pool where the measured DE is
low, fall back to the read-level path where it is not. DE is now measurable per
region, so this is a switch, not a research programme.

### What the correlation structure actually looks like

Pooled correlation vs. CpG–CpG distance, HG00146 chr15 (call-error corrected):

| bp | 3 | 6 | 11 | 19 | 32 | 49 | 73 | 112 | **167** | 245 | 561 | 1225 | 1775 |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| corr | 0.68 | **0.76** | 0.67 | 0.35 | 0.24 | 0.15 | 0.12 | 0.11 | **0.19** | 0.12 | 0.10 | 0.08 | 0.09 |

- **A single exponential does not fit** (r² = 0.25–0.40, against >0.95 on
  simulation). The structure is two-component: a steep decay over ~20–25 bp plus
  a flat long-range floor near 0.07–0.09 that does not decay within 2 kb.
  `simulate_patterns.molecules_markov` assumes a single exponential and is
  therefore misspecified against real data.
- The bump at **167 bp** sits at the nucleosome repeat length. One donor, one
  bin — suggestive, not a finding.
- Call error matters: ε ≈ 0.043 even after filtering at `CONF_THRESH = 0.80`,
  attenuating every correlation by (1−2ε)² = 0.836.

---

## Priorities, in order

### 1. Get the haplotagged BAM and the region tables for chr15 — everything ASM-specific is blocked without them

The extract files alone were enough for the design-effect and decay
measurements, because co-methylation along a molecule **is not a haplotype
quantity** (`P08_decay_from_extract.py` exists to exploit exactly that). Nothing
else works without HP tags:

| available now | blocked on HP tags |
|---|---|
| `decay_bp`, design effect, call-error rate | every ASM test |
| co-methylation structure, stratified | `MatchedNull`, `T_MML/T_NME/T_PDM` |
| the DE-vs-region-scale curve | imprinted-DMR positive controls |

Needed: `bam/haplotagged/by_chrom/{sample}_chr15_haplotagged.bam` and
`tables/dasm_betabinom/{sample}/regions_chr15.tsv`. **This is the single
highest-value transfer.** One chromosome for 2–3 donors is enough; see
`2026-09-07_cpel_review_and_pattern_tests.md` §7.

### 2. Rerun `P05` through the corrected call logic

`P05_design_effect_real_data.py` still uses `mod_qual >= 0.5` on the `m` row
alone, which scores a `p_h = 0.686 / p_m = 0.314` position as confidently
*unmethylated* and ignores strand. Its published numbers (HG00146 1.39,
NA18508 1.13) are attenuated lower bounds. Cheap to fix — `P08` already has the
corrected logic.

### 3. Resolve the `implied_design_effect` vs `design_effect` disagreement

At 500 bp they differ ~2× (2.77 vs 1.29); at 2 kb they agree (0.0733 vs 0.0728).
The likely cause is a summary mismatch — the curve is pooled with variance
weighting, the DE is a median over windows — but it is not proven. **These two
estimators disagreeing is the tripwire that says the decay measurement is not
trustworthy**, so it should not be left standing. Until then, use the direct
per-region `design_effect` for calibration decisions.

### 4. Replace the single-exponential decay model

Both `decay.comethylation_decay`'s fit and `simulate_patterns.molecules_markov`
assume one exponential; real data has two components. The simulator being
misspecified against real correlation structure is the same class of error as the
one that made `simulate_reads.py` unusable for the CPEL comparison. A
two-component fit (short exponential + floor) is a small change.

### 5. Decide whether the entropy axis stays

`benchmarks/compare_pattern_methods.py`: maximum power on equal-mean disorder
imbalance was **0.034** across every method and scenario; `t_nme` scored 0.000 on
that axis in all three. On present evidence it is real, invisible to every
existing test, and **not detectable** at 10–30 reads per haplotype. The one use
that survives is diagnostic rather than discovery — see item 6.

### 6. The two unexplained-overdispersion donors (HG00344, NA21144)

Documented in `scripts/phasing_qc/README.md` and handoff item 3; ρ̂ = 0.097/0.116,
5–10× the cohort, ~zero significant regions. **Not a bug the method must fix — a
property of the inputs that the new paths make decidable**, because they have no
ρ parameter at all. HG00344's chr15 extract is already downloading, so this is
close at hand:

1. `design_effect` on HG00344 vs HG00146 (done above, 1.29). If HG00344's DE is
   genuinely 5–10× higher, the elevated ρ̂ is **correct** and the caller is
   behaving properly on noisier data.
2. `test_region_reads` + `MatchedNull` on it. Normal ⇒ dispersion-estimation
   artifact. Still anomalous ⇒ the haplotype labels are the problem.
3. Switch errors make a mislabelled allele look **bistable**, so they should
   inflate `T_NME` while attenuating `T_MML`. That signature on exactly these two
   donors would be near-positive identification — and it is a genome-wide
   aggregate, so it does not need the per-region power the entropy axis lacks.

### 7. The DE sweep the pattern benchmark is missing

All three benchmark scenarios landed at DE 2.1–3.5, so the comparison never
probed the regime that motivated the read-level path. Less urgent now that real
DE is known to be ~1.3, but it is why "pooled_bb wins everywhere" in that table
should not be quoted.

---

## Resolved / demoted

- **CPEL evaluation** — closed, see `2026-09-07_cpel_review_and_pattern_tests.md`.
  Not adopted; statistics ported, machinery and null not.
- **Permutation p-value floor** — fixed (`null.py`).
- **`decay_bp` unmeasured** — measured (this document). It was priority 1 for two
  days; it is now an input.
