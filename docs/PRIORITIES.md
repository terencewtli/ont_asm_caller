# Priorities

**Living document. Updated 2026-09-07** after the first real read-level
measurement on this project's data (HG00146 chr15). Supersedes the ordering in
`2026-09-07_next_steps_for_discovery_handoff.md`, which was written before that
measurement existed.

---

## UPDATE 2026-09-15: measured on 11 donors across all of chr15, and it is smaller still

`P09_design_effect_cohort.py`, every 10th 500 bp window across the **whole** of chr15 (~16k
windows per donor, identical loci in every donor), same call logic as below
(`benchmarks/results/2026-09-15_chr15_design_effect_12donors.txt`):

| class | donors | median DE, per donor | % windows DE > 5 |
|---|---|---|---|
| sane | 6 | 1.09–1.15 | 0.5–0.8% |
| flagged (high ASM rate) | 3 | 1.03–1.18 | 0.8–0.9% |
| zero-hit (ρ̂ 5–10× cohort) | 2 | 1.11–1.17 | 0.7–0.8% |

NA18508 is excluded: its local extract is truncated (`gzip: unexpected end of file`).

- **Pooling is defensible at the caller's scale in every donor measured**, not just HG00146.
- **The design effect does not distinguish sane, flagged or zero-hit donors.** It cannot
  explain either anomaly, so the switch-error / haplotype-label hypothesis stands, and the
  dipcall re-call is its direct test. This closes the handoff note's item 1 for chr15, and item 6
  step 1 below. Parent project `md/20260907_outstanding_items.md` #10 (a shared root cause for DE
  and the per-donor bulk-%mCG gradient) is also unsupported: DE is flat across donors.
- **The 1.29 below was an overestimate from an unrepresentative stretch.** HG00146's
  whole-chromosome median is **1.09**, against 1.22 on its first 12M rows (the pericentromeric /
  15q11 end). Call error is uniform too: ε = 0.040–0.047 in every donor.
- Not measured: haplotype-split DE, other chromosomes, or the caller's own region boundaries rather
  than fixed windows.

The section below is kept as the record of the first measurement.

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

### 3. ~~Resolve the `implied_design_effect` vs `design_effect` disagreement~~ — CLOSED 2026-09-15

**Not a summary mismatch; a wrong input.** `P10_estimator_reconciliation.py` walks from
one estimator to the other one assumption at a time on the same 9,049 HG00146 windows:

| step | DE |
|---|---|
| median `design_effect`, read halves (the reported number) | 1.22 |
| mean instead of median (right skew) | 1.44 |
| pooled variance-weighted icc, de-attenuated for call error | 1.64 |
| `implied_design_effect` as called by P08 (`n_cpgs` = 13, gap 21 bp) | **2.81** / 3.17 de-att. |
| `implied_design_effect` with C = positions seen by ≥50% of reads (5) | **1.56** |
| pooled curve read at each read's own CpG pairs | 1.52 raw / 1.62 de-att. |

P08 handed it `n_cpgs` = the **union** of positions in a window. Median per-position read
support is **0.02**; 63% of positions are seen by <20% of a window's reads, and the median
read observes **4.3** CpGs. Given the same inputs the two estimators agree within ~0.1. The
remaining 1.22 vs ~1.5 is median vs mean of a right-skewed distribution.

Consequences: the "median 13 CpGs per window" in the headline table counts sparse positions,
not CpGs each molecule carries; use ~5. `design_effect` itself is robust to this (per-read
fractions are unaffected). And **the steep short-range decay is real**: restricting to
well-supported positions leaves the curve unchanged (0.52 → 0.53 at 10 bp), so item 4 is not
chasing a sparse-position artefact.

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

### 8. Ancestry-driven heterozygosity density will confound any cross-population detection-rate comparison

Flagged 2026-09-16 from the sibling `asm_lr_hprc2` project — see `JOURNAL.md`'s
2026-09-16 entry for the full reasoning. Blocked on the same haplotagged-BAM
gap as item 1 (can't measure this against real ASM calls until those exist),
but the shape of the fix affects how the caller should be built, so it belongs
here now rather than after the fact.

**The problem in one line:** African-ancestry genomes carry more heterozygous
sites genome-wide (larger historical effective population size — a standard
population-genetics fact, confirmed in `asm_lr_hprc2`'s own 221-donor panel via
1000G-panel AF divergence and LD-decay measurements this session). Any test
that needs a phasing-informative het site to even attempt a call will therefore
have systematically higher yield in higher-heterozygosity populations, for
reasons with nothing to do with real ASM biology. Compare raw detection counts
across ancestries without accounting for this and you're measuring
heterozygosity density, not biology.

**Do not fix it by loosening the het-site requirement per ancestry** — that
trades a quantifiable confound for an ancestry-conditional analytical choice,
which is worse (see JOURNAL entry for the ascertainment-bias-array analogy).

**Fix direction:** extend this package's existing measured-nuisance-parameter
posture (DE, dispersion, call-error rate are all measured quantities that
propagate into calibration, not assumed constants — see item 3) to per-locus
het-site informativeness. A locus with few phasing-informative reads should
come out *underpowered* (wide interval, explicit low-confidence flag), not
silently folded into a binary "not detected." That's a caller-level output
change (extend `region_statistics`/`test_regions` to report and propagate
informative-read count as an explicit power indicator), not a filtering change
upstream of the caller.

**Measured, same day (see `JOURNAL.md`):** AFR 1470.4 het sites/Mb vs. EAS
1073.4 (~37% higher), ANOVA R²=0.917 across 5 superpopulations, p=1.6×10⁻¹⁰²
— superpopulation alone accounts for essentially all of the variance in
het-site density. The confound is confirmed, not just predicted. **Not yet
done:** the actual detection-rate effect on real ASM calls (blocked on item
1's haplotagged-BAM gap) and deciding the exact form of the power/confidence
output. No code changed in this repo for this item yet.

---

## Resolved / demoted

- **CPEL evaluation** — closed, see `2026-09-07_cpel_review_and_pattern_tests.md`.
  Not adopted; statistics ported, machinery and null not.
- **Permutation p-value floor** — fixed (`null.py`).
- **`decay_bp` unmeasured** — measured (this document). It was priority 1 for two
  days; it is now an input.
