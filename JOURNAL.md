# JOURNAL

A running record of what has been done in this repo, what each step established,
and what is still open. Newest entries at the top of the log.

Companion documents, and the division of labour between them:

| doc | answers |
|---|---|
| `README.md` | what the method is and why it exists |
| `docs/paper/method.pdf` | the models and their derivations |
| `docs/PRIORITIES.md` | **what to do next, in order** (living document) |
| `JOURNAL.md` (this file) | what has happened so far, and the open threads |

`PRIORITIES.md` is the authority on ordering. This file does not restate it; the
"Still to do" section below points into it.

---

## Cross-project flag — 2026-09-16: ancestry-driven heterozygosity density confounds any cross-population detection-rate comparison

Not caller work done in this repo — a design concern raised while working the sibling
`asm_lr_hprc2` project (same population-scale ASM effort, different codebase) that belongs on
record here because it's a property of what any long-read ASM caller has to consume as input,
not a property of that project's specific pipeline. Filed as `PRIORITIES.md` item 8 (see there
for the actionable version); this entry is the reasoning trail.

**The observation.** `asm_lr_hprc2`'s 2026-09-16 session measured real, ancestry-structured
population-genetics divergence in its 221-donor HPRC2 panel (1000G high-coverage panel subset,
chr1-3): among sites common (MAF≥5%) in the pooled panel, EAS donors show 14.6% become
monomorphic and 25.1% become rare within-population vs. AFR's 3.1%/15.4%; LD decay differs
sharply (AFR r²=0.40 at 0-10kb vs. EAS 0.78). Both are well-established population-genetics
facts (African populations' larger historical effective size, no severe out-of-Africa
bottleneck) — the contribution here is confirming they hold in this exact donor panel, not
discovering them.

**Why this is a caller problem, not just a data-prep problem.** Any long-read ASM approach
needs a read to span a het/phasing-informative site to be usable at all (in `asm_lr_hprc2`,
this is H02's `k>=1` het-site filter; this repo's own tests need haplotype tags for the same
reason). Because African-ancestry genomes carry systematically more heterozygous sites
genome-wide, the *fraction of reads that even become testable* will differ by ancestry for
purely mechanical, population-genetic reasons — independent of any real difference in ASM
biology. A population with lower baseline heterozygosity will show a lower apparent ASM
detection rate even when the true underlying regulation is identical, simply because fewer of
its reads clear the informativeness bar. Compare raw "detected loci" counts across ancestries
without accounting for this, and you're measuring heterozygosity density, not biology.

**Why the fix doesn't belong in the filter/threshold.** The tempting fix — loosen the
informativeness requirement for lower-heterozygosity populations — trades a quantifiable
confound for an unquantifiable, ancestry-conditional analytical choice, which is a worse
problem: now the pipeline treats populations differently before looking at the outcome, and
that's much harder to defend or reason about than a uniform filter with a known, measurable
power differential. Historical precedent for the wrong instinct here: genotyping-array
ascertainment bias against non-European variation was corrected downstream (imputation-quality
metrics, ancestry-matched panels), not by redesigning arrays per population — same shape of
problem.

**Where the fix belongs.** In this repo's terms: the design-effect framework already treats
per-locus read support (`n_cpgs`, informative-read count) as a measured quantity that
propagates into calibration (see `docs/PRIORITIES.md` item 3's `implied_design_effect`
reconciliation) rather than a hard pass/fail gate. The same posture should extend to
**het-site informativeness**: a locus with few phasing-informative reads should come out of the
caller as *underpowered* (wide interval, explicit low-confidence flag) rather than silently
folded into a binary "not detected." That converts an invisible ancestry confound into an
explicit, quantified power difference — the kind of thing precision-weighting or a
common-well-powered-subset restriction can correct for in a downstream cross-ancestry
comparison, the same way this package already treats DE and dispersion as measured nuisance
parameters rather than assumed constants.

**Status: flagged, not measured in this repo.** `asm_lr_hprc2` was mid-computation on actual
per-superpopulation het-site density (from its own H01 output) at the time of this entry — not
yet confirmed that African donors in that specific panel show higher het-site density, only
that population genetics strongly predicts they should and that the AF/LD divergence numbers
above are consistent with it. Real numbers, once they land, belong in `PRIORITIES.md` item 8,
not here.

---

## Where things stand — 2026-09-15

**Branch:** `cpel-evaluation-and-pattern-path`, **merged into `master` 2026-09-15** (fast-forward).

**Version:** 0.2.0. **Tests:** 73 passing (`python3 -m pytest`, ~85 s).

**Status of the science:** the package is *no longer* simulation-only — the
design effect, the co-methylation decay curve and the basecall error rate have
all been measured on real HG00146 chr15 data. But **every ASM test in the
package is still simulation-only**, because all of them need haplotype tags and
the haplotagged BAM has not been transferred. That split is the single most
important thing to know about the current state:

| measured on real data | still simulation-only |
|---|---|
| `readlevel.design_effect` | `model.test_locus`, `region.test_regions` |
| `decay.comethylation_decay` | `readlevel.test_region_reads` / `test_region_perm` |
| basecall error rate (ε ≈ 0.043) | `pattern.py` (MML / NME / JSD / KS) |
| region-scale DE curve | `null.MatchedNull` |

**Updated 2026-09-15:** across 11 donors and the whole of chr15 the median design effect is
**1.03–1.18** per donor (see the log entry). **The original headline result of this branch:** the
measured design effect at the caller's own 500 bp region scale is **median 1.29** (IQR 1.05–1.72, 1% of
regions above 5). The calibration critique that motivated the read-level and
pattern paths measured its 2.9×/11× type-I inflation at DE ≈ 5.5 — the 99th
percentile of real regions, not the typical one. Region pooling is therefore
broadly defensible at the scale actually in use, and the right architecture is
**DE-aware** (pool where measured DE is low, fall back to the read-level path
where it is not) rather than a wholesale replacement. DE is measurable per
region, so this is a switch, not a research programme.

---

## Log

### 2026-09-15 — provenance resolved, DE across 11 donors, estimators reconciled, a scoop read (`5cc5840`, `f06d100`, this commit)

- **Merged** `cpel-evaluation-and-pattern-path` into `master` (fast-forward).
- **Production-table provenance: resolved, no re-run needed.** All 140 cluster tables have max span
  exactly 1000 bp (post-cap); `ccabe0f` never touched the code production calls. Two real findings
  instead: production still uses the **global** dispersion estimator, and
  `G02_betabinom_region_chr1.py` is **not in version control** anywhere. → `docs/PROVENANCE.md`
- **The dipcall re-call never ran.** P01/P02 log directories are empty, there is no job accounting,
  and none of the six sane donors has a `_dip.vcf.gz`. Every locus number is still built on
  panel phasing.
- **Design effect, 11 donors × whole chr15** (`P09`): median 1.03–1.18 per donor, 0.5–0.9% of
  windows above 5, **no separation between sane, flagged and zero-hit donors**. DE explains neither
  anomaly. HG00146's whole-chromosome median is 1.09; the earlier 1.22–1.29 came from the first
  few Mb. → `benchmarks/results/2026-09-15_chr15_design_effect_12donors.txt`
- **Estimator disagreement closed** (`P10`): `implied_design_effect` was given the window's union of
  positions (13) as C; reads carry ~4.3. With the right C the estimators agree within ~0.1. The steep
  short-range decay survives filtering to well-supported positions, so it is real.
- **P09 defect found and fixed:** a truncated `.gz` exited 0. NA18508's local extract is truncated.
- **Scoop assessment:** Meredith et al. 2026 (ASM-LR, NIH CARD brain cohort) is a *phased mQTL*,
  not per-individual ASM calling. It takes long-read phased mQTL, SVs as drivers and the name; it
  cannot reach imprinted or non-genetic ASM, per-individual calls or molecule structure. →
  `asm_lr/md/20260915_meredith_2026_asmlr_comparison.md`


### 2026-09-08 — literature survey: what everyone else actually does (`5cf73b2`)

Read the methods sections rather than the abstracts, because the answer turns on
details abstracts omit.

- **Stefansson 2024** (deCODE, 7,179 ONT genomes): the per-individual ASM step
  is a bare threshold rule — ≥3 CpG units each below 0.15 within 500 bp, no
  test, no p-value. All rigor is downstream, in a cohort-scale genotype
  regression at Bonferroni 1e-10.
- **Akbari 2022** (eLife, 12 LCLs) is the closest published design to this
  project, and routes ASM through NanoMethPhase `dma`, a DSS wrapper — i.e. it
  hits the same replicate/haplotype mismatch this repo found independently. No
  per-sample calibration, no permutation, no formal FDR; cross-sample recurrence
  substitutes.
- **Voronina 2025** states no ASM criterion anywhere.

**The fair framing** (recorded in the doc, and worth reusing verbatim in any
manuscript): the field is not careless — the two regimes it has solved are
*n* very large (population statistics carry the inference, so per-individual
calling can be crude) and *n* = 1 (cross-donor comparability never arises). This
project sits between them, at 18 donors with per-individual calls and 8×–18×
coverage variation, where neither escape route is available. Say that, not
"nobody has thought hard about this".

Also: Voronina's headline finding (ASM regions carry 4–10× local SNV density) is
this project's artifact class un-partitioned — `VALIDATION_PLAN.md` D1/D2 are
designed to split it.

→ `docs/2026-09-08_what_do_others_do.md`

### 2026-09-08 — provenance stamping, and the unrecorded cluster runs (`6128282`)

Nothing recorded which commit produced an output, and `pyproject.toml` sat at
0.1.0 across six commits, three of which changed results (`513766e` region
pooling, `6df55b3` the `max_span` cap, `ccabe0f` dispersion calibration — all
landed 2026-09-05, with the production tables written somewhere inside that
window). Not hypothetical: `6df55b3`'s own message describes uncapped behaviour
confirmed on real chr1 data, so it is describing production output that already
existed. P05 reuses those same region boundaries, and DE grows with the number
of CpGs pooled — wrong geometry would mean every DE and FDR number computed on
it describes the wrong regions.

- `provenance.resolve_region_table_version()` settles the geometry half from an
  existing table: one region wider than `max_span` is proof it predates the cap.
  The converse is weaker and the verdict says `consistent-with-post-cap` rather
  than claiming proof.
- The dispersion half leaves **no structural fingerprint**. `PROVENANCE.md` says
  so and recommends re-running rather than reasoning backwards from
  plausibility. **This is still outstanding.**
- `provenance.write_stamp()` drops `<output>.provenance.json` beside any file
  (commit, dirty flag, branch, versions, argv, analysis parameters). Unavailable
  commit is recorded as `null`, not omitted — unknown provenance must be visible
  as unknown.

Version → 0.2.0. → `docs/PROVENANCE.md`, `ont_asm_caller/provenance.py`

### 2026-09-07 — LaTeX methods and theory writeup (`d0ac8fb`)

14 pages covering the three testing paths and the theory behind them.
Derivations included rather than asserted: beta-binomial variance inflation; why
dispersion is pooled across loci rather than replicates (with the ASE
precedent); the boundary failure of the Pearson-χ² dispersion estimator and why
the original benchmark could not have caught it; and `DE = 1 + (C−1)·r̄`, which
makes explicit that the design effect has exactly two ingredients.

Two results derived here for the first time:

- **The basecall-error attenuation law.** For symmetric error ε the correlation
  scales by (1−2ε)² exactly at μ=1/2 and slightly more strongly away from it —
  so de-attenuating by (1−2ε)² under-corrects, and is therefore safe.
- **The BH floor argument, quantified.** A p-value floored at 1/(L+1) needs
  k ≥ n/(α(L+1)) tests tied at the floor before anything is callable: 14,289 at
  CPEL's scale.

The writeup is explicit about where the simulators are *not* high-fidelity, and
records that `simulate_patterns.molecules_markov`'s single-exponential
assumption is refuted by the real data in §9.

→ `docs/paper/method.tex`, `method.pdf`, `make` to rebuild (tectonic, no sudo)

### 2026-09-07 — first real read-level measurement; PRIORITIES.md (`3a69861`)

`readlevel.design_effect` on **HG00146 chr15**, 12M `modkit extract` rows
(~2.5 Mb, 9,049 windows, median 68 reads / 13 CpGs per window), strand-normalised
and confidently called:

| window | median DE | IQR | 90th pct | % DE>2 | % DE>5 |
|---|---|---|---|---|---|
| **500 bp** (the caller's own scale) | **1.29** | 1.05–1.72 | 2.20 | 15% | **1%** |
| 1000 bp (`max_span`) | 1.48 | 1.17–2.06 | 2.78 | 27% | 1% |

This **demotes a large part of the work built in response to the calibration
critique** — the read-level path, pattern path and matched null are all still
correct and still worth having, but the case for making any of them the
*default* is much weaker.

Second finding: **a single exponential does not fit the real decay** (r² =
0.25–0.40, against >0.95 on simulation). The structure is two-component — a
steep ~20–25 bp decay plus a flat long-range floor near 0.07–0.09 that does not
decay within 2 kb — with a bump at **167 bp**, the nucleosome repeat length (one
donor, one bin — suggestive, not a finding). Call error matters: ε ≈ 0.043 even
after filtering at `CONF_THRESH = 0.80`, attenuating every correlation by 0.836.

`implied_design_effect` now integrates the empirical curve rather than the
exponential fit; that narrowed its disagreement with `design_effect` (3.19 →
2.77 at 500 bp) but did not close it against a measured 1.29. **Recorded as
unresolved, not claimed as agreement.**

Adds `P08_decay_from_extract.py`, which measures decay and DE from an extract
file alone — no haplotagged BAM, no region table — because co-methylation along
a molecule is not a haplotype quantity.

→ `docs/PRIORITIES.md`, `benchmarks/results/2026-09-07_HG00146_chr15_decay_500bp.txt`

### 2026-09-07 — two defects in P06, found against real extract output (`c0710f8`)

- **Strand.** `modkit extract` is per-strand and does not combine: the same
  physical CpG is reported at `p` from plus-strand reads and `p+1` from
  minus-strand reads. Keying columns on raw `ref_position` split every CpG into
  two columns 1 bp apart sharing zero reads — which `design_effect` survives but
  which costs the decay curve its short-distance bins and doubles `n_cpgs`, a
  null stratification variable. Now normalised, with the plus/minus split
  reported.
- **Call confidence.** `extract` emits one row per (read, position, mod_code).
  A real row from this cohort reads `p_h=0.686, p_m=0.314`; the old
  `mod_qual >= 0.5` rule on the `m` row alone scores that a confident
  *unmethylated* call. At a true DE of 6.3, ε of 0.02/0.05/0.10/0.20 gives
  5.9/5.3/4.4/2.9 — the whole range from pooling-is-fine to
  pooling-cannot-control-FDR. P06 now takes the argmax over {C, m, h} and
  requires `CONF_THRESH = 0.80`, emitting a no-call otherwise. A no-call costs
  data; a wrong call costs calibration.

**Consequence for an existing number:** the handoff design effects (1.39, 1.13)
came from P05 under the same rule, so they are attenuated **lower bounds**, not
estimates. Direction makes the qualitative conclusion safer, not weaker — but
P05 should be rerun before they are quoted. *(Still outstanding; PRIORITIES #2.)*

### 2026-09-07 — HG00344 doc contradiction; two-donor overdispersion diagnostic (`35cae3c`)

`VALIDATION_PLAN` recommends HG00344 as a first-pass clean donor while
`phasing_qc/README` flags it as one of the two anomalous zero-hit donors
(ρ̂ = 0.097, 5–10× cohort). Both cannot be right. Quickstart now pairs a clean
donor with the anomalous one, so baseline and diagnostic come from the same
extraction.

The diagnostic itself: a global ρ is the only place a pooled-count model can put
within-molecule correlation, so genuinely elevated DE and scrambled haplotype
labels produce **the same symptom** and cannot be separated. The read-level and
pattern paths have no ρ at all, which makes the question decidable — including
the sharp version: switch errors make a mislabelled allele look *bistable*, so
they should inflate `T_NME` while attenuating `T_MML`.

### 2026-09-07 — pattern-method benchmark across three generative models (`3d44493`)

CPEL's statistics do not win anywhere, **including on the Ising-generated
scenario where CPEL is correctly specified** (best `t_pdm` power 0.131 against
the pooled beta-binomial's 0.295). CPEL's own null protocol (Cmin=5, stratified
on CpG count alone) called **zero** regions in all three scenarios, floored at
p ≈ 1e-3, while the identical statistic against the depth-matched null called
14/32/5 — **the null, not the statistic, is the binding constraint.**

The correlation axis is not reachable at this depth: maximum power on equal-mean
disorder imbalance is **0.034** across every method and scenario, and `t_nme`
scored 0.000 on that axis in all three while keeping power on the mean axis —
i.e. at 10–30 reads, NME largely tracks the mean.

Two limitations recorded rather than papered over: all three scenarios sit at
DE 2.1–3.5, the benign end of the curve where the pooled test's FDR goes 0.040 →
0.070 → **0.695** across DE 1.0 → 2.4 → 5.5, so the apparent winner is the
method whose failure mode is not being tested; and the epiallele scenario, built
to favour neither method, turns out to have infinite decay length — still
exchangeable. **"pooled_bb wins everywhere" should not be quoted from this
table.** *(PRIORITIES #7.)*

→ `benchmarks/results/2026-09-07_pattern_methods.txt`

### 2026-09-07 — CPEL evaluated, statistics ported, permutation floor fixed (`1804498`)

Closes handoff item 4. **Verdict: not adopted** — statistics ported, machinery
not. (The paper is **Abante, Fang, Feinberg & Goutsias, Nat Commun 2020**, not
"Jiang et al." as the README and handoff note both had it. Corrected.)

Why not adopt: the Ising model exists to solve a problem long reads do not have
— a ~100 bp bisulfite read observes a fragment of a haplotype, so the joint
pattern must be inferred by marginalising over everything each read missed
(48 h on 20 CPUs per sample). A 35 kb ONT read observes every CpG directly,
collapsing the fit to a small concave complete-data problem (~5 ms). Three
further blockers: null draws at global `Cmin=5` regardless of real depth,
stratified on CpG count alone; p-values floored at 1/(L+1); and testing
restricted to SNP clusters, where long reads phase whole chromosome arms.

What *is* worth taking is the **entropy axis** — `readlevel.py` reduces each
molecule to a fraction, which handles the design effect correctly but discards
the within-read pattern.

- `simulate_patterns.py` — the old simulator could not adjudicate this:
  `simulate_reads._molecules` copies from one molecule-level latent, so
  corr(xᵢ,xⱼ) = q² at *every* lag. Under exchangeability the per-read fraction
  is sufficient, so a pattern method is **guaranteed** to find nothing extra —
  benchmarking CPEL there would be rigged (and simulating from an Ising chain
  rigs it equally hard the other way). Replaced with two orthogonal ASM axes and
  three generative models: epiallele (neither method's — the headline), markov
  (CPEL correctly specified), exchangeable (this package favoured).
- `pattern.py` — MML / NME / JSD via transfer matrix; every recursion checked
  against brute-force enumeration of the 2^N state space to machine precision,
  analytic gradient against finite differences to 1e-9. Plus `read_fraction_ks`.
- `null.py` — **worth more than the CPEL port.** Pools permutation draws across
  regions within (n_cpgs, depth, methylation) strata and fits a generalised
  Pareto tail, removing the floor that makes `test_region_perm` unusable
  genome-wide. Its own first run found a defect in it, recorded because of how
  it hid: unbounded GPD extrapolation returned p = 1e-300 from a 60,000-draw
  null *including on H0 data*, while type I error at α=0.05 and 0.01 still
  looked correct — only the 1e-3 row gave it away. Fixed by rejecting
  pathological shapes, requiring 200 exceedances, and capping extrapolation at
  100× the empirical floor. `calibration_report` now prints the 1e-4 row by
  default so it cannot hide again.
- `decay.py`, `P06`/`P07` — measure `decay_bp`; export read matrices once (tens
  of MB) instead of restreaming 4–5 GB per analysis. `load_read_matrices` yields
  objects duck-typed to `PatternRegion`, so real and simulated data share one
  code path.

Tests 27 → 62. → `docs/2026-09-07_cpel_review_and_pattern_tests.md`

### 2026-09-05 → 2026-09-07 — on `master`, before this branch

`513766e` region-level pooling (recovers power lost to genome-wide BH-FDR) ·
`b449543` README · `6df55b3` `max_span` cap, fixing daisy-chaining through
CpG-dense stretches · `b1570d6` CPEL / field-maturity context · `ccabe0f` the
calibration critique's fixes: region pooling, dispersion calibration, and the
read-level path · `55cce16` the first real-data design-effect measurement plus
the handoff note to the methods-dev arm.

---

## Still to do

`docs/PRIORITIES.md` is the ordered list and the authority. In brief:

1. **Get the haplotagged BAM + chr15 region tables** — `bam/haplotagged/by_chrom/{sample}_chr15_haplotagged.bam`
   and `tables/dasm_betabinom/{sample}/regions_chr15.tsv`, one chromosome for
   2–3 donors. **Every ASM-specific test in the package is blocked on this** and
   nothing else in the list unblocks it. Single highest-value transfer.
2. **Rerun P05 through the corrected call logic** — its published numbers are
   attenuated lower bounds. P08 already has the corrected logic; cheap.
3. **Resolve `implied_design_effect` vs `design_effect`** — ~2× apart at 500 bp,
   agreeing at 2 kb. Two independent estimators of one quantity disagreeing is
   the tripwire that says the decay measurement is not trustworthy.
4. **Replace the single-exponential decay model** in both `decay.comethylation_decay`
   and `simulate_patterns.molecules_markov` with short-exponential + floor.
5. **Decide whether the entropy axis stays** — undetectable at 10–30 reads per
   haplotype; the surviving use is diagnostic (item 6), not discovery.
6. **The two unexplained-overdispersion donors (HG00344, NA21144)** — the
   `T_NME`-inflated / `T_MML`-attenuated signature would be near-positive
   identification of switch error, and it is a genome-wide aggregate, so it does
   not need the per-region power the entropy axis lacks.
7. **The DE sweep the pattern benchmark is missing** — all scenarios landed at
   DE 2.1–3.5, so the regime that motivated the read-level path was never
   probed.

Carried separately from PRIORITIES:

- **Re-run the production tables under a stamped commit.** The dispersion half
  of the provenance question leaves no structural fingerprint and cannot be
  settled by inspection — see `docs/PROVENANCE.md`.
- **Merge decision for this branch.** It is 8 commits ahead of `master` with
  nothing in conflict; it has been left unmerged rather than blocked.
- **Reconcile the entropy axis with the parent project's
  `tables/W12_epiallele_metrics.csv`** (304,514 loci, `within_hp_entropy` /
  `bimodality_coef`) before either arm builds further on it — two arms measuring
  the same quantity independently, which has already happened once on this
  project.
- **Correct the Abante-vs-Jiang misattribution** anywhere it survives outside
  this repo, before it reaches a manuscript.

## Standing caveats

- The measured DE (1.29) is **one donor, one chromosome, one window**. It is the
  best number available and it is a real measurement, but it is not a cohort
  result, and the handoff note's own warning about over-reading a single window
  applies to it too.
- The 167 bp correlation bump is one donor, one bin. Suggestive, not a finding.
- `decay.comethylation_decay` cannot resolve a decay longer than the region
  span: 400 bp regions generated at `decay_bp = 500` return 88. On the caller's
  own regions (`max_span = 1000`), a fitted value near the mean span is a
  **lower bound**.
- No simulator here is high-fidelity, and `method.pdf` says where each one
  fails. Simulation results in this repo bound behaviour; they do not establish
  it.
