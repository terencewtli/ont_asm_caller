# "What do other people do?" — how long-read ASM is actually called

**2026-09-08.** The first question any reviewer will ask of this repo. Written
from a targeted read of the primary sources rather than from abstracts, because
the answer turns on methods-section details that abstracts omit.

**Short answer: almost nobody performs a calibrated per-individual statistical
test for ASM from long reads. The field has three escape routes from the problem
and takes one of them; none of the three is available to a mid-size cohort that
wants per-donor calls.**

---

## The survey

| Study | Design | How an ASM call is actually made | Depth handling | Calibration |
|---|---|---|---|---|
| **Stefansson 2024** (deCODE, *Nat Genet*) | 7,179 ONT genomes, mean 20.6× (range 10–108×) | **A bare threshold rule.** "Three or more CpG units each with 5-mCpG rate <0.15, but located no more than 500 bp apart." No test, no p-value. | **Exclusion.** Restricted to methylomes >20× (n=2,648) and CpG units with ≥10 reads per haplotype. | None per individual. All rigor is downstream: cis-variant regression across the cohort, Bonferroni ~10⁻¹⁰. |
| **Akbari 2022** (*eLife*) | **12 LCLs**, ONT — the closest published design to this project | NanoMethPhase `dma`, which is **a wrapper around DSS**. Fixed cutoffs: p<0.001, \|Δ\|>0.20, recurrence in ≥4 cell lines. | CpGs with ≥5 reads. | **No per-sample calibration, no permutation, no formal FDR.** Cross-sample recurrence substitutes. |
| **Voronina 2025** (*IJMS*) | n = 1, 112× ONT, trio-binned | **Not stated anywhere.** "58,470,893 methylation sites in total. However, only 308,673 (0.5%) *were determined to be* ASM sites." Methods describe Megalodon, Clair3, WhatsHap and NanoMethPhase `phase` — but no ASM criterion, threshold, test or tool. | Subsampling at 12/26/79/112× is used to assess **phasing N50**, not ASM calling. | None. |
| **Tian 2026** (nanoASM, preprint) | prostate tissues | **Does not phase.** Reads grouped by the allele they carry at one variant, pooled across individuals, then `metilene` DMRs between REF and ALT groups. | — | Cross-individual pooling substitutes; a hybrid of ASM and mQTL. |
| **Abante 2020** (CPEL, *Nat Commun*) | WGBS, 1 individual, 10 tissues | The one genuinely per-individual statistical method: Ising model + empirical bootstrap null. | **Null generated at C_min = 5 reads regardless of a region's real depth.** | Null stratified on CpG count alone — not depth, not methylation level. p floored at ~10⁻³. |
| **pycoMeth** (*Genome Biol* 2023) | tool | Supports ASM within one sample via read-group tags; offers a menu of tests chosen by hypothesis and sample count. | — | — |

## The three escape routes

1. **Go to cohort scale, so per-individual power stops mattering.** deCODE's
   per-individual step is a threshold rule that would be indefensible on its own;
   it does not have to be defensible, because 7,179 genomes and a genotype
   regression carry the inference. Tian 2026 does the same thing differently, by
   pooling reads across individuals before testing.
2. **Reuse a group-comparison tool.** Akbari 2022 — the design closest to this
   project — routes ASM through DSS. **DSS estimates overdispersion from
   between-replicate variance, and a single individual's two haplotypes are not
   replicates.** With exactly one "replicate" per group everywhere in the genome
   that quantity is undefined, not merely noisy. This is the same mismatch this
   repo documented independently, and on this project's real data DSS's
   dispersion step ran indefinitely without output.
3. **Do not state a method.** Voronina 2025 reports an ASM rate to six
   significant figures without saying how a site qualified.

## What nobody does

**Nobody calibrates ASM calling across donors of different sequencing depth.**
Searching the coverage-benchmarking literature returns work on *basecaller*
accuracy — how many reads before Nanopolish/Megalodon/DeepSignal agree — and
recommendations like "10× for stable methylation estimates." That is about
whether a methylation *rate* is well estimated. It is not about whether a
*significance threshold* means the same thing in an 8× donor and an 18× donor.
The one paper that explicitly varies coverage (Voronina, at 12/26/79/112×) uses
it to evaluate **phasing N50**, not ASM calibration.

The >50,000× swing in candidate rate this project measured between an 8× and an
18× donor at identical thresholds is therefore an empirically demonstrated
problem that the existing literature does not frame, name, or control for.

## The fair version of the claim

The gap is real, but it is **not** evidence that the field is careless. It is a
gap created by a specific design choice. The two regimes the field has solved are:

- **n very large** (deCODE): per-individual calling can be crude because
  population statistics carry the inference.
- **n = 1** (Voronina, CPEL): cross-donor comparability never arises.

This project sits in between — **a mid-size cohort (18 donors), per-individual
calls, cross-donor comparison, and coverage varying 8×–18× per haplotype by
accident of data availability.** That is exactly where neither escape route is
available: too few donors for a QTL regression to rescue a crude caller, too
many for cross-donor comparability to be ignorable. Akbari 2022 is the one paper
in that regime, and it handles it with fixed thresholds plus a recurrence filter
rather than calibration.

**State it that way.** "Nobody has thought hard about this" is attackable and
slightly unfair. "The published approaches solve the problem by moving to a
scale where it disappears, and that route is not available at 18 donors" is
accurate, specific, and much harder to argue with.

## Two more things worth knowing

**Voronina's headline finding is this project's artifact class.** They report
ASM regions carry a "4- to 10-fold increase in local SNV density" and that
"polymorphisms directly at CpG sites account for 14% of all detected ASM" —
without separating CpG-destroying SNVs, basecaller sequence-context artifacts,
and genuine TF-mediated regulation. That partition is exactly what
`VALIDATION_PLAN.md` D1 (CpG-destroying SNVs) and D2 (1 bp-resolution
distance-to-variant profile) are designed to do. Their number is a ceiling on
how much of the enrichment is mechanistic; this project can say which part.

**Phasing method is a live confound in this literature, not a settled detail.**
deCODE used *sequence-based* (population-panel) phasing on 7,179 genomes —
the same statistical-phasing route this project currently uses, and the leading
suspect for its four flagged donors. Voronina found that without trio binning,
increasing ONT coverage above 25× does not improve phasing N50 at all, and that
CpGs called differentially methylated only without trio binning are phasing
errors. That is independent support for treating switch error as a first-class
source of false ASM — and it is measurable here (`W11_phasing_comparison.sh`)
without needing trios.

## Sources

Stefansson et al., *Nat Genet* 56:1624–1631 (2024), `reference/stefansson_2024/`.
Voronina et al., *Int J Mol Sci* 26:9641 (2025), `reference/voronina_2025/`.
Akbari et al., *eLife* 11:e77898 (2022). Abante et al., *Nat Commun* 11:5238
(2020). Tian et al., bioRxiv 2026.06.17.732357, `reference/tian_2026/`.
Snajder et al., *Genome Biol* 24:83 (2023) (pycoMeth). Akbari et al.,
*Genome Biol* 22:68 (2021) (NanoMethPhase).
