# Real-data validation plan

**2026-09-05.** Everything in this package is simulation-only. This is the
ordered plan for validating it on the parent project's HPRC cohort, written
before the data exists so the analyses are not chosen after seeing results.

Ordering is deliberate: **A blocks everything else.**

---

## Required input

**`modkit extract`, not `modkit pileup`.** This is the one thing that matters
and the two are easy to confuse:

| | what it gives | usable here? |
|---|---|---|
| `modkit pileup` | per-CpG aggregated counts (bedMethyl) | per-CpG test only; **cannot** measure the design effect |
| `modkit extract` | one row per (read, position) | **yes** — required |

`pileup` has already aggregated away molecule identity, which is precisely the
information every conclusion in `docs/2026-09-05_calibration_critique.md`
depends on. Minimum columns needed:

- `read_id`
- `chrom`, `ref_position`
- modification call (`mod_code` + `mod_qual`, so a probability cutoff can be applied)
- **haplotype** — if `modkit extract` does not emit `HP` directly, join `read_id`
  against `samtools view -d HP:1` / `HP:2` from the haplotagged BAM

Reads must come from the **haplotagged** BAM, and untagged reads should be kept
and labelled as such — the untagged fraction is itself a QC metric (see §2 of
`md/20260903_qc_review.md`).

### Which chromosome and which donors

**chr15.** `tables/zink2018_PofO_DMR.bed` holds 229 imprinted DMRs and **98 are
on chr15** (SNRPN/UBE3A); the next-best chromosome has 16. One extraction gives
both the design-effect measurement and a 98-region positive-control set.

**Stage the extraction — do not start with all 18 donors.** chr15 at 80× is on
the order of 10⁸ read×position rows per donor.

1. **Two donors first:** `HG00344` (56×) and `NA19682` (87×). Complete
   downloads, not on the problem-donor list, and a ~1.6× depth contrast.
   Sufficient for analyses A, and for D on one sample.
2. **Then all 18**, but for the depth-invariance figure (B1) you can restrict
   to Zink DMR regions ±5 kb plus a matched random background set — far smaller
   than whole-chromosome and enough for every analysis below except A.

Avoid `NA20762, HG00126, NA19776, HG00253, HG03784, NA20870` — confirmed
incomplete downloads. Avoid `NA18959` (no S3 data) and treat `NA21093` with
caution (ntsm outlier 0.278).

---

## A. Prerequisite: measure the design effect

**This decides the architecture. Nothing else should be run first.**

```python
from ont_asm_caller import design_effect
de = design_effect(m1, m2)   # per region; summarise the distribution
```

- **DE ≈ 1** → CpGs behave independently; `cluster_cpgs` + `test_regions` is
  valid *and more powerful*. Use it.
- **DE ≫ 1** → pooled region testing cannot control FDR (measured: FDR 0.695 at
  DE 5.5). Use the read-level path.

Report DE as a distribution, not one number — it will vary with CpG density and
region width, and the calibration argument needs the upper tail, not the mean.

**Same measurement, second use:** DE as a function of CpG–CpG distance *is* the
co-methylation decay curve of `md/20260903_qc_review.md` §7.3. One extraction
settles both analyses, and they must agree.

## B. Method validation — no external truth set needed

### B1. Depth invariance across donors — the primary figure

Calls per donor against median per-haplotype depth, one line per method.

The original Fisher pipeline produced 113,119 calls in one sample and 2 in
another at matched thresholds. That is the defect this package exists to fix,
so the direct test is whether call count still tracks depth. **Fisher should
show a steep slope; the new caller should be flat.** If it isn't flat, that is
a finding and it invalidates cross-donor and cross-ancestry comparison —
better learned here than after writing the paper.

Confounder to control: coverage correlates with ancestry in this cohort, so
plot ancestry as a covariate, not just depth.

### B2. Cross-donor replication

An ASM locus driven by a common cis variant should recur in donors carrying
that variant and be absent in donors that don't. Stratify replication rate by
carrier status. Requires no annotation and is closer to the paper's actual
claim than any overlap statistic.

### B3. Permutation null on real data

Permute haplotype labels of whole reads within each region and re-run the full
pipeline. The resulting p-value distribution should be uniform. This is the
strongest available calibration check because it uses the real co-methylation,
depth and sequence-context structure rather than a simulator's version of them.

## C. Biological validation — external truth sets

### C1. Imprinted DMR recovery at MATCHED CALL COUNTS

**Do not compare "how many loci did each method call."** Fisher calls more and
most are false; a raw count table makes the better method look worse and hands
a reviewer the wrong summary statistic.

Instead: take the top-*N* ranked calls from each method **for the same *N***,
and measure recovery of the 98 chr15 Zink DMRs. Sweep *N* to get a
precision-recall curve. This removes the count confound and turns the
comparison into a statement about ranking quality.

### C2. Rosenski 2025 atlas

Second truth set — 34k regions where local variation segregates with
methylation, across 39 cell types. Restrict to the lymphoid/LCL-relevant
subset; the atlas's own headline finding is that ASM is frequently cell-type
restricted, so a whole-atlas overlap would understate performance.

### C3. Not available: chrX / XCI

The BAMs are aligned to `GRCh38.autosome.fa` — 22 contigs, no chrX. The XCI
positive control in `W10_control_loci.sh` cannot run without realignment. Noted
so it is not proposed again.

## D. Negative controls and artifact partitioning

These connect to `md/20260903_qc_review.md` §6 and should be run on the *same*
call set, not separately.

### D1. CpG-destroying SNVs (§6.1)

What fraction of calls sit at CpGs destroyed by a heterozygous SNV? These are
real allelic differences but not allele-specific regulation. `W13_cpg_snv_audit.py`
already classifies them. Report the fraction and exclude before any mechanistic
interpretation.

### D2. Distance-to-variant profile (§6.3 test 5)

Effect size against distance from the phasing variant, at **1 bp resolution**
over 0–50 bp. Three mechanisms predict different shapes: CpG-SNP degeneracy at
0–1 bp; basecaller context artifact confined to the sensing window (~±5–10 bp,
decaying); TF-binding disruption at motif scale and beyond. Voronina 2025
reported a 4× SNV enrichment near ASM sites without separating these — this
profile is the separation.

### D3. Strand symmetry

CpG methylation is symmetric under DNMT1 maintenance, so independently called
forward- and reverse-strand cytosines should agree. Requires
`modkit --no-combine-strands`. A sequence-context artifact will not generally
affect both strands equally.

---

## Deferred

- **Cross-platform HiFi replication** (§6.3 test 2) — 5 HPRC samples have both.
  Polymerase kinetics has no k-mer sensing window, so an ONT context artifact
  should not reproduce. The most durable control, since it survives basecaller
  model churn.
- **Switch-error modelling.** `simulate_read_regions` can generate it but no
  test accounts for it. `W11_phasing_comparison.sh` measures the real rate;
  wire that estimate into the test rather than assuming perfect phasing.
- **Read-pattern (epiallele) test.** The current read-level test reduces each
  molecule to a fraction, discarding the within-read pattern. A distributional
  test on read-level methylation *vectors* would use it, and would unify this
  caller with the parent project's epiallele analysis instead of merely being
  compatible with it.

## What must not be reported

- **DSS as "0 loci called."** It hung; that is an implementation failure, not a
  statistical result, and presenting it as a benchmark row would be misleading.
  The correct treatment is a Methods paragraph: DSS estimates dispersion from
  between-replicate variance, a single individual's two haplotypes are not
  replicates, and that quantity is undefined here.
- **Any FDR number from `benchmarks/compare_methods_v2.py`** before analysis A
  has established where the real design effect sits.
