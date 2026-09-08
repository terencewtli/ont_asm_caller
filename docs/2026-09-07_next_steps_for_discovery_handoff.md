# Next steps for the methods-development arm — handoff note

**Written 2026-09-07**, by a session that was doing discovery/QC work in the parent project
(`asm_lr/`) and stumbled onto `2026-09-05_calibration_critique.md` and the `readlevel.py` /
`estimate_dispersion_trend` fix via a `git fetch` that surfaced a commit from a different,
parallel session. Two sessions were independently working the same underlying problem without
knowing it — see the parent repo's `CLAUDE.md` for the guardrail added because of this. This
doc exists so the methods-development session doesn't have to reconstruct that context, and so
the discovery-focused session (this one, going forward) doesn't block on methods work it
doesn't need to do itself.

**Explicit priority framing from the project owner (2026-09-07):** the discovery arm of the
project takes priority right now. This doc's job is to make sure nothing here gets lost or
silently blocks that arm — not to argue the methods work should happen first.

## What's confirmed so far

Implemented `readlevel.design_effect()` against real `modkit extract` data for the first time
(`scripts/phasing_qc/P05_design_effect_real_data.py`, mirrored from the parent project — HP
tags recovered separately from the haplotagged BAM since modkit extract isn't itself
haplotype-split). Pilot: chr1:1-5,000,000, 300 of the production caller's own regions, two
donors:

| Donor | median DE | mean DE | % DE>2 | % DE>5 |
|---|---|---|---|---|
| HG00146 (sane donor) | 1.39 | 2.09 | 32.0% | 5.7% |
| NA18508 (flagged donor, 8.7% pooled-test "significant" on chr1) | 1.13 | 1.52 | 18.3% | 1.0% |

**Design effect is real and substantial on actual data — this is the first non-simulated
confirmation of Defect 1.** But note the direction: in this one window, the *sane* donor has
the higher design effect, not the flagged one. Don't over-read this — it's a single 5Mb
subtelomeric window with quite different HP-tagged read depth between the two donors in that
specific region (18,585 vs 7,079 reads) and may not be representative of the rest of the
chromosome, where the actual candidate-rate divergence between these donors lives. **This is an
open question, not a closed one**: does design effect actually correlate with which donors get
flagged, once measured properly (multiple windows spread across a whole chromosome, more
donors, ideally weighted toward the regions that are actually called significant)? That
measurement hasn't been done yet.

## Concrete next steps, in rough priority order

1. **Broader design-effect measurement.** Extend `P05_design_effect_real_data.py` (or rewrite
   for efficiency — the current version streams a multi-GB `modkit extract` file per sample per
   window, which won't scale to genome-wide) to sample multiple windows across all 22
   autosomes, for at least the 4 flagged donors (NA18508, HG01167, NA21093, NA18620), the 2
   zero-hit donors (HG00344, NA21144), and several of the 6 sane donors. The goal: does DE
   actually separate flagged from sane donors, or is it uniformly elevated cohort-wide (in which
   case it's a general calibration problem, not a donor-heterogeneity explanation)?

2. **Wire `readlevel.test_region_reads`/`test_region_perm` in as a confirmatory pass**, not a
   full re-call. The discovery arm's more immediate need is a defensible locus set for the 6
   sane/good donors' *already-called* significant regions (`tables/dasm_betabinom/{sample}/regions_chr1.tsv`
   with `significant==True`, plus the genome-wide runs) — run the read-level test on just those
   candidates, not genome-wide from scratch. This is a bounded, tractable task (thousands of
   candidate regions per donor, not billions of loci) and gives the discovery arm a
   correctly-calibrated-if-conservative subset to build on without waiting on a full pooled-test
   replacement.

3. **Root-cause the 2 zero-hit donors (HG00344, NA21144)**, still unexplained. Their rho_hat
   (0.097, 0.116) is 5-10x every other donor's, and Defect 2 (boundary bias) predicts rho too
   *low*, not too high — so whatever's happening to these 2 specifically doesn't obviously
   follow from either documented defect and needs its own look. One lead worth checking first:
   HG00344 has the lowest coverage in the cohort (51x manifest / 56.35x measured, vs. 67-110x for
   everyone else) — NA21144 doesn't share that (79x/87.2x), so coverage alone won't fully
   explain it, but it's a place to start.

4. ~~**Evaluate CPEL**~~ — **DONE 2026-09-07**, see
   [`2026-09-07_cpel_review_and_pattern_tests.md`](2026-09-07_cpel_review_and_pattern_tests.md).
   Verdict: not adopted; statistics ported, machinery not. Note the paper is
   **Abante et al.**, not Jiang et al. — that misattribution is in this note and
   in the README's literature section, and both should be corrected before it
   reaches a manuscript. The item as originally written follows.

   **Evaluate CPEL (Nat Commun 2020)** as a possible off-the-shelf alternative
   rather than continuing to patch a custom pooled-count model. It's an Ising-model approach
   that jointly models methylation level *and* within-haplotype correlation/entropy structure —
   i.e. it explicitly models the exact co-methylation phenomenon that's breaking the pooled
   region test here, rather than assuming it away. Already flagged in the parent project's
   literature review as the closest published single-sample-ASM precedent; worth checking
   whether it's usable directly before investing further in an in-house fix.

## What NOT to block on

The discovery arm does not need a fully general, publication-grade single-sample ASM caller to
proceed with the 6 sane/good donors' positive-control analyses (imprinted-DMR enrichment,
cross-donor replicability) — those already show a strong, real signal (83% of hit imprinted DMRs
replicate in ≥2/6 donors vs. 29.9% background) that isn't obviously threatened by design effect
at the aggregate level, even if individual loci carry some inflated false-positive rate. Item 2
above is enough to unblock a defensible discovery-arm locus set; items 1, 3, and 4 are the
actual methods-development research questions and can proceed at whatever pace makes sense
independently.


---

## 2026-09-07 addendum: item 4 closed, and what it produced

Item 4 was worked in a separate session and is closed. Two of its outputs are
**not** CPEL-specific and change the priority of the items above:

- **`ont_asm_caller/null.py` removes the `1/(n_perm+1)` p-value floor.** It
  pools permutation draws across regions within a (n_cpgs, depth, μ) stratum and
  fits a generalised-Pareto tail, so a permutation-style test reaches p-values
  small enough to survive genome-wide BH. This applies to
  `readlevel.test_region_perm` unchanged. That promotes item 2 above from "run
  the read-level test on already-called candidates" to "the read-level test can
  now be run genome-wide", which is a larger change than item 2 assumed.
- **There is a second ASM axis** — haplotypes differing in within-read disorder
  at equal mean — that every test in this package is blind to by construction.
  It overlaps the parent project's `tables/W12_epiallele_metrics.csv`
  (304,514 loci, `within_hp_entropy` / `bimodality_coef`), which is a second
  instance of two arms measuring the same quantity independently. **Reconcile
  before either arm builds further on it.**

Items 1 and 3 are unchanged and still need real data. The new work's binding
constraint is the same measurement item 1 needs — see §6 of the review doc for
what would actually be sufficient (one chromosome, 2–3 donors, and *not* as raw
`modkit extract`).
