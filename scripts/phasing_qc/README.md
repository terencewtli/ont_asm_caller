# phasing_qc — donor-specific ASM-caller diagnostics

**Why this directory exists.** The beta-binomial region-level ASM caller
(`ont_asm_caller`, separate repo) estimates a single global overdispersion
parameter `rho` per (donor, chromosome) and applies it uniformly to every
tested region's likelihood-ratio test (see `ont_asm_caller/dispersion.py` /
`model.py`). Running this on real chr1 data across the 18-donor cohort
surfaced two donor-level anomalies that neither a coverage explanation nor a
switch to this new caller resolved:

- **4 donors with an implausibly high "significant ASM" rate**
  (NA18508 8.7%, HG01167 5.4%, NA21093 1.8%, NA18620 1.9% — vs. deCODE's
  validated 0.51% benchmark and 0.12-0.53% for the cohort's confirmed-sane
  donors), despite **unremarkable, sometimes even lower, rho_hat** than
  several sane donors (0.005-0.01 vs. 0.001-0.018) — so this isn't obviously
  a dispersion-estimation failure.
- **2 different donors (HG00344, NA21144) with a near-opposite failure**:
  rho_hat 5-10x higher than every other donor (0.097, 0.116) and essentially
  **zero significant regions**, i.e. so much estimated noise that the test
  has no power left.

Neither pattern tracks ancestry, sequencing batch, coverage, or run count —
all already tested and ruled out (see the project's `md/2026090*` progress
docs). The leading open hypothesis is that population-panel-only statistical
phasing (no read-based phasing is in production yet) is producing switch
errors that scramble which reads get called "haplotype 1" vs "haplotype 2",
and that this shows up as one of two different signatures in the caller
depending on whether the resulting cross-haplotype noise is diffuse (inflates
`rho`, kills power — the HG00344/NA21144 pattern) or concentrated/systematic
in a subset of regions (survives the test as false-positive "ASM" — the
NA18508/HG01167/NA21093/NA18620 pattern).

## `P04_permutation_null.py`

A donor-specific parametric-bootstrap null, meant to distinguish "the
caller's own model/dispersion fit is miscalibrated for this donor's depth
profile" from "the input haplotype labels feeding the caller are the
problem."

**What it does:** for a given (sample, chromosome), take the REAL per-region
depths (`n1`, `n2`) and the donor's already-fitted `rho_hat` from the
production run (`tables/dasm_betabinom/{sample}/regions_{chrom}.tsv` +
`summary_{chrom}.json`), and simulate new counts under the strict null
(`x1*, x2* ~ BetaBinomial(n, p_hat, rho_hat)`, independently per haplotype,
where `p_hat` is the region's real shared rate estimate). Rerun the exact
same LRT + BH-FDR + delta-floor significance rule the production caller
uses, and report what fraction of regions come out "significant" by chance
alone, averaged over several bootstrap replicates.

Because this null is built from each donor's own fitted `rho_hat` and real
depth distribution — not an external population benchmark or the bare
asymptotic chi-square(1) assumption — it functions as a **learned,
donor-specific filtering baseline**: a donor whose bootstrap null already
predicts an elevated "significant by chance" rate has a caller/dispersion
problem specific to its own depth/coverage profile; a donor whose bootstrap
null stays sane while the real data doesn't implicates the haplotype labels
(phasing) instead, since the model's own assumptions, applied to that
donor's own real depths, don't produce the excess on their own.

**Usage:**
```
python3 P04_permutation_null.py <SAMPLE> <CHROM> [n_boot]
# e.g. python3 P04_permutation_null.py NA18508 chr1 10
```
Reads `tables/dasm_betabinom/{SAMPLE}/regions_{CHROM}.tsv` +
`summary_{CHROM}.json` (must already exist), writes
`tables/phasing_qc/permutation_null_{SAMPLE}_{CHROM}.json` with the real vs.
bootstrap-null `pct_significant` and their ratio.

**Caveat:** this resamples from *counts* (the per-region `x1/n1/x2/n2`
already aggregated by the pileup step), not from raw reads — it tests
whether the caller's own beta-binomial model, at this donor's fitted `rho`
and real depths, would manufacture the observed excess on its own. It does
**not** directly test read-level phenomena (e.g. reference-mapping bias
correlated with a switch point) that a true read-level permutation
(reshuffling HP tags in the BAM and re-running the pileup step) would catch.
Treat a "sane bootstrap null + elevated real data" result as evidence
pointing at phasing/inputs, not as proof — the read-level version is the
follow-up if this doesn't fully settle it.

---

## `P06_export_read_matrices.py` / `P07_comethylation_decay.py`

Added 2026-09-07 alongside the CPEL evaluation
(`github/ont_asm_caller/docs/2026-09-07_cpel_review_and_pattern_tests.md`).

**`P06` exists to stop re-streaming multi-GB files.** `P05` streams a 4–5 GB
gzipped `modkit extract` per sample per window and its own docstring notes this
will not scale genome-wide. Every read-level analysis in `ont_asm_caller` wants
the same thing out of it — a (read × CpG) 0/1 matrix per region per haplotype —
so `P06` extracts that once into a packed `.npz`
(`tables/read_matrices/{sample}_{chrom}.npz`). At 10⁵ regions × ~40 reads × ~8
CpGs that is tens of MB rather than gigabytes, which is what makes the methods
work iterable off-cluster.

`load_read_matrices()` yields objects duck-typed to
`simulate_patterns.PatternRegion`, so `decay.py`, `pattern.py` and `null.py`
consume real and simulated data through the same path with no separate real-data
code path to keep in sync.

**`P07` measures `decay_bp`, and should be run first.** Three analyses are
conditional on that single number — whether `region.cluster_cpgs` pooling is
legal at all, every result in `benchmarks/compare_pattern_methods.py`, and the
parent project's co-methylation-decay curve (`md/20260903_qc_review.md` §7.3) —
and none of them can be quoted before it exists. `P07` prints the fitted decay
length beside `readlevel.design_effect` on the same regions; those are two
independently written estimators of the same quantity, and **if they disagree,
that is the finding**, not a nuisance.

**Resolution caveat.** The curve only sees CpG pairs inside a supplied region, so
it cannot resolve a decay longer than the region span. The caller caps regions at
`max_span = 1000` bp, so run on the caller's own regions `P07` measures decay
only out to ~1 kb; a fitted value near the mean span is a *lower bound*.
Measured on simulation: 400 bp regions generated at `decay_bp = 500` return 88.
