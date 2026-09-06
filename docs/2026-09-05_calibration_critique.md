# Calibration critique of v1, and what changed

**2026-09-05.** Review of the package as of the first release, with measured
numbers. Two defects were found; one invalidates the v1 region-level results.
Everything below is reproducible with `benchmarks/compare_methods_v2.py`.

The v1 diagnosis (Fisher is depth-confounded and unsafe) is correct and is not
in question. The problems are in what replaced it.

---

## Defect 1 — region pooling assumes CpGs are independent observations

`region.cluster_cpgs` sums raw counts across CpGs in a window:

```python
x1 = sum(c.x1 for c in group);  n1 = sum(c.n1 for c in group)
```

With 8 CpGs and 15 reads that is `N = 120`. But a 35 kb ONT read covers every
CpG in a 500 bp window, so there are **15 independent molecules, measured 8
times each**. The likelihood behaves as though it had ~8× more information
than it does.

**The v1 simulator hard-codes the assumption it needs to be tested against.**
In `simulate_spatial_mixture`, each CpG draws its own depth *and* its own
count inside the per-CpG loop. Its docstring says this simulates per-CpG noise
"honestly rather than assuming away" — it is the reverse: it assumes away the
read-sharing. No benchmark built on that simulator can detect this defect.

Measured, on null regions (μ₁ = μ₂ everywhere) with the same molecules
covering every CpG, re-estimating ρ on pooled counts exactly as the v1 README
prescribes:

| within-read co-methylation | type I @ α=0.05 | type I @ α=0.001 |
|---|---|---|
| 0.0 (v1's assumption) | 0.052 ✓ | 0.0015 ✓ |
| 0.5 | 0.060 | 0.0013 |
| **1.0** | **0.144** (2.9×) | **0.0112** (11×) |

It also worsens with region size — at co-methylation 1.0, type I rises from
0.118 (3–4 CpGs) to 0.164 (11–12 CpGs), because one global ρ cannot absorb a
design effect that scales with the number of CpGs pooled.

**The inflation is worst in the extreme tail**, which is exactly where BH-FDR
across millions of regions operates. An 11× inflation at p<0.001 is the regime
that sets the genome-wide false discovery rate.

### Why this matters beyond the package

The parent project exists to measure co-methylation along single molecules
(`md/20260903_qc_review.md` §7 — processivity, epialleles, methylation-LD
decay). The v1 region step assumes that phenomenon has magnitude **zero**.
The two analyses cannot both be right. They are the same parameter viewed from
two directions, and the ASM caller's calibration is downstream of the
co-methylation measurement.

## Defect 2 — the dispersion estimator only works at μ = 0.5

`estimate_dispersion` assumes `E[χ²ᵢ] = 1` under binomial sampling. That is
asymptotic, and it fails at the boundaries: when x₁=n₁ and x₂=n₂ (both
haplotypes fully methylated — the most common configuration in a real bimodal
methylome), p̂ → 1, the variance term collapses, χ²ᵢ → 0, and the locus
contributes `y = −1`, i.e. positive evidence for ρ = 0.

Recovering a known ρ:

| μ distribution | true ρ | ρ̂ | ratio |
|---|---|---|---|
| all μ = 0.5 (**v1 benchmark**) | 0.05 | 0.0493 | 0.99 ✓ |
| uniform(0.2, 0.8) | 0.05 | 0.0489 | 0.98 ✓ |
| **realistic bimodal** | 0.05 | 0.0273 | **0.55** |
| all μ = 0.9 | 0.05 | 0.0354 | 0.71 |
| all μ = 0.95 | 0.05 | 0.0180 | **0.36** |

At ρ=0.15 it is worse (0.28–0.49). `simulate_mixture` defaults
`baseline_mu=0.5` and `compare_methods.py` never overrides it, so **every
number in the v1 Results table comes from the one regime where this estimator
is unbiased**. Underestimating ρ makes the test anti-conservative.

Separately, ρ should be **mean-dependent**. Beta-binomial dispersion is
structurally bounded near μ=0 and μ=1 and largest near μ=0.5. A single global ρ
is simultaneously too conservative at the boundaries — where imprinting-style
ASM lives — and too liberal in the middle. DSS, MOABS, edgeR and DESeq2 all fit
a mean–dispersion trend.

## Defect 3 (minor) — ρ absorbs true signal when ASM is common

The estimator pools all loci, assuming ≪1% are truly ASM. Run on all-ASM data
it returns ρ̂ = 0.32, which then annihilates power. The v1 claim of negligible
bias holds at π=0.005 but not if real prevalence is higher — and the parent
project's deep sample showed 9.3% of CpGs at |Δ|≥0.5. Estimate ρ on a held-out
set excluding strong candidates, or iterate.

---

## What was added

### `simulate_reads.py` — molecules, not independent CpG draws

Reads are drawn **once per region** and span every CpG in it. Two parameters
the old simulators lacked:

- `co_methylation` ∈ [0,1] — probability a CpG's state is copied from the
  molecule's latent state. 0 reproduces v1; 1 is perfect co-methylation.
- `baseline_mu="bimodal"` — realistic 70% high / 20% low / 10% intermediate.

Also fixes a quiet bias: v1 **clipped** μ₁, μ₂ to [0.001, 0.999] after adding
δ/2, so at extreme baselines the realised effect silently differed from
`true_delta`. v2 shrinks δ to fit instead.

### `readlevel.py` — the molecule is the unit of observation

Summarise each read to one number (its methylation fraction over the region)
and compare the haplotypes' read-level distributions. Nothing needs to model
co-methylation: aggregating within a read absorbs it, and the between-read
variance driving the test is estimated from the reads themselves. **No
dispersion parameter is required at all.**

- `test_region_reads` — Welch t and Mann-Whitney on per-read fractions, taking
  the more conservative; analytic, so it produces the small p-values genome-wide
  BH needs.
- **Complete separation is handled exactly.** Every read on HP1 above every
  read on HP2 is the strongest possible ASM signal (an imprinted DMR), and a
  t-test cannot express it — both variances are 0. The exact randomisation
  p-value `2/C(n₁+n₂, n₁)` is used instead. An earlier draft floored it at
  `1/(n₁+n₂)` ≈ 0.03, which never survives BH; that bug cost essentially all
  power at high co-methylation.
- `test_region_perm` — exact randomisation, permuting **whole reads** so
  within-molecule structure is preserved under the null. Correctly calibrated
  at every co-methylation level tested (type I 0.042 / 0.033 / 0.009 against
  nominal 0.05). p is floored at 1/(n_perm+1), so it is a **confirmation step
  for screened candidates, not a genome-wide test**.
- `design_effect` — measures how much pooling inflates effective N, from real
  data. Recovers 1.01 / 2.40 / 5.49 at co-methylation 0 / 0.6 / 1.0.

### `dispersion.estimate_dispersion_trend` — mean-dependent, simulation-calibrated

Does not assume `E[χ²]=1`; estimates it by simulating binomial data at the
observed (n₁, n₂, p̂), which reproduces the same boundary degeneracy, then reads
ρ off the excess. Fitted in bins of p̂ and returned as a callable `ρ(μ)` with
`.knots` exposed — **always plot the curve before trusting a run**.

Better than the global estimator but not unbiased: it overshoots at the
boundaries (≈2–3× at μ=0.05/0.95). Overshoot is conservative, so this is a safe
failure mode, but it costs power exactly where imprinting-style ASM lives.
The read-level path avoids the problem entirely by not needing ρ.

**It must be estimated on the same counts it is applied to.** A per-CpG ρ does
not describe pooled region counts; an earlier version of `compare_methods_v2.py`
made that mistake and produced FDR 0.67.

---

## Benchmark v2 results

40,000 regions, π=0.01, BH α=0.05, bimodal methylome, shared reads.

| co-meth (design effect) | method | called | FDR | power |
|---|---|---|---|---|
| **0.0** (DE 1.01) | pooled read×CpG, global ρ | 226 | 0.040 ✓ | 0.565 |
| | pooled read×CpG, ρ(μ) trend | 220 | 0.041 ✓ | 0.549 |
| | read-level | 108 | 0.000 | 0.281 |
| | permutation (type I) | — | 0.042 ✓ | — |
| **0.6** (DE 2.40) | pooled read×CpG, global ρ | 129 | 0.070 | 0.312 |
| | pooled read×CpG, ρ(μ) trend | 64 | 0.031 ✓ | 0.161 |
| | read-level | 7 | 0.000 | 0.018 |
| | permutation (type I) | — | 0.033 ✓ | — |
| **1.0** (DE 5.49) | pooled read×CpG, global ρ | 105 | **0.695** ✗ | 0.083 |
| | pooled read×CpG, ρ(μ) trend | 4 | 0.250 ✗ | 0.008 |
| | read-level | 0 | — | 0.000 |
| | permutation (type I) | — | 0.009 ✓ | — |

Reading these honestly:

- **When CpGs really are independent, pooling is correct and more powerful.**
  v1's approach is not wrong in its own simulation — it is wrong about long-read
  data.
- **The pooled test degrades catastrophically as the design effect rises**
  (FDR 0.04 → 0.07 → 0.70). The ρ(μ) trend delays the failure but does not
  prevent it.
- **The read-level test is conservative at every level and never breaks.** Its
  low power at high co-methylation is *correct*: at DE 5.5 there genuinely are
  ~5× fewer independent observations. The pooled test's apparent power there is
  not power, it is a 70% false discovery rate.
- **The permutation test is correctly calibrated throughout** and is the
  reference standard.

## The number that decides everything

**The real value of the design effect is an empirical question, and it is
measurable.** Run `readlevel.design_effect()` on real `modkit extract` output.

- DE ≈ 1 → pooled region testing is fine; use it, it is more powerful.
- DE ≫ 1 → pooled region testing cannot control FDR; use read-level.

No FDR number from any benchmark here should be quoted in a paper before that
measurement exists. It is also the same quantity as the co-methylation decay
curve in `md/20260903_qc_review.md` §7.3, so one measurement settles both.

## Still open

- ρ(μ) trend overshoots at boundary μ. Fixable with a better within-bin
  estimator, or sidestepped by using the read-level path.
- Read-level test discards within-read *pattern* (it reduces a molecule to a
  fraction). That pattern is the epiallele signal — a distributional test on
  read-level methylation vectors would use it, and would unify this caller with
  the parent project's epiallele analysis rather than merely being compatible
  with it.
- Switch-error attenuation is simulated but still not modelled in any test.
- No real-data validation yet.
