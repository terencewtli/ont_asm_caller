# Which commit produced the cluster's ASM loci?

> ### RESOLVED 2026-09-15: the existing production tables are post-cap, and `ccabe0f` cannot have changed them
>
> Checked directly against the cluster outputs and scripts rather than inferred:
>
> - **Geometry: post-cap, proven.** All **140** `tables/dasm_betabinom/*/regions_*.tsv`
>   have a maximum `end − start` of exactly **1000 bp**; none exceeds `max_span`.
>   A pre-cap table has regions up to 335 kb. Independently, the production script
>   `scripts/wg/good_donors/G02_betabinom_region_chr1.py` passes
>   `max_span=MAX_SPAN` to `cluster_cpgs`, a keyword that does not exist before
>   `6df55b3` (`b449543` signature: `max_gap, min_cpgs` only), so on older code it
>   would have raised rather than written a table.
> - **Dispersion: unaffected, whichever commit ran.** `ccabe0f` did **not** change
>   how `rho` is estimated or applied. It *added* `estimate_dispersion_trend`,
>   `readlevel.py` and `simulate_reads.py`. `git diff 6df55b3 55cce16` shows zero
>   changes to `region.py` or `model.py` and zero removed lines in `dispersion.py`.
>   `G02` calls the original global `estimate_dispersion` and `test_regions`, so
>   p-values are identical before and after `ccabe0f`.
> - **Timing, consistent with both:** cap committed 09-05 21:15, `G02` edited to
>   pass `max_span` 21:16, the first table (HG00146 chr1) written 21:33, dispersion
>   commit 21:41, and every other table 09-06 15:15–17:14.
>
> **No re-run is needed for provenance.** Two things that remain true and matter
> more:
>
> 1. **Production uses the global dispersion estimator**, the one
>    `docs/2026-09-05_calibration_critique.md` measured as 2–3.5× too low on a
>    bimodal methylome. The trend estimator was never wired into `G02`. That is a
>    calibration question, not a provenance one.
> 2. **`G02_betabinom_region_chr1.py` is not under version control** in either
>    repo (only `G03`/`G04` `.sh` are tracked in `asm_lr`). The production
>    calling logic lives only on the cluster. Add `write_stamp()` to it before the
>    dipcall re-call, which will regenerate every table.
>
> The original analysis follows unchanged, except for the dispersion paragraph,
> corrected below.

**Unrecorded until 2026-09-15; see the resolution above.** This note says what
was ambiguous, how to resolve it from the outputs themselves, and what has been
added so the question stops arising.

## The ambiguity

Nothing in this package recorded its own version until `v0.2.0`.
`pyproject.toml` sat at `0.1.0` across six commits, **three of which changed
results**:

| commit | date | effect on output |
|---|---|---|
| `513766e` | 2026-09-05 | region-level pooling introduced (`cluster_cpgs`) |
| `6df55b3` | 2026-09-05 | **`max_span` cap added — changes region boundaries** |
| `ccabe0f` | 2026-09-05 | ~~changes p-values~~ adds `estimate_dispersion_trend` + read-level path; **production path unchanged** (verified 2026-09-15) |

All three landed on the same day. The production caller was run on the cluster
somewhere in that window, and
`tables/dasm_betabinom/{sample}/regions_{chrom}.tsv` carries no record of which
commit produced it.

**This is not hypothetical.** The commit message for `6df55b3` describes the
uncapped behaviour as *"confirmed on real chr1 data: median region size was
already 20 CpGs … the tail reached 21,163 CpGs spanning 335 kb"* — it is
describing production output that already existed. If those tables were never
regenerated after the cap landed, they contain daisy-chained regions.

That propagates. `P05_design_effect_real_data.py` deliberately reuses *"the SAME
region boundaries the production caller already computed"*, and the design
effect grows with the number of CpGs pooled — `DE = 1 + (C−1)·r̄`. Measured on
HG00146 chr15, DE goes 1.29 → 1.48 → 1.89 as median CpGs per region goes
13 → 26 → 52. A 21,163-CpG region is off that scale entirely. So if the
geometry is wrong, every DE and FDR statement computed on it is about the wrong
regions.

## Resolving it for the tables that already exist

One command, on the cluster:

```bash
python3 -c "
from ont_asm_caller import resolve_region_table_version as r
import json; print(json.dumps(r('tables/dasm_betabinom/HG00146/regions_chr1.tsv'), indent=2))"
```

A region wider than `max_span` **cannot** have been produced by the capped
`cluster_cpgs`, so a single such row is proof the table predates `6df55b3`. The
converse is weaker and the function says so: a table with no wide region is
*consistent with* the later version, not proof of it.

- verdict `pre-cap` → **regenerate the tables**, and treat every downstream
  number computed on them (including the `P05` design effects in the handoff
  note) as measured on the wrong geometry.
- verdict `consistent-with-post-cap` → the geometry is fine; the dispersion
  question below is still open.

### The dispersion half is harder

*(Corrected 2026-09-15: `ccabe0f` added a new estimator but did not change the
one production calls, so for the existing tables this concern is moot. The general
point stands for any future change to `estimate_dispersion` itself.)* A change to
how `rho` is estimated or applied leaves no structural fingerprint in the region
table. If
`summary_{chrom}.json` records `rho_hat`, comparing it against a re-run on the
same inputs will settle it; if it does not, **the honest answer is that it is
not recoverable and the run should be repeated**. Do not reason backwards from
whether the numbers look plausible.

## So that this stops arising

`ont_asm_caller.write_stamp()` drops a `<output>.provenance.json` beside any
file, recording the commit, dirty-tree flag, branch, package version, Python /
NumPy / SciPy versions, `argv`, and the analysis parameters actually used:

```python
from ont_asm_caller import write_stamp
regions.to_csv(path, sep="\t", index=False)
write_stamp(path, max_gap=500, max_span=1000, rho=float(rho_hat),
            fdr=0.05, delta_floor=0.10, input_bam=bam_path)
```

Two deliberate choices. The stamp sits **next to the file it describes**, so it
survives being copied around, where a line in a log does not. And when git is
unavailable or the package was installed from a wheel, the commit is recorded as
`null` rather than omitted — an unknown provenance must be visible as unknown.

The commit alone does not pin a run: pass the actual parameters too, since the
same commit invoked with a different `max_span` or `rho` produces different
tables.

## Two-session working agreement

Changes to this repo are made **only** from the local session; the cluster
session consumes it read-only. That keeps history linear, but it means the
cluster can be running an older checkout than the one the local session is
reasoning about. Two habits close that gap:

1. On the cluster, record `git -C <repo> rev-parse HEAD` in any job script's log
   *before* the run starts, and `git pull` deliberately rather than incidentally.
2. Any result quoted back to the local session should carry its stamp. A number
   without a commit cannot be checked against the code that is being edited here.
