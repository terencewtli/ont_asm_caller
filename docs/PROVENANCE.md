# Which commit produced the cluster's ASM loci?

**Currently unrecorded, and it matters.** This note says what is ambiguous, how
to resolve it from the outputs themselves, and what has been added so the
question stops arising.

## The ambiguity

Nothing in this package recorded its own version until `v0.2.0`.
`pyproject.toml` sat at `0.1.0` across six commits, **three of which changed
results**:

| commit | date | effect on output |
|---|---|---|
| `513766e` | 2026-09-05 | region-level pooling introduced (`cluster_cpgs`) |
| `6df55b3` | 2026-09-05 | **`max_span` cap added — changes region boundaries** |
| `ccabe0f` | 2026-09-05 | **region-pooling + dispersion calibration fixed — changes p-values** |

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

`ccabe0f` changed how `rho` is estimated and applied. Unlike region width, that
leaves no structural fingerprint in the region table. If
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
