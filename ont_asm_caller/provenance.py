"""
Record which version of this code produced a given output.

WHY THIS EXISTS
---------------
Nothing in this package recorded its own version until now, and
`pyproject.toml` sat at 0.1.0 across six commits, three of which changed
results:

  513766e  region-level pooling introduced (cluster_cpgs)
  6df55b3  max_span cap added -- CHANGES REGION BOUNDARIES
  ccabe0f  region-pooling and dispersion calibration fixed -- CHANGES p-values

All three landed on 2026-09-05. The production caller was run on the cluster
somewhere in that window and its outputs
(`tables/dasm_betabinom/{sample}/regions_{chrom}.tsv`) carry no record of which
commit produced them. Everything downstream inherits that ambiguity, including
the design-effect measurement in `P05`, which deliberately reuses "the SAME
region boundaries the production caller already computed".

This is not hypothetical. The commit message for 6df55b3 describes the
uncapped behaviour as "confirmed on real chr1 data: median region size was
already 20 CpGs ... the tail reached 21,163 CpGs spanning 335kb" -- i.e. it is
describing production output that already existed. If those tables were never
regenerated after the cap landed, they contain daisy-chained regions, and since
the design effect grows with the number of CpGs pooled (DE = 1 + (C-1)*rbar),
any DE or FDR statement computed on them is about the wrong geometry.

HOW TO RESOLVE IT FOR AN EXISTING OUTPUT
----------------------------------------
`resolve_region_table_version()` below answers it from the table itself: a
region wider than `max_span` or holding hundreds of CpGs can only have been
produced before 6df55b3.

HOW TO AVOID IT IN FUTURE
-------------------------
Call `write_stamp()` next to every output a cluster run produces.
"""
from __future__ import annotations

import json
import os
import subprocess
import sys
from datetime import datetime, timezone

_PKG_DIR = os.path.dirname(os.path.abspath(__file__))


def git_commit(repo_dir: str | None = None) -> dict:
    """Commit of the checkout this module was imported from, plus whether the
    working tree was dirty. Returns `{"commit": None, ...}` rather than raising
    when git is unavailable or the package was installed from a wheel -- an
    unknown provenance must be recorded as unknown, never omitted."""
    d = repo_dir or _PKG_DIR
    out = {"commit": None, "dirty": None, "branch": None, "repo": d}
    try:
        def _git(*args):
            return subprocess.run(("git", "-C", d) + args, capture_output=True,
                                  text=True, timeout=10).stdout.strip()
        commit = _git("rev-parse", "HEAD")
        if not commit:
            return out
        out["commit"] = commit
        out["branch"] = _git("rev-parse", "--abbrev-ref", "HEAD") or None
        out["dirty"] = bool(_git("status", "--porcelain"))
    except Exception:
        pass
    return out


def stamp(**params) -> dict:
    """Everything needed to reproduce a run. `params` should carry the actual
    analysis parameters (rho, max_gap, max_span, thresholds, input paths), since
    the commit alone does not pin what the caller was invoked with."""
    try:
        from . import __version__ as ver
    except Exception:
        ver = None
    if ver is None:
        try:
            from importlib.metadata import version
            ver = version("ont-asm-caller")
        except Exception:
            ver = "unknown"
    import numpy, scipy
    return {
        "package": "ont_asm_caller",
        "version": ver,
        "git": git_commit(),
        "generated_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "python": sys.version.split()[0],
        "numpy": numpy.__version__,
        "scipy": scipy.__version__,
        "argv": sys.argv,
        "params": params,
    }


def write_stamp(output_path: str, **params) -> str:
    """Write `<output_path>.provenance.json` beside an output file.

    Call this from every cluster script that writes a table. A stamp that sits
    next to the file it describes survives being copied around; a line in a log
    does not."""
    path = str(output_path) + ".provenance.json"
    with open(path, "w") as fh:
        json.dump(stamp(**params), fh, indent=2)
    return path


def resolve_region_table_version(regions_tsv: str, max_span: int = 1000) -> dict:
    """Decide, from an existing region table alone, whether it predates the
    `max_span` cap (6df55b3).

    A region wider than `max_span` cannot have been produced by the capped
    `cluster_cpgs`, so a single such row is proof the table is older. The
    converse is weaker -- a table with no wide region is *consistent* with
    either version -- so the verdict is reported as one of
    "pre-cap" / "consistent-with-post-cap" / "unknown".
    """
    import csv
    n = 0
    max_cpgs = 0
    max_width = 0
    with open(regions_tsv) as fh:
        rd = csv.DictReader(fh, delimiter="\t")
        cols = rd.fieldnames or []
        if not {"start", "end"} <= set(cols):
            return {"verdict": "unknown", "reason": f"columns are {cols}"}
        for row in rd:
            n += 1
            w = int(row["end"]) - int(row["start"])
            max_width = max(max_width, w)
            if "n_cpgs" in row and row["n_cpgs"]:
                max_cpgs = max(max_cpgs, int(row["n_cpgs"]))
    verdict = ("pre-cap" if max_width > max_span
               else "consistent-with-post-cap")
    return {
        "verdict": verdict,
        "n_regions": n,
        "max_region_width_bp": max_width,
        "max_cpgs_in_a_region": max_cpgs or None,
        "max_span_tested": max_span,
        "note": ("A region wider than max_span proves the table predates "
                 "6df55b3. Absence of one is consistent with either version, "
                 "not proof of the later one."),
    }
