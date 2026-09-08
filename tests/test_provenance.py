"""Provenance stamping, and the diagnostic that recovers it after the fact.

The package ran on the cluster for six commits -- three of which changed
results -- while `pyproject.toml` sat at 0.1.0 and nothing recorded a commit.
These tests cover both halves of the fix: stamping future outputs, and deciding
from an existing region table which side of the `max_span` cap it came from."""
import json

import pytest

from ont_asm_caller import provenance
from ont_asm_caller.region import CpG, cluster_cpgs


def test_stamp_records_version_git_and_parameters():
    s = provenance.stamp(max_gap=500, max_span=1000, rho=0.0123)
    assert s["package"] == "ont_asm_caller"
    assert s["version"] and s["version"] != "unknown"
    assert s["params"]["max_span"] == 1000
    for key in ("git", "generated_utc", "python", "numpy", "scipy", "argv"):
        assert key in s
    # unknown provenance must be RECORDED as unknown, never silently omitted
    assert set(s["git"]) >= {"commit", "dirty", "branch"}


def test_git_commit_never_raises_outside_a_repo(tmp_path):
    out = provenance.git_commit(str(tmp_path))
    assert out["commit"] is None          # not an exception, not a guess


def test_write_stamp_lands_beside_the_output(tmp_path):
    out = tmp_path / "regions_chr1.tsv"
    out.write_text("chrom\tstart\tend\n")
    path = provenance.write_stamp(str(out), max_span=1000)
    assert path == str(out) + ".provenance.json"
    assert json.loads(open(path).read())["params"]["max_span"] == 1000


def _write_regions(path, spans):
    with open(path, "w") as fh:
        fh.write("chrom\tstart\tend\tn_cpgs\n")
        for i, (s, w, n) in enumerate(spans):
            fh.write(f"chr1\t{s}\t{s+w}\t{n}\n")


def test_wide_region_proves_the_table_predates_the_cap(tmp_path):
    """6df75b3 capped region width at max_span. A region wider than that cannot
    have come from the capped code, so one such row is proof."""
    p = tmp_path / "old.tsv"
    _write_regions(p, [(1000, 400, 8), (5000, 335_000, 21_163)])
    r = provenance.resolve_region_table_version(str(p), max_span=1000)
    assert r["verdict"] == "pre-cap"
    assert r["max_region_width_bp"] == 335_000
    assert r["max_cpgs_in_a_region"] == 21_163


def test_narrow_table_is_only_CONSISTENT_with_the_capped_version(tmp_path):
    """The converse is weaker and the verdict must say so: a table with no wide
    region could have come from either version."""
    p = tmp_path / "new.tsv"
    _write_regions(p, [(1000, 400, 8), (5000, 900, 20)])
    r = provenance.resolve_region_table_version(str(p), max_span=1000)
    assert r["verdict"] == "consistent-with-post-cap"
    assert "not proof" in r["note"]


def test_diagnostic_agrees_with_what_cluster_cpgs_actually_produces():
    """Tie the diagnostic to the real code path rather than to a hand-made
    fixture: regions from the capped clusterer must never trip the pre-cap
    verdict."""
    cpgs = [CpG(pos=p, x1=5, n1=10, x2=5, n2=10) for p in range(1000, 6000, 40)]
    regions = cluster_cpgs("chr1", cpgs, max_gap=500, max_span=1000)
    assert max(r.end - r.start for r in regions) <= 1000
