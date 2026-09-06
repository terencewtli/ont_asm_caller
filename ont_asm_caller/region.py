"""
Region-level pooling: cluster neighboring CpGs by genomic distance and pool
their read counts into one effective (X, N) pair per haplotype, then run the
same beta-binomial test on the pooled counts.

Why this fixes the power problem found in benchmarks/compare_methods.py:
that benchmark showed the per-CpG test, while correctly calibrated, has
near-zero power once BH-FDR is applied genome-wide at a realistic (~0.5%)
true-ASM prevalence -- with millions of tests and few true positives, a true
positive's raw p-value must be extremely small to survive correction, and
per-CpG depth/dispersion often can't produce p-values that small.

Region pooling attacks BOTH halves of that problem, not just one:
1. Fewer tests (regions, not individual CpGs) directly loosens the BH
   penalty -- the (k/n)*alpha threshold is less punishing when n is smaller.
2. Pooling reads across multiple CpGs that share a true rate genuinely
   reduces per-test noise, not just the correction: each CpG's beta-binomial
   draw carries its own site-specific extra-binomial noise (rho), and
   combining several roughly-independent draws averages some of that out --
   this is a real power gain, not merely fewer chances to be penalized. It
   relies on the (biologically reasonable, for a real ASM-driving locus)
   assumption that CpGs within ~500bp of each other are governed by the same
   underlying haplotype-specific rate, which is why clustering is
   distance-based (matching this project's existing GAP_MERGE_BP=500 and
   DSS's own smoothing.span=500) rather than an arbitrary fixed-width
   sliding window -- a real DMR can span a window boundary or be shorter
   than a fixed tile, and overlapping tiles reintroduce their own multiple-
   testing double-counting problem that distance-based clustering avoids.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from .model import ASMResult, test_locus


@dataclass
class CpG:
    pos: int
    x1: int
    n1: int
    x2: int
    n2: int


@dataclass
class Region:
    chrom: str
    start: int
    end: int
    n_cpgs: int
    x1: int
    n1: int
    x2: int
    n2: int


def cluster_cpgs(chrom: str, cpgs: list[CpG], max_gap: int = 500,
                  min_cpgs: int = 1, max_span: int = 1000) -> list[Region]:
    """Merge CpGs into regions wherever consecutive positions (after sorting)
    are within `max_gap` bp of each other. `cpgs` need not be pre-sorted.
    `min_cpgs` filters out regions smaller than this after clustering (set
    to e.g. 3 to match DSS callDMR's default minCG, or 1 to keep singletons
    as their own one-CpG "regions" -- equivalent to the plain per-CpG test).

    `max_span` caps each region's total width (region.end - region.start).
    Without this, a simple adjacent-gap rule can DAISY-CHAIN through any
    CpG-dense stretch (e.g. a long CpG island, or just a gene-rich region)
    into one absurdly large "region" -- confirmed on real chr1 data: with no
    cap, median region size was already 20 CpGs (bigger than the 4-8 the
    benchmark in benchmarks/compare_methods.py validated), and the tail
    reached 21,163 CpGs spanning 335kb in a single region. That's not a
    coherent regulatory locus sharing one true rate -- it's an artifact of
    the merge rule having no sense of overall scale. `max_span=1000`
    (matching a plausible single-regulatory-element/CpG-island scale) forces
    a split once a region would exceed it, even if every individual gap
    within it is `<= max_gap`.
    """
    if not cpgs:
        return []

    ordered = sorted(cpgs, key=lambda c: c.pos)
    regions = []
    cur = [ordered[0]]

    def _flush(group: list[CpG]) -> Region:
        return Region(
            chrom=chrom, start=group[0].pos, end=group[-1].pos,
            n_cpgs=len(group),
            x1=sum(c.x1 for c in group), n1=sum(c.n1 for c in group),
            x2=sum(c.x2 for c in group), n2=sum(c.n2 for c in group),
        )

    for c in ordered[1:]:
        would_span = c.pos - cur[0].pos
        if c.pos - cur[-1].pos <= max_gap and would_span <= max_span:
            cur.append(c)
        else:
            regions.append(_flush(cur))
            cur = [c]
    regions.append(_flush(cur))

    return [r for r in regions if r.n_cpgs >= min_cpgs]


def test_regions(regions: list[Region], rho: float) -> list[ASMResult]:
    """Run the beta-binomial LRT on each region's pooled counts. Uses
    `region.start` as the reported position (the ASMResult.pos field) --
    callers that need the full span should keep the Region objects around
    alongside the results (e.g. `dict(zip(regions, results))`)."""
    return [
        test_locus(r.chrom, r.start, r.x1, r.n1, r.x2, r.n2, rho)
        for r in regions
    ]
