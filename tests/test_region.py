"""
Tests for distance-based CpG clustering and region-level testing: the
clustering logic itself (does it merge/split at the right boundaries, pool
counts correctly), and -- the actual point of building this -- whether
pooling recovers power that the per-CpG test loses to multiple testing.

Note: `cluster_cpgs`/`test_regions`/`test_locus` are imported under aliases
below purely to avoid pytest auto-collecting `test_locus`/`test_regions` as
test items themselves (see tests/test_model.py for the same issue).
"""
import numpy as np

from ont_asm_caller.region import CpG
from ont_asm_caller.region import cluster_cpgs as cluster
from ont_asm_caller.region import test_regions as run_region_tests
from ont_asm_caller.model import test_locus as run_asm_test
from ont_asm_caller.dispersion import estimate_dispersion
from ont_asm_caller.simulate_mixture import simulate_spatial_mixture


def test_clusters_nearby_cpgs_together():
    cpgs = [
        CpG(pos=100, x1=5, n1=10, x2=5, n2=10),
        CpG(pos=200, x1=5, n1=10, x2=5, n2=10),
        CpG(pos=350, x1=5, n1=10, x2=5, n2=10),
    ]
    regions = cluster("chr1", cpgs, max_gap=500)
    assert len(regions) == 1
    r = regions[0]
    assert r.n_cpgs == 3
    assert r.x1 == 15 and r.n1 == 30
    assert r.x2 == 15 and r.n2 == 30
    assert r.start == 100 and r.end == 350


def test_splits_distant_cpgs_into_separate_regions():
    cpgs = [
        CpG(pos=100, x1=1, n1=10, x2=1, n2=10),
        CpG(pos=200, x1=1, n1=10, x2=1, n2=10),
        CpG(pos=5000, x1=1, n1=10, x2=1, n2=10),  # far from the first two
    ]
    regions = cluster("chr1", cpgs, max_gap=500)
    assert len(regions) == 2
    assert regions[0].n_cpgs == 2
    assert regions[1].n_cpgs == 1


def test_min_cpgs_filter():
    cpgs = [
        CpG(pos=100, x1=1, n1=10, x2=1, n2=10),
        CpG(pos=5000, x1=1, n1=10, x2=1, n2=10),
        CpG(pos=5100, x1=1, n1=10, x2=1, n2=10),
        CpG(pos=5200, x1=1, n1=10, x2=1, n2=10),
    ]
    regions = cluster("chr1", cpgs, max_gap=500, min_cpgs=3)
    assert len(regions) == 1
    assert regions[0].n_cpgs == 3


def test_handles_unsorted_input():
    cpgs = [
        CpG(pos=350, x1=5, n1=10, x2=5, n2=10),
        CpG(pos=100, x1=5, n1=10, x2=5, n2=10),
        CpG(pos=200, x1=5, n1=10, x2=5, n2=10),
    ]
    regions = cluster("chr1", cpgs, max_gap=500)
    assert len(regions) == 1
    assert regions[0].start == 100 and regions[0].end == 350


def test_empty_input():
    assert cluster("chr1", [], max_gap=500) == []


def test_pooled_dispersion_is_lower_than_per_cpg_dispersion():
    """The key statistical fact this module relies on: pooling several
    independent per-CpG beta-binomial draws into one combined count reduces
    the EFFECTIVE dispersion of that combined count (averaging out
    site-specific noise) -- it does NOT leave the same rho valid at the
    pooled level. If this weren't true, reusing the per-CpG rho for pooled
    counts would be fine; it isn't, which is why region.py's docstring
    warns against it and callers must re-estimate rho on pooled data."""
    rho_true = 0.05
    cpgs = simulate_spatial_mixture(
        n_regions=4000, pi=0.0, cpgs_per_region=(6, 6),  # fixed 6 CpGs/region, all null
        depth=(15, 15), rho=rho_true, seed=5,
    )
    x1 = np.array([c.x1 for c in cpgs]); n1 = np.array([c.n1 for c in cpgs])
    x2 = np.array([c.x2 for c in cpgs]); n2 = np.array([c.n2 for c in cpgs])
    rho_percpg = estimate_dispersion(x1, n1, x2, n2)

    region_cpgs = [CpG(pos=c.pos, x1=c.x1, n1=c.n1, x2=c.x2, n2=c.n2) for c in cpgs]
    regions = cluster("chr1", region_cpgs, max_gap=400, min_cpgs=1)
    rx1 = np.array([r.x1 for r in regions]); rn1 = np.array([r.n1 for r in regions])
    rx2 = np.array([r.x2 for r in regions]); rn2 = np.array([r.n2 for r in regions])
    rho_region = estimate_dispersion(rx1, rn1, rx2, rn2)

    assert rho_region < rho_percpg, (
        f"expected pooled dispersion < per-CpG dispersion, "
        f"got rho_percpg={rho_percpg:.4f} rho_region={rho_region:.4f}"
    )


def test_region_pooling_with_correctly_reestimated_dispersion_improves_power():
    """The actual point of this module, done correctly: cluster CpGs into
    regions, estimate dispersion SEPARATELY at the region level (not reused
    from the per-CpG estimate -- see the test above for why that would be
    wrong), then compare region-level BH-FDR power against per-CpG BH-FDR
    power on the same underlying data."""
    rho_true = 0.05
    cpgs = simulate_spatial_mixture(
        n_regions=4000, pi=0.02, cpgs_per_region=(4, 8),
        depth=(8, 20), rho=rho_true, effect_size_dist="uniform_0.2_0.5",
        seed=11,
    )

    # per-CpG BH-FDR, using the correct per-CpG dispersion estimate
    x1 = np.array([c.x1 for c in cpgs]); n1 = np.array([c.n1 for c in cpgs])
    x2 = np.array([c.x2 for c in cpgs]); n2 = np.array([c.n2 for c in cpgs])
    rho_cpg_hat = estimate_dispersion(x1, n1, x2, n2)

    pvals_cpg = np.array([
        run_asm_test("chr1", c.pos, c.x1, c.n1, c.x2, c.n2, rho=rho_cpg_hat).pval
        for c in cpgs
    ])
    is_asm_cpg = np.array([c.is_true_asm for c in cpgs])
    n = len(pvals_cpg)
    order = np.argsort(pvals_cpg)
    q = np.minimum.accumulate((pvals_cpg[order] * n / (np.arange(n) + 1))[::-1])[::-1]
    q_cpg = np.empty(n); q_cpg[order] = np.clip(q, 0, 1)
    power_cpg = ((q_cpg < 0.05) & is_asm_cpg).sum() / is_asm_cpg.sum()

    # region-level BH-FDR, with dispersion RE-estimated on pooled counts
    region_cpgs = [CpG(pos=c.pos, x1=c.x1, n1=c.n1, x2=c.x2, n2=c.n2) for c in cpgs]
    regions = cluster("chr1", region_cpgs, max_gap=400, min_cpgs=1)
    rx1 = np.array([r.x1 for r in regions]); rn1 = np.array([r.n1 for r in regions])
    rx2 = np.array([r.x2 for r in regions]); rn2 = np.array([r.n2 for r in regions])
    rho_region_hat = estimate_dispersion(rx1, rn1, rx2, rn2)

    results = run_region_tests(regions, rho=rho_region_hat)
    pvals_region = np.array([r.pval for r in results])

    pos_to_asm = {c.pos: c.is_true_asm for c in cpgs}
    region_truth = np.array([
        any(pos_to_asm.get(p, False) for p in range(reg.start, reg.end + 1))
        for reg in regions
    ])

    nr = len(pvals_region)
    order_r = np.argsort(pvals_region)
    qr = np.minimum.accumulate((pvals_region[order_r] * nr / (np.arange(nr) + 1))[::-1])[::-1]
    q_region = np.empty(nr); q_region[order_r] = np.clip(qr, 0, 1)
    power_region = ((q_region < 0.05) & region_truth).sum() / region_truth.sum()

    assert power_region > power_cpg, (
        f"region pooling (with correctly re-estimated dispersion) did not "
        f"improve power: per-CpG={power_cpg:.3f}, region={power_region:.3f}"
    )
