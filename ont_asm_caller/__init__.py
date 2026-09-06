from .model import ASMResult, test_locus, fit_mu, bb_logpmf
from .dispersion import estimate_dispersion, estimate_dispersion_trend
from .simulate import simulate_loci, SimulatedLocus
from .simulate_mixture import (
    simulate_mixture, MixtureLocus,
    simulate_spatial_mixture, SpatialCpG,
)
from .region import CpG, Region, cluster_cpgs, test_regions
from .simulate_reads import (
    ReadRegion, simulate_read_regions, per_cpg_counts, draw_baseline_mu,
)
from .readlevel import (
    ReadLevelResult, test_region_reads, test_region_perm,
    read_fractions, design_effect,
)

__all__ = [
    # per-CpG / pooled-region path (see docs/2026-09-05_calibration_critique.md
    # before using cluster_cpgs + test_regions: valid only when design_effect ~ 1)
    "ASMResult", "test_locus", "fit_mu", "bb_logpmf",
    "estimate_dispersion", "estimate_dispersion_trend",
    "CpG", "Region", "cluster_cpgs", "test_regions",
    # read-level path (molecule = unit of observation; no dispersion needed)
    "ReadLevelResult", "test_region_reads", "test_region_perm",
    "read_fractions", "design_effect",
    # simulation
    "simulate_loci", "SimulatedLocus",
    "simulate_mixture", "MixtureLocus",
    "simulate_spatial_mixture", "SpatialCpG",
    "ReadRegion", "simulate_read_regions", "per_cpg_counts", "draw_baseline_mu",
]
