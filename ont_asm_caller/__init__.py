__version__ = "0.2.0"

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
# read-PATTERN path -- CPEL's statistics (Abante et al., Nat Commun 2020) in the
# form long reads make available, plus the cheap non-parametric alternative.
# See docs/2026-09-07_cpel_review_and_pattern_tests.md.
from .pattern import (
    fit_chain, chain_mml, chain_entropy_bits, chain_logZ, chain_marginals,
    jsd_chains, region_statistics, read_fraction_ks,
)
from .null import MatchedNull, calibration_report, split_reads
from .decay import (DecayCurve, comethylation_decay, implied_design_effect,
                    deattenuate)
from .provenance import stamp, write_stamp, git_commit, resolve_region_table_version
from .simulate_patterns import (
    PatternRegion, simulate_pattern_regions,
    molecules_exchangeable, molecules_markov, molecules_epiallele,
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
    # read-pattern path: MML / NME / PDM + the depth-matched, tail-extrapolated
    # null that lets a permutation-style test survive genome-wide BH
    "fit_chain", "chain_mml", "chain_entropy_bits", "chain_logZ",
    "chain_marginals", "jsd_chains", "region_statistics", "read_fraction_ks",
    "MatchedNull", "calibration_report", "split_reads",
    "DecayCurve", "comethylation_decay", "implied_design_effect",
    "deattenuate",
    # provenance: stamp every cluster output with the commit that made it
    "stamp", "write_stamp", "git_commit", "resolve_region_table_version",
    "__version__",
    "PatternRegion", "simulate_pattern_regions",
    "molecules_exchangeable", "molecules_markov", "molecules_epiallele",
]
