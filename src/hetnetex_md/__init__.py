"""
HetNetEX-MD
===========
Exact and asymptotic inference for Multi-DWPC on heterogeneous knowledge graphs.

Two entry points cover the two null models Multi-DWPC uses:

    hetnetex_md_resampling(...)   # N1: gene-set resampling  -> EXACT
    hetnetex_md_network(...)      # N2: network permutation  -> asymptotic

plus an exact null for the median aggregate:

    exact_median_pvalue(...)

See README.md for a worked example and the mapping from functions to the
theorems they implement.
"""

from .core import (  # noqa: F401
    # data structures
    EdgeType,
    HetNet,
    # observed scores
    dwpc_matrix,
    # N1: exact resampling null
    exact_resampling_moments,
    edgeworth_upper_tail,
    hetnetex_md_resampling,
    exact_median_pvalue,
    # N2: asymptotic network null
    network_null_moments,
    aggregate_network_null,
    hetnetex_md_network,
    correction_factors,
    fit_soft_cm,
    soft_cm_ratio,
    admissible_self_kernels,
)

__version__ = "0.1.0"

__all__ = [
    "EdgeType",
    "HetNet",
    "dwpc_matrix",
    "exact_resampling_moments",
    "edgeworth_upper_tail",
    "hetnetex_md_resampling",
    "exact_median_pvalue",
    "network_null_moments",
    "aggregate_network_null",
    "hetnetex_md_network",
    "correction_factors",
    "fit_soft_cm",
    "soft_cm_ratio",
    "admissible_self_kernels",
    "__version__",
]
