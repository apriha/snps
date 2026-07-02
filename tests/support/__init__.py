"""Shared, committed test support: a fixture-backed resource provider and canned data."""

from tests.support.data import (
    GRCh37_GRCh38,
    GRCh37_GRCh38_PAR,
    GRCh37_NCBI36,
    NCBI36_GRCh37,
    chip_clusters_df,
    get_test_assembly_mapping_data,
    low_quality_snps_df,
    standard_assembly_mappings,
)
from tests.support.fake_resources import FakeResources

__all__ = [
    "FakeResources",
    "GRCh37_GRCh38",
    "GRCh37_GRCh38_PAR",
    "GRCh37_NCBI36",
    "NCBI36_GRCh37",
    "chip_clusters_df",
    "get_test_assembly_mapping_data",
    "low_quality_snps_df",
    "standard_assembly_mappings",
]
