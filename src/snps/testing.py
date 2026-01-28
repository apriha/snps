"""Shared test utilities for snps."""

from __future__ import annotations

import os
from typing import Any

import numpy as np
import pandas as pd
from pandas.api.types import is_object_dtype, is_string_dtype

# Standard dtypes for normalized SNP DataFrames
NORMALIZED_DTYPES = {
    "rsid": object,
    "chrom": object,
    "pos": np.uint32,
    "genotype": object,
}


def get_complement(base: str) -> str:
    """Get the complement of a DNA base.

    Parameters
    ----------
    base : str
        A single DNA base (A, C, G, or T)

    Returns
    -------
    str
        The complementary base (A<->T, C<->G), or the original if not a valid base
    """
    complements = {"A": "T", "T": "A", "C": "G", "G": "C"}
    return complements.get(base, base)


def complement_genotype(genotype: str) -> str:
    """Get the complement of a genotype (both alleles).

    Parameters
    ----------
    genotype : str
        A two-character genotype string (e.g., "AT", "CG")

    Returns
    -------
    str
        The complemented genotype, or np.nan if input is null
    """
    if pd.isnull(genotype):
        return np.nan
    return "".join(get_complement(base) for base in genotype)


def complement_one_allele(genotype: str) -> str:
    """Get the complement of only the first allele of a genotype.

    The second allele is preserved unchanged. This is useful for simulating
    partial strand complementation in test data.

    Parameters
    ----------
    genotype : str
        A two-character genotype string (e.g., "AT", "CG")

    Returns
    -------
    str
        Genotype with first allele complemented, or np.nan if input is null
    """
    if pd.isnull(genotype):
        return np.nan
    return get_complement(genotype[0]) + genotype[1]


def create_snp_df(
    rsid: list[str],
    chrom: list[str],
    pos: list[int],
    genotype: list[str],
) -> pd.DataFrame:
    """Create a normalized SNP DataFrame.

    Parameters
    ----------
    rsid : list of str
        SNP identifiers (becomes the index)
    chrom : list of str
        Chromosome values
    pos : list of int
        Position values
    genotype : list of str
        Genotype values

    Returns
    -------
    ~pandas.DataFrame
        DataFrame with rsid index and chrom, pos, genotype columns
    """
    df = pd.DataFrame(
        {"rsid": rsid, "chrom": chrom, "pos": pos, "genotype": genotype},
        columns=["rsid", "chrom", "pos", "genotype"],
    )
    df = df.astype(NORMALIZED_DTYPES)
    df = df.set_index("rsid")
    return df


def create_simulated_snp_df(
    chrom: str = "1",
    pos_start: int = 1,
    pos_max: int = 248140902,
    pos_step: int = 100,
    pos_dtype: type = np.uint32,
    genotype: str = "AA",
    insert_nulls: bool = True,
    null_snp_step: int = 101,
    complement_genotype_one_allele: bool = False,
    complement_genotype_two_alleles: bool = False,
    complement_snp_step: int = 50,
) -> pd.DataFrame:
    """Create a simulated SNP DataFrame for testing.

    This is the core logic for creating simulated SNP data. Each project
    can wrap this to assign to their specific object types.

    Parameters
    ----------
    chrom : str
        Chromosome value for all SNPs (default: "1")
    pos_start : int
        Starting position (default: 1)
    pos_max : int
        Maximum position (default: 248140902)
    pos_step : int
        Step between positions (default: 100)
    pos_dtype : type
        Numpy dtype for positions (default: np.uint32)
    genotype : str
        Default genotype for all SNPs (default: "AA")
    insert_nulls : bool
        Whether to insert null genotypes (default: True)
    null_snp_step : int
        Insert null every N SNPs (default: 101)
    complement_genotype_one_allele : bool
        Complement first allele at intervals (default: False)
    complement_genotype_two_alleles : bool
        Complement both alleles at intervals (default: False)
    complement_snp_step : int
        Apply complement every N SNPs (default: 50)

    Returns
    -------
    ~pandas.DataFrame
        DataFrame with rsid index and chrom, pos, genotype columns
    """
    positions = np.arange(pos_start, pos_max, pos_step, dtype=pos_dtype)
    snps = pd.DataFrame(
        {"chrom": chrom},
        index=pd.Index([f"rs{x + 1}" for x in range(len(positions))], name="rsid"),
    )
    snps["pos"] = positions
    snps["genotype"] = genotype

    if insert_nulls:
        snps.loc[snps.iloc[0::null_snp_step, :].index, "genotype"] = np.nan

    indices = snps.iloc[0::complement_snp_step, :].index
    if complement_genotype_two_alleles:
        snps.loc[indices, "genotype"] = snps.loc[indices, "genotype"].apply(
            complement_genotype
        )
    elif complement_genotype_one_allele:
        snps.loc[indices, "genotype"] = snps.loc[indices, "genotype"].apply(
            complement_one_allele
        )

    return snps


def assert_series_equal_with_string_dtype(
    left: pd.Series,
    right: pd.Series,
    test_case: Any = None,
    **kwargs,
) -> None:
    """Assert Series are equal, accepting both object and StringDtype for string data.

    In Python 3.14+, pandas infers StringDtype for string data instead of object.
    This function compares Series without strict dtype matching for string data.

    Parameters
    ----------
    left : ~pandas.Series
        First Series to compare
    right : ~pandas.Series
        Second Series to compare
    test_case : object, optional
        Object with assertTrue method for assertions (uses assert if None)
    **kwargs : dict
        Additional arguments passed to pd.testing.assert_series_equal
    """
    import pandas as pd

    # Verify string series have string or object dtypes
    if is_string_dtype(left.dtype) or is_object_dtype(left.dtype):
        right_is_string = is_string_dtype(right.dtype) or is_object_dtype(right.dtype)
        if test_case:
            test_case.assertTrue(
                right_is_string,
                f"Right series dtype {right.dtype} should be string/object type",
            )
        else:
            assert right_is_string, (
                f"Right series dtype {right.dtype} should be string/object type"
            )
    # Compare Series without strict dtype matching
    pd.testing.assert_series_equal(left, right, check_dtype=False, **kwargs)


def assert_frame_equal_with_string_index(
    left: pd.DataFrame,
    right: pd.DataFrame,
    test_case: Any = None,
    **kwargs,
) -> None:
    """Assert DataFrames are equal, accepting both object and StringDtype for string columns.

    In Python 3.14+, pandas infers StringDtype for string columns/indices instead of object.
    This function validates that string columns have string types, then compares the
    DataFrames without strict dtype matching for object/string columns.

    Parameters
    ----------
    left : ~pandas.DataFrame
        First DataFrame to compare
    right : ~pandas.DataFrame
        Second DataFrame to compare
    test_case : object, optional
        Object with assertTrue method for assertions (uses assert if None)
    **kwargs : dict
        Additional arguments passed to pd.testing.assert_frame_equal
    """
    import pandas as pd

    def _assert(condition: bool, message: str) -> None:
        if test_case:
            test_case.assertTrue(condition, message)
        else:
            assert condition, message

    # Verify index dtypes are string types if they're named 'rsid'
    if left.index.name == "rsid":
        _assert(
            is_string_dtype(left.index.dtype),
            f"Left index dtype {left.index.dtype} is not a string type",
        )
    if right.index.name == "rsid":
        _assert(
            is_string_dtype(right.index.dtype),
            f"Right index dtype {right.index.dtype} is not a string type",
        )

    # Verify string columns (chrom, genotype) have string dtypes
    for col in ["chrom", "genotype"]:
        if col in left.columns:
            _assert(
                is_string_dtype(left[col].dtype) or is_object_dtype(left[col].dtype),
                f"Left column '{col}' dtype {left[col].dtype} is not a string/object type",
            )
        if col in right.columns:
            _assert(
                is_string_dtype(right[col].dtype) or is_object_dtype(right[col].dtype),
                f"Right column '{col}' dtype {right[col].dtype} is not a string/object type",
            )

    # Compare DataFrames without strict dtype matching for string columns
    pd.testing.assert_frame_equal(
        left, right, check_index_type=False, check_dtype=False, **kwargs
    )


class SNPsTestMixin:
    """Mixin class providing common test assertions and utilities for SNP DataFrames.

    This mixin can be combined with unittest.TestCase to add convenient
    assertion methods for comparing SNP DataFrames with flexible string dtype handling,
    plus common test utilities like creating test DataFrames.

    Example
    -------
    >>> class MyTestCase(SNPsTestMixin, TestCase):
    ...     def test_something(self):
    ...         df = self.generic_snps()
    ...         self.assert_frame_equal_with_string_index(df, expected_df)
    """

    @property
    def downloads_enabled(self) -> bool:
        """Check if external downloads are enabled for tests.

        Only download from external resources when an environment variable named
        "DOWNLOADS_ENABLED" is set to "true".

        Returns
        -------
        bool
        """
        return os.getenv("DOWNLOADS_ENABLED") == "true"

    @staticmethod
    def get_complement(base: str) -> str:
        """Get the complement of a DNA base.

        See :func:`get_complement` for details.
        """
        return get_complement(base)

    def complement_genotype(self, genotype: str) -> str:
        """Get the complement of a genotype (both alleles).

        See :func:`complement_genotype` for details.
        """
        return complement_genotype(genotype)

    def complement_one_allele(self, genotype: str) -> str:
        """Get the complement of only the first allele of a genotype.

        See :func:`complement_one_allele` for details.
        """
        return complement_one_allele(genotype)

    @staticmethod
    def create_snp_df(
        rsid: list[str],
        chrom: list[str],
        pos: list[int],
        genotype: list[str],
    ) -> pd.DataFrame:
        """Create a normalized SNP DataFrame.

        See :func:`create_snp_df` for details.
        """
        return create_snp_df(rsid, chrom, pos, genotype)

    def generic_snps(self) -> pd.DataFrame:
        """Create a generic SNP DataFrame for testing.

        Returns
        -------
        ~pandas.DataFrame
            DataFrame with 8 SNPs (rs1-rs8) on chromosome 1
        """
        return create_snp_df(
            rsid=[f"rs{i}" for i in range(1, 9)],
            chrom=["1"] * 8,
            pos=list(range(101, 109)),
            genotype=["AA", "CC", "GG", "TT", np.nan, "GC", "TC", "AT"],
        )

    def assert_series_equal_with_string_dtype(self, left, right, **kwargs):
        """Assert Series are equal, accepting both object and StringDtype for string data.

        See :func:`assert_series_equal_with_string_dtype` for details.
        """
        assert_series_equal_with_string_dtype(left, right, test_case=self, **kwargs)

    def assert_frame_equal_with_string_index(self, left, right, **kwargs):
        """Assert DataFrames are equal, accepting both object and StringDtype for string columns.

        See :func:`assert_frame_equal_with_string_index` for details.
        """
        assert_frame_equal_with_string_index(left, right, test_case=self, **kwargs)
