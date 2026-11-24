"""AncestryDNA raw data files."""

from __future__ import annotations

import numpy as np
import pandas as pd

from snps.io.readers._common import TWO_ALLELE_DTYPES
from snps.io.readers._registry import Header, ParsedData, ParseInput, ReaderFormat


def _detect(header: Header) -> bool:
    return "AncestryDNA" in header.first_line


def _parse(parse_input: ParseInput) -> ParsedData:
    """Read and parse Ancestry.com file.

    http://www.ancestry.com

    Parameters
    ----------
    parse_input : ParseInput
        file to parse

    Returns
    -------
    ParsedData
        parsed SNPs
    """
    with parse_input.source.open() as file:
        df = pd.read_csv(  # type: ignore[call-overload]  # pandas-stubs rejects a scalar na_values
            file,
            comment="#",
            header=0,
            engine="c",
            sep=r"\s+",
            # delim_whitespace=True,  # https://stackoverflow.com/a/15026839
            na_values=0,
            names=["rsid", "chrom", "pos", "allele1", "allele2"],
            index_col=0,
            dtype=TWO_ALLELE_DTYPES,
        )

    # create genotype column from allele columns
    df["genotype"] = df["allele1"] + df["allele2"]

    # delete allele columns
    # http://stackoverflow.com/a/13485766
    del df["allele1"]
    del df["allele2"]

    # https://redd.it/5y90un
    df.iloc[np.where(df["chrom"] == "23")[0], 0] = "X"
    df.iloc[np.where(df["chrom"] == "24")[0], 0] = "Y"
    df.iloc[np.where(df["chrom"] == "25")[0], 0] = "PAR"
    df.iloc[np.where(df["chrom"] == "26")[0], 0] = "MT"

    return ParsedData(df)


FORMAT = ReaderFormat(name="AncestryDNA", detect=_detect, parse=_parse)
