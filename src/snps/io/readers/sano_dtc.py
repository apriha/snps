"""Sano Genetics direct-to-consumer raw data files."""

from __future__ import annotations

import re

import pandas as pd

from snps.io.readers._common import TWO_ALLELE_DTYPES
from snps.io.readers._registry import Header, ParsedData, ParseInput, ReaderFormat


def _detect(header: Header) -> bool:
    return (
        "# This file was derived from the corresponding VCF" in header.comments
        or re.match(
            r"^\s*rsid\s+chromosome\s+position\s+allele_1\s+allele_2\s*$",
            header.first_line,
        )
        is not None
    )


def _source(header: Header) -> str:
    return "Sano"


def _parse(parse_input: ParseInput) -> ParsedData:
    """Read and parse Sano Genetics DTC file.

    https://sanogenetics.com

    Parameters
    ----------
    parse_input : ParseInput
        file to parse

    Returns
    -------
    ParsedData
        parsed SNPs
    """
    file = parse_input.file
    compression = parse_input.compression

    df = pd.read_csv(
        file,
        comment="#",
        header=0,
        engine="c",
        sep=r"\s+",
        na_values="-",
        names=["rsid", "chrom", "pos", "allele1", "allele2"],
        index_col=0,
        dtype=TWO_ALLELE_DTYPES,
        compression=compression,
    )

    # create genotype column from allele columns and keep only relevant columns
    df = df.assign(genotype=df["allele1"] + df["allele2"].fillna(""))[
        ["chrom", "pos", "genotype"]
    ]

    return ParsedData(df)


FORMAT = ReaderFormat(name="sano_dtc", detect=_detect, parse=_parse, source=_source)
