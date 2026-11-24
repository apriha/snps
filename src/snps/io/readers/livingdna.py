"""LivingDNA raw data files."""

from __future__ import annotations

import pandas as pd

from snps.constants import NORMALIZED_DTYPES
from snps.io.readers._registry import Header, ParsedData, ParseInput, ReaderFormat


def _detect(header: Header) -> bool:
    return "Living DNA" in header.first_line


def _parse(parse_input: ParseInput) -> ParsedData:
    """Read and parse LivingDNA file.

    https://livingdna.com/

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
        return ParsedData(
            pd.read_csv(
                file,
                comment="#",
                sep="\t",
                na_values="--",
                names=["rsid", "chrom", "pos", "genotype"],
                index_col=0,
                dtype=NORMALIZED_DTYPES,
            )
        )


FORMAT = ReaderFormat(name="LivingDNA", detect=_detect, parse=_parse)
