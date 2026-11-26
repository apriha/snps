"""AncestryDNA raw data files.

https://www.ancestry.com
"""

from __future__ import annotations

import numpy as np

from snps.io.readers._registry import Header, ParsedData, ParseInput, ReaderFormat
from snps.io.readers._table import (
    TWO_ALLELE_COLUMNS,
    TWO_ALLELE_DTYPES,
    TableLayout,
    read_table,
)

_LAYOUT = TableLayout(
    columns=TWO_ALLELE_COLUMNS,
    separator=r"\s+",
    na_values=("0",),
    header=0,
    dtypes=TWO_ALLELE_DTYPES,
    alleles="both",
)


def _detect(header: Header) -> bool:
    return "AncestryDNA" in header.first_line


def _parse(parse_input: ParseInput) -> ParsedData:
    df = read_table(parse_input.source, _LAYOUT)

    # https://redd.it/5y90un
    df.iloc[np.where(df["chrom"] == "23")[0], 0] = "X"
    df.iloc[np.where(df["chrom"] == "24")[0], 0] = "Y"
    df.iloc[np.where(df["chrom"] == "25")[0], 0] = "PAR"
    df.iloc[np.where(df["chrom"] == "26")[0], 0] = "MT"

    return ParsedData(df)


FORMAT = ReaderFormat(name="AncestryDNA", detect=_detect, parse=_parse)
