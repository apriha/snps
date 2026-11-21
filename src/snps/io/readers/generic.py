"""Generic CSV or TSV files."""

from __future__ import annotations

import io
import re

import pandas as pd

from snps.constants import NORMALIZED_DTYPES
from snps.io.readers._registry import Header, ParsedData, ParseInput, ReaderFormat

_HEADER = "^#*[ \t]*rsid[, \t]*chr"
_NO_HEADER = "^rs[0-9]*[, \t]{1}[1]"


def _detect(header: Header) -> bool:
    return (
        re.match(_HEADER, header.first_line) is not None
        or re.match(_NO_HEADER, header.first_line) is not None
    )


def _parse(parse_input: ParseInput) -> ParsedData:
    """Read and parse generic CSV or TSV file.

    Notes
    -----
    Assumes columns are 'rsid', 'chrom' / 'chromosome', 'pos' / 'position', and 'genotype';
    values are comma separated; unreported genotypes are indicated by '--'; and one header row
    precedes data. For example:

        rsid,chromosome,position,genotype
        rs1,1,1,AA
        rs2,1,2,CC
        rs3,1,3,--

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
    skip = 1 if re.match(_HEADER, parse_input.header.first_line) else 0

    def parse(sep: str) -> pd.DataFrame:
        return pd.read_csv(
            file,
            sep=sep,
            skiprows=skip,
            na_values="--",
            names=["rsid", "chrom", "pos", "genotype"],
            index_col=0,
            dtype=NORMALIZED_DTYPES,
            compression=compression,
        )

    try:
        df = parse(",")
    except ValueError:
        try:
            if isinstance(file, io.BufferedIOBase):
                file.seek(0)

            df = parse("\t")
        except ValueError:
            if isinstance(file, io.BufferedIOBase):
                file.seek(0)

            df = pd.read_csv(
                file,
                sep=None,
                na_values="--",
                skiprows=skip,
                engine="python",
                names=["rsid", "chrom", "pos", "genotype"],
                usecols=[0, 1, 2, 3],
                index_col=0,
                dtype=NORMALIZED_DTYPES,
                compression=compression,
            )
    return ParsedData(df)


FORMAT = ReaderFormat(name="generic", detect=_detect, parse=_parse)
