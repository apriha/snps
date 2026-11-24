"""Family Tree DNA (FTDNA) raw data files, including "famfinder" files."""

from __future__ import annotations

import pandas as pd

from snps.constants import NORMALIZED_DTYPES
from snps.exceptions import ParseError
from snps.io.readers._common import TWO_ALLELE_DTYPES, extract_comments
from snps.io.readers._registry import Header, ParsedData, ParseInput, ReaderFormat


def _detect_ftdna(header: Header) -> bool:
    return header.first_line.startswith("RSID")


def _detect_famfinder(header: Header) -> bool:
    return "famfinder" in header.first_line


def _source_famfinder(header: Header) -> str:
    return "FTDNA"


def _parse_ftdna(parse_input: ParseInput) -> ParsedData:
    """Read and parse Family Tree DNA (FTDNA) file.

    https://www.familytreedna.com

    Parameters
    ----------
    parse_input : ParseInput
        file to parse

    Returns
    -------
    ParsedData
        parsed SNPs
    """
    source = parse_input.source

    try:
        with source.open() as file:
            df = pd.read_csv(
                file,
                skiprows=1,
                na_values="--",
                names=["rsid", "chrom", "pos", "genotype"],
                index_col=0,
                dtype=NORMALIZED_DTYPES,
            )
    except ParseError:
        raise
    except ValueError:
        # read files with second header for concatenated data
        with source.open() as file:
            parts = extract_comments(file, include_data=True)
        # reconstruct file content from `extract_comments` results
        text = parts[0] + parts[2]
        lines = [line.strip() for line in text.split("\n")]
        # find index of second header
        second_header_idx = lines.index("RSID,CHROMOSOME,POSITION,RESULT", 1)

        with source.open() as file:
            df = pd.read_csv(
                file,
                skiprows=[0, second_header_idx],
                na_values="--",
                names=["rsid", "chrom", "pos", "genotype"],
                index_col=0,
                dtype=NORMALIZED_DTYPES,
            )

    return ParsedData(df)


def _parse_famfinder(parse_input: ParseInput) -> ParsedData:
    """Read and parse Family Tree DNA (FTDNA) "famfinder" file.

    https://www.familytreedna.com

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
        df = pd.read_csv(
            file,
            comment="#",
            na_values="-",
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

    return ParsedData(df)


FTDNA = ReaderFormat(name="FTDNA", detect=_detect_ftdna, parse=_parse_ftdna)

FAMFINDER = ReaderFormat(
    name="ftdna_famfinder",
    detect=_detect_famfinder,
    parse=_parse_famfinder,
    source=_source_famfinder,
)
