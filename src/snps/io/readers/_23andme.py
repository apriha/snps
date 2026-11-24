"""23andMe raw data files."""

from __future__ import annotations

import numpy as np
import pandas as pd

from snps.constants import NORMALIZED_DTYPES
from snps.io.readers._common import TWO_ALLELE_DTYPES
from snps.io.readers._registry import Header, ParsedData, ParseInput, ReaderFormat


def _joined(header: Header) -> bool | None:
    """Return whether the genotype is one column, or None for an unsupported layout."""
    comments = header.comments
    # some 23andMe files have separate alleles
    if comments.endswith(
        "# rsid\tchromosome\tposition\tallele1\tallele2\n"
    ) or comments.endswith("# rsid\tchromosome\tposition\tallele1\tallele2\r\n"):
        return False
    # some 23andMe files have a combined genotype
    elif comments.endswith(
        "# rsid\tchromosome\tposition\tgenotype\n"
    ) or comments.endswith("# rsid\tchromosome\tposition\tgenotype\r\n"):
        return True
    # something we havent seen before and can't handle
    else:
        return None


def _detect(header: Header) -> bool:
    return "23andMe" in header.first_line


def _source(header: Header) -> str | None:
    return None if _joined(header) is None else "23andMe"


def _parse(parse_input: ParseInput) -> ParsedData | None:
    """Read and parse 23andMe file.

    https://www.23andme.com

    Parameters
    ----------
    parse_input : ParseInput
        file to parse

    Returns
    -------
    ParsedData or None
        parsed SNPs, or None if the column layout is not supported
    """
    joined = _joined(parse_input.header)
    if joined is None:
        return None
    if joined:
        columnnames = ["rsid", "chrom", "pos", "genotype"]
        dtype = NORMALIZED_DTYPES.copy()
    else:
        columnnames = ["rsid", "chrom", "pos", "allele1", "allele2"]
        dtype = TWO_ALLELE_DTYPES.copy()

    # Temporarily use nullable UInt32 for 'pos' column
    dtype["pos"] = pd.UInt32Dtype()  # type: ignore[assignment]  # read_csv takes dtype instances

    with parse_input.source.open() as file:
        df = pd.read_csv(
            file,
            comment="#",
            sep="\t",
            na_values=["--", "-"],
            names=columnnames,
            dtype=dtype,
        )

    # Drop rows with NaN values in 'pos' column
    df = df.dropna(subset=["pos"])

    # Convert 'pos' column to np.uint32
    df["pos"] = df["pos"].astype(np.uint32)

    if not joined:
        # stick separate alleles together
        df["genotype"] = df["allele1"] + df["allele2"]
        del df["allele1"]
        del df["allele2"]
    df = df.dropna(subset=["rsid", "chrom", "pos"])
    df = df.astype(dtype=NORMALIZED_DTYPES)
    df = df.set_index("rsid")
    return ParsedData(df)


FORMAT = ReaderFormat(name="23andMe", detect=_detect, parse=_parse, source=_source)
