"""tellmeGen raw data files."""

from __future__ import annotations

import pandas as pd

from snps.constants import NORMALIZED_DTYPES
from snps.io.readers._registry import Header, ParsedData, ParseInput, ReaderFormat


def _detect(header: Header) -> bool:
    return "rsid\tChromosome\tposition\tgenotype" == header.first_line.strip()


def _parse(parse_input: ParseInput) -> ParsedData:
    """Read and parse tellmeGen files.

    https://www.tellmegen.com/

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
    resources = parse_input.resources

    df = pd.read_csv(
        file,
        sep="\t",
        skiprows=1,
        na_values="--",
        names=["rsid", "chrom", "pos", "genotype"],
        dtype=NORMALIZED_DTYPES,
        compression=compression,
    )

    # use the external mapping to turn
    # SNP names into rsids where possible
    df = df.merge(
        resources.get_gsa_rsid(),  # type: ignore[no-untyped-call]  # snps.resources is untyped
        how="left",  # left-hand join, gsa rsids may be NA
        left_on="rsid",
        right_on="gsaname_rsid",
        suffixes=(None, "_gsa_rsid"),
    )
    df.loc[~pd.isna(df["gsaname_rsid"]), "rsid"] = df.loc[
        ~pd.isna(df["gsaname_rsid"]), "gsarsid"
    ]

    # keep only the columns we want
    df = df.filter(items=("rsid", "chrom", "pos", "genotype"), axis=1)

    # reindex for the new identifiers
    df = df.set_index(["rsid"])

    return ParsedData(df)


FORMAT = ReaderFormat(name="tellmeGen", detect=_detect, parse=_parse)
