"""tellmeGen raw data files.

https://www.tellmegen.com
"""

from __future__ import annotations

import pandas as pd

from snps.io.readers._registry import Header, ParsedData, ParseInput, ReaderFormat
from snps.io.readers._table import TableLayout, read_table

_LAYOUT = TableLayout(comment=None, skiprows=1, index=None)


def _detect(header: Header) -> bool:
    return "rsid\tChromosome\tposition\tgenotype" == header.first_line.strip()


def _parse(parse_input: ParseInput) -> ParsedData:
    resources = parse_input.resources

    df = read_table(parse_input.source, _LAYOUT)

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
