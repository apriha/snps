"""Mapmygenome raw data files.

https://mapmygenome.in
"""

from __future__ import annotations

import numpy as np

from snps.constants import NORMALIZED_DTYPES
from snps.io.readers._registry import Header, ParsedData, ParseInput, ReaderFormat
from snps.io.readers._table import TableLayout, read_table

_GENOTYPE_LAYOUT = TableLayout(
    columns=None,
    header=0,
    dtypes={
        "rsid": object,
        "chromosome": object,
        "position": np.uint32,
        "genotype": object,
    },
    index=None,
)


def _detect(header: Header) -> bool:
    return (
        "SNP Name\trsID" in header.first_line
        or "SNP.Name\tSample.ID" in header.first_line
        or "Mapmygenome" in header.comments
    )


def _allele_layout(rsid_col_name: str, rsid_col_idx: int) -> TableLayout:
    """Return the layout of a file with allele columns, indexed by its rsid column."""
    return TableLayout(
        columns=None,
        header=0,
        dtypes={
            rsid_col_name: object,
            "Chr": object,
            "Position": np.uint32,
            "Allele1...Top": object,
            "Allele2...Top": object,
        },
        index=rsid_col_idx,
    )


def _parse(parse_input: ParseInput) -> ParsedData:
    source = parse_input.source
    header = parse_input.header.first_line
    comments = parse_input.header.comments

    if "rsid\tchromosome\tposition\tgenotype" in header or "Mapmygenome" in comments:
        df = read_table(source, _GENOTYPE_LAYOUT)
        df = df.rename(columns={"chromosome": "chrom", "position": "pos"})
        df = df.astype(dtype=NORMALIZED_DTYPES)
        df = df.set_index("rsid")
        return ParsedData(df)
    else:
        if "rsID" in header:
            df = read_table(source, _allele_layout("rsID", 1))
        else:
            df = read_table(source, _allele_layout("SNP.Name", 0))

        df["genotype"] = df["Allele1...Plus"] + df["Allele2...Plus"]
        df = df.rename(columns={"Chr": "chrom", "Position": "pos"})
        df.index.name = "rsid"
        df = df[["chrom", "pos", "genotype"]]

        return ParsedData(df)


FORMAT = ReaderFormat(name="Mapmygenome", detect=_detect, parse=_parse)
