"""23andMe raw data files.

https://www.23andme.com
"""

from __future__ import annotations

import pandas as pd

from snps.constants import NORMALIZED_DTYPES
from snps.io.readers._registry import Header, ParsedData, ParseInput, ReaderFormat
from snps.io.readers._table import (
    TWO_ALLELE_COLUMNS,
    TWO_ALLELE_DTYPES,
    TableLayout,
    read_table,
)

_GENOTYPE_LAYOUT = TableLayout(
    na_values=("--", "-"),
    dtypes={**NORMALIZED_DTYPES, "pos": pd.UInt32Dtype()},
    drop_incomplete=True,
)

_ALLELES_LAYOUT = TableLayout(
    columns=TWO_ALLELE_COLUMNS,
    na_values=("--", "-"),
    dtypes={**TWO_ALLELE_DTYPES, "pos": pd.UInt32Dtype()},
    alleles="both",
    drop_incomplete=True,
)


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
    joined = _joined(parse_input.header)
    if joined is None:
        return None
    layout = _GENOTYPE_LAYOUT if joined else _ALLELES_LAYOUT
    return ParsedData(read_table(parse_input.source, layout))


FORMAT = ReaderFormat(name="23andMe", detect=_detect, parse=_parse, source=_source)
