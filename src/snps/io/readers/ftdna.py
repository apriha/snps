"""Family Tree DNA (FTDNA) raw data files, including "famfinder" files.

https://www.familytreedna.com
"""

from __future__ import annotations

from dataclasses import replace
from functools import partial

from snps.io.readers._registry import Header, ParsedData, ParseInput, ReaderFormat
from snps.io.readers._table import (
    TWO_ALLELE_COLUMNS,
    TWO_ALLELE_DTYPES,
    TableLayout,
    parse_table,
    read_table,
)

_FTDNA_LAYOUT = TableLayout(separator=",", comment=None, skiprows=1)

_CONCATENATED_LAYOUT = replace(
    _FTDNA_LAYOUT, repeated_header="RSID,CHROMOSOME,POSITION,RESULT"
)

_FAMFINDER_LAYOUT = TableLayout(
    columns=TWO_ALLELE_COLUMNS,
    separator=",",
    na_values=("-",),
    dtypes=TWO_ALLELE_DTYPES,
    alleles="both",
)


def _detect_ftdna(header: Header) -> bool:
    return header.first_line.startswith("RSID")


def _detect_famfinder(header: Header) -> bool:
    return "famfinder" in header.first_line


def _source_famfinder(header: Header) -> str:
    return "FTDNA"


def _parse_ftdna(parse_input: ParseInput) -> ParsedData:
    frame = read_table(parse_input.source, _FTDNA_LAYOUT, _CONCATENATED_LAYOUT)
    return ParsedData(frame)


FTDNA = ReaderFormat(name="FTDNA", detect=_detect_ftdna, parse=_parse_ftdna)

FAMFINDER = ReaderFormat(
    name="ftdna_famfinder",
    detect=_detect_famfinder,
    parse=partial(parse_table, _FAMFINDER_LAYOUT),
    source=_source_famfinder,
)
