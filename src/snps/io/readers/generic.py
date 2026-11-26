"""Generic CSV or TSV files."""

from __future__ import annotations

import re
from dataclasses import replace

from snps.constants import NORMALIZED_DTYPES
from snps.io.readers._registry import Header, ParsedData, ParseInput, ReaderFormat
from snps.io.readers._table import TableLayout, read_table

_HEADER = "^#*[ \t]*rsid[, \t]*chr"
_NO_HEADER = "^rs[0-9]*[, \t]{1}[1]"

_COMMA_LAYOUT = TableLayout(separator=",", comment=None)

_LAYOUTS = (
    _COMMA_LAYOUT,
    replace(_COMMA_LAYOUT, separator="\t"),
    replace(
        _COMMA_LAYOUT,
        separator=None,
        usecols=(0, 1, 2, 3),
        dtypes={**NORMALIZED_DTYPES, "pos": str},
    ),
)


def _detect(header: Header) -> bool:
    return (
        re.match(_HEADER, header.first_line) is not None
        or re.match(_NO_HEADER, header.first_line) is not None
    )


def _parse(parse_input: ParseInput) -> ParsedData:
    skip = 1 if re.match(_HEADER, parse_input.header.first_line) else 0
    layouts = [replace(layout, skiprows=skip) for layout in _LAYOUTS]
    return ParsedData(read_table(parse_input.source, *layouts))


FORMAT = ReaderFormat(name="generic", detect=_detect, parse=_parse)
