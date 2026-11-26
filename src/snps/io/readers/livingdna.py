"""LivingDNA raw data files.

https://livingdna.com
"""

from __future__ import annotations

from functools import partial

from snps.io.readers._registry import Header, ReaderFormat
from snps.io.readers._table import TableLayout, parse_table

_LAYOUT = TableLayout()


def _detect(header: Header) -> bool:
    return "Living DNA" in header.first_line


FORMAT = ReaderFormat(
    name="LivingDNA", detect=_detect, parse=partial(parse_table, _LAYOUT)
)
