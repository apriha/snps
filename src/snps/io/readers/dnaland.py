"""DNA.Land raw data files.

https://dna.land
"""

from __future__ import annotations

from functools import partial

from snps.io.readers._registry import Header, ReaderFormat
from snps.io.readers._table import TableLayout, parse_table

_LAYOUT = TableLayout()


def _detect(header: Header) -> bool:
    return "DNA.Land" in header.comments


FORMAT = ReaderFormat(
    name="DNA.Land", detect=_detect, parse=partial(parse_table, _LAYOUT)
)
