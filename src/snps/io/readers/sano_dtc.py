"""Sano Genetics direct-to-consumer raw data files.

https://sanogenetics.com
"""

from __future__ import annotations

import re
from functools import partial

from snps.io.readers._registry import Header, ReaderFormat
from snps.io.readers._table import (
    TWO_ALLELE_COLUMNS,
    TWO_ALLELE_DTYPES,
    TableLayout,
    parse_table,
)

_LAYOUT = TableLayout(
    columns=TWO_ALLELE_COLUMNS,
    separator=r"\s+",
    na_values=("-",),
    header=0,
    dtypes=TWO_ALLELE_DTYPES,
    alleles="allele2_optional",
)


def _detect(header: Header) -> bool:
    return (
        "# This file was derived from the corresponding VCF" in header.comments
        or re.match(
            r"^\s*rsid\s+chromosome\s+position\s+allele_1\s+allele_2\s*$",
            header.first_line,
        )
        is not None
    )


def _source(header: Header) -> str:
    return "Sano"


FORMAT = ReaderFormat(
    name="sano_dtc",
    detect=_detect,
    parse=partial(parse_table, _LAYOUT),
    source=_source,
)
