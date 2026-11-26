"""SelfDecode raw data files.

https://selfdecode.com
"""

from __future__ import annotations

from functools import partial

import pandas as pd

from snps.constants import NORMALIZED_DTYPES
from snps.io.readers._registry import Header, ReaderFormat
from snps.io.readers._table import TableLayout, parse_table

_LAYOUT = TableLayout(
    dtypes={**NORMALIZED_DTYPES, "pos": pd.UInt32Dtype()}, drop_incomplete=True
)


def _detect(header: Header) -> bool:
    return "SelfDecode" in header.comments


FORMAT = ReaderFormat(
    name="SelfDecode", detect=_detect, parse=partial(parse_table, _LAYOUT)
)
