"""23Mofang raw data files.

https://www.23mofang.com
"""

from __future__ import annotations

from functools import partial

import pandas as pd

from snps.constants import NORMALIZED_DTYPES
from snps.io.readers._registry import Header, ReaderFormat
from snps.io.readers._table import TableLayout, parse_table

_LAYOUT = TableLayout(
    na_values=("--", "-"),
    dtypes={**NORMALIZED_DTYPES, "pos": pd.UInt32Dtype()},
    drop_incomplete=True,
)


def _detect(header: Header) -> bool:
    return "23Mofang" in header.first_line or "23Mofang" in header.comments


FORMAT = ReaderFormat(
    name="23Mofang", detect=_detect, parse=partial(parse_table, _LAYOUT)
)
