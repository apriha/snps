"""PLINK text files.

https://www.cog-genomics.org/plink/
"""

from __future__ import annotations

from collections.abc import Mapping
from functools import partial
from types import MappingProxyType

from snps.io.readers._registry import Header, ReaderFormat
from snps.io.readers._table import TableLayout, parse_table

# integer keys match codes pandas reads as numbers (floats when a code is
# missing); a row whose code is not listed is dropped
_CHROMOSOMES: Mapping[str | int, str] = MappingProxyType(
    {code: str(n) for n in range(1, 23) for code in (str(n), n)}
    | {"X": "X", "Y": "Y", "XY": "Y", "MT": "MT"}
)

_LAYOUT = TableLayout(
    na_values=("--", "-"),
    dtypes={"pos": str},
    chromosomes=_CHROMOSOMES,
    drop_incomplete=True,
)


def _detect(header: Header) -> bool:
    return "# Below is a text version of your data." in header.comments


FORMAT = ReaderFormat(name="PLINK", detect=_detect, parse=partial(parse_table, _LAYOUT))
