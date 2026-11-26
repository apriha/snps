"""Genes for Good raw data files.

https://genesforgood.sph.umich.edu

References
----------
1. Genes for Good, University of Michigan. Genes for Good Genotypes Download
   Version 1.2 README. 2016. https://genesforgood.sph.umich.edu/readme/readme1.2.txt
"""

from __future__ import annotations

from functools import partial

from snps.io.readers._registry import Header, ReaderFormat
from snps.io.readers._table import TableLayout, parse_table

_LAYOUT = TableLayout()


def _detect(header: Header) -> bool:
    return "Genes for Good" in header.comments or "PLINK" in header.comments


FORMAT = ReaderFormat(
    name="GenesForGood", detect=_detect, parse=partial(parse_table, _LAYOUT)
)
