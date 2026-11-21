"""Readers for raw data file formats, in the order detection tries them."""

from __future__ import annotations

from snps.io.readers import (
    _23andme,
    _23mofang,
    ancestry,
    circledna,
    dnaland,
    ftdna,
    generic,
    genesforgood,
    gsa,
    livingdna,
    mapmygenome,
    myheritage,
    plink,
    sano_dtc,
    selfdecode,
    snps_csv,
    tellmegen,
    vcf,
)
from snps.io.readers._registry import ReaderFormat

# detection uses the first format that claims a header, and several claims overlap
_FORMATS: tuple[ReaderFormat, ...] = (
    _23andme.FORMAT,
    ancestry.FORMAT,
    ftdna.FTDNA,
    ftdna.FAMFINDER,
    myheritage.FORMAT,
    livingdna.FORMAT,
    mapmygenome.FORMAT,
    snps_csv.FORMAT,
    tellmegen.FORMAT,
    _23mofang.FORMAT,
    sano_dtc.FORMAT,
    generic.FORMAT,
    vcf.FORMAT,
    genesforgood.FORMAT,
    dnaland.FORMAT,
    selfdecode.FORMAT,
    gsa.SANO,
    gsa.CODIGO46,
    gsa.GSA,
    circledna.FORMAT,
    plink.FORMAT,
)


def formats() -> tuple[ReaderFormat, ...]:
    """Return the known formats, in the order detection tries them."""
    return _FORMATS
