"""The shared reader of raw data files whose data is one delimited table."""

from __future__ import annotations

from collections.abc import Mapping
from dataclasses import dataclass, field
from types import MappingProxyType
from typing import TYPE_CHECKING, Literal

import numpy as np
import pandas as pd

from snps.constants import NORMALIZED_DTYPES
from snps.exceptions import ParseError
from snps.io.readers._common import parse_positions
from snps.io.readers._registry import ParsedData

if TYPE_CHECKING:
    from pandas.api.extensions import ExtensionDtype

    from snps.io.readers._registry import ParseInput
    from snps.io.readers._source import InputSource

TWO_ALLELE_COLUMNS = ("rsid", "chrom", "pos", "allele1", "allele2")
TWO_ALLELE_DTYPES: Mapping[str, type] = MappingProxyType(
    {**dict.fromkeys(TWO_ALLELE_COLUMNS, object), "pos": np.uint32}
)


@dataclass(frozen=True, slots=True)
class TableLayout:
    """How the table of a tabular raw data file format is laid out.

    Attributes
    ----------
    columns : tuple of str or None
        names given to the columns, or None to keep the names in the header row
    separator : str or None
        field separator; None sniffs it with pandas' python engine
    na_values : tuple of str
        strings read as missing, besides pandas' defaults
    comment : str or None
        character that starts a comment, or None
    skiprows : int
        number of lines skipped at the start of the file
    header : int or None
        row of column names, replaced by ``columns`` if given; None for none
    dtypes : mapping of str to dtype
        dtype of each named column as it is read; a ``pos`` column read as
        ``str`` is then converted with ``parse_positions``
    index : int or None
        column read as the index; None for none
    usecols : tuple of int or None
        fields read from each row, by position; None for all
    alleles : {"both", "allele2_optional"} or None
        how columns ``allele1`` and ``allele2`` join into column ``genotype``:
        ``"both"`` gives a genotype only if both alleles are given,
        ``"allele2_optional"`` also gives ``allele1`` alone if ``allele2`` is
        missing; None if the file has a genotype column
    chromosomes : mapping or None
        chromosome for each code in column ``chrom``, a code not in it being
        missing; None keeps the column as read
    drop_incomplete : bool
        whether rows with no rsid, chromosome or position are dropped and the
        rest indexed by rsid, with normalized dtypes; ``index`` is then ignored
    repeated_header : str or None
        header line that the file repeats, as files of concatenated tables do;
        the read also skips its first repeat after the first line
    """

    columns: tuple[str, ...] | None = ("rsid", "chrom", "pos", "genotype")
    separator: str | None = "\t"
    na_values: tuple[str, ...] = ("--",)
    comment: str | None = "#"
    skiprows: int = 0
    header: int | None = None
    dtypes: Mapping[str, type | ExtensionDtype] = field(
        default_factory=NORMALIZED_DTYPES.copy
    )
    index: int | None = 0
    usecols: tuple[int, ...] | None = None
    alleles: Literal["both", "allele2_optional"] | None = None
    chromosomes: Mapping[str | int, str] | None = None
    drop_incomplete: bool = False
    repeated_header: str | None = None

    def __post_init__(self) -> None:
        # the layout owns read-only copies of the mappings it is given
        object.__setattr__(self, "dtypes", MappingProxyType(dict(self.dtypes)))
        if self.chromosomes is not None:
            chromosomes = MappingProxyType(dict(self.chromosomes))
            object.__setattr__(self, "chromosomes", chromosomes)


def read_table(
    source: InputSource, layout: TableLayout, *fallbacks: TableLayout
) -> pd.DataFrame:
    """Read the table of a raw data file.

    Parameters
    ----------
    source : InputSource
        the file
    layout : TableLayout
        how its table is laid out
    *fallbacks : TableLayout
        layouts tried in turn while reading with the one before fails with
        ValueError

    Returns
    -------
    pandas.DataFrame
        the table

    Raises
    ------
    ParseError
        if the compressed data is corrupt or truncated
    ValueError
        if the data fits none of the layouts
    """
    try:
        frame = _read_csv(source, layout)
    except ParseError:
        # corrupt compressed data is no reason to try another layout
        raise
    except ValueError:
        if not fallbacks:
            raise
        # inside this handler, so that a later error carries this one
        return read_table(source, *fallbacks)
    if layout.chromosomes is not None:
        frame["chrom"] = frame["chrom"].map(dict(layout.chromosomes))
    if layout.drop_incomplete:
        # positions first, so the text columns are checked on fewer rows
        frame = frame.dropna(subset=["pos"]).dropna(subset=["rsid", "chrom"])
    if layout.alleles is not None:
        allele2 = frame["allele2"]
        if layout.alleles == "allele2_optional":
            allele2 = allele2.fillna("")
        frame["genotype"] = frame["allele1"] + allele2
        frame = frame.drop(columns=["allele1", "allele2"])
    if layout.dtypes.get("pos") is str:
        frame["pos"] = parse_positions(frame["pos"])
    if layout.drop_incomplete:
        frame = frame.astype(NORMALIZED_DTYPES).set_index("rsid")
    return frame


def parse_table(layout: TableLayout, parse_input: ParseInput) -> ParsedData:
    """Parse a file in a format that its table layout describes alone.

    Parameters
    ----------
    layout : TableLayout
        how the table is laid out
    parse_input : ParseInput
        file to parse

    Returns
    -------
    ParsedData
        parsed SNPs
    """
    return ParsedData(read_table(parse_input.source, layout))


def _read_csv(source: InputSource, layout: TableLayout) -> pd.DataFrame:
    """Read the table as pandas parses it."""
    skiprows: int | list[int] = layout.skiprows
    if layout.repeated_header is not None:
        repeat = _find_repeat(source, layout.repeated_header)
        skiprows = [*range(layout.skiprows), repeat]
    with source.open() as file:
        return pd.read_csv(
            file,
            sep=layout.separator,
            # only the python engine sniffs a separator
            engine="python" if layout.separator is None else None,
            header=layout.header,
            names=None if layout.columns is None else list(layout.columns),
            index_col=None if layout.drop_incomplete else layout.index,
            usecols=None if layout.usecols is None else list(layout.usecols),
            dtype=dict(layout.dtypes),
            skiprows=skiprows,
            na_values=list(layout.na_values),
            comment=layout.comment,
        )


def _find_repeat(source: InputSource, header: str) -> int:
    """Return the number of the first line after the first that is ``header``."""
    with source.open() as file:
        text = file.read().decode()
    lines = [line.strip() for line in text.split("\n")]
    return lines.index(header, 1)
