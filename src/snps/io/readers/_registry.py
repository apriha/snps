"""Format descriptors for raw data files and the detection that picks one."""

from __future__ import annotations

import io
from collections.abc import Callable, Sequence
from dataclasses import dataclass
from typing import TYPE_CHECKING, Literal

import pandas as pd

if TYPE_CHECKING:
    from snps.resources import ResourceProvider


@dataclass(frozen=True, slots=True)
class Header:
    """Start of a raw data file, which is all that format detection looks at.

    Attributes
    ----------
    first_line : str
        first line of the file
    comments : str
        leading comment lines (or ``[Header]`` section) of the file
    """

    first_line: str
    comments: str


@dataclass(frozen=True, slots=True)
class ParseInput:
    """Everything a parser receives.

    Attributes
    ----------
    file : str or io.BytesIO
        path to the file, or a buffer holding its bytes
    compression : str or None
        compression of the file: ``"zip"``, ``"gzip"`` or None
    header : Header
        start of the file
    resources : ResourceProvider or None
        provider of external resources
    rsids : frozenset of str
        rsids to extract from a VCF file; empty for all
    """

    file: str | io.BytesIO
    compression: Literal["zip", "gzip"] | None
    header: Header
    resources: ResourceProvider | None
    rsids: frozenset[str]


@dataclass(frozen=True, slots=True, eq=False)
class ParsedData:
    """What a parser returns.

    Attributes
    ----------
    frame : pandas.DataFrame
        parsed SNPs; may carry a boolean ``phased`` column
    phased : bool
        whether the whole file is phased
    build : int
        build found in the file itself, or 0
    """

    frame: pd.DataFrame
    phased: bool = False
    build: int = 0


@dataclass(frozen=True, slots=True)
class ReaderFormat:
    """Description of one raw data file format.

    Attributes
    ----------
    name : str
        unique key of the format, reported as the source unless ``source`` is set
    detect : callable
        returns whether a header belongs to this format
    parse : callable
        returns the parsed file, or None if its layout is not supported
    source : callable or None
        returns the source to report for a header: ``""`` if the file does not
        name one, or None if its layout is not supported
    """

    name: str
    detect: Callable[[Header], bool]
    parse: Callable[[ParseInput], ParsedData | None]
    source: Callable[[Header], str | None] | None = None


def detect(header: Header, formats: Sequence[ReaderFormat]) -> ReaderFormat | None:
    """Find the format of a file from its header.

    Parameters
    ----------
    header : Header
        start of the file
    formats : sequence of ReaderFormat
        formats to try, in order

    Returns
    -------
    ReaderFormat or None
        first format that claims the header, or None if none does
    """
    for fmt in formats:
        if fmt.detect(header):
            return fmt
    return None
