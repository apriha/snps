"""CSV files written by ``snps``."""

from __future__ import annotations

from dataclasses import replace

import pandas as pd

from snps.io.readers._registry import Header, ParsedData, ParseInput, ReaderFormat
from snps.io.readers._table import TableLayout, read_table

_COMMA_LAYOUT = TableLayout(separator=",", header=0)

_TAB_LAYOUT = replace(_COMMA_LAYOUT, separator="\t")


def _detect(header: Header) -> bool:
    return "lineage" in header.first_line or "snps" in header.first_line


def _read_comments(comments: str) -> tuple[str, bool, int]:
    """Read the source, phase and build from the comments of a ``snps`` CSV file.

    Parameters
    ----------
    comments : str
        comments at beginning of file

    Returns
    -------
    tuple
        source (``""`` if not given), whether the file is phased, and build (0
        unless it was detected)

    Raises
    ------
    ValueError
        if the build is not an integer
    """
    source = ""
    phased = False
    build = 0

    comment_lines = comments.split("\n")
    for comment1 in comment_lines:
        if "Source(s):" in comment1:
            source = comment1.split("Source(s):")[1].strip()
        elif "Phased:" in comment1:
            if comment1.split("Phased:")[1].strip() == "True":
                phased = True
        elif "Build:" in comment1:
            temp = int(comment1.split("Build:")[1].strip())
            for comment2 in comment_lines:
                if "Build Detected:" in comment2:
                    # only assign build if it was detected
                    if comment2.split("Build Detected:")[1].strip() == "True":
                        build = temp
                    break

    return source, phased, build


def _source(header: Header) -> str:
    source, _, _ = _read_comments(header.comments)
    return source


def _parse(parse_input: ParseInput) -> ParsedData:
    _, phased, build = _read_comments(parse_input.header.comments)
    try:
        return ParsedData(read_table(parse_input.source, _COMMA_LAYOUT), phased, build)
    except pd.errors.ParserError:
        return ParsedData(read_table(parse_input.source, _TAB_LAYOUT), phased, build)


FORMAT = ReaderFormat(name="snps_csv", detect=_detect, parse=_parse, source=_source)
