"""Readers for raw data file formats, in the order detection tries them."""

from __future__ import annotations

import errno
import gzip
import io
import os
import zipfile
from collections.abc import Collection
from typing import Literal

import pandas as pd

from snps.data import SNPData
from snps.exceptions import ParseError, SNPsError, UnsupportedFormatError
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
from snps.io.readers._common import (
    detect_build_from_comments,
    extract_comments,
    handle_bytes_data,
)
from snps.io.readers._registry import (
    Header,
    ParseInput,
    ReaderFormat,
    detect,
)
from snps.resources import ResourceProvider, get_default_provider

__all__ = ["formats", "read_file"]

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


def _read_header(
    file: str | bytes,
) -> tuple[str | io.BytesIO, Literal["zip", "gzip"] | None, Header]:
    """Read the start of an existing file at a path, or of a file held in bytes.

    Parameters
    ----------
    file : str or bytes
        path to an existing file, or the contents of a file

    Returns
    -------
    tuple
        the path, or a buffer holding the bytes; the compression of the file
        (``"zip"``, ``"gzip"`` or None); and the start of the file
    """
    buffer: str | io.BytesIO
    compression: Literal["zip", "gzip"] | None
    if isinstance(file, str):
        buffer = file
        if ".zip" in file:
            with zipfile.ZipFile(file) as z:
                with z.open(z.namelist()[0], "r") as f:
                    first_line, comments, _ = extract_comments(f, decode=True)
                compression = "zip"
        elif ".gz" in file:
            with gzip.open(file, "rt") as f:
                first_line, comments, _ = extract_comments(f)
            compression = "gzip"
        else:
            with open(file, "rb") as f:
                first_line, comments, _, compression = handle_bytes_data(f.read())
    else:
        first_line, comments, _, compression = handle_bytes_data(file)
        buffer = io.BytesIO(file)
    return buffer, compression, Header(first_line, comments)


def _problem(exc: Exception) -> str:
    """Describe what a parser found wrong with a file."""
    if isinstance(exc, UnicodeDecodeError):
        return "the file could not be decoded as text"
    if isinstance(exc, pd.errors.ParserError):
        return f"rows do not have the expected fields ({exc})"
    if isinstance(exc, KeyError):
        return f"expected column or field {exc} is missing"
    if isinstance(exc, IndexError):
        return f"a row or field is shorter than expected ({exc})"
    return f"unexpected value ({exc})"


def _parse_error(
    name: str, format_name: str, problem: str, *, advise: bool = True
) -> ParseError:
    """Return the error for a file in a recognized format that cannot be parsed."""
    message = f"Cannot read {name} ({format_name} format): {problem.rstrip('.')}"
    if advise:
        message += "; check that the file is complete and unmodified"
    return ParseError(message + ".")


def _unsupported_layout(name: str, format_name: str) -> UnsupportedFormatError:
    """Return the error for a file in a recognized format with an unknown layout."""
    return UnsupportedFormatError(
        f"{name} appears to be in {format_name} format, "
        "but its column layout is not supported."
    )


def _raised_by(exc: BaseException, resources: ResourceProvider) -> bool:
    """Return whether a method of ``resources`` is in the traceback of ``exc``."""
    tb = exc.__traceback__
    while tb is not None:
        if tb.tb_frame.f_locals.get("self") is resources:
            return True
        tb = tb.tb_next
    return False


def read_file(
    file: str | bytes,
    *,
    resources: ResourceProvider | None = None,
    rsids: Collection[str] = (),
) -> SNPData:
    """Read SNPs from a raw data file.

    Parameters
    ----------
    file : str or bytes
        path to the file, or the contents of the file; either may be zip or
        gzip compressed
    resources : ResourceProvider, optional
        provider of external resources; defaults to the default provider
    rsids : collection of str, optional
        rsids to extract from a VCF file; empty for all

    Returns
    -------
    SNPData
        SNPs in the order of the file and not cleaned, with the source of the
        file and the build named by the file (0 if it names none)

    Raises
    ------
    TypeError
        if ``file`` is not a str or bytes
    FileNotFoundError
        if no file exists at the path
    IsADirectoryError
        if the path is a directory
    UnsupportedFormatError
        if the file is in no supported format or layout
    ParseError
        if the file is in a supported format but its content cannot be parsed

    Errors raised in a method of ``resources`` pass through unchanged.
    """
    if isinstance(file, str):
        if not os.path.exists(file):
            raise FileNotFoundError(errno.ENOENT, os.strerror(errno.ENOENT), file)
        if os.path.isdir(file):
            raise IsADirectoryError(
                errno.EISDIR, "Is a directory; pass the path of one raw data file", file
            )
        name = f"'{file}'"
    elif isinstance(file, bytes):
        name = "bytes input"
    else:
        raise TypeError(
            f"file must be a str path or the bytes of a file, not {type(file).__name__}"
        )

    try:
        buffer, compression, header = _read_header(file)
    except UnicodeDecodeError as exc:
        raise UnsupportedFormatError(
            f"{name} is not a text raw data file in a supported format."
        ) from exc

    fmt = detect(header, formats())
    if fmt is None:
        raise UnsupportedFormatError(
            f"{name} is not in a supported raw data format; "
            "check that it is a genotype raw data file."
        )

    file_source = fmt.source(header) if fmt.source is not None else fmt.name
    if file_source is None:
        raise _unsupported_layout(name, fmt.name)

    if resources is None:
        resources = get_default_provider()  # type: ignore[no-untyped-call]  # snps.resources is untyped

    try:
        parsed = fmt.parse(
            ParseInput(buffer, compression, header, resources, frozenset(rsids))
        )
        if parsed is None:
            raise _unsupported_layout(name, fmt.name)
        build = parsed.build or detect_build_from_comments(header.comments, file_source)
    except ParseError as exc:
        raise _parse_error(name, fmt.name, str(exc), advise=False) from exc
    except SNPsError:
        raise
    except (
        ValueError,
        KeyError,
        IndexError,
        TypeError,
        pd.errors.ParserError,
        UnicodeDecodeError,
    ) as exc:
        # a failure inside the provider is not a fault of the file
        if _raised_by(exc, resources):
            raise
        raise _parse_error(name, fmt.name, _problem(exc)) from exc

    frame = parsed.frame
    # a parser that reports the whole file as phased yields phased calls
    if parsed.phased and "phased" not in frame.columns:
        frame = frame.assign(phased=frame["genotype"].notna())

    return SNPData(
        frame,
        build=build,
        build_detected=bool(build),
        build_original=build,
        source=tuple(file_source.split(", ")),
    )
