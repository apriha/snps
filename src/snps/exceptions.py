"""Exceptions raised by ``snps``."""

from __future__ import annotations

__all__ = ["ParseError", "SNPsError", "UnsupportedFormatError"]


class SNPsError(Exception):
    """Base class for errors raised by ``snps``."""


class UnsupportedFormatError(SNPsError, ValueError):
    """A file is not in any raw data format that ``snps`` can read."""


class ParseError(SNPsError, ValueError):
    """A file in a recognized format has content that cannot be parsed."""
