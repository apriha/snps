"""Helpers shared by the raw data file readers."""

from __future__ import annotations

import binascii
import gzip
import io
import zipfile
from copy import deepcopy
from typing import IO, Literal

import numpy as np

from snps.build_constants import CHROM_SIZES
from snps.exceptions import ParseError

TWO_ALLELE_DTYPES: dict[str, type] = {
    "rsid": object,
    "chrom": object,
    "pos": np.uint32,
    "allele1": object,
    "allele2": object,
}


def extract_comments(
    f: IO[bytes] | IO[str] | io.BufferedIOBase,
    decode: bool = False,
    include_data: bool = False,
) -> tuple[str, str, str]:
    """Read the first line, the leading comments and optionally the data of a file.

    Parameters
    ----------
    f : file object
        file to read; rewound afterwards unless it is a zip member
    decode : bool
        decode the bytes read from ``f`` as UTF-8
    include_data : bool
        also read the rest of the file

    Returns
    -------
    tuple of str
        first line, comments, and data (empty unless ``include_data``)
    """
    line = read_line(f, decode)

    first_line = line
    comments = ""
    data = ""

    if first_line.startswith("#"):
        while line.startswith("#"):
            comments += line
            line = read_line(f, decode)
        if include_data:
            while line:
                data += line
                line = read_line(f, decode)

    elif first_line.startswith("[Header]"):
        while not line.startswith("[Data]"):
            comments += line
            line = read_line(f, decode)
        # Ignore the [Data] row
        line = read_line(f, decode)
        if include_data:
            while line:
                data += line
                line = read_line(f, decode)
    if not data and include_data:
        rest = f.read()
        if decode:
            assert isinstance(rest, bytes)
            data = rest.decode()
        else:
            assert isinstance(rest, str)
            data = rest
    if not isinstance(f, zipfile.ZipExtFile):
        f.seek(0)
    return first_line, comments, data


def detect_build_from_comments(comments: str, source: str) -> int:
    """Detect the build of a file from its comments.

    Parameters
    ----------
    comments : str
        leading comments of the file
    source : str
        source reported for the file; ``"vcf"`` selects the VCF header rules

    Returns
    -------
    int
        detected build, or 0

    Raises
    ------
    ParseError
        if a VCF header line is malformed
    """
    # if its a VCF parse these properly
    if source == "vcf":
        for line in comments.split("\n"):
            line = line.strip()
            if not line:
                # skip blanks
                continue
            if not line.startswith("#"):
                raise ParseError(f"Header line does not start with '#': {line}")
            if not line.startswith("##"):
                # skip comments but not preamble
                continue
            if "=" not in line:
                # skip lines without key/value in
                continue
            header_line = line
            line = line[2:].strip()
            key = line[: line.index("=")]
            value = line[line.index("=") + 1 :]
            if key.lower() == "contig":
                if not (value.startswith("<") and value.endswith(">")):
                    raise ParseError(
                        f"Contig line is not of the form ##contig=<...>: {header_line}"
                    )
                parts = value[1:-1].split(",")
                for part in parts:
                    if "=" not in part:
                        raise ParseError(
                            f"Contig line has a field without '=': {header_line}"
                        )
                    part_key = part[: part.index("=")]
                    part_value = part[part.index("=") + 1 :]
                    if part_key.lower() == "assembly":
                        if "36" in part_value:
                            return 36
                        elif "37" in part_value or "hg19" in part_value:
                            return 37
                        elif "38" in part_value:
                            return 38
                    elif part_key.lower() == "length":
                        # Check chromosome 1 length against known builds
                        if part_value == str(CHROM_SIZES[37]["1"]):
                            return 37
                        elif part_value == str(CHROM_SIZES[38]["1"]):
                            return 38
                        elif part_value == str(CHROM_SIZES[36]["1"]):
                            return 36
        # couldn't find anything
        return 0
    else:
        # not a vcf
        if "build 36" in comments.lower():
            return 36
        elif "build 37" in comments.lower():
            return 37
        elif "build 38" in comments.lower():
            return 38
        # these can cause false positives
        # elif "b37" in comments.lower():
        #    return 37
        # elif "b38" in comments.lower():
        #    return 38
        # elif "hg19" in comments.lower():
        #    return 37
        # elif "hg38" in comments.lower():
        #    return 38
        elif "grch38" in comments.lower():
            return 38
        elif "grch37" in comments.lower():
            return 37
        # Check for chromosome 1 length in comments
        elif str(CHROM_SIZES[37]["1"]) in comments.lower():
            return 37
        elif str(CHROM_SIZES[38]["1"]) in comments.lower():
            return 38
        elif str(CHROM_SIZES[36]["1"]) in comments.lower():
            return 36
        else:
            return 0


def handle_bytes_data(
    file: bytes, include_data: bool = False
) -> tuple[str, str, str, Literal["zip", "gzip"] | None]:
    """Read the first line, comments and optionally data of a file held in bytes.

    Parameters
    ----------
    file : bytes
        contents of the file, which may be zip or gzip compressed
    include_data : bool
        also read the rest of the file

    Returns
    -------
    tuple
        first line, comments, data (empty unless ``include_data``), and the
        compression of the file (``"zip"``, ``"gzip"`` or None)
    """
    compression: Literal["zip", "gzip"] | None
    if is_zip(file):
        compression = "zip"
        with zipfile.ZipFile(io.BytesIO(file)) as z:
            namelist = z.namelist()
            key = "GFG_filtered_unphased_genotypes_23andMe.txt"
            key_search = [key in name for name in namelist]

            if any(key_search):
                filename = namelist[key_search.index(True)]
            else:
                filename = namelist[0]

            with z.open(filename, "r") as f:
                first_line, comments, data = extract_comments(
                    f, decode=True, include_data=include_data
                )

    elif is_gzip(file):
        compression = "gzip"

        with gzip.open(io.BytesIO(file), "rb") as f:
            first_line, comments, data = extract_comments(
                f, decode=True, include_data=include_data
            )

    else:
        compression = None
        buffer = io.BytesIO(file)
        first_line, comments, data = extract_comments(
            deepcopy(buffer), decode=True, include_data=include_data
        )
        buffer.seek(0)
    return first_line, comments, data, compression


def is_zip(bytes_data: bytes) -> bool:
    """Check whether or not a bytes_data file is a valid Zip file."""
    return zipfile.is_zipfile(io.BytesIO(bytes_data))


def is_gzip(bytes_data: bytes) -> bool:
    """Check whether or not a bytes_data file is a valid gzip file."""
    return binascii.hexlify(bytes_data[:2]) == b"1f8b"


def read_line(f: IO[bytes] | IO[str] | io.BufferedIOBase, decode: bool) -> str:
    """Read a line from a file, decoding it from UTF-8 bytes if ``decode``."""
    line = f.readline()
    if decode:
        assert isinstance(line, bytes)
        # https://stackoverflow.com/a/606199
        return line.decode("utf-8")
    else:
        assert isinstance(line, str)
        return line
