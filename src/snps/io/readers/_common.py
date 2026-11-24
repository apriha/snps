"""Helpers shared by the raw data file readers."""

from __future__ import annotations

from typing import IO

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


def extract_comments(f: IO[bytes], include_data: bool = False) -> tuple[str, str, str]:
    """Read the first line, the leading comments and optionally the data of a file.

    Parameters
    ----------
    f : file object
        binary file to read, decoded as UTF-8
    include_data : bool
        also read the rest of the file

    Returns
    -------
    tuple of str
        first line, comments, and data (empty unless ``include_data``)
    """
    line = read_line(f)

    first_line = line
    comments = ""
    data = ""

    if first_line.startswith("#"):
        while line.startswith("#"):
            comments += line
            line = read_line(f)
        if include_data:
            while line:
                data += line
                line = read_line(f)

    elif first_line.startswith("[Header]"):
        while line and not line.startswith("[Data]"):
            comments += line
            line = read_line(f)
        # Ignore the [Data] row
        line = read_line(f)
        if include_data:
            while line:
                data += line
                line = read_line(f)
    if not data and include_data:
        data = f.read().decode()
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


def read_line(f: IO[bytes]) -> str:
    """Read a line from a binary file, decoding it from UTF-8."""
    line = f.readline()
    # https://stackoverflow.com/a/606199
    return line.decode("utf-8")
