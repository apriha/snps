"""Family Tree DNA (FTDNA) raw data files, including "famfinder" files."""

from __future__ import annotations

import io
import zlib

import pandas as pd

from snps.constants import NORMALIZED_DTYPES
from snps.io.readers._common import TWO_ALLELE_DTYPES, handle_bytes_data
from snps.io.readers._registry import Header, ParsedData, ParseInput, ReaderFormat


def _detect_ftdna(header: Header) -> bool:
    return header.first_line.startswith("RSID")


def _detect_famfinder(header: Header) -> bool:
    return "famfinder" in header.first_line


def _source_famfinder(header: Header) -> str:
    return "FTDNA"


def _parse_ftdna(parse_input: ParseInput) -> ParsedData:
    """Read and parse Family Tree DNA (FTDNA) file.

    https://www.familytreedna.com

    Parameters
    ----------
    parse_input : ParseInput
        file to parse

    Returns
    -------
    ParsedData
        parsed SNPs
    """
    file = parse_input.file
    compression = parse_input.compression

    try:
        df = pd.read_csv(
            file,
            skiprows=1,
            na_values="--",
            names=["rsid", "chrom", "pos", "genotype"],
            index_col=0,
            dtype=NORMALIZED_DTYPES,
            compression=compression,
        )
    except ValueError:
        # read files with second header for concatenated data
        if isinstance(file, io.BytesIO):
            file.seek(0)
            parts = handle_bytes_data(file.read(), include_data=True)
            file.seek(0)
        else:
            with open(file, "rb") as f:
                parts = handle_bytes_data(f.read(), include_data=True)
        # reconstruct file content from `handle_bytes_data` results
        text = parts[0] + parts[2]
        lines = [line.strip() for line in text.split("\n")]
        # find index of second header
        second_header_idx = lines.index("RSID,CHROMOSOME,POSITION,RESULT", 1)

        df = pd.read_csv(
            file,
            skiprows=[0, second_header_idx],
            na_values="--",
            names=["rsid", "chrom", "pos", "genotype"],
            index_col=0,
            dtype=NORMALIZED_DTYPES,
            compression=compression,
        )
    except OSError:
        # read concatenated gzip files with extra data
        if isinstance(file, io.BytesIO):
            file.seek(0)
            data: memoryview | bytes = file.getbuffer()
        else:
            with open(file, "rb") as f:
                data = f.read()

        # https://stackoverflow.com/q/4928560
        # https://stackoverflow.com/a/37042747
        decompressor = zlib.decompressobj(31)

        # decompress data from first concatenated gzip file
        data = decompressor.decompress(data)

        # decompress data from second concatenated gzip file
        additional_data = zlib.decompress(decompressor.unused_data, 31)
        data += additional_data[33:]  # skip over second header

        new_file = io.BytesIO(data)

        df = pd.read_csv(
            new_file,
            skiprows=1,
            na_values="--",
            names=["rsid", "chrom", "pos", "genotype"],
            index_col=0,
            dtype=NORMALIZED_DTYPES,
            compression=None,  # already decompressed
        )

    return ParsedData(df)


def _parse_famfinder(parse_input: ParseInput) -> ParsedData:
    """Read and parse Family Tree DNA (FTDNA) "famfinder" file.

    https://www.familytreedna.com

    Parameters
    ----------
    parse_input : ParseInput
        file to parse

    Returns
    -------
    ParsedData
        parsed SNPs
    """
    file = parse_input.file
    compression = parse_input.compression

    df = pd.read_csv(
        file,
        comment="#",
        na_values="-",
        names=["rsid", "chrom", "pos", "allele1", "allele2"],
        index_col=0,
        dtype=TWO_ALLELE_DTYPES,
        compression=compression,
    )

    # create genotype column from allele columns
    df["genotype"] = df["allele1"] + df["allele2"]

    # delete allele columns
    # http://stackoverflow.com/a/13485766
    del df["allele1"]
    del df["allele2"]

    return ParsedData(df)


FTDNA = ReaderFormat(name="FTDNA", detect=_detect_ftdna, parse=_parse_ftdna)

FAMFINDER = ReaderFormat(
    name="ftdna_famfinder",
    detect=_detect_famfinder,
    parse=_parse_famfinder,
    source=_source_famfinder,
)
