"""23Mofang raw data files."""

from __future__ import annotations

import numpy as np
import pandas as pd

from snps.constants import NORMALIZED_DTYPES
from snps.io.readers._registry import Header, ParsedData, ParseInput, ReaderFormat


def _detect(header: Header) -> bool:
    return "23Mofang" in header.first_line or "23Mofang" in header.comments


def _parse(parse_input: ParseInput) -> ParsedData:
    """Read and parse 23Mofang file.

    https://www.23mofang.com/

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

    columnnames = ["rsid", "chrom", "pos", "genotype"]
    dtype = NORMALIZED_DTYPES.copy()

    # Temporarily use nullable UInt32 for 'pos' column
    dtype["pos"] = pd.UInt32Dtype()  # type: ignore[assignment]  # read_csv takes dtype instances

    df = pd.read_csv(
        file,
        comment="#",
        sep="\t",
        na_values=["--", "-"],
        names=columnnames,
        compression=compression,
        dtype=dtype,
    )

    # Drop rows with NaN values in 'pos' column
    df = df.dropna(subset=["pos"])

    # Convert 'pos' column to np.uint32
    df["pos"] = df["pos"].astype(np.uint32)

    df = df.dropna(subset=["rsid", "chrom", "pos"])
    df = df.astype(dtype=NORMALIZED_DTYPES)
    df = df.set_index("rsid")
    return ParsedData(df)


FORMAT = ReaderFormat(name="23Mofang", detect=_detect, parse=_parse)
