"""CircleDNA raw data files."""

from __future__ import annotations

import pandas as pd

from snps.constants import NORMALIZED_DTYPES
from snps.io.readers._common import parse_positions
from snps.io.readers._registry import Header, ParsedData, ParseInput, ReaderFormat


def _detect(header: Header) -> bool:
    return "Circle" in header.first_line


def _parse(parse_input: ParseInput) -> ParsedData:
    """Read and parse CircleDNA file.

    https://circledna.com/

    Notes
    -----
    This function attempts to read and parse a whole exome file, optionally compressed
    with gzip or zip. Some assumptions are made throughout this process:

        * SNPs that are not annotated with an RSID are skipped
        * Insertions and deletions are skipped

    Parameters
    ----------
    parse_input : ParseInput
        file to parse

    Returns
    -------
    ParsedData
        parsed SNPs
    """
    rs_chunks = []
    with (
        parse_input.source.open() as file,
        pd.read_csv(
            file,
            comment="#",
            sep="\t",
            chunksize=10000,
            names=["rsid", "chrom", "pos", "genotype"],
        ) as reader,
    ):
        for chunk in reader:
            # filter for SNPs with rsids
            tmp = chunk.loc[
                (chunk.rsid.str.startswith("rs")) & (chunk.genotype.str.len() == 3)
            ]
            if len(tmp) > 0:
                rs_chunks.append(tmp)

    df = pd.concat(rs_chunks)
    df["chrom"] = df.chrom.str[3:]
    df["genotype"] = df.genotype.apply(lambda x: "".join(x.split("/")))
    df["pos"] = parse_positions(df["pos"])
    df = df.astype(NORMALIZED_DTYPES)
    df = df.set_index("rsid")

    return ParsedData(df)


FORMAT = ReaderFormat(name="CircleDNA", detect=_detect, parse=_parse)
