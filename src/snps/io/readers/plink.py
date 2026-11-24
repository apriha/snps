"""PLINK text files."""

from __future__ import annotations

import pandas as pd

from snps.constants import NORMALIZED_DTYPES
from snps.io.readers._registry import Header, ParsedData, ParseInput, ReaderFormat


def _detect(header: Header) -> bool:
    return "# Below is a text version of your data." in header.comments


def _parse(parse_input: ParseInput) -> ParsedData:
    """Read and parse plink file.

    Parameters
    ----------
    parse_input : ParseInput
        file to parse

    Returns
    -------
    ParsedData
        parsed SNPs
    """
    columnnames = ["rsid", "chrom", "pos", "genotype"]
    with parse_input.source.open() as file:
        df = pd.read_csv(
            file,
            comment="#",
            sep="\t",
            na_values=["--", "-"],
            names=columnnames,
        )
    # turn number numbers into string numbers
    df["chrom"] = df["chrom"].map(
        {
            "1": "1",
            "2": "2",
            "3": "3",
            "4": "4",
            "5": "5",
            "6": "6",
            "7": "7",
            "8": "8",
            "9": "9",
            "10": "10",
            "11": "11",
            "12": "12",
            "13": "13",
            "14": "14",
            "15": "15",
            "16": "16",
            "17": "17",
            "18": "18",
            "19": "19",
            "20": "20",
            "21": "21",
            "22": "22",
            1: "1",
            2: "2",
            3: "3",
            4: "4",
            5: "5",
            6: "6",
            7: "7",
            8: "8",
            9: "9",
            10: "10",
            11: "11",
            12: "12",
            13: "13",
            14: "14",
            15: "15",
            16: "16",
            17: "17",
            18: "18",
            19: "19",
            20: "20",
            21: "21",
            22: "22",
            "X": "X",
            "Y": "Y",
            "XY": "Y",
            "MT": "MT",
        }
    )
    df = df.dropna(subset=["rsid", "chrom", "pos"])
    df = df.astype(dtype=NORMALIZED_DTYPES)
    df = df.set_index("rsid")
    return ParsedData(df)


FORMAT = ReaderFormat(name="PLINK", detect=_detect, parse=_parse)
