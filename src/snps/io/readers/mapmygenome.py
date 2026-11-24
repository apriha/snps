"""Mapmygenome raw data files."""

from __future__ import annotations

import numpy as np
import pandas as pd

from snps.constants import NORMALIZED_DTYPES
from snps.io.readers._registry import Header, ParsedData, ParseInput, ReaderFormat


def _detect(header: Header) -> bool:
    return (
        "SNP Name\trsID" in header.first_line
        or "SNP.Name\tSample.ID" in header.first_line
        or "Mapmygenome" in header.comments
    )


def _parse(parse_input: ParseInput) -> ParsedData:
    """Read and parse Mapmygenome file.

    https://mapmygenome.in

    Parameters
    ----------
    parse_input : ParseInput
        file to parse

    Returns
    -------
    ParsedData
        parsed SNPs
    """
    source = parse_input.source
    header = parse_input.header.first_line
    comments = parse_input.header.comments

    if "rsid\tchromosome\tposition\tgenotype" in header or "Mapmygenome" in comments:
        with source.open() as file:
            df = pd.read_csv(
                file,
                comment="#",
                sep="\t",
                na_values="--",
                header=0,
                dtype={
                    "rsid": object,
                    "chromosome": object,
                    "position": np.uint32,
                    "genotype": object,
                },
            )
        df = df.rename(columns={"chromosome": "chrom", "position": "pos"})
        df = df.astype(dtype=NORMALIZED_DTYPES)
        df = df.set_index("rsid")
        return ParsedData(df)
    else:

        def parse(rsid_col_name: str, rsid_col_idx: int) -> pd.DataFrame:
            with source.open() as file:
                return pd.read_csv(
                    file,
                    comment="#",
                    sep="\t",
                    na_values="--",
                    header=0,
                    index_col=rsid_col_idx,
                    dtype={
                        rsid_col_name: object,
                        "Chr": object,
                        "Position": np.uint32,
                        "Allele1...Top": object,
                        "Allele2...Top": object,
                    },
                )

        if "rsID" in header:
            df = parse("rsID", 1)
        else:
            df = parse("SNP.Name", 0)

        df["genotype"] = df["Allele1...Plus"] + df["Allele2...Plus"]
        df = df.rename(columns={"Chr": "chrom", "Position": "pos"})
        df.index.name = "rsid"
        df = df[["chrom", "pos", "genotype"]]

        return ParsedData(df)


FORMAT = ReaderFormat(name="Mapmygenome", detect=_detect, parse=_parse)
