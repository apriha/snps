"""Illumina Global Screening Array files, including Sano and Codigo46 files."""

from __future__ import annotations

import io

import numpy as np
import pandas as pd

from snps.constants import NORMALIZED_DTYPES
from snps.exceptions import ParseError
from snps.io.readers._common import extract_comments
from snps.io.readers._registry import Header, ParsedData, ParseInput, ReaderFormat


def _detect_gsa(header: Header) -> bool:
    return header.first_line.startswith("[Header]")


def _detect_sano(header: Header) -> bool:
    return _detect_gsa(header) and "SANO" in header.comments


def _detect_codigo46(header: Header) -> bool:
    return _detect_gsa(header) and "CODIGO46" in header.comments


def _parse(parse_input: ParseInput) -> ParsedData:
    """Read and parse Illumina Global Screening Array files.

    Parameters
    ----------
    parse_input : ParseInput
        file to parse

    Returns
    -------
    ParsedData
        parsed SNPs

    Raises
    ------
    ParseError
        if the data section has a column named rsid, chrom, pos or genotype
    ValueError
        if the file has no supported allele columns
    """
    resources = parse_input.resources

    # read the comments so we get to the actual data
    with parse_input.source.open() as file:
        _, _, data = extract_comments(file, include_data=True)

    # turn the data into a pandas dataframe for manipulation
    df = pd.read_csv(
        io.StringIO(data),
        sep="\t",
        engine="c",
        dtype={
            "Position": NORMALIZED_DTYPES["pos"],
            "Chr": NORMALIZED_DTYPES["chrom"],
        },
    )

    # reserve columns we want out
    reserved = [c for c in ("rsid", "chrom", "pos", "genotype") if c in df.columns]
    if reserved:
        raise ParseError(
            f"The data section has columns reserved by snps: {', '.join(reserved)}; "
            "rename or remove those columns."
        )

    # prefer the specified chromosome and position, in prsent
    # this uses the gsa names, so needs to be done before rsid
    if "Chr" in df.columns and "Position" in df.columns:
        # put the chromosome in the right column with the right type
        df["chrom"] = df["Chr"]
        # put the position in the right column with the right type
        df["pos"] = df["Position"]

    else:
        # use an external source to map snp names to chromosome and position
        df = df.merge(
            resources.get_gsa_chrpos(),  # type: ignore[no-untyped-call]  # snps.resources is untyped
            how="inner",  # inner join, everything must have a position
            left_on="SNP Name",
            right_on="gsaname_chrpos",
            suffixes=(None, "_gsa"),
        )
        # make sure its the right types
        df["chrom"] = df["gsachr"].apply(str).astype(NORMALIZED_DTYPES["chrom"])
        df["pos"] = df["gsapos"].astype(NORMALIZED_DTYPES["pos"])

    # use the given rsid when avaliable, SNP Name when unavaliable
    df["rsid"] = df["SNP Name"].astype(NORMALIZED_DTYPES["rsid"])
    if "RsID" in df.columns:
        df.loc[df["RsID"] != ".", "rsid"] = df.loc[df["rsid"] != ".", "RsID"]
    else:
        # if given RSIDs are not avaliable, then use the external mapping to turn
        # SNP names into rsids where possible
        df = df.merge(
            resources.get_gsa_rsid(),  # type: ignore[no-untyped-call]  # snps.resources is untyped
            how="left",  # left-hand join, gsa rsids may be NA
            left_on="SNP Name",
            right_on="gsaname_rsid",
            suffixes=(None, "_gsa_rsid"),
        )
        df.loc[~pd.isna(df["gsaname_rsid"]), "rsid"] = df.loc[
            ~pd.isna(df["gsaname_rsid"]), "gsarsid"
        ]

    # combine the alleles into genotype
    # prefer Plus strand as that is forward reference
    if "Allele1 - Plus" in df.columns and "Allele2 - Plus" in df.columns:
        df["genotype"] = (df["Allele1 - Plus"] + df["Allele2 - Plus"]).astype(
            NORMALIZED_DTYPES["genotype"]
        )
    elif "Allele1 - Forward" in df.columns and "Allele2 - Forward" in df.columns:
        # if strand is forward, need to take reverse complement of *some* rsids
        # this is because it is Illumina forward, which is dbSNP strand, which
        # is reverse reference for some RSIDs before dbSNP 151.

        # load list of reversable rsids
        dbsnp151 = resources.get_dbsnp_151_37_reverse()  # type: ignore[no-untyped-call]  # snps.resources is untyped
        # keep only the rsids
        dbsnp151 = dbsnp151.filter(items=("dbsnp151revrsid",), axis=1)

        # add it as an extra column
        df = df.merge(
            dbsnp151,
            how="left",
            left_on="rsid",
            right_on="dbsnp151revrsid",
            suffixes=(None, "_dbsnp151rev"),
        )

        # create plus strand columns from the forward alleles and flip them if appropriate
        for i in (1, 2):
            df[f"Allele{i} - Plus"] = df[f"Allele{i} - Forward"]
            df.loc[
                (df[f"Allele{i} - Forward"] == "A") & (~pd.isna(df["dbsnp151revrsid"])),
                f"Allele{i} - Plus",
            ] = "T"
            df.loc[
                (df[f"Allele{i} - Forward"] == "T") & (~pd.isna(df["dbsnp151revrsid"])),
                f"Allele{i} - Plus",
            ] = "A"
            df.loc[
                (df[f"Allele{i} - Forward"] == "C") & (~pd.isna(df["dbsnp151revrsid"])),
                f"Allele{i} - Plus",
            ] = "G"
            df.loc[
                (df[f"Allele{i} - Forward"] == "G") & (~pd.isna(df["dbsnp151revrsid"])),
                f"Allele{i} - Plus",
            ] = "C"

        # create a genotype by combining the new plus columns
        df["genotype"] = (df["Allele1 - Plus"] + df["Allele2 - Plus"]).astype(
            NORMALIZED_DTYPES["genotype"]
        )
    else:
        raise ValueError("No supported allele column")

    # mark -- genotype as na
    df.loc[df["genotype"] == "--", "genotype"] = np.nan

    # keep only the columns we want
    df = df.filter(items=("rsid", "chrom", "pos", "genotype"), axis=1)

    # discard rows without values
    df = df.dropna(subset=["rsid", "chrom", "pos"])

    # reindex for the new identifiers
    df = df.set_index(["rsid"])

    return ParsedData(df)


SANO = ReaderFormat(name="Sano", detect=_detect_sano, parse=_parse)

CODIGO46 = ReaderFormat(name="Codigo46", detect=_detect_codigo46, parse=_parse)

GSA = ReaderFormat(name="GSA", detect=_detect_gsa, parse=_parse)
