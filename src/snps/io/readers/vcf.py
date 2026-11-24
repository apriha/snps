"""VCF files."""

from __future__ import annotations

import io
import logging
from collections.abc import Set
from typing import IO

import numpy as np
import pandas as pd

from snps.constants import NORMALIZED_DTYPES
from snps.io.readers._registry import Header, ParsedData, ParseInput, ReaderFormat

logger = logging.getLogger(__name__)


def _detect(header: Header) -> bool:
    return "vcf" in header.comments.lower() or "##contig" in header.comments.lower()


def _source(header: Header) -> str:
    """Return the company named by a ``##detectedCompany`` line, else ``"vcf"``."""
    provider = "vcf"
    header_lines = header.comments.replace("\r\n", "\n").split("\n")

    detected_company = ""
    for line in header_lines:
        if line.startswith("##detectedCompany="):
            detected_company = line.split("=")[1].strip('"')
            break

    if detected_company:
        provider = detected_company

    return provider


def _parse(parse_input: ParseInput) -> ParsedData:
    """Read and parse VCF file.

    Notes
    -----
    This function attempts to read and parse a VCF file or buffer, optionally
    compressed with gzip or zip. Some assumptions are made throughout this process:

        * SNPs that are not annotated with an RSID are skipped
        * If the VCF contains multiple samples, only the first sample is used to
          lookup the genotype
        * Precise insertions and deletions are skipped
        * If a sample allele is not specified, the genotype is reported as NaN
        * If a sample allele refers to a REF or ALT allele that is not specified,
          the genotype is reported as NaN

    Parameters
    ----------
    parse_input : ParseInput
        file to parse, and the rsids to extract (all if empty)

    Returns
    -------
    ParsedData
        parsed SNPs, with a ``phased`` column
    """
    rsids = parse_input.rsids

    with parse_input.source.open() as f:
        df, phased = _parse_vcf(f, rsids)

    return ParsedData(df, phased)


def _parse_vcf(buffer: IO[bytes], rsids: Set[str]) -> tuple[pd.DataFrame, bool]:
    """Parse the records of a VCF file."""
    rows = []
    phased = True

    logged_multi_sample = False

    with io.TextIOWrapper(buffer) as file:
        for line in file:
            line_strip = line.strip("\n")

            # skip blank lines
            if not line_strip:
                continue

            # skip comment lines
            if line_strip.startswith("#"):
                continue

            rsid = line_strip.split("\t")[2]

            # skip SNPs with missing rsIDs.
            if rsid == ".":
                continue
            if rsids:
                if rsid not in rsids:
                    continue

            line_split = line_strip.split("\t")

            # snps does not yet support multi-sample vcf.
            if not logged_multi_sample and len(line_split) > 10:
                logger.info("Multiple samples detected in the vcf file")
                logged_multi_sample = True

            ref = line_split[3]
            alt = line_split[4]
            if len(alt.split(",")) > 1 and alt.split(",")[1] == "<NON_REF>":
                alt = alt.split(",")[0]

            # Handle <INS> and <DEL> alleles (imprecise insertions and deletions)
            ref_alt = [ref] + [
                "I" if allele == "<INS>" else "D" if allele == "<DEL>" else allele
                for allele in alt.split(",")
            ]

            # skip precise insertions and deletions
            if sum(map(len, ref_alt)) > len(ref_alt):
                continue

            # GT (genotype) is usually the first sample-specific field
            # | = diploid phased
            # / = diploid unphased
            # or haploid e.g. male sex chromosome
            genotype: str | float = ""
            zygote = line_split[9]
            zygote = zygote.split(":")[0]
            for z in zygote.replace("|", "/").split("/"):
                if z == ".":
                    # missing genotype
                    genotype = np.nan
                    break
                allele_index = int(z)
                if allele_index >= len(ref_alt):
                    # invalid genotype number
                    genotype = np.nan
                    break
                elif ref_alt[allele_index] == ".":
                    # missing genotype in ref or alt
                    genotype = np.nan
                    break
                else:
                    assert isinstance(genotype, str)
                    genotype = genotype + ref_alt[allele_index]

            if "/" in zygote and pd.notna(genotype):
                phased = False

            record_array = [
                rsid,
                f"{line_split[0]}".strip("chr"),
                line_split[1],
                genotype,
                # phased: a call with "|" or a single allele
                "/" not in zygote and pd.notna(genotype),
            ]
            rows.append(record_array)

        if len(rows) == 0:
            phased = False

        df = pd.DataFrame(rows, columns=["rsid", "chrom", "pos", "genotype", "phased"])
        df = df.astype({**NORMALIZED_DTYPES, "phased": bool})

        df = df.set_index("rsid", drop=True)

    return (df, phased)


FORMAT = ReaderFormat(name="vcf", detect=_detect, parse=_parse, source=_source)
