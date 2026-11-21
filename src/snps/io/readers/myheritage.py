"""MyHeritage raw data files."""

from __future__ import annotations

import io

import pandas as pd

from snps.constants import NORMALIZED_DTYPES
from snps.io.readers._common import handle_bytes_data
from snps.io.readers._registry import Header, ParsedData, ParseInput, ReaderFormat


def _detect(header: Header) -> bool:
    return "MyHeritage" in header.first_line


def _parse(parse_input: ParseInput) -> ParsedData:
    """Read and parse MyHeritage file.

    https://www.myheritage.com

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

    if isinstance(file, str):
        with open(file, "rb") as f:
            first_line, comments, data, comrpession = handle_bytes_data(
                f.read(), include_data=True
            )
    else:
        first_line, comments, data, compression = handle_bytes_data(
            file.read(), include_data=True
        )

    file_string_in = io.StringIO(data)
    file_string_out = io.StringIO()
    for line in file_string_in:
        # user the number of quotes in a line to tell old from new
        if line.count('"') == 14:
            # extra-quoted new variant file
            # can all be stripped so pandas can read it normally
            line = line.replace('"', "")
            # take it apart and put it back together so it looks
            # like the older MyHeritage files
            line = '"' + '","'.join(line.strip().split(",")) + '"\n'
        file_string_out.write(line)

    return ParsedData(
        pd.read_csv(
            io.StringIO(file_string_out.getvalue()),
            comment="#",
            header=0,
            na_values="--",
            names=["rsid", "chrom", "pos", "genotype"],
            index_col=0,
            dtype=NORMALIZED_DTYPES,
        )
    )


FORMAT = ReaderFormat(name="MyHeritage", detect=_detect, parse=_parse)
