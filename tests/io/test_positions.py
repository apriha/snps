"""Tests for how ``read_file`` reads positions.

The VCF, PLINK, CircleDNA and generic readers accept a position whose value is a
whole number from 0 to 4294967295, such as ``102`` or ``102.0``, and raise
``ParseError`` for any other position.
"""

import re

import numpy as np
import pandas as pd
import pytest

from snps.exceptions import ParseError
from snps.io.readers import read_file
from tests.support import FakeResources

# Each layout reads rs1 at 101 and rs2 at the position under test.
LAYOUTS = {
    "vcf": (
        ("vcf",),
        "##fileformat=VCFv4.2\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\n"
        "1\t101\trs1\tA\tG\t.\t.\t.\tGT\t0/1\n"
        "1\t{pos}\trs2\tC\tT\t.\t.\t.\tGT\t1/1\n",
    ),
    "plink": (
        ("PLINK",),
        "# Below is a text version of your data.\n"
        "rsid\tchromosome\tposition\tgenotype\n"
        "rs1\t1\t101\tAA\n"
        "rs2\t1\t{pos}\tCC\n",
    ),
    "circledna": (
        ("CircleDNA",),
        "# Circle\n#\n#\n# MARKERNAME\tCHROM\tPOS\tGT\n"
        "rs1\tchr1\t101\tA/A\n"
        "rs2\tchr1\t{pos}\tC/C\n",
    ),
    "generic-csv": (
        ("generic",),
        "rsid,chromosome,position,genotype\nrs1,1,101,AA\nrs2,1,{pos},CC\n",
    ),
    "generic-tsv": (
        ("generic",),
        "rsid\tchromosome\tposition\tgenotype\nrs1\t1\t101\tAA\nrs2\t1\t{pos}\tCC\n",
    ),
}

# position text -> whole number it reads as
ACCEPTED = {
    "102": 102,
    "102.0": 102,  # 102.0 = 102
    "102.000": 102,  # 102.000 = 102
    "1.02e2": 102,  # 1.02 * 100
    " 102 ": 102,  # spaces around the number are ignored
    "0": 0,  # lowest uint32
    "4294967295": 4294967295,  # 2**32 - 1, largest uint32
    "4294967295.0": 4294967295,  # 2**32 - 1
}

REJECTED = [
    "102.5",  # has a fractional part
    "4294967295.5",  # largest uint32 plus a fraction
    "-1",  # below 0
    "4294967296",  # 2**32, one past the largest uint32
    "abc",  # not a number
]

MESSAGE = "is not a whole number from 0 to 4294967295"


def read(layout, pos):
    _, template = LAYOUTS[layout]
    content = template.format(pos=pos).encode("utf-8")
    return read_file(content, resources=FakeResources())


@pytest.mark.parametrize("text", ACCEPTED, ids=repr)
@pytest.mark.parametrize("layout", LAYOUTS)
def test_accepted_position(layout, text):
    source, _ = LAYOUTS[layout]
    expected = ACCEPTED[text]

    data = read(layout, text)

    assert data.source == source
    assert data.frame.loc["rs2", "pos"] == expected
    assert data.frame["pos"].dtype == np.uint32
    pd.testing.assert_frame_equal(data.frame, read(layout, str(expected)).frame)


@pytest.mark.parametrize("text", REJECTED, ids=repr)
@pytest.mark.parametrize("layout", LAYOUTS)
def test_rejected_position(layout, text):
    with pytest.raises(ParseError, match=re.escape(MESSAGE)):
        read(layout, text)
