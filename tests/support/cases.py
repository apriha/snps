"""Raw data files of each vendor format, with what reading each one gives."""

from collections.abc import Callable
from dataclasses import dataclass
from functools import partial
from pathlib import Path

import pandas as pd
from pandas.api.types import is_object_dtype, is_string_dtype, is_unsigned_integer_dtype

from snps.testing import (
    SNPsTestMixin,
    assert_frame_equal_with_string_index,
    create_snp_df,
)

INPUT = Path(__file__).resolve().parent.parent / "input"


@dataclass(frozen=True, eq=False)
class ReadCase:
    """A raw data file and what reading it gives.

    Attributes
    ----------
    id : str
        name of the case in test ids
    file : str or callable
        name of a fixture in ``tests/input``, or a function returning the content
        of a generated file
    source : str
        expected source
    expected : callable
        function returning the expected SNPs
    phased : bool
        whether the SNPs are expected to be phased
    build : int
        expected build
    build_detected : bool
        whether the build is expected to be detected
    """

    id: str
    file: str | Callable[[], bytes]
    source: str
    expected: Callable[[], pd.DataFrame]
    phased: bool = False
    build: int = 37
    build_detected: bool = False

    def path(self, directory: Path) -> Path:
        """Return the path of the file, writing a generated file to ``directory``."""
        if isinstance(self.file, str):
            return INPUT / self.file
        path = directory / f"{self.id}.txt"
        path.write_bytes(self.file())
        return path


def generic_snps() -> pd.DataFrame:
    """SNPs of most fixtures: rs1 to rs8 on chromosome 1, rs5 not called."""
    return SNPsTestMixin().generic_snps()


def homozygous_snps(count: int) -> pd.DataFrame:
    """SNPs rs1 to ``rs<count>`` on chromosome 1 from position 101, all AA."""
    return create_snp_df(
        rsid=[f"rs{1 + i}" for i in range(count)],
        chrom=["1"] * count,
        pos=[101 + i for i in range(count)],
        genotype=["AA"] * count,
    )


def _called_generic_snps() -> pd.DataFrame:
    return generic_snps().drop("rs5")


def _with_comment(fixture: str, comment: bytes) -> bytes:
    """Return a fixture with a comment line inserted after its second line."""
    lines = (INPUT / fixture).read_bytes().splitlines(keepends=True)
    return b"".join([*lines[:2], comment, *lines[2:]])


def _ancestry_extra_tab() -> bytes:
    """AncestryDNA rows, the first of which has an extra tab separator."""
    s = "#AncestryDNA\r\n"
    s += "rsid\tchromosome\tposition\tallele1\tallele2\r\n"
    s += "rs1\t1\t101\t\tA\tA\r\n"
    for i in range(1, 100):
        s += f"rs{1 + i}\t1\t{101 + i}\tA\tA\r\n"
    return s.encode()


def _ftdna_second_header() -> bytes:
    """FTDNA rows with a second header line, as in concatenated files."""
    s = "RSID,CHROMOSOME,POSITION,RESULT\n"
    for i in range(0, 100):
        s += f'"rs{1 + i}","1","{101 + i}","AA"\n'
    s += "RSID,CHROMOSOME,POSITION,RESULT\n"
    for i in range(0, 10):
        s += f'"rs{100 + 1 + i}","1","{100 + 101 + i}","AA"\n'
    return s.encode()


CASES = (
    ReadCase("23andme", "23andme.txt", "23andMe", generic_snps),
    # some files have split allele columns
    ReadCase("23andme_allele", "23andme_allele.txt", "23andMe", generic_snps),
    # windows line endings
    ReadCase("23andme_win", "23andme_win.txt", "23andMe", generic_snps),
    ReadCase(
        "23andme_build36",
        partial(_with_comment, "23andme.txt", b"# build 36\n"),
        "23andMe",
        generic_snps,
        build=36,
        build_detected=True,
    ),
    ReadCase(
        "23andme_build37",
        partial(_with_comment, "23andme.txt", b"# build 37\n"),
        "23andMe",
        generic_snps,
        build=37,
        build_detected=True,
    ),
    ReadCase("ancestry", "ancestry.txt", "AncestryDNA", generic_snps),
    ReadCase(
        "ancestry_extra_tab",
        _ancestry_extra_tab,
        "AncestryDNA",
        partial(homozygous_snps, 100),
    ),
    ReadCase(
        "ancestry_multi_sep", "ancestry_multi_sep.txt", "AncestryDNA", generic_snps
    ),
    ReadCase("codigo46", "codigo46.txt", "Codigo46", generic_snps),
    ReadCase("tellmeGen", "tellmeGen.txt", "tellmeGen", generic_snps),
    ReadCase("DNALand", "DNALand.txt", "DNA.Land", generic_snps),
    # only called genotypes
    ReadCase("circledna", "circledna.txt", "CircleDNA", _called_generic_snps),
    ReadCase("plink", "plink.txt", "PLINK", generic_snps),
    ReadCase("SelfDecode", "selfdecode.txt", "SelfDecode", generic_snps),
    ReadCase("23Mofang", "23mofang.txt", "23Mofang", generic_snps),
    ReadCase("ftdna", "ftdna.csv", "FTDNA", generic_snps),
    ReadCase("ftdna_famfinder", "ftdna_famfinder.csv", "FTDNA", generic_snps),
    ReadCase(
        "ftdna_second_header",
        _ftdna_second_header,
        "FTDNA",
        partial(homozygous_snps, 110),
    ),
    ReadCase("generic_csv", "generic.csv", "generic", generic_snps),
    ReadCase("generic_tsv", "generic.tsv", "generic", generic_snps),
    ReadCase(
        "generic_extra_column", "generic_extra_column.tsv", "generic", generic_snps
    ),
    ReadCase(
        "generic_header_comment", "generic_header_comment.tsv", "generic", generic_snps
    ),
    ReadCase("generic_multi_rsid", "generic_multi_rsid.tsv", "generic", generic_snps),
    ReadCase("generic_no_header", "generic_no_header.tsv", "generic", generic_snps),
    ReadCase(
        "generic_non_standard_columns",
        "generic_non_standard_columns.tsv",
        "generic",
        generic_snps,
    ),
    ReadCase("genes_for_good", "genesforgood.txt", "GenesForGood", generic_snps),
    ReadCase("livingdna", "livingdna.csv", "LivingDNA", generic_snps),
    ReadCase("mapmygenome", "mapmygenome.txt", "Mapmygenome", generic_snps),
    ReadCase(
        "mapmygenome_alt_header",
        "mapmygenome_alt_header.txt",
        "Mapmygenome",
        generic_snps,
    ),
    ReadCase(
        "mapmygenome_new_format",
        "mapmygenome_new_format.txt",
        "Mapmygenome",
        generic_snps,
        build_detected=True,
    ),
    ReadCase("myheritage", "myheritage.csv", "MyHeritage", generic_snps),
    ReadCase(
        "myheritage_extra_quotes",
        "myheritage_extra_quotes.csv",
        "MyHeritage",
        generic_snps,
    ),
    ReadCase("sano", "sano.txt", "Sano", generic_snps),
    ReadCase("sano_dtc", "sano_dtc.txt", "Sano", generic_snps),
    ReadCase("sano_dtc_no_comments", "sano_dtc_no_comments.txt", "Sano", generic_snps),
)


def assert_normalized_frame(df: pd.DataFrame) -> None:
    """Assert that a frame has the index and column dtypes of normalized SNPs."""
    assert df.index.name == "rsid"
    # Accept both object dtype and StringDtype (used in Python 3.14+)
    assert is_string_dtype(df.index.dtype)
    assert is_string_dtype(df.chrom.dtype) or is_object_dtype(df.chrom.dtype)
    assert is_unsigned_integer_dtype(df.pos.dtype)
    assert is_string_dtype(df.genotype.dtype) or is_object_dtype(df.genotype.dtype)


def assert_matches_case(snps, case: ReadCase) -> None:
    """Assert that the ``SNPs`` read from the file of a case are what it expects."""
    assert snps.source == case.source
    assert_frame_equal_with_string_index(snps.snps, case.expected(), check_exact=True)
    assert snps.phased is case.phased
    assert snps.build == case.build
    assert snps.build_detected is case.build_detected
    assert_normalized_frame(snps.snps)
