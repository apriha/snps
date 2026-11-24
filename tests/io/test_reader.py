import os
import tempfile
from functools import partial

import pytest
from atomicwrites import atomic_write

from snps import SNPs
from snps.io.readers import read_file
from snps.utils import gzip_file
from tests import BaseSNPsTestCase
from tests.support.cases import CASES, ReadCase, assert_matches_case, homozygous_snps
from tests.support.inputs import (
    INPUT_IDS,
    INPUTS,
    assert_same_data,
    given_as,
    package,
)


@pytest.mark.parametrize("case", CASES, ids=lambda case: case.id)
@pytest.mark.parametrize(("kind", "packaging"), INPUTS, ids=INPUT_IDS)
def test_read(tmp_path, case, kind, packaging):
    plain = case.path(tmp_path)
    packed = package(plain, packaging, tmp_path / "packed")
    with given_as(packed, kind) as file:
        assert_matches_case(SNPs(file), case)
    with given_as(packed, kind) as file:
        data = read_file(file)
    assert_same_data(data, read_file(plain))


@pytest.mark.parametrize("kind", ["str", "bytes", "file"])
def test_read_ftdna_concat_gzip_extra_data(tmp_path, kind):
    # https://www.familytreedna.com

    total_snps1 = 10
    total_snps2 = 10
    # generate content of first file
    s1 = "RSID,CHROMOSOME,POSITION,RESULT\r\n"
    for i in range(0, total_snps1):
        s1 += f'"rs{1 + i}","1","{101 + i}","AA"\r\n'

    # generate content of second file
    s2 = "RSID,CHROMOSOME,POSITION,RESULT\r\n"
    for i in range(0, total_snps2):
        s2 += f'"rs{total_snps1 + 1 + i}","1","{total_snps1 + 101 + i}","AA"\r\n'

    file1 = tmp_path / "ftdna_concat_gzip1.csv"
    file2 = tmp_path / "ftdna_concat_gzip2.csv"
    path = tmp_path / "ftdna_concat_gzip.csv.gz"

    # write and compress individual files
    file1.write_bytes(s1.encode())
    file2.write_bytes(s2.encode())
    gzip_file(str(file1), f"{file1}.gz")
    gzip_file(str(file2), f"{file2}.gz")

    # concatenate gzips and add extra data
    data = (tmp_path / f"{file1.name}.gz").read_bytes()
    data += (tmp_path / f"{file2.name}.gz").read_bytes()
    data += b"extra data"
    path.write_bytes(data)

    case = ReadCase(
        "ftdna_concat_gzip_extra_data",
        path.read_bytes,
        "FTDNA",
        partial(homozygous_snps, total_snps1 + total_snps2),
    )
    with given_as(path, kind) as file:
        assert_matches_case(SNPs(file), case)


class TestReader(BaseSNPsTestCase):
    def run_build_detection_test(
        self,
        run_parsing_tests_func,
        build_str,
        build_int,
        file="tests/input/testvcf.vcf",
        source="vcf",
        comment_str="##{}\n",
        insertion_line=1,
    ):
        with tempfile.TemporaryDirectory() as tmpdir:
            s = ""
            with open(file, "r") as f:
                for i, line in enumerate(f.readlines()):
                    s += line
                    # insert comment from which to detect build
                    if i == insertion_line:
                        s += comment_str.format(build_str)

            file_build_comment = os.path.join(tmpdir, os.path.basename(file))
            with atomic_write(file_build_comment, mode="w") as f:
                f.write(s)

            run_parsing_tests_func(
                file_build_comment, source, build=build_int, build_detected=True
            )

    def test_read_vcf(self):
        self.run_parsing_tests_vcf("tests/input/testvcf.vcf")

    def test_read_vcf_NCBI36(self):
        self.run_build_detection_test(
            self.run_parsing_tests_vcf,
            "NCBI36",
            36,
            comment_str="##contig=<ID=1,assembly={}>\n",
        )

    def test_read_vcf_b37(self):
        self.run_build_detection_test(
            self.run_parsing_tests_vcf,
            "b37",
            37,
            comment_str="##contig=<ID=1,assembly={}>\n",
        )

    def test_read_vcf_hg19(self):
        self.run_build_detection_test(
            self.run_parsing_tests_vcf,
            "hg19",
            37,
            comment_str="##contig=<ID=1,assembly={}>\n",
        )

    def test_read_vcf_hg38(self):
        self.run_build_detection_test(
            self.run_parsing_tests_vcf,
            "hg38",
            38,
            comment_str="##contig=<ID=1,assembly={}>\n",
        )

    def test_read_vcf_GRCh38(self):
        self.run_build_detection_test(
            self.run_parsing_tests_vcf,
            "GRCh38",
            38,
            comment_str="##contig=<ID=1,assembly={}>\n",
        )

    def test_read_vcf_build38(self):
        self.run_build_detection_test(
            self.run_parsing_tests_vcf,
            "Build 38",
            38,
            comment_str="##contig=<ID=1,assembly={}>\n",
        )

    def test_read_vcf_b38(self):
        self.run_build_detection_test(
            self.run_parsing_tests_vcf,
            "b38",
            38,
            comment_str="##contig=<ID=1,assembly={}>\n",
        )

    def test_read_vcf_length38(self):
        self.run_build_detection_test(
            self.run_parsing_tests_vcf,
            "",
            38,
            comment_str="##contig=<ID=1,length=248956422>\n",
        )

    def test_read_vcf_length37(self):
        self.run_build_detection_test(
            self.run_parsing_tests_vcf,
            "",
            37,
            comment_str="##contig=<ID=1,length=249250621>\n",
        )

    def test_read_vcf_multi_sample(self):
        self.run_parsing_tests_vcf("tests/input/testvcf_multi_sample.vcf")

    def test_read_vcf_phased(self):
        self.run_parsing_tests_vcf("tests/input/testvcf_phased.vcf", phased=True)

    def test_read_vcf_mixed_phase(self):
        expected = self.mixed_phase_snps()
        self.run_parsing_tests_vcf(
            "tests/input/testvcf_mixed_phase.vcf",
            snps_df=expected.drop(columns="phased"),
        )
        s = SNPs("tests/input/testvcf_mixed_phase.vcf")
        self.assert_frame_equal_with_string_index(s._snps, expected, check_exact=True)
        self.assertEqual(s._snps["phased"].dtype, bool)

    def test_read_no_calls_unphased(self):
        vcf = SNPs(
            b"##fileformat=VCFv4.1\n"
            b"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE\n"
            b"1\t101\trs1\tA\tT\t.\t.\t.\tGT\t0|1\n"
            b"1\t102\trs2\tA\tT\t.\t.\t.\tGT\t.|.\n"
            b"1\t103\trs3\tA\tT\t.\t.\t.\tGT\t./.\n"
            b"1\t104\trs4\tA\tT\t.\t.\t.\tGT\t.\n"
        )
        csv = SNPs(
            b"# Generated by snps\n# Phased: True\n"
            b"rsid,chromosome,position,genotype\nrs1,1,101,AT\nrs2,1,102,--\n"
        )
        self.assertListEqual(vcf._snps["phased"].tolist(), [True, False, False, False])
        self.assertListEqual(csv._snps["phased"].tolist(), [True, False])

    def test_read_vcf_rsids(self):
        self.run_parsing_tests_vcf("tests/input/testvcf.vcf", rsids=["rs1", "rs2"])

    def test_read_unannotated_vcf(self):
        self.run_parsing_tests_vcf(
            "tests/input/unannotated_testvcf.vcf", "vcf", unannotated=True, build=0
        )

    def test_read_vcf_chr_prefix(self):
        self.run_parsing_tests_vcf("tests/input/testvcf_chr_prefix.vcf")
