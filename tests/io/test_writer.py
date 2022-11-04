"""
BSD 3-Clause License

Copyright (c) 2019, Andrew Riha
All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this
   list of conditions and the following disclaimer.

2. Redistributions in binary form must reproduce the above copyright notice,
   this list of conditions and the following disclaimer in the documentation
   and/or other materials provided with the distribution.

3. Neither the name of the copyright holder nor the names of its
   contributors may be used to endorse or promote products derived from
   this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

"""

import os
import tempfile

import numpy as np
import pandas as pd

from snps import SNPs
from snps.resources import Resources, ReferenceSequence
from snps.utils import gzip_file
from tests import BaseSNPsTestCase


class TestWriter(BaseSNPsTestCase):
    def run_writer_test(
        self, func_str, filename="", output_file="", expected_output="", **kwargs
    ):
        if func_str == "to_vcf":
            with tempfile.TemporaryDirectory() as tmpdir1:
                s = SNPs("tests/input/testvcf.vcf", output_dir=tmpdir1)

                r = Resources()
                r._reference_sequences["GRCh37"] = {}

                output = os.path.join(tmpdir1, output_file)
                with tempfile.TemporaryDirectory() as tmpdir2:
                    dest = os.path.join(tmpdir2, "generic.fa.gz")
                    gzip_file("tests/input/generic.fa", dest)

                    seq = ReferenceSequence(ID="1", path=dest)

                    r._reference_sequences["GRCh37"]["1"] = seq

                    if not filename:
                        result = s.to_vcf(**kwargs)
                    else:
                        result = s.to_vcf(filename, **kwargs)

                    self.assertEqual(result, output)

                    if expected_output:
                        # read result
                        with open(output, "r") as f:
                            actual = f.read()

                        # read expected result
                        with open(expected_output, "r") as f:
                            expected = f.read()

                        self.assertIn(expected, actual)

                self.run_parsing_tests_vcf(output)
        else:
            with tempfile.TemporaryDirectory() as tmpdir:
                snps = SNPs("tests/input/generic.csv", output_dir=tmpdir)
                dest = os.path.join(tmpdir, output_file)
                # https://stackoverflow.com/a/3071
                if not filename:
                    self.assertEqual(getattr(snps, func_str)(), dest)
                else:
                    self.assertEqual(getattr(snps, func_str)(filename), dest)
                self.run_parsing_tests(dest, "generic")

    def test_to_csv(self):
        self.run_writer_test("to_csv", output_file="generic_GRCh37.csv")

    def test_to_csv_filename(self):
        self.run_writer_test(
            "to_csv", filename="generic.csv", output_file="generic.csv"
        )

    def test_to_tsv(self):
        self.run_writer_test("to_tsv", output_file="generic_GRCh37.txt")

    def test_to_tsv_filename(self):
        self.run_writer_test(
            "to_tsv", filename="generic.txt", output_file="generic.txt"
        )

    def test_to_vcf(self):
        self.run_writer_test(
            "to_vcf",
            output_file="vcf_GRCh37.vcf",
            expected_output="tests/output/vcf_generic.vcf",
        )

    def test_to_vcf_filename(self):
        self.run_writer_test("to_vcf", filename="vcf.vcf", output_file="vcf.vcf")

    def test_to_vcf_chrom_prefix(self):
        self.run_writer_test(
            "to_vcf",
            output_file="vcf_GRCh37.vcf",
            expected_output="tests/output/vcf_chrom_prefix_chr.vcf",
            chrom_prefix="chr",
        )

    def test_save_snps_false_positive_build(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            snps = SNPs("tests/input/generic.csv", output_dir=tmpdir)
            output = os.path.join(tmpdir, "generic_GRCh37.txt")
            self.assertEqual(snps.to_tsv(), output)

            s = ""
            with open(output, "r") as f:
                for line in f.readlines():
                    if "snps v" in line:
                        s += "# Generated by snps v1.2.3.post85.dev0+gb386302, https://pypi.org/project/snps/\n"
                    else:
                        s += line

            with open(output, "w") as f:
                f.write(s)

            self.run_parsing_tests(output, "generic")

    def test_save_snps_vcf_false_positive_build(self):
        with tempfile.TemporaryDirectory() as tmpdir1:
            snps = SNPs("tests/input/testvcf.vcf", output_dir=tmpdir1)

            r = Resources()
            r._reference_sequences["GRCh37"] = {}

            output = os.path.join(tmpdir1, "vcf_GRCh37.vcf")
            with tempfile.TemporaryDirectory() as tmpdir2:
                dest = os.path.join(tmpdir2, "generic.fa.gz")
                gzip_file("tests/input/generic.fa", dest)

                seq = ReferenceSequence(ID="1", path=dest)

                r._reference_sequences["GRCh37"]["1"] = seq

                self.assertEqual(snps.to_vcf(), output)

                s = ""
                with open(output, "r") as f:
                    for line in f.readlines():
                        if "snps v" in line:
                            s += '##source="vcf; snps v1.2.3.post85.dev0+gb386302; https://pypi.org/project/snps/"\n'
                        else:
                            s += line

                with open(output, "w") as f:
                    f.write(s)

            self.run_parsing_tests_vcf(output)

    def test_save_snps_vcf_discrepant_pos(self):
        with tempfile.TemporaryDirectory() as tmpdir1:
            s = SNPs("tests/input/testvcf.vcf", output_dir=tmpdir1)

            r = Resources()
            r._reference_sequences["GRCh37"] = {}

            output = os.path.join(tmpdir1, "vcf_GRCh37.vcf")
            with tempfile.TemporaryDirectory() as tmpdir2:
                dest = os.path.join(tmpdir2, "generic.fa.gz")
                gzip_file("tests/input/generic.fa", dest)

                seq = ReferenceSequence(ID="1", path=dest)

                r._reference_sequences["GRCh37"]["1"] = seq

                # create discrepant SNPs by setting positions outside reference sequence
                s._snps.loc["rs1", "pos"] = 0
                s._snps.loc["rs17", "pos"] = 118

                # esnure this is the right type after manual tweaking
                s._snps = s._snps.astype({"pos": np.uint32})

                self.assertEqual(s.to_vcf(), output)

            pd.testing.assert_frame_equal(
                s.discrepant_vcf_position,
                self.create_snp_df(
                    rsid=["rs1", "rs17"],
                    chrom=["1", "1"],
                    pos=[0, 118],
                    genotype=["AA", np.nan],
                ),
                check_exact=True,
            )

            expected = self.generic_snps_vcf().drop(["rs1", "rs17"])
            self.run_parsing_tests_vcf(output, snps_df=expected)

    def test_save_snps_vcf_phased(self):
        with tempfile.TemporaryDirectory() as tmpdir1:
            # read phased data
            s = SNPs("tests/input/testvcf_phased.vcf", output_dir=tmpdir1)

            # setup resource to use test FASTA reference sequence
            r = Resources()
            r._reference_sequences["GRCh37"] = {}

            output = os.path.join(tmpdir1, "vcf_GRCh37.vcf")
            with tempfile.TemporaryDirectory() as tmpdir2:
                dest = os.path.join(tmpdir2, "generic.fa.gz")
                gzip_file("tests/input/generic.fa", dest)

                seq = ReferenceSequence(ID="1", path=dest)

                r._reference_sequences["GRCh37"]["1"] = seq

                # save phased data to VCF
                self.assertEqual(s.to_vcf(), output)

            # read saved VCF
            self.run_parsing_tests_vcf(output, phased=True)

    def test_save_snps_phased(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            # read phased data
            s = SNPs("tests/input/testvcf_phased.vcf", output_dir=tmpdir)
            dest = os.path.join(tmpdir, "vcf_GRCh37.txt")
            # save phased data to TSV
            self.assertEqual(s.to_tsv(), dest)
            # read saved TSV
            self.run_parsing_tests_vcf(dest, phased=True)

    def run_vcf_qc_test(
        self, expected_output, vcf_qc_only, vcf_qc_filter, cluster="c1"
    ):
        def f():
            with tempfile.TemporaryDirectory() as tmpdir1:
                s = SNPs("tests/input/generic.csv", output_dir=tmpdir1)

                # setup resource to use test FASTA reference sequence
                r = Resources()
                r._reference_sequences["GRCh37"] = {}

                output = os.path.join(tmpdir1, "generic_GRCh37.vcf")
                with tempfile.TemporaryDirectory() as tmpdir2:
                    dest = os.path.join(tmpdir2, "generic.fa.gz")
                    gzip_file("tests/input/generic.fa", dest)

                    seq = ReferenceSequence(ID="1", path=dest)

                    r._reference_sequences["GRCh37"]["1"] = seq

                    # save phased data to VCF
                    self.assertEqual(
                        s.to_vcf(
                            qc_only=vcf_qc_only,
                            qc_filter=vcf_qc_filter,
                        ),
                        output,
                    )

                    # read result
                    with open(output, "r") as f:
                        actual = f.read()

                    # read expected result
                    with open(expected_output, "r") as f:
                        expected = f.read()

                    self.assertIn(expected, actual)

                    if not vcf_qc_filter or not cluster:
                        self.assertNotIn("##FILTER=<ID=lq", actual)

        self.run_low_quality_snps_test(f, self.get_low_quality_snps(), cluster=cluster)

    def test_save_vcf_qc_only_F_qc_filter_F(self):
        self.run_vcf_qc_test(
            "tests/output/vcf_qc/qc_only_F_qc_filter_F.vcf",
            vcf_qc_only=False,
            vcf_qc_filter=False,
        )

    def test_save_vcf_qc_only_F_qc_filter_T(self):
        self.run_vcf_qc_test(
            "tests/output/vcf_qc/qc_only_F_qc_filter_T.vcf",
            vcf_qc_only=False,
            vcf_qc_filter=True,
        )

    def test_save_vcf_qc_only_T_qc_filter_F(self):
        self.run_vcf_qc_test(
            "tests/output/vcf_qc/qc_only_T_qc_filter_F.vcf",
            vcf_qc_only=True,
            vcf_qc_filter=False,
        )

    def test_save_vcf_qc_only_T_qc_filter_T(self):
        self.run_vcf_qc_test(
            "tests/output/vcf_qc/qc_only_T_qc_filter_T.vcf",
            vcf_qc_only=True,
            vcf_qc_filter=True,
        )

    def test_save_vcf_no_cluster(self):
        self.run_vcf_qc_test(
            "tests/output/vcf_qc/qc_only_F_qc_filter_F.vcf",
            vcf_qc_only=False,
            vcf_qc_filter=False,
            cluster="",
        )
        self.run_vcf_qc_test(
            "tests/output/vcf_qc/qc_only_F_qc_filter_F.vcf",
            vcf_qc_only=False,
            vcf_qc_filter=True,
            cluster="",
        )
        self.run_vcf_qc_test(
            "tests/output/vcf_qc/qc_only_F_qc_filter_F.vcf",
            vcf_qc_only=True,
            vcf_qc_filter=False,
            cluster="",
        )
        self.run_vcf_qc_test(
            "tests/output/vcf_qc/qc_only_F_qc_filter_F.vcf",
            vcf_qc_only=True,
            vcf_qc_filter=True,
            cluster="",
        )
