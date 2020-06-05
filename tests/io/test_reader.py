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
import warnings

import pandas as pd

from snps.resources import Resources
from snps.utils import gzip_file
from tests import BaseSNPsTestCase


class TestReader(BaseSNPsTestCase):
    def setUp(self):
        # set warnings filter (test runner overrides statement in `io/reader.py`)
        warnings.filterwarnings("error", category=pd.errors.DtypeWarning)
        super().setUp()

    @staticmethod
    def _setup_gsa_test(resources_dir):
        # reset resource if already loaded
        r = Resources()
        r._resources_dir = resources_dir
        r._gsa_resources = {}

        gzip_file(
            "tests/resources/gsa_rsid_map.txt",
            os.path.join(resources_dir, "gsa_rsid_map.txt.gz"),
        )
        gzip_file(
            "tests/resources/gsa_chrpos_map.txt",
            os.path.join(resources_dir, "gsa_chrpos_map.txt.gz"),
        )

    @staticmethod
    def _teardown_gsa_test():
        r = Resources()
        r._resources_dir = "resources"
        r._gsa_resources = {}

    def test_read_23andme(self):
        # https://www.23andme.com
        self.run_parsing_tests("tests/input/23andme.txt", "23andMe")

    def test_read_23andme_build36(self):
        self.run_parsing_tests(
            "tests/input/23andme_build36.txt", "23andMe", build=36, build_detected=True
        )

    def test_read_23andme_build37(self):
        self.run_parsing_tests(
            "tests/input/23andme_build37.txt", "23andMe", build=37, build_detected=True
        )

    def test_read_ancestry(self):
        # https://www.ancestry.com
        self.run_parsing_tests("tests/input/ancestry.txt", "AncestryDNA")

    def test_read_ancestry_extra_tab(self):
        # https://www.ancestry.com

        # we need a large file to generate the `pd.errors.DtypeWarning`
        total_snps = 500000
        s = "#Ancestry\r\n"
        s += "rsid\tchromosome\tposition\tallele1\tallele2\r\n"
        # add extra tab separator in first line
        s += "rs1\t1\t101\t\tA\tA\r\n"
        # generate remainder of lines
        for i in range(1, total_snps):
            s += "rs{}\t1\t{}\tA\tA\r\n".format(1 + i, 101 + i)

        snps_df = self.create_snp_df(
            rsid=["rs{}".format(1 + i) for i in range(0, total_snps)],
            chrom="1",
            pos=[101 + i for i in range(0, total_snps)],
            genotype="AA",
        )

        with tempfile.TemporaryDirectory() as tmpdir:
            path = os.path.join(tmpdir, "ancestry_extra_tab.txt")
            with open(path, "w") as f:
                f.write(s)

            self.run_parsing_tests(path, "AncestryDNA", snps_df=snps_df)

    def test_read_ancestry_multi_sep(self):
        # https://www.ancestry.com
        self.run_parsing_tests("tests/input/ancestry_multi_sep.txt", "AncestryDNA")

    def test_read_codigo46(self):
        # https://codigo46.com.mx
        with tempfile.TemporaryDirectory() as tmpdir:
            self._setup_gsa_test(tmpdir)
            self.run_parsing_tests("tests/input/codigo46.txt", "Codigo46")
            self._teardown_gsa_test()

    def test_read_DNALand(self):
        # https://dna.land/
        self.run_parsing_tests("tests/input/DNALand.txt", "DNA.Land")

    def test_read_ftdna(self):
        # https://www.familytreedna.com
        self.run_parsing_tests("tests/input/ftdna.csv", "FTDNA")

    def test_read_ftdna_concat_gzip_extra_data(self):
        # https://www.familytreedna.com

        total_snps1 = 10
        total_snps2 = 10
        # generate content of first file
        s1 = "RSID,CHROMOSOME,POSITION,RESULT\r\n"
        for i in range(0, total_snps1):
            s1 += '"rs{}","1","{}","AA"\r\n'.format(1 + i, 101 + i)

        # generate content of second file
        s2 = "RSID,CHROMOSOME,POSITION,RESULT\r\n"
        for i in range(0, total_snps2):
            s2 += '"rs{}","1","{}","AA"\r\n'.format(
                total_snps1 + 1 + i, total_snps1 + 101 + i
            )

        snps_df = self.create_snp_df(
            rsid=["rs{}".format(1 + i) for i in range(0, total_snps1 + total_snps2)],
            chrom="1",
            pos=[101 + i for i in range(0, total_snps1 + total_snps2)],
            genotype="AA",
        )

        with tempfile.TemporaryDirectory() as tmpdir:
            file1 = os.path.join(tmpdir, "ftdna_concat_gzip1.csv")
            file1_gz = "{}.gz".format(file1)
            file2 = os.path.join(tmpdir, "ftdna_concat_gzip2.csv")
            file2_gz = "{}.gz".format(file2)
            path = os.path.join(tmpdir, "ftdna_concat_gzip.csv.gz")

            # write individual files
            with open(file1, "w") as f:
                f.write(s1)
            with open(file2, "w") as f:
                f.write(s2)

            # compress files
            gzip_file(file1, file1_gz)
            gzip_file(file2, file2_gz)

            # concatenate gzips
            with open(file1_gz, "rb") as f:
                data = f.read()
            with open(file2_gz, "rb") as f:
                data += f.read()

            # add extra data
            data += b"extra data"

            # write file with concatenated gzips and extra data
            with open(path, "wb") as f:
                f.write(data)

            self.make_parsing_assertions(
                self.parse_file(path), "FTDNA", False, 37, False, snps_df
            )
            self.make_parsing_assertions(
                self.parse_bytes(path), "FTDNA", False, 37, False, snps_df
            )

    def test_read_ftdna_famfinder(self):
        # https://www.familytreedna.com
        self.run_parsing_tests("tests/input/ftdna_famfinder.csv", "FTDNA")

    def test_read_ftdna_second_header(self):
        # https://www.familytreedna.com

        # we need a large file to generate the `pd.errors.DtypeWarning`
        total_snps1 = 500000
        total_snps2 = 10000
        s = "RSID,CHROMOSOME,POSITION,RESULT\n"
        # generate first chunk of lines
        for i in range(0, total_snps1):
            s += '"rs{}","1","{}","AA"\n'.format(1 + i, 101 + i)
        # add second header
        s += "RSID,CHROMOSOME,POSITION,RESULT\n"
        # generate second chunk of lines
        for i in range(0, total_snps2):
            s += '"rs{}","1","{}","AA"\n'.format(
                total_snps1 + 1 + i, total_snps1 + 101 + i
            )

        snps_df = self.create_snp_df(
            rsid=["rs{}".format(1 + i) for i in range(0, total_snps1 + total_snps2)],
            chrom="1",
            pos=[101 + i for i in range(0, total_snps1 + total_snps2)],
            genotype="AA",
        )

        with tempfile.TemporaryDirectory() as tmpdir:
            path = os.path.join(tmpdir, "ftdna_second_header.txt")
            with open(path, "w") as f:
                f.write(s)

            self.run_parsing_tests(path, "FTDNA", snps_df=snps_df)

    def test_read_generic_csv(self):
        self.run_parsing_tests("tests/input/generic.csv", "generic")

    def test_read_generic_tsv(self):
        self.run_parsing_tests("tests/input/generic.tsv", "generic")

    def test_read_generic_extra_column_tsv(self):
        self.run_parsing_tests("tests/input/generic_extra_column.tsv", "generic")

    def test_read_generic_header_comment(self):
        self.run_parsing_tests("tests/input/generic_header_comment.tsv", "generic")

    def test_read_generic_multi_rsid_tsv(self):
        self.run_parsing_tests("tests/input/generic_multi_rsid.tsv", "generic")

    def test_read_generic_no_header(self):
        self.run_parsing_tests("tests/input/generic_no_header.tsv", "generic")

    def test_read_generic_non_standard_columns(self):
        self.run_parsing_tests(
            "tests/input/generic_non_standard_columns.tsv", "generic"
        )

    def test_read_genes_for_good(self):
        # https://genesforgood.sph.umich.edu/
        self.run_parsing_tests("tests/input/genesforgood.txt", "GenesForGood")

    def test_read_livingdna(self):
        # https://livingdna.com
        self.run_parsing_tests("tests/input/livingdna.csv", "LivingDNA")

    def test_read_mapmygenome(self):
        # https://mapmygenome.in
        self.run_parsing_tests("tests/input/mapmygenome.txt", "Mapmygenome")

    def test_read_myheritage(self):
        # https://www.myheritage.com
        self.run_parsing_tests("tests/input/myheritage.csv", "MyHeritage")

    def test_read_sano(self):
        # https://sanogenetics.com
        with tempfile.TemporaryDirectory() as tmpdir:
            self._setup_gsa_test(tmpdir)
            self.run_parsing_tests("tests/input/sano.txt", "Sano")
            self._teardown_gsa_test()

    def test_read_vcf(self):
        self.run_parsing_tests_vcf("tests/input/testvcf.vcf")

    def test_read_vcf_b37(self):
        self.run_parsing_tests_vcf(
            "tests/input/testvcf_b37.vcf", build=37, build_detected=True
        )

    def test_read_vcf_hg19(self):
        self.run_parsing_tests_vcf(
            "tests/input/testvcf_hg19.vcf", build=37, build_detected=True
        )

    def test_read_vcf_multi_sample(self):
        self.run_parsing_tests_vcf("tests/input/testvcf_multi_sample.vcf")

    def test_read_vcf_phased(self):
        self.run_parsing_tests_vcf("tests/input/testvcf_phased.vcf", phased=True)

    def test_read_vcf_rsids(self):
        self.run_parsing_tests_vcf("tests/input/testvcf.vcf", rsids=["rs1", "rs2"])

    def test_read_unannotated_vcf(self):
        self.run_parsing_tests_vcf(
            "tests/input/unannotated_testvcf.vcf", unannotated=True, build=0
        )
