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

    def test_read_ftdna_famfinder(self):
        # https://www.familytreedna.com
        self.run_parsing_tests("tests/input/ftdna_famfinder.csv", "FTDNA")

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
