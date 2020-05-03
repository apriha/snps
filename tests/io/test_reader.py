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
import pandas as pd

from snps import SNPs
from snps.utils import create_dir, zip_file, gzip_file
from tests import BaseSNPsTestCase


class TestReader(BaseSNPsTestCase):
    def test_read_DNALand(self):
        # https://dna.land/
        self.run_parsing_tests("tests/input/DNALand.txt", "DNA.Land")

    def test_read_23andme(self):
        # https://www.23andme.com
        self.run_parsing_tests("tests/input/23andme.txt", "23andMe")

    def test_read_23andme_zip(self):
        zip_file(
            "tests/input/23andme.txt", "tests/input/23andme.txt.zip", "23andme.txt"
        )
        self.run_parsing_tests("tests/input/23andme.txt.zip", "23andMe")

    def test_read_ftdna(self):
        # https://www.familytreedna.com
        self.run_parsing_tests("tests/input/ftdna.csv", "FTDNA")

    def test_read_ftdna_gzip(self):
        gzip_file("tests/input/ftdna.csv", "tests/input/ftdna.csv.gz")
        self.run_parsing_tests("tests/input/ftdna.csv.gz", "FTDNA")

    def test_read_ftdna_famfinder(self):
        # https://www.familytreedna.com
        self.run_parsing_tests("tests/input/ftdna_famfinder.csv", "FTDNA")

    def test_read_ancestry(self):
        # https://www.ancestry.com
        self.run_parsing_tests("tests/input/ancestry.txt", "AncestryDNA")

    def test_read_genes_for_good(self):
        # https://genesforgood.sph.umich.edu/
        self.run_parsing_tests("tests/input/genesforgood.txt", "GenesForGood")

    def test_read_myheritage(self):
        # https://www.myheritage.com
        self.run_parsing_tests("tests/input/myheritage.csv", "MyHeritage")

    @staticmethod
    def _setup_gsa_test():
        # reset resource if already loaded
        temp = SNPs()
        temp._resources._gsa_resources = {}

        create_dir("../resources")

        gzip_file("tests/resources/gsa_rsid_map.txt", "resources/gsa_rsid_map.txt.gz")
        gzip_file(
            "tests/resources/gsa_chrpos_map.txt", "resources/gsa_chrpos_map.txt.gz"
        )

    @staticmethod
    def _teardown_gsa_test():
        os.remove("resources/gsa_rsid_map.txt.gz")
        os.remove("resources/gsa_chrpos_map.txt.gz")

    def test_read_codigo46(self):
        # https://codigo46.com.mx
        self._setup_gsa_test()
        self.run_parsing_tests("tests/input/codigo46.txt", "Codigo46")
        self._teardown_gsa_test()

    def test_read_sano(self):
        # https://sanogenetics.com
        self._setup_gsa_test()
        self.run_parsing_tests("tests/input/sano.txt", "Sano")
        self._teardown_gsa_test()

    def test_read_livingdna(self):
        # https://livingdna.com
        self.run_parsing_tests("tests/input/livingdna.csv", "LivingDNA")

    def test_read_mapmygenome(self):
        # https://mapmygenome.in
        self.run_parsing_tests("tests/input/mapmygenome.txt", "Mapmygenome")

    def test_read_vcf(self):
        # https://samtools.github.io/hts-specs/VCFv4.2.pdf
        # this tests for homozygous snps, heterozygous snps, multiallelic snps,
        # phased snps, and snps with missing rsID
        s = SNPs("tests/input/testvcf.vcf")
        assert s.source == "vcf"
        assert not s.unannotated_vcf
        assert not s.phased
        pd.testing.assert_frame_equal(s.snps, self.generic_snps_vcf())

    def test_read_vcf_phased(self):
        # https://samtools.github.io/hts-specs/VCFv4.2.pdf
        # this tests for homozygous snps, heterozygous snps, multiallelic snps,
        # phased snps, and snps with missing rsID
        s = SNPs("tests/input/testvcf_phased.vcf")
        assert s.source == "vcf"
        assert not s.unannotated_vcf
        assert s.phased
        pd.testing.assert_frame_equal(s.snps, self.generic_snps_vcf())

    def test_read_vcf_rsids(self):
        # https://samtools.github.io/hts-specs/VCFv4.2.pdf
        # this tests for homozygous snps, heterozygous snps, multiallelic snps,
        # phased snps, and snps with missing rsID
        rsids = ["rs1", "rs2"]
        s = SNPs("tests/input/testvcf.vcf", rsids=rsids)
        assert s.source == "vcf"
        assert not s.unannotated_vcf
        pd.testing.assert_frame_equal(s.snps, self.generic_snps_vcf().loc[rsids])

    def test_read_vcf_gz(self):
        # https://samtools.github.io/hts-specs/VCFv4.2.pdf
        # this tests for homozygous snps, heterozygous snps, multiallelic snps,
        # phased snps, and snps with missing rsID
        gzip_file("tests/input/testvcf.vcf", "tests/input/testvcf.vcf.gz")
        s = SNPs("tests/input/testvcf.vcf.gz")
        assert s.source == "vcf"
        pd.testing.assert_frame_equal(s.snps, self.generic_snps_vcf())

    def test_read_vcf_gz_rsids(self):
        # https://samtools.github.io/hts-specs/VCFv4.2.pdf
        # this tests for homozygous snps, heterozygous snps, multiallelic snps,
        # phased snps, and snps with missing rsID
        gzip_file("tests/input/testvcf.vcf", "tests/input/testvcf.vcf.gz")
        rsids = ["rs1", "rs2"]
        s = SNPs("tests/input/testvcf.vcf.gz", rsids=rsids)
        assert s.source == "vcf"
        pd.testing.assert_frame_equal(s.snps, self.generic_snps_vcf().loc[rsids])

    def test_read_unannotated_vcf(self):
        # https://samtools.github.io/hts-specs/VCFv4.2.pdf
        # this tests for homozygous snps, heterozygous snps, multiallelic snps,
        # phased snps, and snps with missing rsID
        s = SNPs("tests/input/unannotated_testvcf.vcf")
        assert s.source == "vcf"
        assert s.unannotated_vcf

    def test_read_not_phased(self):
        s = SNPs("tests/input/generic.csv")
        assert not s.phased

    def test_read_generic_csv(self):
        self.run_parsing_tests("tests/input/generic.csv", "generic")

    def test_read_generic_tsv(self):
        self.run_parsing_tests("tests/input/generic.tsv", "generic")

    def test_read_generic_non_standard_columns(self):
        self.run_parsing_tests(
            "tests/input/generic_non_standard_columns.tsv", "generic"
        )

    def test_read_generic_multi_rsid_tsv(self):
        self.run_parsing_tests("tests/input/generic_multi_rsid.tsv", "generic")

    def test_read_generic_extra_column_tsv(self):
        self.run_parsing_tests("tests/input/generic_extra_column.tsv", "generic")

    def test_read_vcf_buffer(self):
        with open("tests/input/testvcf.vcf", "r") as f:
            snps_vcf_buffer = SNPs(f.read().encode("utf-8"))
        # https://samtools.github.io/hts-specs/VCFv4.2.pdf
        # this tests for homozygous snps, heterozygous snps, multiallelic snps,
        # phased snps, and snps with missing rsID
        assert snps_vcf_buffer.source == "vcf"
        pd.testing.assert_frame_equal(snps_vcf_buffer.snps, self.generic_snps_vcf())

    def test_read_vcf_buffer_rsids(self):
        with open("tests/input/testvcf.vcf", "r") as f:
            rsids = ["rs1", "rs2"]
            df = SNPs(f.read().encode("utf-8"), rsids=rsids)
        # https://samtools.github.io/hts-specs/VCFv4.2.pdf
        # this tests for homozygous snps, heterozygous snps, multiallelic snps,
        # phased snps, and snps with missing rsID
        assert df.source == "vcf"
        pd.testing.assert_frame_equal(df.snps, self.generic_snps_vcf().loc[rsids])

    def test_read_vcf_buffer_gz(self):
        gzip_file("tests/input/testvcf.vcf", "tests/input/testvcf.vcf.gz")
        with open("tests/input/testvcf.vcf.gz", "rb") as f:
            data = f.read()
            s = SNPs(data)
        os.remove("tests/input/testvcf.vcf.gz")
        # https://samtools.github.io/hts-specs/VCFv4.2.pdf
        # this tests for homozygous snps, heterozygous snps, multiallelic snps,
        # phased snps, and snps with missing rsID
        assert s.source == "vcf"
        pd.testing.assert_frame_equal(s.snps, self.generic_snps_vcf())

    def test_read_vcf_buffer_gz_rsids(self):
        gzip_file("tests/input/testvcf.vcf", "tests/input/testvcf.vcf.gz")
        with open("tests/input/testvcf.vcf.gz", "rb") as f:
            rsids = ["rs1", "rs2"]
            data = f.read()
            s = SNPs(data, rsids=rsids)
        os.remove("tests/input/testvcf.vcf.gz")
        # https://samtools.github.io/hts-specs/VCFv4.2.pdf
        # this tests for homozygous snps, heterozygous snps, multiallelic snps,
        # phased snps, and snps with missing rsID
        assert s.source == "vcf"
        pd.testing.assert_frame_equal(s.snps, self.generic_snps_vcf().loc[rsids])

    def test_source_generic(self):
        s = SNPs("tests/input/NCBI36.csv")
        assert s.source == "generic"
