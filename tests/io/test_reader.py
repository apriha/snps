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

from atomicwrites import atomic_write

from snps.resources import Resources
from snps.utils import gzip_file
from tests import BaseSNPsTestCase


class TestReader(BaseSNPsTestCase):
    @staticmethod
    def _setup_gsa_test(resources_dir):
        # reset resource if already loaded
        r = Resources()
        r._resources_dir = resources_dir
        r._init_resource_attributes()

        gzip_file(
            "tests/resources/gsa_rsid_map.txt",
            os.path.join(resources_dir, "gsa_rsid_map.txt.gz"),
        )
        gzip_file(
            "tests/resources/gsa_chrpos_map.txt",
            os.path.join(resources_dir, "gsa_chrpos_map.txt.gz"),
        )
        gzip_file(
            "tests/resources/dbsnp_151_37_reverse.txt",
            os.path.join(resources_dir, "dbsnp_151_37_reverse.txt.gz"),
        )

    @staticmethod
    def _teardown_gsa_test():
        r = Resources()
        r._resources_dir = "resources"
        r._init_resource_attributes()

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

    def test_read_23andme(self):
        # https://www.23andme.com
        self.run_parsing_tests("tests/input/23andme.txt", "23andMe")

    def test_read_23andme_build36(self):
        self.run_build_detection_test(
            self.run_parsing_tests,
            "build 36",
            36,
            file="tests/input/23andme.txt",
            source="23andMe",
            comment_str="# {}\n",
        )

    def test_read_23andme_build37(self):
        self.run_build_detection_test(
            self.run_parsing_tests,
            "build 37",
            37,
            file="tests/input/23andme.txt",
            source="23andMe",
            comment_str="# {}\n",
        )

    def test_read_ancestry(self):
        # https://www.ancestry.com
        self.run_parsing_tests("tests/input/ancestry.txt", "AncestryDNA")

    def test_read_ancestry_extra_tab(self):
        # https://www.ancestry.com

        total_snps = 100
        s = "#Ancestry\r\n"
        s += "rsid\tchromosome\tposition\tallele1\tallele2\r\n"
        # add extra tab separator in first line
        s += "rs1\t1\t101\t\tA\tA\r\n"
        # generate remainder of lines
        for i in range(1, total_snps):
            s += f"rs{1 + i}\t1\t{101 + i}\tA\tA\r\n"

        snps_df = self.create_snp_df(
            rsid=[f"rs{1 + i}" for i in range(0, total_snps)],
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

    def test_read_tellmeGen(self):
        # https://www.tellmegen.com/
        with tempfile.TemporaryDirectory() as tmpdir:
            self._setup_gsa_test(tmpdir)
            self.run_parsing_tests("tests/input/tellmeGen.txt", "tellmeGen")
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
            s1 += f'"rs{1 + i}","1","{101 + i}","AA"\r\n'

        # generate content of second file
        s2 = "RSID,CHROMOSOME,POSITION,RESULT\r\n"
        for i in range(0, total_snps2):
            s2 += f'"rs{total_snps1 + 1 + i}","1","{ total_snps1 + 101 + i}","AA"\r\n'
        snps_df = self.create_snp_df(
            rsid=[f"rs{1 + i}" for i in range(0, total_snps1 + total_snps2)],
            chrom="1",
            pos=[101 + i for i in range(0, total_snps1 + total_snps2)],
            genotype="AA",
        )

        with tempfile.TemporaryDirectory() as tmpdir:
            file1 = os.path.join(tmpdir, "ftdna_concat_gzip1.csv")
            file1_gz = f"{file1}.gz"
            file2 = os.path.join(tmpdir, "ftdna_concat_gzip2.csv")
            file2_gz = f"{file2}.gz"
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

        total_snps1 = 100
        total_snps2 = 10
        s = "RSID,CHROMOSOME,POSITION,RESULT\n"
        # generate first chunk of lines
        for i in range(0, total_snps1):
            s += f'"rs{1 + i}","1","{101 + i}","AA"\n'
        # add second header
        s += "RSID,CHROMOSOME,POSITION,RESULT\n"
        # generate second chunk of lines
        for i in range(0, total_snps2):
            s += f'"rs{total_snps1 + 1 + i}","1","{total_snps1 + 101 + i}","AA"\n'

        snps_df = self.create_snp_df(
            rsid=[f"rs{1 + i}" for i in range(0, total_snps1 + total_snps2)],
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

    def test_read_mapmygenome_alt_header(self):
        # https://mapmygenome.in
        self.run_parsing_tests("tests/input/mapmygenome_alt_header.txt", "Mapmygenome")

    def test_read_myheritage(self):
        # https://www.myheritage.com
        self.run_parsing_tests("tests/input/myheritage.csv", "MyHeritage")

    def test_read_myheritage_extra_quotes(self):
        # https://www.myheritage.com
        self.run_parsing_tests("tests/input/myheritage_extra_quotes.csv", "MyHeritage")

    def test_read_sano(self):
        # https://sanogenetics.com
        with tempfile.TemporaryDirectory() as tmpdir:
            self._setup_gsa_test(tmpdir)
            self.run_parsing_tests("tests/input/sano.txt", "Sano")
            self._teardown_gsa_test()

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

    def test_read_vcf_rsids(self):
        self.run_parsing_tests_vcf("tests/input/testvcf.vcf", rsids=["rs1", "rs2"])

    def test_read_unannotated_vcf(self):
        self.run_parsing_tests_vcf(
            "tests/input/unannotated_testvcf.vcf", "vcf", unannotated=True, build=0
        )
