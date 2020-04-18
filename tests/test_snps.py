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
import io
from atomicwrites import atomic_write
import zipfile
import gzip
import pandas as pd
import shutil

import numpy as np

from snps import SNPs
from snps.utils import create_dir
from snps.resources import Resources, ReferenceSequence
from tests import BaseSNPsTestCase


class TestSnps(BaseSNPsTestCase):
    def setUp(self):
        self.snps_GRCh38 = SNPs("tests/input/GRCh38.csv")
        self.snps = SNPs("tests/input/chromosomes.csv")
        self.snps_only_detect_source = SNPs(
            "tests/input/chromosomes.csv", only_detect_source=True
        )
        self.snps_none = SNPs(None)

        with open("tests/input/chromosomes.csv", "r") as f:
            self.snps_buffer = SNPs(f.read().encode("utf-8"))

        with atomic_write(
            "tests/input/chromosomes.csv.zip", mode="wb", overwrite=True
        ) as f:
            with zipfile.ZipFile(f, "w") as f_zip:
                f_zip.write("tests/input/chromosomes.csv", arcname="chromosomes.csv")

        with open("tests/input/chromosomes.csv.zip", "rb") as f:
            data = f.read()
            self.snps_buffer_zip = SNPs(data)
        os.remove("tests/input/chromosomes.csv.zip")

        with open("tests/input/chromosomes.csv", "rb") as f_in:
            with atomic_write(
                "tests/input/chromosomes.csv.gz", mode="wb", overwrite=True
            ) as f_out:
                with gzip.open(f_out, "wb") as f_gzip:
                    shutil.copyfileobj(f_in, f_gzip)

        with open("tests/input/chromosomes.csv.gz", "rb") as f:
            data = f.read()
            self.snps_buffer_gz = SNPs(data)
        os.remove("tests/input/chromosomes.csv.gz")

    def snps_discrepant_pos(self):
        return self.create_snp_df(
            rsid=["rs3094315"], chrom=["1"], pos=[1], genotype=["AA"]
        )

    def test_assembly(self):
        assert self.snps_GRCh38.assembly == "GRCh38"

    def test_assembly_no_snps(self):
        assert self.snps_none.assembly == ""

    def test_snp_buffer_zip(self):
        assert self.snps_buffer_zip.snp_count == 6

    def test_snp_buffer_gz(self):

        assert self.snps_buffer_gz.snp_count == 6

    def test_snp_buffer(self):
        assert self.snps_buffer.snp_count == 6

    def test_build_no_snps(self):
        assert not self.snps_none.build

    def test_build_detected_no_snps(self):
        assert not self.snps_none.build_detected

    def test_build_detected_PAR_snps(self):
        if (
            not os.getenv("DOWNLOADS_ENABLED")
            or os.getenv("DOWNLOADS_ENABLED") == "true"
        ):
            snps = SNPs("tests/input/GRCh37_PAR.csv")
            assert snps.build == 37
            assert snps.build_detected

    def test_sex_no_snps(self):
        assert self.snps_none.sex == ""

    def test_get_summary(self):
        assert self.snps_GRCh38.get_summary() == {
            "source": "generic",
            "assembly": "GRCh38",
            "build": 38,
            "build_detected": True,
            "snp_count": 4,
            "chromosomes": "1, 3",
            "sex": "",
        }

    def test_get_summary_no_snps(self):
        assert not self.snps_none.get_summary()

    def test_is_valid_True(self):
        assert self.snps_GRCh38.is_valid()

    def test_is_valid_False(self):
        assert not self.snps_none.is_valid()

    def test__read_raw_data(self):
        assert self.snps_none.snps.empty
        assert self.snps_none.source == ""

    def test__lookup_build_with_snp_pos_None(self):
        snps = SNPs()
        snps._snps = self.snps_discrepant_pos()
        assert not snps.detect_build()

    def test_get_assembly_None(self):
        snps = SNPs()
        assert snps.get_assembly() is ""

    def test_save_snps_source(self):
        assert (
            os.path.relpath(self.snps_GRCh38.save_snps()) == "output/generic_GRCh38.csv"
        )
        snps = SNPs("output/generic_GRCh38.csv")
        pd.testing.assert_frame_equal(snps.snps, self.snps_GRCh38.snps)

    def test_save_snps_buffer(self):
        out = io.StringIO()
        self.snps.save_snps(out)
        assert out.read().startswith("# Generated by snps")

    def test_snps_only_detect_source(self):
        assert self.snps_only_detect_source.source == "generic"

    def test_duplicate_rsids(self):
        snps = SNPs("tests/input/duplicate_rsids.csv")
        result = self.create_snp_df(
            rsid=["rs1"], chrom=["1"], pos=[101], genotype=["AA"]
        )
        duplicate_snps = self.create_snp_df(
            rsid=["rs1", "rs1"], chrom=["1", "1"], pos=[102, 103], genotype=["CC", "GG"]
        )
        pd.testing.assert_frame_equal(snps.snps, result)
        pd.testing.assert_frame_equal(snps.duplicate_snps, duplicate_snps)

    def test_deduplicate_false(self):
        snps = SNPs("tests/input/duplicate_rsids.csv", deduplicate=False)

        result = self.create_snp_df(
            rsid=["rs1", "rs1", "rs1"],
            chrom=["1", "1", "1"],
            pos=[101, 102, 103],
            genotype=["AA", "CC", "GG"],
        )
        pd.testing.assert_frame_equal(snps.snps, result)

    def test_empty_dataframe(self):
        s = SNPs()
        assert s.snp_count == 0
        assert list(s.snps.columns.values) == ["chrom", "pos", "genotype"]
        assert s.snps.index.name == "rsid"

    def test_empty_dataframe_file(self):
        with atomic_write("tests/input/empty.txt", mode="w", overwrite=True):
            pass
        s = SNPs("tests/input/empty.txt")
        assert s.snp_count == 0
        assert list(s.snps.columns.values) == ["chrom", "pos", "genotype"]
        assert s.snps.index.name == "rsid"

    def generic_snps(self):
        return self.create_snp_df(
            rsid=["rs" + str(i) for i in range(1, 9)],
            chrom=["1"] * 8,
            pos=list(range(101, 109)),
            genotype=["AA", "CC", "GG", "TT", np.nan, "GC", "TC", "AT"],
        )

    def generic_snps_vcf(self):
        df = self.generic_snps()
        return df.append(
            self.create_snp_df(
                rsid=["rs" + str(i) for i in range(12, 18)],
                chrom=["1"] * 6,
                pos=list(range(112, 118)),
                genotype=[np.nan] * 6,
            )
        )

    def snps_GRCh38_func(self):
        return self.create_snp_df(
            rsid=["rs3094315", "rsIndelTest", "rs2500347", "rs11928389"],
            chrom=["1", "1", "1", "3"],
            pos=[817186, 148946168, 148946169, 50889578],
            genotype=["AA", "ID", np.nan, "TC"],
        )

    def snps_GRCh38_PAR(self):
        return self.create_snp_df(
            rsid=["rs28736870", "rs113378274", "rs113313554"],
            chrom=["X", "X", "Y"],
            pos=[304103, 93431058, 624523],
            genotype=["AA", "AA", "AA"],
        )

    def run_parsing_tests(self, file, source):
        self.make_parsing_assertions(self.parse_file(file), source)
        self.make_parsing_assertions(self.parse_bytes(file), source)

    def parse_file(self, file):
        return SNPs(file)

    def parse_bytes(self, file):
        with open(file, "rb") as f:
            return SNPs(f.read())

    def make_parsing_assertions(self, snps, source):
        self.assertEqual(snps.source, source)
        pd.testing.assert_frame_equal(snps.snps, self.generic_snps())

    def test_snps_DNALand(self):
        # https://dna.land/
        self.run_parsing_tests("tests/input/DNALand.txt", "DNA.Land")

    def test_snps_23andme(self):
        # https://www.23andme.com
        self.run_parsing_tests("tests/input/23andme.txt", "23andMe")

    def test_snps_23andme_zip(self):
        with atomic_write(
            "tests/input/23andme.txt.zip", mode="wb", overwrite=True
        ) as f:
            with zipfile.ZipFile(f, "w") as f_zip:
                # https://stackoverflow.com/a/16104667
                f_zip.write("tests/input/23andme.txt", arcname="23andme.txt")
        self.run_parsing_tests("tests/input/23andme.txt.zip", "23andMe")

    def test_snps_ftdna(self):
        # https://www.familytreedna.com
        self.run_parsing_tests("tests/input/ftdna.csv", "FTDNA")

    def test_snps_ftdna_gzip(self):
        with open("tests/input/ftdna.csv", "rb") as f_in:
            with atomic_write(
                "tests/input/ftdna.csv.gz", mode="wb", overwrite=True
            ) as f_out:
                with gzip.open(f_out, "wb") as f_gzip:
                    shutil.copyfileobj(f_in, f_gzip)
        self.run_parsing_tests("tests/input/ftdna.csv.gz", "FTDNA")

    def test_snps_ftdna_famfinder(self):
        # https://www.familytreedna.com
        self.run_parsing_tests("tests/input/ftdna_famfinder.csv", "FTDNA")

    def test_snps_ancestry(self):
        # https://www.ancestry.com
        self.run_parsing_tests("tests/input/ancestry.txt", "AncestryDNA")

    def test_snps_genes_for_good(self):
        # https://genesforgood.sph.umich.edu/
        self.run_parsing_tests("tests/input/genesforgood.txt", "GenesForGood")

    def test_snps_myheritage(self):
        # https://www.myheritage.com
        self.run_parsing_tests("tests/input/myheritage.csv", "MyHeritage")

    @staticmethod
    def _setup_gsa_test():
        # reset resource if already loaded
        temp = SNPs()
        temp._resources._gsa_resources = {}

        create_dir("resources")

        with open("tests/resources/gsa_rsid_map.txt", "rb") as f_in:
            with atomic_write(
                "resources/gsa_rsid_map.txt.gz", mode="wb", overwrite=True
            ) as f_out:
                with gzip.open(f_out, "wb") as f_gzip:
                    shutil.copyfileobj(f_in, f_gzip)

        with open("tests/resources/gsa_chrpos_map.txt", "rb") as f_in:
            with atomic_write(
                "resources/gsa_chrpos_map.txt.gz", mode="wb", overwrite=True
            ) as f_out:
                with gzip.open(f_out, "wb") as f_gzip:
                    shutil.copyfileobj(f_in, f_gzip)

    @staticmethod
    def _teardown_gsa_test():
        os.remove("resources/gsa_rsid_map.txt.gz")
        os.remove("resources/gsa_chrpos_map.txt.gz")

    def test_snps_codigo46(self):
        # https://codigo46.com.mx
        self._setup_gsa_test()
        self.run_parsing_tests("tests/input/codigo46.txt", "Codigo46")
        self._teardown_gsa_test()

    def test_snps_sano(self):
        # https://sanogenetics.com
        self._setup_gsa_test()
        self.run_parsing_tests("tests/input/sano.txt", "Sano")
        self._teardown_gsa_test()

    def test_snps_livingdna(self):
        # https://livingdna.com
        self.run_parsing_tests("tests/input/livingdna.csv", "LivingDNA")

    def test_snps_mapmygenome(self):
        # https://mapmygenome.in
        self.run_parsing_tests("tests/input/mapmygenome.txt", "Mapmygenome")

    def test_snps_vcf(self):
        # https://samtools.github.io/hts-specs/VCFv4.2.pdf
        # this tests for homozygous snps, heterozygous snps, multiallelic snps,
        # phased snps, and snps with missing rsID
        s = SNPs("tests/input/testvcf.vcf")
        assert s.source == "vcf"
        assert not s.unannotated_vcf
        assert not s.phased
        pd.testing.assert_frame_equal(s.snps, self.generic_snps_vcf())

    def test_snps_vcf_phased(self):
        # https://samtools.github.io/hts-specs/VCFv4.2.pdf
        # this tests for homozygous snps, heterozygous snps, multiallelic snps,
        # phased snps, and snps with missing rsID
        s = SNPs("tests/input/testvcf_phased.vcf")
        assert s.source == "vcf"
        assert not s.unannotated_vcf
        assert s.phased
        pd.testing.assert_frame_equal(s.snps, self.generic_snps_vcf())

    def test_snps_vcf_rsids(self):
        # https://samtools.github.io/hts-specs/VCFv4.2.pdf
        # this tests for homozygous snps, heterozygous snps, multiallelic snps,
        # phased snps, and snps with missing rsID
        rsids = ["rs1", "rs2"]
        s = SNPs("tests/input/testvcf.vcf", rsids=rsids)
        assert s.source == "vcf"
        assert not s.unannotated_vcf
        pd.testing.assert_frame_equal(s.snps, self.generic_snps_vcf().loc[rsids])

    def test_snps_vcf_gz(self):
        # https://samtools.github.io/hts-specs/VCFv4.2.pdf
        # this tests for homozygous snps, heterozygous snps, multiallelic snps,
        # phased snps, and snps with missing rsID
        with open("tests/input/testvcf.vcf", "rb") as f_in:
            with atomic_write(
                "tests/input/testvcf.vcf.gz", mode="wb", overwrite=True
            ) as f_out:
                with gzip.open(f_out, "wb") as f_gzip:
                    shutil.copyfileobj(f_in, f_gzip)

        s = SNPs("tests/input/testvcf.vcf.gz")
        assert s.source == "vcf"
        pd.testing.assert_frame_equal(s.snps, self.generic_snps_vcf())

    def test_snps_vcf_gz_rsids(self):
        # https://samtools.github.io/hts-specs/VCFv4.2.pdf
        # this tests for homozygous snps, heterozygous snps, multiallelic snps,
        # phased snps, and snps with missing rsID
        with open("tests/input/testvcf.vcf", "rb") as f_in:
            with atomic_write(
                "tests/input/testvcf.vcf.gz", mode="wb", overwrite=True
            ) as f_out:
                with gzip.open(f_out, "wb") as f_gzip:
                    shutil.copyfileobj(f_in, f_gzip)

        rsids = ["rs1", "rs2"]
        s = SNPs("tests/input/testvcf.vcf.gz", rsids=rsids)
        assert s.source == "vcf"
        pd.testing.assert_frame_equal(s.snps, self.generic_snps_vcf().loc[rsids])

    def test_snps_unannotated_vcf(self):
        # https://samtools.github.io/hts-specs/VCFv4.2.pdf
        # this tests for homozygous snps, heterozygous snps, multiallelic snps,
        # phased snps, and snps with missing rsID
        s = SNPs("tests/input/unannotated_testvcf.vcf")
        assert s.source == "vcf"
        assert s.unannotated_vcf

    def test_snps_not_phased(self):
        s = SNPs("tests/input/generic.csv")
        assert s.source == "generic"
        assert not s.phased

    def test_snps_vcf_buffer(self):
        with open("tests/input/testvcf.vcf", "r") as f:
            snps_vcf_buffer = SNPs(f.read().encode("utf-8"))
        # https://samtools.github.io/hts-specs/VCFv4.2.pdf
        # this tests for homozygous snps, heterozygous snps, multiallelic snps,
        # phased snps, and snps with missing rsID
        assert snps_vcf_buffer.source == "vcf"
        pd.testing.assert_frame_equal(snps_vcf_buffer.snps, self.generic_snps_vcf())

    def test_snps_vcf_buffer_rsids(self):
        with open("tests/input/testvcf.vcf", "r") as f:
            rsids = ["rs1", "rs2"]
            df = SNPs(f.read().encode("utf-8"), rsids=rsids)
        # https://samtools.github.io/hts-specs/VCFv4.2.pdf
        # this tests for homozygous snps, heterozygous snps, multiallelic snps,
        # phased snps, and snps with missing rsID
        assert df.source == "vcf"
        pd.testing.assert_frame_equal(df.snps, self.generic_snps_vcf().loc[rsids])

    def test_snps_vcf_buffer_gz(self):
        with open("tests/input/testvcf.vcf", "rb") as f_in:
            with atomic_write(
                "tests/input/testvcf.vcf.gz", mode="wb", overwrite=True
            ) as f_out:
                with gzip.open(f_out, "wb") as f_gzip:
                    shutil.copyfileobj(f_in, f_gzip)

        with open("tests/input/testvcf.vcf.gz", "rb") as f:
            data = f.read()
            s = SNPs(data)
        os.remove("tests/input/testvcf.vcf.gz")
        # https://samtools.github.io/hts-specs/VCFv4.2.pdf
        # this tests for homozygous snps, heterozygous snps, multiallelic snps,
        # phased snps, and snps with missing rsID
        assert s.source == "vcf"
        pd.testing.assert_frame_equal(s.snps, self.generic_snps_vcf())

    def test_snps_vcf_buffer_gz_rsids(self):
        with open("tests/input/testvcf.vcf", "rb") as f_in:
            with atomic_write(
                "tests/input/testvcf.vcf.gz", mode="wb", overwrite=True
            ) as f_out:
                with gzip.open(f_out, "wb") as f_gzip:
                    shutil.copyfileobj(f_in, f_gzip)

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

    def test_snps_no_snps(self):
        s = SNPs()
        assert s.snps.empty

    def test_snp_count(self):
        s = SNPs("tests/input/NCBI36.csv")
        assert s.snp_count == 4

    def test_snp_count_no_snps(self):
        s = SNPs()
        assert s.snp_count == 0

    def test_chromosomes(self):
        s = SNPs("tests/input/chromosomes.csv")
        assert s.chromosomes == ["1", "2", "3", "5", "PAR", "MT"]

    def test_chromosomes_no_snps(self):
        s = SNPs()
        assert s.chromosomes == []

    def test_chromosomes_summary(self):
        s = SNPs("tests/input/chromosomes.csv")
        assert s.chromosomes_summary == "1-3, 5, PAR, MT"

    def test_chromosomes_summary_no_snps(self):
        s = SNPs()
        assert s.chromosomes_summary == ""

    def test_build(self):
        s = SNPs("tests/input/NCBI36.csv")
        assert s.build == 36
        assert s.assembly == "NCBI36"

    def test_sex_Male_Y_chrom(self):
        s = self.simulate_snps(chrom="Y", pos_start=1, pos_max=59373566, pos_step=10000)
        assert s.sex == "Male"

    def test_sex_Female_Y_chrom(self):
        s = self.simulate_snps(
            chrom="Y", pos_start=1, pos_max=59373566, pos_step=10000, null_snp_step=1
        )
        assert s.sex == "Female"

    def test_sex_Female_X_chrom(self):
        s = self.simulate_snps(
            chrom="X", pos_start=1, pos_max=155270560, pos_step=10000, genotype="AC"
        )
        assert s.sex == "Female"

    def test_sex_Male_X_chrom(self):
        s = self.simulate_snps(
            chrom="X", pos_start=1, pos_max=155270560, pos_step=10000, genotype="AA"
        )
        assert s.snp_count == 15528
        s._deduplicate_XY_chrom()
        assert s.snp_count == 15528
        assert len(s.discrepant_XY_snps) == 0
        assert s.sex == "Male"

    def test_sex_Male_X_chrom_discrepant_XY_snps(self):
        s = self.simulate_snps(
            chrom="X", pos_start=1, pos_max=155270560, pos_step=10000, genotype="AA"
        )
        assert s.snp_count == 15528
        s._snps.loc["rs8001", "genotype"] = "AC"
        s._deduplicate_XY_chrom()
        assert s.snp_count == 15527
        result = self.create_snp_df(
            rsid=["rs8001"], chrom=["X"], pos=[80000001], genotype=["AC"]
        )
        pd.testing.assert_frame_equal(s.discrepant_XY_snps, result)
        assert s.sex == "Male"

    def test_sex_not_determined(self):
        s = self.simulate_snps(
            chrom="1", pos_start=1, pos_max=249250621, pos_step=10000
        )
        assert s.sex == ""

    def test_save_snps(self):
        snps = SNPs("tests/input/GRCh37.csv")
        assert os.path.relpath(snps.save_snps()) == "output/generic_GRCh37.csv"
        s_saved = SNPs("output/generic_GRCh37.csv")
        pd.testing.assert_frame_equal(s_saved.snps, self.snps_GRCh37())

    def test_save_snps_bytes(self):
        snps = SNPs("tests/input/GRCh37.csv")
        assert os.path.relpath(snps.save_snps()) == "output/generic_GRCh37.csv"
        with open("output/generic_GRCh37.csv", "rb") as f:
            s_saved = SNPs(f.read())
        pd.testing.assert_frame_equal(s_saved.snps, self.snps_GRCh37())

    def test_save_snps_tsv(self):
        snps = SNPs("tests/input/generic.csv")
        assert (
            os.path.relpath(snps.save_snps("generic.tsv", sep="\t"))
            == "output/generic.tsv"
        )
        self.run_parsing_tests("output/generic.tsv", "generic")

    def test_save_snps_vcf(self):
        s = SNPs("tests/input/testvcf.vcf")

        r = Resources()
        r._reference_sequences["GRCh37"] = {}
        with open("tests/input/generic.fa", "rb") as f_in:
            with atomic_write(
                "tests/input/generic.fa.gz", mode="wb", overwrite=True
            ) as f_out:
                with gzip.open(f_out, "wb") as f_gzip:
                    shutil.copyfileobj(f_in, f_gzip)

        seq = ReferenceSequence(ID="1", path="tests/input/generic.fa.gz")

        r._reference_sequences["GRCh37"]["1"] = seq

        assert os.path.relpath(s.save_snps(vcf=True)) == "output/vcf_GRCh37.vcf"
        s = SNPs("output/vcf_GRCh37.vcf")
        assert not s.phased
        pd.testing.assert_frame_equal(s.snps, self.generic_snps_vcf())

    def test_save_snps_vcf_phased(self):
        # read phased data
        s = SNPs("tests/input/testvcf_phased.vcf")

        # setup resource to use test FASTA reference sequence
        r = Resources()
        r._reference_sequences["GRCh37"] = {}
        with open("tests/input/generic.fa", "rb") as f_in:
            with atomic_write(
                "tests/input/generic.fa.gz", mode="wb", overwrite=True
            ) as f_out:
                with gzip.open(f_out, "wb") as f_gzip:
                    shutil.copyfileobj(f_in, f_gzip)

        seq = ReferenceSequence(ID="1", path="tests/input/generic.fa.gz")

        r._reference_sequences["GRCh37"]["1"] = seq

        # save phased data to VCF
        assert os.path.relpath(s.save_snps(vcf=True)) == "output/vcf_GRCh37.vcf"
        # read saved VCF
        s = SNPs("output/vcf_GRCh37.vcf")
        assert s.phased
        pd.testing.assert_frame_equal(s.snps, self.generic_snps_vcf())

    def test_save_snps_csv_phased(self):
        # read phased data
        s = SNPs("tests/input/testvcf_phased.vcf")
        # save phased data to CSV
        assert os.path.relpath(s.save_snps()) == "output/vcf_GRCh37.csv"
        # read saved CSV
        s = SNPs("output/vcf_GRCh37.csv")
        assert s.phased
        pd.testing.assert_frame_equal(s.snps, self.generic_snps_vcf())

    def test_save_snps_specify_file(self):
        s = SNPs("tests/input/GRCh37.csv")
        assert os.path.relpath(s.save_snps("snps.csv")) == "output/snps.csv"
        s_saved = SNPs("output/snps.csv")
        pd.testing.assert_frame_equal(s_saved.snps, self.snps_GRCh37())

    def test_save_snps_no_snps(self):
        s = SNPs()
        assert not s.save_snps()

    def test_save_snps_no_snps_vcf(self):
        s = SNPs()
        assert not s.save_snps(vcf=True)

    def test_remap_snps_36_to_37(self):
        s = SNPs("tests/input/NCBI36.csv")
        chromosomes_remapped, chromosomes_not_remapped = s.remap_snps(37)
        assert s.build == 37
        assert s.assembly == "GRCh37"
        assert len(chromosomes_remapped) == 2
        assert len(chromosomes_not_remapped) == 0
        pd.testing.assert_frame_equal(s.snps, self.snps_GRCh37())

    def test_remap_snps_36_to_37_multiprocessing(self):
        s = SNPs("tests/input/NCBI36.csv", parallelize=True)
        chromosomes_remapped, chromosomes_not_remapped = s.remap_snps(37)
        assert s.build == 37
        assert s.assembly == "GRCh37"
        assert len(chromosomes_remapped) == 2
        assert len(chromosomes_not_remapped) == 0
        pd.testing.assert_frame_equal(s.snps, self.snps_GRCh37())

    def test_remap_snps_37_to_36(self):
        s = SNPs("tests/input/GRCh37.csv")
        chromosomes_remapped, chromosomes_not_remapped = s.remap_snps(36)
        assert s.build == 36
        assert s.assembly == "NCBI36"
        assert len(chromosomes_remapped) == 2
        assert len(chromosomes_not_remapped) == 0
        pd.testing.assert_frame_equal(s.snps, self.snps_NCBI36())

    def test_remap_snps_37_to_38(self):
        s = SNPs("tests/input/GRCh37.csv")
        chromosomes_remapped, chromosomes_not_remapped = s.remap_snps(38)
        assert s.build == 38
        assert s.assembly == "GRCh38"
        assert len(chromosomes_remapped) == 2
        assert len(chromosomes_not_remapped) == 0
        pd.testing.assert_frame_equal(s.snps, self.snps_GRCh38_func())

    def test_remap_snps_37_to_38_with_PAR_SNP(self):
        if (
            not os.getenv("DOWNLOADS_ENABLED")
            or os.getenv("DOWNLOADS_ENABLED") == "true"
        ):
            s = SNPs("tests/input/GRCh37_PAR.csv")
            assert s.snp_count == 4
            chromosomes_remapped, chromosomes_not_remapped = s.remap_snps(38)
            assert s.build == 38
            assert s.assembly == "GRCh38"
            assert len(chromosomes_remapped) == 2
            assert len(chromosomes_not_remapped) == 1
            assert s.snp_count == 3
            pd.testing.assert_frame_equal(s.snps, self.snps_GRCh38_PAR())

    def test_remap_snps_37_to_37(self):
        s = SNPs("tests/input/GRCh37.csv")
        chromosomes_remapped, chromosomes_not_remapped = s.remap_snps(37)
        assert s.build == 37
        assert s.assembly == "GRCh37"
        assert len(chromosomes_remapped) == 0
        assert len(chromosomes_not_remapped) == 2
        pd.testing.assert_frame_equal(s.snps, self.snps_GRCh37())

    def test_remap_snps_no_snps(self):
        s = SNPs()
        chromosomes_remapped, chromosomes_not_remapped = s.remap_snps(38)
        assert not s.build
        assert len(chromosomes_remapped) == 0
        assert len(chromosomes_not_remapped) == 0

    def test_remap_snps_invalid_assembly(self):
        s = SNPs("tests/input/GRCh37.csv")
        chromosomes_remapped, chromosomes_not_remapped = s.remap_snps(-1)
        assert s.build == 37
        assert s.assembly == "GRCh37"
        assert len(chromosomes_remapped) == 0
        assert len(chromosomes_not_remapped) == 2

    def test___repr__snps(self):
        s = SNPs("tests/input/GRCh37.csv")
        assert "SNPs('tests/input/GRCh37.csv')" == s.__repr__()

    def test_load_opensnp_datadump_file(self):
        # temporarily set resources dir to tests
        r = Resources()
        r._resources_dir = "tests/resources"

        # write test openSNP datadump zip
        with atomic_write(
            "tests/resources/opensnp_datadump.current.zip", mode="wb", overwrite=True
        ) as f:
            with zipfile.ZipFile(f, "w") as f_zip:
                f_zip.write("tests/input/generic.csv", arcname="generic1.csv")
                f_zip.write("tests/input/generic.csv", arcname="generic2.csv")

        snps1 = SNPs(r.load_opensnp_datadump_file("generic1.csv"))
        snps2 = SNPs(r.load_opensnp_datadump_file("generic2.csv"))

        pd.testing.assert_frame_equal(snps1.snps, self.generic_snps())
        pd.testing.assert_frame_equal(snps2.snps, self.generic_snps())

        r._resources_dir = "resources"

    def test_heterozygous_snps(self):
        s = SNPs("tests/input/generic.csv")
        pd.testing.assert_frame_equal(
            s.heterozygous_snps(),
            self.create_snp_df(
                rsid=["rs6", "rs7", "rs8"],
                chrom=["1", "1", "1"],
                pos=[106, 107, 108],
                genotype=["GC", "TC", "AT"],
            ),
        )

    def test_not_null_snps(self):
        s = SNPs("tests/input/generic.csv")
        snps = self.generic_snps()
        snps.drop("rs5", inplace=True)
        pd.testing.assert_frame_equal(s.not_null_snps(), snps)
