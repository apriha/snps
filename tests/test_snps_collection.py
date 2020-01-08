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

import gzip
import os
import shutil
import zipfile

from atomicwrites import atomic_write
import numpy as np
import pandas as pd

from snps import SNPs, SNPsCollection
from snps.resources import Resources, ReferenceSequence
from snps.utils import create_dir
from tests import BaseSNPsTestCase


class TestSNPsCollection(BaseSNPsTestCase):
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

    def snps_NCBI36(self):
        return self.create_snp_df(
            rsid=["rs3094315", "rs2500347", "rsIndelTest", "rs11928389"],
            chrom=["1", "1", "1", "3"],
            pos=[742429, 143649677, 143649678, 50908372],
            genotype=["AA", np.nan, "ID", "AG"],
        )

    def snps_NCBI36_discrepant_snps(self):
        return self.create_snp_df(
            rsid=["rs3094315", "rs2500347", "rsIndelTest", "rs11928389"],
            chrom=["1", "1", "1", "3"],
            pos=[742429, 143649677, 143649678, 50908372],
            genotype=["AA", np.nan, "ID", np.nan],
        )

    def snps_GRCh37(self):
        return self.create_snp_df(
            rsid=["rs3094315", "rs2500347", "rsIndelTest", "rs11928389"],
            chrom=["1", "1", "1", "3"],
            pos=[752566, 144938320, 144938321, 50927009],
            genotype=["AA", np.nan, "ID", "TC"],
        )

    def snps_GRCh38(self):
        return self.create_snp_df(
            rsid=["rs3094315", "rsIndelTest", "rs2500347", "rs11928389"],
            chrom=["1", "1", "1", "3"],
            pos=[817186, 148946168, 148946169, 50889578],
            genotype=["AA", "ID", np.nan, "TC"],
        )

    def snps_GRCh38_PAR(self):
        return self.create_snp_df(
            rsid=["rs28736870", "rs113313554"],
            chrom=["X", "Y"],
            pos=[304103, 624523],
            genotype=["AA", "AA"],
        )

    def test_snps_23andme(self):
        # https://www.23andme.com
        s = SNPs("tests/input/23andme.txt")
        assert s.source == "23andMe"
        pd.testing.assert_frame_equal(s.snps, self.generic_snps())

    def test_snps_23andme_zip(self):
        with atomic_write(
            "tests/input/23andme.txt.zip", mode="wb", overwrite=True
        ) as f:
            with zipfile.ZipFile(f, "w") as f_zip:
                # https://stackoverflow.com/a/16104667
                f_zip.write("tests/input/23andme.txt", arcname="23andme.txt")
        s = SNPs("tests/input/23andme.txt.zip")
        assert s.source == "23andMe"
        pd.testing.assert_frame_equal(s.snps, self.generic_snps())

    def test_snps_ftdna(self):
        # https://www.familytreedna.com
        s = SNPs("tests/input/ftdna.csv")
        assert s.source == "FTDNA"
        pd.testing.assert_frame_equal(s.snps, self.generic_snps())

    def test_snps_ftdna_gzip(self):
        with open("tests/input/ftdna.csv", "rb") as f_in:
            with atomic_write(
                "tests/input/ftdna.csv.gz", mode="wb", overwrite=True
            ) as f_out:
                with gzip.open(f_out, "wb") as f_gzip:
                    shutil.copyfileobj(f_in, f_gzip)
        s = SNPs("tests/input/ftdna.csv.gz")
        assert s.source == "FTDNA"
        pd.testing.assert_frame_equal(s.snps, self.generic_snps())

    def test_snps_ftdna_famfinder(self):
        # https://www.familytreedna.com
        s = SNPs("tests/input/ftdna_famfinder.csv")
        assert s.source == "FTDNA"
        pd.testing.assert_frame_equal(s.snps, self.generic_snps())

    def test_snps_ancestry(self):
        # https://www.ancestry.com
        s = SNPs("tests/input/ancestry.txt")
        assert s.source == "AncestryDNA"
        pd.testing.assert_frame_equal(s.snps, self.generic_snps())

    def test_snps_genes_for_good(self):
        # https://genesforgood.sph.umich.edu/
        s = SNPs("tests/input/genesforgood.txt")
        assert s.source == "GenesForGood"
        pd.testing.assert_frame_equal(s.snps, self.generic_snps())

    def test_snps_myheritage(self):
        # https://www.myheritage.com
        s = SNPs("tests/input/myheritage.csv")
        assert s.source == "MyHeritage"
        pd.testing.assert_frame_equal(s.snps, self.generic_snps())

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

    def test_snps_codigo46_bytes(self):
        # https://codigo46.com.mx

        self._setup_gsa_test()

        with open("tests/input/codigo46.txt", "rb") as f:
            s = SNPs(f.read())
        assert s.source == "Codigo46"
        pd.testing.assert_frame_equal(s.snps, self.generic_snps())

        self._teardown_gsa_test()

    def test_snps_codigo46(self):
        # https://codigo46.com.mx

        self._setup_gsa_test()

        s = SNPs("tests/input/codigo46.txt")
        assert s.source == "Codigo46"
        pd.testing.assert_frame_equal(s.snps, self.generic_snps())

        self._teardown_gsa_test()

    def test_snps_sano_bytes(self):
        # https://sanogenetics.com

        self._setup_gsa_test()

        with open("tests/input/sano.txt", "rb") as f:
            s = SNPs(f.read())
        assert s.source == "Sano"
        pd.testing.assert_frame_equal(s.snps, self.generic_snps())

        self._teardown_gsa_test()

    def test_snps_sano(self):
        # https://sanogenetics.com

        self._setup_gsa_test()

        s = SNPs("tests/input/sano.txt")
        assert s.source == "Sano"
        pd.testing.assert_frame_equal(s.snps, self.generic_snps())

        self._teardown_gsa_test()

    def test_snps_livingdna(self):
        # https://livingdna.com
        s = SNPs("tests/input/livingdna.csv")
        assert s.source == "LivingDNA"
        pd.testing.assert_frame_equal(s.snps, self.generic_snps())

    def test_snps_mapmygenome(self):
        # https://mapmygenome.in
        s = SNPs("tests/input/mapmygenome.txt")
        assert s.source == "Mapmygenome"
        pd.testing.assert_frame_equal(s.snps, self.generic_snps())

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

    def test_source_lineage_file(self):
        sc = SNPsCollection("tests/input/GRCh37.csv")
        assert sc.source == "generic"
        sc.load_snps("tests/input/23andme.txt")
        assert sc.source == "generic, 23andMe"
        file = sc.save_snps()
        s = SNPs(file)
        assert s.source == "generic, 23andMe"
        pd.testing.assert_frame_equal(sc.snps, s.snps)

    def test_source_lineage_file_gzip(self):
        sc = SNPsCollection("tests/input/GRCh37.csv")
        assert sc.source == "generic"
        sc.load_snps("tests/input/23andme.txt")
        assert sc.source == "generic, 23andMe"
        file = sc.save_snps()
        with open(file, "rb") as f_in:
            with atomic_write(file + ".gz", mode="wb", overwrite=True) as f_out:
                with gzip.open(f_out, "wb") as f_gzip:
                    shutil.copyfileobj(f_in, f_gzip)
        s = SNPs(file + ".gz")
        assert s.source == "generic, 23andMe"
        pd.testing.assert_frame_equal(sc.snps, s.snps)

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

    def test_load_snps_list(self):
        sc = SNPsCollection()
        sc.load_snps(["tests/input/GRCh37.csv", "tests/input/GRCh37.csv"])
        pd.testing.assert_frame_equal(sc.snps, self.snps_GRCh37())
        assert sc.source == "generic, generic"

    def test_load_snps_None(self):
        sc = SNPsCollection()
        with self.assertRaises(TypeError):
            sc.load_snps(None)

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

    def test_discrepant_positions(self):
        sc = SNPsCollection()
        sc.load_snps(["tests/input/NCBI36.csv", "tests/input/GRCh37.csv"])
        assert len(sc.discrepant_positions) == 4

    def test_discrepant_genotypes(self):
        sc = SNPsCollection()
        sc.load_snps(["tests/input/NCBI36.csv", "tests/input/GRCh37.csv"])
        assert len(sc.discrepant_genotypes) == 1

    def test_discrepant_snps(self):
        sc = SNPsCollection()
        sc.load_snps(["tests/input/NCBI36.csv", "tests/input/GRCh37.csv"])
        assert len(sc.discrepant_snps) == 4

    def test_load_snps_non_existent_file(self):
        sc = SNPsCollection()
        sc.load_snps(["tests/input/GRCh37.csv", "tests/input/non_existent_file.csv"])
        pd.testing.assert_frame_equal(sc.snps, self.snps_GRCh37())

    def test_load_snps_invalid_file(self):
        sc = SNPsCollection()
        with atomic_write("tests/input/empty.txt", mode="w", overwrite=True):
            pass
        sc.load_snps(["tests/input/GRCh37.csv", "tests/input/empty.txt"])
        pd.testing.assert_frame_equal(sc.snps, self.snps_GRCh37())

    def test_load_snps_assembly_mismatch(self):
        sc = SNPsCollection()
        sc.load_snps(["tests/input/NCBI36.csv", "tests/input/GRCh37.csv"])
        assert not os.path.exists("output/ind_discrepant_positions_1.csv")
        assert not os.path.exists("output/ind_discrepant_genotypes_1.csv")
        assert len(sc.discrepant_positions) == 4
        assert len(sc.discrepant_genotypes) == 1
        pd.testing.assert_frame_equal(sc.snps, self.snps_NCBI36_discrepant_snps())

    def test_load_snps_assembly_mismatch_save_output(self):
        sc = SNPsCollection()
        sc.load_snps(
            ["tests/input/NCBI36.csv", "tests/input/GRCh37.csv"], save_output=True
        )
        assert os.path.exists("output/discrepant_positions_1.csv")
        assert os.path.exists("output/discrepant_genotypes_1.csv")
        assert len(sc.discrepant_positions) == 4
        assert len(sc.discrepant_genotypes) == 1
        pd.testing.assert_frame_equal(sc.snps, self.snps_NCBI36_discrepant_snps())

    def test_load_snps_assembly_mismatch_exceed_discrepant_positions_threshold(self):
        sc = SNPsCollection()
        sc.load_snps(
            ["tests/input/NCBI36.csv", "tests/input/GRCh37.csv"],
            discrepant_snp_positions_threshold=0,
        )
        assert not os.path.exists("output/discrepant_positions_1.csv")
        assert not os.path.exists("output/discrepant_genotypes_1.csv")
        assert len(sc.discrepant_positions) == 4
        assert len(sc.discrepant_genotypes) == 0
        pd.testing.assert_frame_equal(sc.snps, self.snps_NCBI36())

    def test_load_snps_assembly_mismatch_exceed_discrepant_genotypes_threshold(self):
        sc = SNPsCollection()
        sc.load_snps(
            ["tests/input/NCBI36.csv", "tests/input/GRCh37.csv"],
            discrepant_genotypes_threshold=0,
        )
        assert not os.path.exists("output/discrepant_positions_1.csv")
        assert not os.path.exists("output/discrepant_genotypes_1.csv")
        assert len(sc.discrepant_positions) == 4
        assert len(sc.discrepant_genotypes) == 1
        pd.testing.assert_frame_equal(sc.snps, self.snps_NCBI36())

    def test_merging_files_discrepant_snps(self):
        df = pd.read_csv(
            "tests/input/discrepant_snps.csv",
            skiprows=1,
            na_values="--",
            names=[
                "rsid",
                "chrom",
                "pos_file1",
                "pos_file2",
                "genotype_file1",
                "genotype_file2",
                "discrepant_position",
                "discrepant_genotype",
                "expected_position",
                "expected_genotype",
            ],
            index_col=0,
            dtype={
                "chrom": object,
                "pos_file1": np.int64,
                "pos_file2": np.int64,
                "discrepant_position": bool,
                "discrepant_genotype": bool,
            },
        )

        df1 = df[["chrom", "pos_file1", "genotype_file1"]]
        df2 = df[["chrom", "pos_file2", "genotype_file2"]]

        df1.to_csv(
            "tests/input/discrepant_snps1.csv",
            na_rep="--",
            header=["chromosome", "position", "genotype"],
        )

        df2.to_csv(
            "tests/input/discrepant_snps2.csv",
            na_rep="--",
            header=["chromosome", "position", "genotype"],
        )

        sc = SNPsCollection(
            ["tests/input/discrepant_snps1.csv", "tests/input/discrepant_snps2.csv"]
        )

        expected = df[
            [
                "chrom",
                "discrepant_position",
                "discrepant_genotype",
                "expected_position",
                "expected_genotype",
            ]
        ]
        expected = expected.rename(
            columns={"expected_position": "pos", "expected_genotype": "genotype"}
        )
        expected_snps = SNPs()
        expected_snps._snps = expected
        expected_snps.sort_snps()
        expected = expected_snps.snps

        pd.testing.assert_index_equal(
            sc.discrepant_positions.index,
            expected.loc[expected["discrepant_position"] == True].index,
        )

        pd.testing.assert_index_equal(
            sc.discrepant_genotypes.index,
            expected.loc[expected["discrepant_genotype"] == True].index,
        )

        pd.testing.assert_series_equal(sc.snps["pos"], expected["pos"])
        pd.testing.assert_series_equal(sc.snps["genotype"], expected["genotype"])

    def test_save_snps(self):
        snps = SNPs("tests/input/GRCh37.csv")
        assert os.path.relpath(snps.save_snps()) == "output/generic_GRCh37.csv"
        s_saved = SNPs("output/generic_GRCh37.csv")
        pd.testing.assert_frame_equal(s_saved.snps, self.snps_GRCh37())

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

    def test_save_discrepant_positions(self):
        sc = SNPsCollection()
        sc.load_snps(["tests/input/NCBI36.csv", "tests/input/GRCh37.csv"])
        assert len(sc.discrepant_positions) == 4
        discrepant_positions_file = sc.save_discrepant_positions()
        assert (
            os.path.relpath(discrepant_positions_file)
            == "output/discrepant_positions.csv"
        )
        assert os.path.exists(discrepant_positions_file)

    def test_save_discrepant_positions_specify_file(self):
        sc = SNPsCollection()
        sc.load_snps(["tests/input/NCBI36.csv", "tests/input/GRCh37.csv"])
        assert len(sc.discrepant_positions) == 4
        discrepant_positions_file = sc.save_discrepant_positions(
            "discrepant_positions.csv"
        )
        assert (
            os.path.relpath(discrepant_positions_file)
            == "output/discrepant_positions.csv"
        )
        assert os.path.exists(discrepant_positions_file)

    def test_save_discrepant_positions_no_discrepant_snps(self):
        sc = SNPsCollection()
        assert len(sc.discrepant_positions) == 0
        assert not sc.save_discrepant_positions()

    def test_save_discrepant_positions_exception(self):
        sc = SNPsCollection()
        sc._discrepant_positions = "invalid"
        assert not sc.save_discrepant_positions()

    def test_save_discrepant_genotypes(self):
        sc = SNPsCollection()
        sc.load_snps(["tests/input/NCBI36.csv", "tests/input/GRCh37.csv"])
        assert len(sc.discrepant_genotypes) == 1
        discrepant_genotypes_file = sc.save_discrepant_genotypes()
        assert (
            os.path.relpath(discrepant_genotypes_file)
            == "output/discrepant_genotypes.csv"
        )
        assert os.path.exists(discrepant_genotypes_file)

    def test_save_discrepant_genotypes_specify_file(self):
        sc = SNPsCollection()
        sc.load_snps(["tests/input/NCBI36.csv", "tests/input/GRCh37.csv"])
        assert len(sc.discrepant_genotypes) == 1
        discrepant_genotypes_file = sc.save_discrepant_genotypes(
            "discrepant_genotypes.csv"
        )
        assert (
            os.path.relpath(discrepant_genotypes_file)
            == "output/discrepant_genotypes.csv"
        )
        assert os.path.exists(discrepant_genotypes_file)

    def test_save_discrepant_genotypes_no_discrepant_snps(self):
        sc = SNPsCollection()
        assert len(sc.discrepant_genotypes) == 0
        assert not sc.save_discrepant_genotypes()

    def test_save_discrepant_genotypes_exception(self):
        sc = SNPsCollection()
        sc._discrepant_genotypes = "invalid"
        assert not sc.save_discrepant_genotypes()

    def test_save_discrepant_snps(self):
        sc = SNPsCollection()
        sc.load_snps(["tests/input/NCBI36.csv", "tests/input/GRCh37.csv"])
        assert len(sc.discrepant_snps) == 4
        discrepant_snps_file = sc.save_discrepant_snps()
        assert os.path.relpath(discrepant_snps_file) == "output/discrepant_snps.csv"
        assert os.path.exists(discrepant_snps_file)

    def test_save_discrepant_snps_specify_file(self):
        sc = SNPsCollection()
        sc.load_snps(["tests/input/NCBI36.csv", "tests/input/GRCh37.csv"])
        assert len(sc.discrepant_snps) == 4
        discrepant_snps_file = sc.save_discrepant_snps("discrepant_snps.csv")
        assert os.path.relpath(discrepant_snps_file) == "output/discrepant_snps.csv"
        assert os.path.exists(discrepant_snps_file)

    def test_save_discrepant_snps_no_discrepant_snps(self):
        sc = SNPsCollection()
        assert len(sc.discrepant_snps) == 0
        assert not sc.save_discrepant_snps()

    def test_save_discrepant_snps_exception(self):
        sc = SNPsCollection()
        sc._discrepant_snps = "invalid"
        assert not sc.save_discrepant_snps()

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
        pd.testing.assert_frame_equal(s.snps, self.snps_GRCh38())

    def test_remap_snps_37_to_38_with_PAR_SNP(self):
        if (
            not os.getenv("DOWNLOADS_ENABLED")
            or os.getenv("DOWNLOADS_ENABLED") == "true"
        ):
            s = SNPs("tests/input/GRCh37_PAR.csv")
            assert s.snp_count == 3
            chromosomes_remapped, chromosomes_not_remapped = s.remap_snps(38)
            assert s.build == 38
            assert s.assembly == "GRCh38"
            assert len(chromosomes_remapped) == 2
            assert len(chromosomes_not_remapped) == 1
            assert s.snp_count == 2
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

    def test___repr__snps_collection(self):
        sc = SNPsCollection()
        assert "SNPsCollection(name='')" == sc.__repr__()

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
