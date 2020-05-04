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

import io
import os
import tempfile

from atomicwrites import atomic_write
import numpy as np
import pandas as pd

from snps import SNPs
from tests import BaseSNPsTestCase


class TestSnps(BaseSNPsTestCase):
    def setUp(self):
        self.snps_GRCh38 = SNPs("tests/input/GRCh38.csv")
        self.snps = SNPs("tests/input/chromosomes.csv")
        self.snps_only_detect_source = SNPs(
            "tests/input/chromosomes.csv", only_detect_source=True
        )
        self.snps_none = SNPs(None)

    def snps_discrepant_pos(self):
        return self.create_snp_df(
            rsid=["rs3094315"], chrom=["1"], pos=[1], genotype=["AA"]
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

    def test___repr__snps(self):
        s = SNPs("tests/input/GRCh37.csv")
        assert "SNPs('tests/input/GRCh37.csv')" == s.__repr__()

    def test__lookup_build_with_snp_pos_None(self):
        snps = SNPs()
        snps._snps = self.snps_discrepant_pos()
        assert not snps.detect_build()

    def test__read_raw_data(self):
        assert self.snps_none.snps.empty
        assert self.snps_none.source == ""

    def test_assembly(self):
        assert self.snps_GRCh38.assembly == "GRCh38"

    def test_assembly_no_snps(self):
        assert self.snps_none.assembly == ""

    def test_build(self):
        s = SNPs("tests/input/NCBI36.csv")
        assert s.build == 36
        assert s.assembly == "NCBI36"

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

    def test_build_no_snps(self):
        assert not self.snps_none.build

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

    def test_deduplicate_false(self):
        snps = SNPs("tests/input/duplicate_rsids.csv", deduplicate=False)

        result = self.create_snp_df(
            rsid=["rs1", "rs1", "rs1"],
            chrom=["1", "1", "1"],
            pos=[101, 102, 103],
            genotype=["AA", "CC", "GG"],
        )
        pd.testing.assert_frame_equal(snps.snps, result)

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

    def test_empty_dataframe(self):
        s = SNPs()
        assert s.snp_count == 0
        assert list(s.snps.columns.values) == ["chrom", "pos", "genotype"]
        assert s.snps.index.name == "rsid"

    def test_empty_dataframe_file(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            dest = os.path.join(tmpdir, "empty.txt")

            with atomic_write(dest, mode="w", overwrite=True):
                pass

            s = SNPs(dest)
            assert s.snp_count == 0
            assert list(s.snps.columns.values) == ["chrom", "pos", "genotype"]
            assert s.snps.index.name == "rsid"

    def test_get_assembly_None(self):
        snps = SNPs()
        assert snps.get_assembly() is ""

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

    def test_is_valid_False(self):
        assert not self.snps_none.is_valid()

    def test_is_valid_True(self):
        assert self.snps_GRCh38.is_valid()

    def test_no_snps(self):
        s = SNPs()
        assert s.snps.empty

    def test_not_null_snps(self):
        s = SNPs("tests/input/generic.csv")
        snps = self.generic_snps()
        snps.drop("rs5", inplace=True)
        pd.testing.assert_frame_equal(s.not_null_snps(), snps)

    def test_only_detect_source(self):
        assert self.snps_only_detect_source.source == "generic"

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

    def test_remap_snps_invalid_assembly(self):
        s = SNPs("tests/input/GRCh37.csv")
        chromosomes_remapped, chromosomes_not_remapped = s.remap_snps(-1)
        assert s.build == 37
        assert s.assembly == "GRCh37"
        assert len(chromosomes_remapped) == 0
        assert len(chromosomes_not_remapped) == 2

    def test_remap_snps_no_snps(self):
        s = SNPs()
        chromosomes_remapped, chromosomes_not_remapped = s.remap_snps(38)
        assert not s.build
        assert len(chromosomes_remapped) == 0
        assert len(chromosomes_not_remapped) == 0

    def test_save_snps_buffer(self):
        out = io.StringIO()
        self.snps.save_snps(out)
        assert out.read().startswith("# Generated by snps")

    def test_save_snps_no_snps(self):
        s = SNPs()
        assert not s.save_snps()

    def test_save_snps_no_snps_vcf(self):
        s = SNPs()
        assert not s.save_snps(vcf=True)

    def test_save_snps_source(self):
        assert (
            os.path.relpath(self.snps_GRCh38.save_snps()) == "output/generic_GRCh38.csv"
        )
        snps = SNPs("output/generic_GRCh38.csv")
        pd.testing.assert_frame_equal(snps.snps, self.snps_GRCh38.snps)

    def test_sex_Female_X_chrom(self):
        s = self.simulate_snps(
            chrom="X", pos_start=1, pos_max=155270560, pos_step=10000, genotype="AC"
        )
        assert s.sex == "Female"

    def test_sex_Female_Y_chrom(self):
        s = self.simulate_snps(
            chrom="Y", pos_start=1, pos_max=59373566, pos_step=10000, null_snp_step=1
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

    def test_sex_Male_Y_chrom(self):
        s = self.simulate_snps(chrom="Y", pos_start=1, pos_max=59373566, pos_step=10000)
        assert s.sex == "Male"

    def test_sex_not_determined(self):
        s = self.simulate_snps(
            chrom="1", pos_start=1, pos_max=249250621, pos_step=10000
        )
        assert s.sex == ""

    def test_sex_no_snps(self):
        assert self.snps_none.sex == ""

    def test_snp_count(self):
        s = SNPs("tests/input/NCBI36.csv")
        assert s.snp_count == 4

    def test_snp_count_no_snps(self):
        s = SNPs()
        assert s.snp_count == 0
