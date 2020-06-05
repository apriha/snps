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
from unittest.mock import Mock, patch

import pandas as pd

from snps import SNPs
from tests import BaseSNPsTestCase


class TestSnps(BaseSNPsTestCase):
    @staticmethod
    def empty_snps():
        return [SNPs(), SNPs(b""), SNPs("tests/input/empty.txt")]

    def test___repr__snps(self):
        s = SNPs("tests/input/GRCh37.csv")
        self.assertEqual("SNPs('tests/input/GRCh37.csv')", s.__repr__())

    def test__lookup_build_with_snp_pos_None(self):
        snps = SNPs()
        snps._snps = self.create_snp_df(
            rsid=["rs3094315"], chrom=["1"], pos=[1], genotype=["AA"]
        )
        self.assertFalse(snps.detect_build())

    def test_assembly(self):
        s = SNPs("tests/input/GRCh38.csv")
        self.assertEqual(s.assembly, "GRCh38")

    def test_assembly_no_snps(self):
        for snps in self.empty_snps():
            self.assertFalse(snps.assembly)

    def test_build(self):
        s = SNPs("tests/input/NCBI36.csv")
        self.assertEqual(s.build, 36)
        self.assertEqual(s.assembly, "NCBI36")

    def test_build_detected_no_snps(self):
        for snps in self.empty_snps():
            self.assertFalse(snps.build_detected)

    def test_build_detected_PAR_snps(self):
        snps = self.load_assign_PAR_SNPs("tests/input/GRCh37_PAR.csv")
        self.assertEqual(snps.build, 37)
        self.assertTrue(snps.build_detected)
        pd.testing.assert_frame_equal(
            snps.snps, self.snps_GRCh37_PAR(), check_exact=True
        )

    def test_build_no_snps(self):
        for snps in self.empty_snps():
            self.assertFalse(snps.build)

    def test_chromosomes(self):
        s = SNPs("tests/input/chromosomes.csv")
        self.assertListEqual(s.chromosomes, ["1", "2", "3", "5", "PAR", "MT"])

    def test_chromosomes_no_snps(self):
        for snps in self.empty_snps():
            self.assertFalse(snps.chromosomes)

    def test_chromosomes_summary(self):
        s = SNPs("tests/input/chromosomes.csv")
        self.assertEqual(s.chromosomes_summary, "1-3, 5, PAR, MT")

    def test_chromosomes_summary_no_snps(self):
        for snps in self.empty_snps():
            self.assertFalse(snps.chromosomes_summary)

    def test_deduplicate_false(self):
        snps = SNPs("tests/input/duplicate_rsids.csv", deduplicate=False)

        result = self.create_snp_df(
            rsid=["rs1", "rs1", "rs1"],
            chrom=["1", "1", "1"],
            pos=[101, 102, 103],
            genotype=["AA", "CC", "GG"],
        )
        pd.testing.assert_frame_equal(snps.snps, result, check_exact=True)

    def test_duplicate_rsids(self):
        snps = SNPs("tests/input/duplicate_rsids.csv")
        result = self.create_snp_df(
            rsid=["rs1"], chrom=["1"], pos=[101], genotype=["AA"]
        )
        duplicate_snps = self.create_snp_df(
            rsid=["rs1", "rs1"], chrom=["1", "1"], pos=[102, 103], genotype=["CC", "GG"]
        )
        pd.testing.assert_frame_equal(snps.snps, result, check_exact=True)
        pd.testing.assert_frame_equal(
            snps.duplicate_snps, duplicate_snps, check_exact=True
        )

    def test_empty_dataframe(self):
        for snps in self.empty_snps():
            self.assertListEqual(
                list(snps.snps.columns.values), ["chrom", "pos", "genotype"]
            )
            self.assertEqual(snps.snps.index.name, "rsid")

    def test_get_assembly_None(self):
        snps = SNPs()
        self.assertFalse(snps.get_assembly())

    def test_get_summary(self):
        s = SNPs("tests/input/GRCh38.csv")
        self.assertDictEqual(
            s.get_summary(),
            {
                "source": "generic",
                "assembly": "GRCh38",
                "build": 38,
                "build_detected": True,
                "snp_count": 4,
                "chromosomes": "1, 3",
                "sex": "",
            },
        )

    def test_get_summary_no_snps(self):
        for snps in self.empty_snps():
            self.assertFalse(snps.get_summary())

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
            check_exact=True,
        )

    def test_is_valid_False(self):
        for snps in self.empty_snps():
            self.assertFalse(snps.is_valid())

    def test_is_valid_True(self):
        s = SNPs("tests/input/generic.csv")
        self.assertTrue(s.is_valid())

    def test_not_null_snps(self):
        s = SNPs("tests/input/generic.csv")
        snps = self.generic_snps()
        snps.drop("rs5", inplace=True)
        pd.testing.assert_frame_equal(s.not_null_snps(), snps, check_exact=True)

    def test_only_detect_source(self):
        s = SNPs("tests/input/generic.csv", only_detect_source=True)
        self.assertEqual(s.source, "generic")
        self.assertEqual(s.snp_count, 0)

    def _run_remap_test(self, f, mappings):
        if self.downloads_enabled:
            f()
        else:
            mock = Mock(return_value=mappings)
            with patch("snps.resources.Resources.get_assembly_mapping_data", mock):
                f()

    def test_remap_snps_36_to_37(self):
        def f():
            s = SNPs("tests/input/NCBI36.csv")
            chromosomes_remapped, chromosomes_not_remapped = s.remap_snps(37)
            self.assertEqual(s.build, 37)
            self.assertEqual(s.assembly, "GRCh37")
            self.assertEqual(len(chromosomes_remapped), 2)
            self.assertEqual(len(chromosomes_not_remapped), 0)
            pd.testing.assert_frame_equal(s.snps, self.snps_GRCh37(), check_exact=True)

        self._run_remap_test(f, self.NCBI36_GRCh37())

    def test_remap_snps_36_to_37_multiprocessing(self):
        def f():
            s = SNPs("tests/input/NCBI36.csv", parallelize=True)
            chromosomes_remapped, chromosomes_not_remapped = s.remap_snps(37)
            self.assertEqual(s.build, 37)
            self.assertEqual(s.assembly, "GRCh37")
            self.assertEqual(len(chromosomes_remapped), 2)
            self.assertEqual(len(chromosomes_not_remapped), 0)
            pd.testing.assert_frame_equal(s.snps, self.snps_GRCh37(), check_exact=True)

        self._run_remap_test(f, self.NCBI36_GRCh37())

    def test_remap_snps_37_to_36(self):
        def f():
            s = SNPs("tests/input/GRCh37.csv")
            chromosomes_remapped, chromosomes_not_remapped = s.remap_snps(36)
            self.assertEqual(s.build, 36)
            self.assertEqual(s.assembly, "NCBI36")
            self.assertEqual(len(chromosomes_remapped), 2)
            self.assertEqual(len(chromosomes_not_remapped), 0)
            pd.testing.assert_frame_equal(s.snps, self.snps_NCBI36(), check_exact=True)

        self._run_remap_test(f, self.GRCh37_NCBI36())

    def test_remap_snps_37_to_38(self):
        def f():
            s = SNPs("tests/input/GRCh37.csv")
            chromosomes_remapped, chromosomes_not_remapped = s.remap_snps(38)
            self.assertEqual(s.build, 38)
            self.assertEqual(s.assembly, "GRCh38")
            self.assertEqual(len(chromosomes_remapped), 2)
            self.assertEqual(len(chromosomes_not_remapped), 0)
            pd.testing.assert_frame_equal(s.snps, self.snps_GRCh38(), check_exact=True)

        self._run_remap_test(f, self.GRCh37_GRCh38())

    def test_remap_snps_37_to_38_with_PAR_SNP(self):
        def f():
            s = self.load_assign_PAR_SNPs("tests/input/GRCh37_PAR.csv")
            self.assertEqual(s.snp_count, 4)
            chromosomes_remapped, chromosomes_not_remapped = s.remap_snps(38)
            self.assertEqual(s.build, 38)
            self.assertEqual(s.assembly, "GRCh38")
            self.assertEqual(len(chromosomes_remapped), 2)
            self.assertEqual(len(chromosomes_not_remapped), 1)
            self.assertEqual(s.snp_count, 3)
            pd.testing.assert_frame_equal(
                s.snps, self.snps_GRCh38_PAR(), check_exact=True
            )

        self._run_remap_test(f, self.GRCh37_GRCh38_PAR())

    def test_remap_snps_37_to_37(self):
        s = SNPs("tests/input/GRCh37.csv")
        chromosomes_remapped, chromosomes_not_remapped = s.remap_snps(37)
        self.assertEqual(s.build, 37)
        self.assertEqual(s.assembly, "GRCh37")
        self.assertEqual(len(chromosomes_remapped), 0)
        self.assertEqual(len(chromosomes_not_remapped), 2)
        pd.testing.assert_frame_equal(s.snps, self.snps_GRCh37(), check_exact=True)

    def test_remap_snps_invalid_assembly(self):
        s = SNPs("tests/input/GRCh37.csv")
        chromosomes_remapped, chromosomes_not_remapped = s.remap_snps(-1)
        self.assertEqual(s.build, 37)
        self.assertEqual(s.assembly, "GRCh37")
        self.assertEqual(len(chromosomes_remapped), 0)
        self.assertEqual(len(chromosomes_not_remapped), 2)

    def test_remap_snps_no_snps(self):
        s = SNPs()
        chromosomes_remapped, chromosomes_not_remapped = s.remap_snps(38)
        self.assertFalse(s.build)
        self.assertEqual(len(chromosomes_remapped), 0)
        self.assertEqual(len(chromosomes_not_remapped), 0)

    def test_save_snps_buffer(self):
        s = SNPs("tests/input/generic.csv")
        out = io.StringIO()
        s.save_snps(out)
        self.assertTrue(out.read().startswith("# Generated by snps"))

    def test_save_snps_no_snps(self):
        s = SNPs()
        self.assertFalse(s.save_snps())

    def test_save_snps_no_snps_vcf(self):
        s = SNPs()
        self.assertFalse(s.save_snps(vcf=True))

    def test_save_snps_source(self):
        s = SNPs("tests/input/GRCh38.csv")
        self.assertEqual(os.path.relpath(s.save_snps()), "output/generic_GRCh38.csv")
        snps = SNPs("output/generic_GRCh38.csv")
        pd.testing.assert_frame_equal(snps.snps, self.snps_GRCh38(), check_exact=True)

    def test_sex_Female_X_chrom(self):
        s = self.simulate_snps(
            chrom="X", pos_start=1, pos_max=155270560, pos_step=10000, genotype="AC"
        )
        self.assertEqual(s.sex, "Female")

    def test_sex_Female_Y_chrom(self):
        s = self.simulate_snps(
            chrom="Y", pos_start=1, pos_max=59373566, pos_step=10000, null_snp_step=1
        )
        self.assertEqual(s.sex, "Female")

    def test_sex_Male_X_chrom(self):
        s = self.simulate_snps(
            chrom="X", pos_start=1, pos_max=155270560, pos_step=10000, genotype="AA"
        )
        self.assertEqual(s.snp_count, 15528)
        s._deduplicate_XY_chrom()
        self.assertEqual(s.snp_count, 15528)
        self.assertEqual(len(s.discrepant_XY_snps), 0)
        self.assertEqual(s.sex, "Male")

    def test_sex_Male_X_chrom_discrepant_XY_snps(self):
        s = self.simulate_snps(
            chrom="X", pos_start=1, pos_max=155270560, pos_step=10000, genotype="AA"
        )
        self.assertEqual(s.snp_count, 15528)
        s._snps.loc["rs8001", "genotype"] = "AC"
        s._deduplicate_XY_chrom()
        self.assertEqual(s.snp_count, 15527)
        result = self.create_snp_df(
            rsid=["rs8001"], chrom=["X"], pos=[80000001], genotype=["AC"]
        )
        pd.testing.assert_frame_equal(s.discrepant_XY_snps, result, check_exact=True)
        self.assertEqual(s.sex, "Male")

    def test_sex_Male_Y_chrom(self):
        s = self.simulate_snps(chrom="Y", pos_start=1, pos_max=59373566, pos_step=10000)
        self.assertEqual(s.sex, "Male")

    def test_sex_not_determined(self):
        s = self.simulate_snps(
            chrom="1", pos_start=1, pos_max=249250621, pos_step=10000
        )
        self.assertEqual(s.sex, "")

    def test_sex_no_snps(self):
        for snps in self.empty_snps():
            self.assertFalse(snps.sex)

    def test_source(self):
        s = SNPs("tests/input/generic.csv")
        self.assertEqual(s.source, "generic")

    def test_source_no_snps(self):
        for snps in self.empty_snps():
            self.assertFalse(snps.source)

    def test_snp_count(self):
        s = SNPs("tests/input/NCBI36.csv")
        self.assertEqual(s.snp_count, 4)

    def test_snp_count_no_snps(self):
        for snps in self.empty_snps():
            self.assertFalse(snps.snp_count)
            self.assertTrue(snps.snps.empty)
