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
from unittest.mock import Mock, patch
import warnings

import numpy as np
import pandas as pd

from snps import SNPs
from snps.io import get_empty_snps_dataframe
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

    def test_deduplicate_MT_chrom(self):
        snps = SNPs("tests/input/ancestry_mt.txt")

        result = self.create_snp_df(
            rsid=["rs1", "rs2"],
            chrom=["MT", "MT"],
            pos=[101, 102],
            genotype=["A", np.nan],
        )
        pd.testing.assert_frame_equal(snps.snps, result, check_exact=True)

        heterozygous_MT_snps = self.create_snp_df(
            rsid=["rs3"], chrom=["MT"], pos=[103], genotype=["GC"]
        )
        pd.testing.assert_frame_equal(
            snps.heterozygous_MT, heterozygous_MT_snps, check_exact=True
        )

    def test_deduplicate_MT_chrom_false(self):
        snps = SNPs("tests/input/ancestry_mt.txt", deduplicate_MT_chrom=False)

        result = self.create_snp_df(
            rsid=["rs1", "rs2", "rs3"],
            chrom=["MT", "MT", "MT"],
            pos=[101, 102, 103],
            genotype=["AA", np.nan, "GC"],
        )
        pd.testing.assert_frame_equal(snps.snps, result, check_exact=True)

    def test_duplicate_rsids(self):
        snps = SNPs("tests/input/duplicate_rsids.csv")
        result = self.create_snp_df(
            rsid=["rs1"], chrom=["1"], pos=[101], genotype=["AA"]
        )
        duplicate = self.create_snp_df(
            rsid=["rs1", "rs1"], chrom=["1", "1"], pos=[102, 103], genotype=["CC", "GG"]
        )
        pd.testing.assert_frame_equal(snps.snps, result, check_exact=True)
        pd.testing.assert_frame_equal(snps.duplicate, duplicate, check_exact=True)

    def test_empty_dataframe(self):
        for snps in self.empty_snps():
            self.assertListEqual(
                list(snps.snps.columns.values), ["chrom", "pos", "genotype"]
            )
            self.assertEqual(snps.snps.index.name, "rsid")

    def test_assembly_None(self):
        snps = SNPs()
        self.assertFalse(snps.assembly)

    def test_summary(self):
        s = SNPs("tests/input/GRCh38.csv")
        self.assertDictEqual(
            s.summary,
            {
                "source": "generic",
                "assembly": "GRCh38",
                "build": 38,
                "build_detected": True,
                "count": 4,
                "chromosomes": "1, 3",
                "sex": "",
            },
        )

    def test_summary_no_snps(self):
        for snps in self.empty_snps():
            self.assertDictEqual(snps.summary, {})

    def test_heterozygous(self):
        s = SNPs("tests/input/generic.csv")
        pd.testing.assert_frame_equal(
            s.heterozygous(),
            self.create_snp_df(
                rsid=["rs6", "rs7", "rs8"],
                chrom=["1", "1", "1"],
                pos=[106, 107, 108],
                genotype=["GC", "TC", "AT"],
            ),
            check_exact=True,
        )

    def test_homozygous(self):
        s = SNPs("tests/input/generic.csv")
        pd.testing.assert_frame_equal(
            s.homozygous(),
            self.create_snp_df(
                rsid=["rs1", "rs2", "rs3", "rs4"],
                chrom=["1", "1", "1", "1"],
                pos=[101, 102, 103, 104],
                genotype=["AA", "CC", "GG", "TT"],
            ),
            check_exact=True,
        )

    def test_homozygous_chrom(self):
        s = SNPs("tests/input/generic.csv")
        pd.testing.assert_frame_equal(
            s.homozygous("1"),
            self.create_snp_df(
                rsid=["rs1", "rs2", "rs3", "rs4"],
                chrom=["1", "1", "1", "1"],
                pos=[101, 102, 103, 104],
                genotype=["AA", "CC", "GG", "TT"],
            ),
            check_exact=True,
        )

    def test_valid_False(self):
        for snps in self.empty_snps():
            self.assertFalse(snps.valid)

    def test_valid_True(self):
        s = SNPs("tests/input/generic.csv")
        self.assertTrue(s.valid)

    def test_notnull(self):
        s = SNPs("tests/input/generic.csv")
        snps = self.generic_snps()
        snps.drop("rs5", inplace=True)
        pd.testing.assert_frame_equal(s.notnull(), snps, check_exact=True)

    def test_only_detect_source(self):
        s = SNPs("tests/input/generic.csv", only_detect_source=True)
        self.assertEqual(s.source, "generic")
        self.assertEqual(s.count, 0)

    def _run_remap_test(self, f, mappings):
        if self.downloads_enabled:
            f()
        else:
            mock = Mock(return_value=mappings)
            with patch("snps.resources.Resources.get_assembly_mapping_data", mock):
                f()

    def test_remap_36_to_37(self):
        def f():
            s = SNPs("tests/input/NCBI36.csv")
            chromosomes_remapped, chromosomes_not_remapped = s.remap(37)
            self.assertEqual(s.build, 37)
            self.assertEqual(s.assembly, "GRCh37")
            self.assertEqual(len(chromosomes_remapped), 2)
            self.assertEqual(len(chromosomes_not_remapped), 0)
            pd.testing.assert_frame_equal(s.snps, self.snps_GRCh37(), check_exact=True)

        self._run_remap_test(f, self.NCBI36_GRCh37())

    def test_remap_36_to_37_multiprocessing(self):
        def f():
            s = SNPs("tests/input/NCBI36.csv", parallelize=True)
            chromosomes_remapped, chromosomes_not_remapped = s.remap(37)
            self.assertEqual(s.build, 37)
            self.assertEqual(s.assembly, "GRCh37")
            self.assertEqual(len(chromosomes_remapped), 2)
            self.assertEqual(len(chromosomes_not_remapped), 0)
            pd.testing.assert_frame_equal(s.snps, self.snps_GRCh37(), check_exact=True)

        self._run_remap_test(f, self.NCBI36_GRCh37())

    def test_remap_37_to_36(self):
        def f():
            s = SNPs("tests/input/GRCh37.csv")
            chromosomes_remapped, chromosomes_not_remapped = s.remap(36)
            self.assertEqual(s.build, 36)
            self.assertEqual(s.assembly, "NCBI36")
            self.assertEqual(len(chromosomes_remapped), 2)
            self.assertEqual(len(chromosomes_not_remapped), 0)
            pd.testing.assert_frame_equal(s.snps, self.snps_NCBI36(), check_exact=True)

        self._run_remap_test(f, self.GRCh37_NCBI36())

    def test_remap_37_to_38(self):
        def f():
            s = SNPs("tests/input/GRCh37.csv")
            chromosomes_remapped, chromosomes_not_remapped = s.remap(38)
            self.assertEqual(s.build, 38)
            self.assertEqual(s.assembly, "GRCh38")
            self.assertEqual(len(chromosomes_remapped), 2)
            self.assertEqual(len(chromosomes_not_remapped), 0)
            pd.testing.assert_frame_equal(s.snps, self.snps_GRCh38(), check_exact=True)

        self._run_remap_test(f, self.GRCh37_GRCh38())

    def test_remap_37_to_38_with_PAR_SNP(self):
        def f():
            s = self.load_assign_PAR_SNPs("tests/input/GRCh37_PAR.csv")
            self.assertEqual(s.count, 4)
            chromosomes_remapped, chromosomes_not_remapped = s.remap(38)
            self.assertEqual(s.build, 38)
            self.assertEqual(s.assembly, "GRCh38")
            self.assertEqual(len(chromosomes_remapped), 2)
            self.assertEqual(len(chromosomes_not_remapped), 1)
            self.assertEqual(s.count, 3)
            pd.testing.assert_frame_equal(
                s.snps, self.snps_GRCh38_PAR(), check_exact=True
            )

        self._run_remap_test(f, self.GRCh37_GRCh38_PAR())

    def test_remap_37_to_37(self):
        s = SNPs("tests/input/GRCh37.csv")
        chromosomes_remapped, chromosomes_not_remapped = s.remap(37)
        self.assertEqual(s.build, 37)
        self.assertEqual(s.assembly, "GRCh37")
        self.assertEqual(len(chromosomes_remapped), 0)
        self.assertEqual(len(chromosomes_not_remapped), 2)
        pd.testing.assert_frame_equal(s.snps, self.snps_GRCh37(), check_exact=True)

    def test_remap_invalid_assembly(self):
        s = SNPs("tests/input/GRCh37.csv")
        chromosomes_remapped, chromosomes_not_remapped = s.remap(-1)
        self.assertEqual(s.build, 37)
        self.assertEqual(s.assembly, "GRCh37")
        self.assertEqual(len(chromosomes_remapped), 0)
        self.assertEqual(len(chromosomes_not_remapped), 2)

    def test_remap_no_snps(self):
        s = SNPs()
        chromosomes_remapped, chromosomes_not_remapped = s.remap(38)
        self.assertFalse(s.build)
        self.assertEqual(len(chromosomes_remapped), 0)
        self.assertEqual(len(chromosomes_not_remapped), 0)

    def test_save_buffer(self):
        s = SNPs("tests/input/generic.csv")
        out = io.StringIO()
        s.save(out)
        self.assertTrue(out.read().startswith("# Generated by snps"))

    def test_save_no_snps(self):
        s = SNPs()
        self.assertFalse(s.save())

    def test_save_no_snps_vcf(self):
        s = SNPs()
        self.assertFalse(s.save(vcf=True))

    def test_save_source(self):
        s = SNPs("tests/input/GRCh38.csv")
        self.assertEqual(os.path.relpath(s.save()), "output/generic_GRCh38.txt")
        snps = SNPs("output/generic_GRCh38.txt")
        self.assertEqual(snps.build, 38)
        self.assertTrue(snps.build_detected)
        self.assertEqual(snps.source, "generic")
        self.assertListEqual(snps._source, ["generic"])
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
        self.assertEqual(s.count, 15528)
        s._deduplicate_XY_chrom()
        self.assertEqual(s.count, 15528)
        self.assertEqual(len(s.discrepant_XY), 0)
        self.assertEqual(s.sex, "Male")

    def test_sex_Male_X_chrom_discrepant_XY(self):
        s = self.simulate_snps(
            chrom="X", pos_start=1, pos_max=155270560, pos_step=10000, genotype="AA"
        )
        self.assertEqual(s.count, 15528)
        s._snps.loc["rs8001", "genotype"] = "AC"
        s._deduplicate_XY_chrom()
        self.assertEqual(s.count, 15527)
        result = self.create_snp_df(
            rsid=["rs8001"], chrom=["X"], pos=[80000001], genotype=["AC"]
        )
        pd.testing.assert_frame_equal(s.discrepant_XY, result, check_exact=True)
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
        self.assertEqual(s._source, ["generic"])

    def test_source_no_snps(self):
        for snps in self.empty_snps():
            self.assertFalse(snps.source)

    def test_count(self):
        s = SNPs("tests/input/NCBI36.csv")
        self.assertEqual(s.count, 4)

    def test_count_no_snps(self):
        for snps in self.empty_snps():
            self.assertEqual(snps.count, 0)
            self.assertTrue(snps.snps.empty)


class TestSNPsMerge(TestSnps):
    def assert_results(self, results, expected_results):
        self.assertEqual(len(results), len(expected_results))

        for i, result in enumerate(results):
            expected_result = expected_results[i]

            self.assertListEqual(
                [
                    "common_rsids",
                    "discrepant_genotype_rsids",
                    "discrepant_position_rsids",
                    "merged",
                ],
                sorted(list(result.keys())),
            )

            if "merged" in expected_result:
                if expected_result["merged"]:
                    self.assertTrue(result["merged"])
                else:
                    self.assertFalse(result["merged"])
            else:
                self.assertFalse(result["merged"])

            for key in [
                "common_rsids",
                "discrepant_position_rsids",
                "discrepant_genotype_rsids",
            ]:
                if key in expected_result:
                    pd.testing.assert_index_equal(
                        result[key],
                        expected_result[key],
                        check_exact=True,
                        check_names=True,
                    )
                else:
                    pd.testing.assert_index_equal(
                        result[key],
                        pd.Index([], name="rsid"),
                        check_exact=True,
                        check_names=True,
                    )

    def test_source_snps(self):
        s = SNPs("tests/input/GRCh37.csv")
        self.assertEqual(s.source, "generic")
        results = s.merge((SNPs("tests/input/23andme.txt"),))
        self.assertEqual(s.source, "generic, 23andMe")
        self.assertListEqual(s._source, ["generic", "23andMe"])
        self.assertEqual(
            os.path.relpath(s.save()), "output/generic__23andMe_GRCh37.txt"
        )
        s = SNPs("output/generic__23andMe_GRCh37.txt")
        self.assertEqual(s.source, "generic, 23andMe")
        self.assertListEqual(s._source, ["generic", "23andMe"])
        pd.testing.assert_frame_equal(s.snps, s.snps, check_exact=True)
        self.assert_results(results, [{"merged": True}])

    def test_merge_list(self):
        s = SNPs()
        results = s.merge(
            [SNPs("tests/input/GRCh37.csv"), SNPs("tests/input/GRCh37.csv")]
        )
        pd.testing.assert_frame_equal(s.snps, self.snps_GRCh37(), check_exact=True)
        self.assertEqual(s.source, "generic, generic")
        self.assertListEqual(s._source, ["generic", "generic"])
        self.assert_results(
            results,
            [
                {"merged": True},
                {
                    "merged": True,
                    "common_rsids": pd.Index(
                        ["rs3094315", "rs2500347", "rsIndelTest", "rs11928389"],
                        name="rsid",
                    ),
                },
            ],
        )

    def test_merge_remapping(self):
        def f():
            s = SNPs("tests/input/NCBI36.csv")
            results = s.merge([SNPs("tests/input/GRCh37.csv")])
            self.assertEqual(len(s.discrepant_merge_positions), 0)
            self.assertEqual(len(s.discrepant_merge_genotypes), 0)
            pd.testing.assert_frame_equal(s.snps, self.snps_NCBI36(), check_exact=True)
            self.assert_results(
                results,
                [
                    {
                        "merged": True,
                        "common_rsids": pd.Index(
                            ["rs3094315", "rs2500347", "rsIndelTest", "rs11928389"],
                            name="rsid",
                        ),
                    }
                ],
            )

        self._run_remap_test(f, self.GRCh37_NCBI36())

    def test_merge_remap_false(self):
        s = SNPs("tests/input/NCBI36.csv")
        results = s.merge([SNPs("tests/input/GRCh37.csv")], remap=False)
        self.assertEqual(len(s.discrepant_merge_positions), 4)
        pd.testing.assert_index_equal(
            s.discrepant_merge_positions.index,
            results[0]["discrepant_position_rsids"],
            check_exact=True,
            check_names=True,
        )
        self.assertEqual(len(s.discrepant_merge_genotypes), 1)
        pd.testing.assert_index_equal(
            s.discrepant_merge_genotypes.index,
            results[0]["discrepant_genotype_rsids"],
            check_exact=True,
            check_names=True,
        )
        self.assertEqual(len(s.discrepant_merge_positions_genotypes), 4)
        pd.testing.assert_index_equal(
            s.discrepant_merge_positions_genotypes.index,
            results[0]["discrepant_position_rsids"],
            check_exact=True,
            check_names=True,
        )
        expected = self.snps_NCBI36()
        expected.loc[
            "rs11928389", "genotype"
        ] = np.nan  # discrepant genotype is set to null / NA
        pd.testing.assert_frame_equal(s.snps, expected, check_exact=True)
        self.assert_results(
            results,
            [
                {
                    "merged": True,
                    "common_rsids": pd.Index(
                        ["rs3094315", "rs2500347", "rsIndelTest", "rs11928389"],
                        name="rsid",
                    ),
                    "discrepant_position_rsids": pd.Index(
                        ["rs3094315", "rs2500347", "rsIndelTest", "rs11928389"],
                        name="rsid",
                    ),
                    "discrepant_genotype_rsids": pd.Index(["rs11928389"], name="rsid"),
                }
            ],
        )

    def test_merge_phased(self):
        s1 = SNPs("tests/input/generic.csv")
        s2 = SNPs("tests/input/generic.csv")
        s1._phased = True
        s2._phased = True

        results = s1.merge([s2])
        self.assertTrue(s1.phased)
        pd.testing.assert_frame_equal(s1.snps, self.generic_snps(), check_exact=True)
        self.assert_results(
            results,
            [
                {
                    "merged": True,
                    "common_rsids": pd.Index(
                        ["rs1", "rs2", "rs3", "rs4", "rs5", "rs6", "rs7", "rs8"],
                        name="rsid",
                    ),
                }
            ],
        )

    def test_merge_unphased(self):
        s1 = SNPs("tests/input/generic.csv")
        s2 = SNPs("tests/input/generic.csv")
        s1._phased = True

        results = s1.merge([s2])
        self.assertFalse(s1.phased)
        pd.testing.assert_frame_equal(s1.snps, self.generic_snps(), check_exact=True)
        self.assert_results(
            results,
            [
                {
                    "merged": True,
                    "common_rsids": pd.Index(
                        ["rs1", "rs2", "rs3", "rs4", "rs5", "rs6", "rs7", "rs8"],
                        name="rsid",
                    ),
                }
            ],
        )

    def test_merge_non_existent_file(self):
        s = SNPs()
        results = s.merge(
            [SNPs("tests/input/non_existent_file.csv"), SNPs("tests/input/GRCh37.csv")]
        )
        pd.testing.assert_frame_equal(s.snps, self.snps_GRCh37(), check_exact=True)
        self.assert_results(results, [{}, {"merged": True}])

    def test_merge_invalid_file(self):
        s = SNPs()
        results = s.merge(
            [SNPs("tests/input/GRCh37.csv"), SNPs("tests/input/empty.txt")]
        )
        pd.testing.assert_frame_equal(s.snps, self.snps_GRCh37(), check_exact=True)
        self.assert_results(results, [{"merged": True}, {}])

    def test_merge_exceed_discrepant_positions_threshold(self):
        s1 = SNPs("tests/input/generic.csv")
        s2 = SNPs("tests/input/generic.csv")
        s2._snps.loc["rs1", "pos"] = 100

        results = s1.merge([s2], discrepant_positions_threshold=0)
        self.assertEqual(len(s1.discrepant_merge_positions), 0)
        self.assertEqual(len(s1.discrepant_merge_genotypes), 0)
        self.assertEqual(len(s1.discrepant_merge_positions_genotypes), 0)
        pd.testing.assert_frame_equal(s1.snps, self.generic_snps(), check_exact=True)
        self.assert_results(results, [{}])

    def test_merge_exceed_discrepant_genotypes_threshold(self):
        s1 = SNPs("tests/input/generic.csv")
        s2 = SNPs("tests/input/generic.csv")
        s2._snps.loc["rs1", "genotype"] = "CC"

        results = s1.merge([s2], discrepant_genotypes_threshold=0)
        self.assertEqual(len(s1.discrepant_merge_positions), 0)
        self.assertEqual(len(s1.discrepant_merge_genotypes), 0)
        self.assertEqual(len(s1.discrepant_merge_positions_genotypes), 0)
        pd.testing.assert_frame_equal(s1.snps, self.generic_snps(), check_exact=True)
        self.assert_results(results, [{}])

    def test_merging_files_discrepant_snps(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            dest1 = os.path.join(tmpdir, "discrepant_snps1.csv")
            dest2 = os.path.join(tmpdir, "discrepant_snps2.csv")

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
                    "pos_file1": np.uint32,
                    "pos_file2": np.uint32,
                    "discrepant_position": bool,
                    "discrepant_genotype": bool,
                    "expected_position": np.uint32,
                },
            )

            df1 = df[["chrom", "pos_file1", "genotype_file1"]]
            df2 = df[["chrom", "pos_file2", "genotype_file2"]]

            df1.to_csv(
                dest1, na_rep="--", header=["chromosome", "position", "genotype"]
            )

            df2.to_csv(
                dest2, na_rep="--", header=["chromosome", "position", "genotype"]
            )

            s = SNPs()
            s.merge([SNPs(dest1), SNPs(dest2)])

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
            expected_snps.sort()
            expected = expected_snps.snps

            pd.testing.assert_index_equal(
                s.discrepant_merge_positions.index,
                expected.loc[expected["discrepant_position"] == True].index,
            )

            pd.testing.assert_index_equal(
                s.discrepant_merge_genotypes.index,
                expected.loc[expected["discrepant_genotype"] == True].index,
            )

            pd.testing.assert_series_equal(s.snps["pos"], expected["pos"])
            pd.testing.assert_series_equal(s.snps["genotype"], expected["genotype"])

    def test_appending_dfs(self):
        s = SNPs()
        s._snps = self.create_snp_df(
            rsid=["rs1"], chrom=["1"], pos=[1], genotype=["AA"]
        )
        s._duplicate = self.create_snp_df(
            rsid=["rs1"], chrom=["1"], pos=[1], genotype=["AA"]
        )
        s._discrepant_XY = self.create_snp_df(
            rsid=["rs1"], chrom=["1"], pos=[1], genotype=["AA"]
        )
        s.merge([s])
        df = self.create_snp_df(
            rsid=["rs1", "rs1"], chrom=["1", "1"], pos=[1, 1], genotype=["AA", "AA"]
        )
        pd.testing.assert_frame_equal(s.duplicate, df, check_exact=True)
        pd.testing.assert_frame_equal(s.discrepant_XY, df, check_exact=True)
        pd.testing.assert_frame_equal(
            s.heterozygous_MT, get_empty_snps_dataframe(), check_exact=True
        )
        pd.testing.assert_frame_equal(
            s.discrepant_vcf_position, get_empty_snps_dataframe(), check_exact=True
        )

    def test_merge_chrom(self):
        s1 = SNPs("tests/input/generic.csv")
        df = s1.snps.append(
            self.create_snp_df(
                rsid=["rs100", "rs101", "rs102", "rs103"],
                chrom=["Y", "Y", "Y", "Y"],
                pos=[100, 101, 102, 103],
                genotype=["A", np.nan, "A", "A"],
            )
        )
        s1._snps = df.copy()
        s2 = SNPs()
        s2._build = 37
        s2._snps = df.copy()

        # set values for chrom that will be ignored (that would otherwise result in
        # identification of discrepant SNPs or updating genotype)
        s2._snps.loc["rs3", "pos"] = 1003  # discrepant position
        s2._snps.loc["rs4", "genotype"] = "AA"  # discrepant genotype
        s2._snps.loc["rs5", "genotype"] = "AA"

        # set values for chrom to be merged
        s2._snps.loc["rs100", "genotype"] = "T"  # discrepant genotype
        s2._snps.loc["rs101", "genotype"] = "A"
        s2._snps.loc["rs102", "pos"] = 1002  # discrepant position

        # set expected values for merge result
        df.loc["rs100", "genotype"] = np.nan  # discrepant genotype sets to np.nan
        df.loc["rs101", "genotype"] = "A"  # updates np.nan

        results = s1.merge([s2], chrom="Y")

        pd.testing.assert_frame_equal(s1.snps, df, check_exact=True)

        self.assert_results(
            results,
            [
                {
                    "merged": True,
                    "common_rsids": pd.Index(
                        ["rs100", "rs101", "rs102", "rs103"], name="rsid"
                    ),
                    "discrepant_position_rsids": pd.Index(["rs102"], name="rsid"),
                    "discrepant_genotype_rsids": pd.Index(["rs100"], name="rsid"),
                }
            ],
        )

        self.assertEqual(len(s1.discrepant_merge_positions), 1)
        self.assertEqual(len(s1.discrepant_merge_genotypes), 1)


class TestDeprecatedMethods(TestSnps):
    def run_deprecated_test(self, f, msg):
        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")
            f()
            self.assertEqual(len(w), 1)
            self.assertTrue(issubclass(w[-1].category, DeprecationWarning))
            self.assertEqual(str(w[-1].message), msg)

    def test_sort_snps(self):
        def f():
            s = SNPs("tests/input/generic.csv")
            s.sort_snps()

        self.run_deprecated_test(f, "This method has been renamed to `sort`.")

    def test_remap_snps(self):
        def f1():
            def f2():
                s = SNPs("tests/input/NCBI36.csv")
                chromosomes_remapped, chromosomes_not_remapped = s.remap_snps(37)
                self.assertEqual(s.build, 37)
                self.assertEqual(s.assembly, "GRCh37")
                self.assertEqual(len(chromosomes_remapped), 2)
                self.assertEqual(len(chromosomes_not_remapped), 0)
                pd.testing.assert_frame_equal(
                    s.snps, self.snps_GRCh37(), check_exact=True
                )

            self._run_remap_test(f2, self.NCBI36_GRCh37())

        self.run_deprecated_test(f1, "This method has been renamed to `remap`.")

    def test_save_snps(self):
        def f():
            snps = SNPs("tests/input/generic.csv")
            self.assertEqual(
                os.path.relpath(snps.save_snps(sep=",")), "output/generic_GRCh37.csv"
            )
            self.run_parsing_tests("output/generic_GRCh37.csv", "generic")

        self.run_deprecated_test(f, "This method has been renamed to `save`.")

    def test_snp_count(self):
        def f():
            s = SNPs("tests/input/NCBI36.csv")
            self.assertEqual(s.snp_count, 4)

        self.run_deprecated_test(f, "This property has been renamed to `count`.")

    def test_get_snp_count(self):
        def f():
            s = SNPs("tests/input/NCBI36.csv")
            self.assertEqual(s.get_snp_count(), 4)

        self.run_deprecated_test(f, "This method has been renamed to `get_count`.")

    def test_not_null_snps(self):
        def f():
            s = SNPs("tests/input/generic.csv")
            snps = self.generic_snps()
            snps.drop("rs5", inplace=True)
            pd.testing.assert_frame_equal(s.not_null_snps(), snps, check_exact=True)

        self.run_deprecated_test(f, "This method has been renamed to `notnull`.")

    def test_get_summary(self):
        def f():
            s = SNPs("tests/input/GRCh38.csv")
            self.assertDictEqual(
                s.get_summary(),
                {
                    "source": "generic",
                    "assembly": "GRCh38",
                    "build": 38,
                    "build_detected": True,
                    "count": 4,
                    "chromosomes": "1, 3",
                    "sex": "",
                },
            )

        self.run_deprecated_test(
            f, "This method has been renamed to `summary` and is now a property."
        )

    def test_get_assembly(self):
        def f():
            s = SNPs("tests/input/GRCh38.csv")
            self.assertEqual(s.get_assembly(), "GRCh38")

        self.run_deprecated_test(f, "See the `assembly` property.")

    def test_get_chromosomes(self):
        def f():
            s = SNPs("tests/input/chromosomes.csv")
            self.assertListEqual(s.get_chromosomes(), ["1", "2", "3", "5", "PAR", "MT"])

        self.run_deprecated_test(f, "See the `chromosomes` property.")

    def test_get_chromosomes_summary(self):
        def f():
            s = SNPs("tests/input/chromosomes.csv")
            self.assertEqual(s.get_chromosomes_summary(), "1-3, 5, PAR, MT")

        self.run_deprecated_test(f, "See the `chromosomes_summary` property.")

    def test_duplicate_snps(self):
        def f():
            snps = SNPs("tests/input/duplicate_rsids.csv")
            result = self.create_snp_df(
                rsid=["rs1"], chrom=["1"], pos=[101], genotype=["AA"]
            )
            duplicate = self.create_snp_df(
                rsid=["rs1", "rs1"],
                chrom=["1", "1"],
                pos=[102, 103],
                genotype=["CC", "GG"],
            )
            pd.testing.assert_frame_equal(snps.snps, result, check_exact=True)
            pd.testing.assert_frame_equal(
                snps.duplicate_snps, duplicate, check_exact=True
            )

        self.run_deprecated_test(f, "This property has been renamed to `duplicate`.")

    def test_discrepant_XY_snps(self):
        def f():
            s = self.simulate_snps(
                chrom="X", pos_start=1, pos_max=155270560, pos_step=10000, genotype="AA"
            )
            self.assertEqual(s.count, 15528)
            s._snps.loc["rs8001", "genotype"] = "AC"
            s._deduplicate_XY_chrom()
            self.assertEqual(s.count, 15527)
            result = self.create_snp_df(
                rsid=["rs8001"], chrom=["X"], pos=[80000001], genotype=["AC"]
            )
            pd.testing.assert_frame_equal(
                s.discrepant_XY_snps, result, check_exact=True
            )
            self.assertEqual(s.sex, "Male")

        self.run_deprecated_test(
            f, "This property has been renamed to `discrepant_XY`."
        )

    def test_heterozygous_MT_snps(self):
        def f():
            snps = SNPs("tests/input/ancestry_mt.txt")

            result = self.create_snp_df(
                rsid=["rs1", "rs2"],
                chrom=["MT", "MT"],
                pos=[101, 102],
                genotype=["A", np.nan],
            )
            pd.testing.assert_frame_equal(snps.snps, result, check_exact=True)

            heterozygous_MT_snps = self.create_snp_df(
                rsid=["rs3"], chrom=["MT"], pos=[103], genotype=["GC"]
            )
            pd.testing.assert_frame_equal(
                snps.heterozygous_MT_snps, heterozygous_MT_snps, check_exact=True
            )

        self.run_deprecated_test(
            f, "This property has been renamed to `heterozygous_MT`."
        )

    def test_heterozygous_snps(self):
        def f():
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

        self.run_deprecated_test(f, "This method has been renamed to `heterozygous`.")

    def test_homozygous_snps(self):
        def f():
            s = SNPs("tests/input/generic.csv")
            pd.testing.assert_frame_equal(
                s.homozygous_snps(),
                self.create_snp_df(
                    rsid=["rs1", "rs2", "rs3", "rs4"],
                    chrom=["1", "1", "1", "1"],
                    pos=[101, 102, 103, 104],
                    genotype=["AA", "CC", "GG", "TT"],
                ),
                check_exact=True,
            )

        self.run_deprecated_test(f, "This method has been renamed to `homozygous`.")

    def test_discrepant_positions(self):
        def f():
            s = SNPs("tests/input/NCBI36.csv")
            results = s.merge([SNPs("tests/input/GRCh37.csv")], remap=False)
            df = s.discrepant_positions
            self.assertEqual(len(df), 4)
            pd.testing.assert_index_equal(
                df.index,
                results[0]["discrepant_position_rsids"],
                check_exact=True,
                check_names=True,
            )

        self.run_deprecated_test(
            f, "This property has been renamed to `discrepant_merge_positions`."
        )

    def test_discrepant_genotypes(self):
        def f():
            s = SNPs("tests/input/NCBI36.csv")
            results = s.merge([SNPs("tests/input/GRCh37.csv")], remap=False)
            df = s.discrepant_genotypes
            self.assertEqual(len(df), 1)
            pd.testing.assert_index_equal(
                df.index,
                results[0]["discrepant_genotype_rsids"],
                check_exact=True,
                check_names=True,
            )

        self.run_deprecated_test(
            f, "This property has been renamed to `discrepant_merge_genotypes`."
        )

    def test_discrepant_snps(self):
        def f():
            s = SNPs("tests/input/NCBI36.csv")
            results = s.merge([SNPs("tests/input/GRCh37.csv")], remap=False)
            df = s.discrepant_snps
            self.assertEqual(len(df), 4)
            pd.testing.assert_index_equal(
                df.index,
                results[0]["discrepant_position_rsids"],
                check_exact=True,
                check_names=True,
            )

        self.run_deprecated_test(
            f,
            "This property has been renamed to `discrepant_merge_positions_genotypes`.",
        )

    def test_is_valid(self):
        def f():
            s = SNPs("tests/input/generic.csv")
            self.assertTrue(s.is_valid())

        self.run_deprecated_test(
            f, "This method has been renamed to `valid` and is now a property."
        )
