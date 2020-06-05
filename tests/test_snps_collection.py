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
import tempfile

from atomicwrites import atomic_write
import numpy as np
import pandas as pd

from snps import SNPs, SNPsCollection

from tests import BaseSNPsTestCase


class TestSNPsCollection(BaseSNPsTestCase):
    def snps_NCBI36_discrepant_snps(self):
        return self.create_snp_df(
            rsid=["rs3094315", "rs2500347", "rsIndelTest", "rs11928389"],
            chrom=["1", "1", "1", "3"],
            pos=[742429, 143649677, 143649678, 50908372],
            genotype=["AA", np.nan, "ID", np.nan],
        )

    def test_source_lineage_file(self):
        sc = SNPsCollection("tests/input/GRCh37.csv")
        self.assertEqual(sc.source, "generic")
        sc.load_snps("tests/input/23andme.txt")
        self.assertEqual(sc.source, "generic, 23andMe")
        file = sc.save_snps()
        s = SNPs(file)
        self.assertEqual(s.source, "generic, 23andMe")
        pd.testing.assert_frame_equal(sc.snps, s.snps, check_exact=True)

    def test_source_lineage_file_gzip(self):
        sc = SNPsCollection("tests/input/GRCh37.csv")
        self.assertEqual(sc.source, "generic")
        sc.load_snps("tests/input/23andme.txt")
        self.assertEqual(sc.source, "generic, 23andMe")
        file = sc.save_snps()
        with open(file, "rb") as f_in:
            with atomic_write(file + ".gz", mode="wb", overwrite=True) as f_out:
                with gzip.open(f_out, "wb") as f_gzip:
                    shutil.copyfileobj(f_in, f_gzip)
        s = SNPs(file + ".gz")
        self.assertEqual(s.source, "generic, 23andMe")
        pd.testing.assert_frame_equal(sc.snps, s.snps, check_exact=True)

    def test_load_snps_list(self):
        sc = SNPsCollection()
        sc.load_snps(["tests/input/GRCh37.csv", "tests/input/GRCh37.csv"])
        pd.testing.assert_frame_equal(sc.snps, self.snps_GRCh37(), check_exact=True)
        self.assertEqual(sc.source, "generic, generic")

    def test_load_snps_None(self):
        sc = SNPsCollection()
        with self.assertRaises(TypeError):
            sc.load_snps(None)

    def test_discrepant_positions(self):
        sc = SNPsCollection()
        sc.load_snps(["tests/input/NCBI36.csv", "tests/input/GRCh37.csv"])
        self.assertEqual(len(sc.discrepant_positions), 4)

    def test_discrepant_genotypes(self):
        sc = SNPsCollection()
        sc.load_snps(["tests/input/NCBI36.csv", "tests/input/GRCh37.csv"])
        self.assertEqual(len(sc.discrepant_genotypes), 1)

    def test_discrepant_snps(self):
        sc = SNPsCollection()
        sc.load_snps(["tests/input/NCBI36.csv", "tests/input/GRCh37.csv"])
        self.assertEqual(len(sc.discrepant_snps), 4)

    def test_load_snps_non_existent_file(self):
        sc = SNPsCollection()
        sc.load_snps(["tests/input/GRCh37.csv", "tests/input/non_existent_file.csv"])
        pd.testing.assert_frame_equal(sc.snps, self.snps_GRCh37(), check_exact=True)

    def test_load_snps_invalid_file(self):
        sc = SNPsCollection()
        sc.load_snps(["tests/input/GRCh37.csv", "tests/input/empty.txt"])
        pd.testing.assert_frame_equal(sc.snps, self.snps_GRCh37(), check_exact=True)

    def test_load_snps_assembly_mismatch(self):
        sc = SNPsCollection()
        sc.load_snps(["tests/input/NCBI36.csv", "tests/input/GRCh37.csv"])
        self.assertFalse(os.path.exists("output/ind_discrepant_positions_1.csv"))
        self.assertFalse(os.path.exists("output/ind_discrepant_genotypes_1.csv"))
        self.assertEqual(len(sc.discrepant_positions), 4)
        self.assertEqual(len(sc.discrepant_genotypes), 1)
        pd.testing.assert_frame_equal(
            sc.snps, self.snps_NCBI36_discrepant_snps(), check_exact=True
        )

    def test_load_snps_assembly_mismatch_save_output(self):
        sc = SNPsCollection()
        sc.load_snps(
            ["tests/input/NCBI36.csv", "tests/input/GRCh37.csv"], save_output=True
        )
        self.assertTrue(os.path.exists("output/discrepant_positions_1.csv"))
        self.assertTrue(os.path.exists("output/discrepant_genotypes_1.csv"))
        self.assertEqual(len(sc.discrepant_positions), 4)
        self.assertEqual(len(sc.discrepant_genotypes), 1)
        pd.testing.assert_frame_equal(
            sc.snps, self.snps_NCBI36_discrepant_snps(), check_exact=True
        )

    def test_load_snps_assembly_mismatch_exceed_discrepant_positions_threshold(self):
        sc = SNPsCollection()
        sc.load_snps(
            ["tests/input/NCBI36.csv", "tests/input/GRCh37.csv"],
            discrepant_snp_positions_threshold=0,
        )
        self.assertFalse(os.path.exists("output/discrepant_positions_1.csv"))
        self.assertFalse(os.path.exists("output/discrepant_genotypes_1.csv"))
        self.assertEqual(len(sc.discrepant_positions), 4)
        self.assertEqual(len(sc.discrepant_genotypes), 0)
        pd.testing.assert_frame_equal(sc.snps, self.snps_NCBI36(), check_exact=True)

    def test_load_snps_assembly_mismatch_exceed_discrepant_genotypes_threshold(self):
        sc = SNPsCollection()
        sc.load_snps(
            ["tests/input/NCBI36.csv", "tests/input/GRCh37.csv"],
            discrepant_genotypes_threshold=0,
        )
        self.assertFalse(os.path.exists("output/discrepant_positions_1.csv"))
        self.assertFalse(os.path.exists("output/discrepant_genotypes_1.csv"))
        self.assertEqual(len(sc.discrepant_positions), 4)
        self.assertEqual(len(sc.discrepant_genotypes), 1)
        pd.testing.assert_frame_equal(sc.snps, self.snps_NCBI36(), check_exact=True)

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
                    "pos_file1": np.int64,
                    "pos_file2": np.int64,
                    "discrepant_position": bool,
                    "discrepant_genotype": bool,
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

            sc = SNPsCollection([dest1, dest2])

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

    def test_save_discrepant_positions(self):
        sc = SNPsCollection()
        sc.load_snps(["tests/input/NCBI36.csv", "tests/input/GRCh37.csv"])
        self.assertEqual(len(sc.discrepant_positions), 4)
        discrepant_positions_file = sc.save_discrepant_positions()
        self.assertEqual(
            os.path.relpath(discrepant_positions_file),
            "output/discrepant_positions.csv",
        )
        self.assertTrue(os.path.exists(discrepant_positions_file))

    def test_save_discrepant_positions_specify_file(self):
        sc = SNPsCollection()
        sc.load_snps(["tests/input/NCBI36.csv", "tests/input/GRCh37.csv"])
        self.assertEqual(len(sc.discrepant_positions), 4)
        discrepant_positions_file = sc.save_discrepant_positions(
            "discrepant_positions.csv"
        )
        self.assertEqual(
            os.path.relpath(discrepant_positions_file),
            "output/discrepant_positions.csv",
        )
        self.assertTrue(os.path.exists(discrepant_positions_file))

    def test_save_discrepant_positions_no_discrepant_snps(self):
        sc = SNPsCollection()
        self.assertEqual(len(sc.discrepant_positions), 0)
        self.assertFalse(sc.save_discrepant_positions())

    def test_save_discrepant_positions_exception(self):
        sc = SNPsCollection()
        sc._discrepant_positions = "invalid"
        self.assertFalse(sc.save_discrepant_positions())

    def test_save_discrepant_genotypes(self):
        sc = SNPsCollection()
        sc.load_snps(["tests/input/NCBI36.csv", "tests/input/GRCh37.csv"])
        self.assertEqual(len(sc.discrepant_genotypes), 1)
        discrepant_genotypes_file = sc.save_discrepant_genotypes()
        self.assertEqual(
            os.path.relpath(discrepant_genotypes_file),
            "output/discrepant_genotypes.csv",
        )
        self.assertTrue(os.path.exists(discrepant_genotypes_file))

    def test_save_discrepant_genotypes_specify_file(self):
        sc = SNPsCollection()
        sc.load_snps(["tests/input/NCBI36.csv", "tests/input/GRCh37.csv"])
        self.assertEqual(len(sc.discrepant_genotypes), 1)
        discrepant_genotypes_file = sc.save_discrepant_genotypes(
            "discrepant_genotypes.csv"
        )
        self.assertEqual(
            os.path.relpath(discrepant_genotypes_file),
            "output/discrepant_genotypes.csv",
        )
        self.assertTrue(os.path.exists(discrepant_genotypes_file))

    def test_save_discrepant_genotypes_no_discrepant_snps(self):
        sc = SNPsCollection()
        self.assertEqual(len(sc.discrepant_genotypes), 0)
        self.assertFalse(sc.save_discrepant_genotypes())

    def test_save_discrepant_genotypes_exception(self):
        sc = SNPsCollection()
        sc._discrepant_genotypes = "invalid"
        self.assertFalse(sc.save_discrepant_genotypes())

    def test_save_discrepant_snps(self):
        sc = SNPsCollection()
        sc.load_snps(["tests/input/NCBI36.csv", "tests/input/GRCh37.csv"])
        self.assertEqual(len(sc.discrepant_snps), 4)
        discrepant_snps_file = sc.save_discrepant_snps()
        self.assertEqual(
            os.path.relpath(discrepant_snps_file), "output/discrepant_snps.csv"
        )
        self.assertTrue(os.path.exists(discrepant_snps_file))

    def test_save_discrepant_snps_specify_file(self):
        sc = SNPsCollection()
        sc.load_snps(["tests/input/NCBI36.csv", "tests/input/GRCh37.csv"])
        self.assertEqual(len(sc.discrepant_snps), 4)
        discrepant_snps_file = sc.save_discrepant_snps("discrepant_snps.csv")
        self.assertEqual(
            os.path.relpath(discrepant_snps_file), "output/discrepant_snps.csv"
        )
        self.assertTrue(os.path.exists(discrepant_snps_file))

    def test_save_discrepant_snps_no_discrepant_snps(self):
        sc = SNPsCollection()
        self.assertEqual(len(sc.discrepant_snps), 0)
        self.assertFalse(sc.save_discrepant_snps())

    def test_save_discrepant_snps_exception(self):
        sc = SNPsCollection()
        sc._discrepant_snps = "invalid"
        self.assertFalse(sc.save_discrepant_snps())

    def test___repr__snps_collection(self):
        sc = SNPsCollection()
        self.assertEqual("SNPsCollection(name='')", sc.__repr__())
