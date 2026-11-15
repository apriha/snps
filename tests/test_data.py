import pandas as pd

from snps.data import SNPData
from tests import BaseSNPsTestCase


class TestSNPData(BaseSNPsTestCase):
    def test_frame(self):
        for data in (SNPData(), SNPData(self.mixed_phase_snps())):
            self.assertListEqual(list(data.frame.columns), ["chrom", "pos", "genotype"])
            pd.testing.assert_index_equal(data.phased.index, data.frame.index)

    def test_phase_default_unphased(self):
        df = self.generic_snps()
        data = SNPData(df)
        self.assertListEqual(data.phased.tolist(), [False] * 8)
        self.assertEqual(data.phased.dtype, bool)
        # the frame passed in is not modified
        self.assertListEqual(list(df.columns), ["chrom", "pos", "genotype"])

    def test_sort_phase(self):
        df = self.create_snp_df(
            rsid=["rs3", "rs1", "rs2", "rs1"],
            chrom=["2", "1", "MT", "1"],
            pos=[1, 2, 1, 1],
            genotype=["AC", "GT", "A", "TA"],
        ).assign(phased=[True, False, True, True])
        data = SNPData(df).sort()
        self.assert_frame_equal_with_string_index(
            data.frame, df.iloc[[3, 1, 0, 2]].drop(columns="phased"), check_exact=True
        )
        self.assertListEqual(data.phased.tolist(), [True, False, True, True])
