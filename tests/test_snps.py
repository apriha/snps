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
from tests import BaseSNPsTestCase


class TestSnps(BaseSNPsTestCase):
    def setUp(self):
        self.snps_GRCh38 = SNPs("tests/input/GRCh38.csv")
        self.snps = SNPs("tests/input/chromosomes.csv")
        self.snps_none = SNPs(None)
        self.del_output_dir_helper()

    def snps_discrepant_pos(self):
        return self.create_snp_df(
            rsid=["rs3094315"], chrom=["1"], pos=[1], genotype=["AA"]
        )

    def test_assembly(self):
        assert self.snps_GRCh38.assembly == "GRCh38"

    def test_assembly_no_snps(self):
        assert self.snps_none.assembly == ""

    def test_snp_count(self):
        assert self.snps.snp_count == 6

    def test_snp_count_no_snps(self):
        assert self.snps_none.snp_count == 0

    def test_chromosomes(self):
        assert self.snps.chromosomes == ["1", "2", "3", "5", "PAR", "MT"]

    def test_chromosomes_no_snps(self):
        assert self.snps_none.chromosomes == []

    def test_chromosomes_summary(self):
        assert self.snps.chromosomes_summary == "1-3, 5, PAR, MT"

    def test_chromosomes_summary_no_snps(self):
        assert self.snps_none.chromosomes_summary == ""

    def test_build_no_snps(self):
        assert not self.snps_none.build

    def test_build_detected_no_snps(self):
        assert not self.snps_none.build_detected

    def test_build_detected_PAR_snps(self):
        if os.getenv("DOWNLOADS_ENABLED"):
            snps = SNPs("tests/input/GRCh37_PAR.csv")
            assert snps.build == 37
            assert snps.build_detected

    def test_sex_no_snps(self):
        assert self.snps_none.sex == ""

    def test_sex_Male_Y_chrom(self):
        s = self.simulate_snps(chrom="Y", pos_start=1, pos_max=59373566, pos_step=10000)
        file = s.save_snps()
        snps = SNPs(file)
        assert snps.sex == "Male"

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
