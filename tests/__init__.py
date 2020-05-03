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
import shutil
from unittest import TestCase

import numpy as np
import pandas as pd

from snps import SNPs


class BaseSNPsTestCase(TestCase):
    def setUp(self):
        self.del_output_dir_helper()

    def tearDown(self):
        self.del_output_dir_helper()

    @staticmethod
    def del_output_dir_helper():
        if os.path.exists("output"):
            shutil.rmtree("output")

    def simulate_snps(
        self,
        chrom="1",
        pos_start=1,
        pos_max=248140902,
        pos_step=100,
        genotype="AA",
        insert_nulls=True,
        null_snp_step=101,
        complement_genotype_one_chrom=False,
        complement_genotype_two_chroms=False,
        complement_snp_step=50,
    ):
        s = SNPs()

        s._build = 37

        positions = np.arange(pos_start, pos_max, pos_step, dtype=np.int64)
        snps = pd.DataFrame(
            {"chrom": chrom},
            index=pd.Index(
                ["rs" + str(x + 1) for x in range(len(positions))], name="rsid"
            ),
        )
        snps["pos"] = positions
        snps["genotype"] = genotype

        if insert_nulls:
            snps.loc[snps.iloc[0::null_snp_step, :].index, "genotype"] = np.nan

        indices = snps.iloc[0::complement_snp_step, :].index
        if complement_genotype_two_chroms:
            snps.loc[indices, "genotype"] = snps.loc[indices, "genotype"].apply(
                self.complement_two_chroms
            )
        elif complement_genotype_one_chrom:
            snps.loc[indices, "genotype"] = snps.loc[indices, "genotype"].apply(
                self.complement_one_chrom
            )

        s._snps = snps

        return s

    @staticmethod
    def get_complement(base):
        if base == "A":
            return "T"
        elif base == "G":
            return "C"
        elif base == "C":
            return "G"
        elif base == "T":
            return "A"
        else:
            return base

    def complement_one_chrom(self, genotype):
        if pd.isnull(genotype):
            return np.nan

        complement = ""

        for base in list(genotype):
            complement += self.get_complement(base)
            complement += genotype[1]
            return complement

    def complement_two_chroms(self, genotype):
        if pd.isnull(genotype):
            return np.nan

        complement = ""

        for base in list(genotype):
            complement += self.get_complement(base)

        return complement

    @staticmethod
    def create_snp_df(rsid, chrom, pos, genotype):
        df = pd.DataFrame(
            {"rsid": rsid, "chrom": chrom, "pos": pos, "genotype": genotype},
            columns=["rsid", "chrom", "pos", "genotype"],
        )
        df = df.set_index("rsid")
        return df

    def snps_NCBI36(self):
        return self.create_snp_df(
            rsid=["rs3094315", "rs2500347", "rsIndelTest", "rs11928389"],
            chrom=["1", "1", "1", "3"],
            pos=[742429, 143649677, 143649678, 50908372],
            genotype=["AA", np.nan, "ID", "AG"],
        )

    def snps_GRCh37(self):
        return self.create_snp_df(
            rsid=["rs3094315", "rs2500347", "rsIndelTest", "rs11928389"],
            chrom=["1", "1", "1", "3"],
            pos=[752566, 144938320, 144938321, 50927009],
            genotype=["AA", np.nan, "ID", "TC"],
        )

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

    def run_parsing_tests(self, file, source):
        self.make_parsing_assertions(self.parse_file(file), source)
        self.make_parsing_assertions(self.parse_bytes(file), source)

    def run_parsing_tests_vcf(
        self, file, source="vcf", phased=False, unannotated=False, rsids=()
    ):
        # https://samtools.github.io/hts-specs/VCFv4.2.pdf
        # this tests for homozygous snps, heterozygous snps, multiallelic snps,
        # phased snps, and snps with missing rsID
        self.make_parsing_assertions_vcf(
            self.parse_file(file, rsids), source, phased, unannotated, rsids
        )
        self.make_parsing_assertions_vcf(
            self.parse_bytes(file, rsids), source, phased, unannotated, rsids
        )

    def parse_file(self, file, rsids=()):
        return SNPs(file, rsids=rsids)

    def parse_bytes(self, file, rsids=()):
        with open(file, "rb") as f:
            return SNPs(f.read(), rsids=rsids)

    def make_parsing_assertions(self, snps, source):
        self.assertEqual(snps.source, source)
        pd.testing.assert_frame_equal(snps.snps, self.generic_snps())

    def make_parsing_assertions_vcf(self, snps, source, phased, unannotated, rsids):
        self.assertEqual(snps.source, source)

        if unannotated:
            self.assertTrue(snps.unannotated_vcf)
            self.assertEqual(0, snps.snp_count)
        else:
            self.assertFalse(snps.unannotated_vcf)
            pd.testing.assert_frame_equal(
                snps.snps, self.generic_snps_vcf().loc[rsids]
            ) if rsids else pd.testing.assert_frame_equal(
                snps.snps, self.generic_snps_vcf()
            )

        self.assertTrue(snps.phased) if phased else self.assertFalse(snps.phased)
