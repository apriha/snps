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
from atomicwrites import atomic_write
import gzip
import pandas as pd
import shutil

from snps import SNPs
from snps.resources import Resources, ReferenceSequence
from tests import BaseSNPsTestCase


class TestWriter(BaseSNPsTestCase):
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
