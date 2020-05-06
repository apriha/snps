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
import tempfile

from snps import SNPs
from snps.resources import Resources, ReferenceSequence
from snps.utils import gzip_file
from tests import BaseSNPsTestCase


class TestWriter(BaseSNPsTestCase):
    def test_save_snps(self):
        snps = SNPs("tests/input/generic.csv")
        self.assertEqual(os.path.relpath(snps.save_snps()), "output/generic_GRCh37.csv")
        self.run_parsing_tests("output/generic_GRCh37.csv", "generic")

    def test_save_snps_tsv(self):
        snps = SNPs("tests/input/generic.csv")
        self.assertEqual(
            os.path.relpath(snps.save_snps("generic.tsv", sep="\t")),
            "output/generic.tsv",
        )
        self.run_parsing_tests("output/generic.tsv", "generic")

    def test_save_snps_vcf(self):
        s = SNPs("tests/input/testvcf.vcf")

        r = Resources()
        r._reference_sequences["GRCh37"] = {}

        with tempfile.TemporaryDirectory() as tmpdir:
            dest = os.path.join(tmpdir, "generic.fa.gz")
            gzip_file("tests/input/generic.fa", dest)

            seq = ReferenceSequence(ID="1", path=dest)

            r._reference_sequences["GRCh37"]["1"] = seq

            self.assertEqual(
                os.path.relpath(s.save_snps(vcf=True)), "output/vcf_GRCh37.vcf"
            )

        self.run_parsing_tests_vcf("output/vcf_GRCh37.vcf")

    def test_save_snps_vcf_phased(self):
        # read phased data
        s = SNPs("tests/input/testvcf_phased.vcf")

        # setup resource to use test FASTA reference sequence
        r = Resources()
        r._reference_sequences["GRCh37"] = {}

        with tempfile.TemporaryDirectory() as tmpdir:
            dest = os.path.join(tmpdir, "generic.fa.gz")
            gzip_file("tests/input/generic.fa", dest)

            seq = ReferenceSequence(ID="1", path=dest)

            r._reference_sequences["GRCh37"]["1"] = seq

            # save phased data to VCF
            self.assertEqual(
                os.path.relpath(s.save_snps(vcf=True)), "output/vcf_GRCh37.vcf"
            )

        # read saved VCF
        self.run_parsing_tests_vcf("output/vcf_GRCh37.vcf", phased=True)

    def test_save_snps_csv_phased(self):
        # read phased data
        s = SNPs("tests/input/testvcf_phased.vcf")
        # save phased data to CSV
        self.assertEqual(os.path.relpath(s.save_snps()), "output/vcf_GRCh37.csv")
        # read saved CSV
        self.run_parsing_tests_vcf("output/vcf_GRCh37.csv", phased=True)

    def test_save_snps_specify_file(self):
        s = SNPs("tests/input/generic.csv")
        self.assertEqual(os.path.relpath(s.save_snps("snps.csv")), "output/snps.csv")
        self.run_parsing_tests("output/snps.csv", "generic")
