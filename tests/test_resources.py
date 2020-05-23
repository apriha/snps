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
import socket
import tempfile
from unittest.mock import Mock, patch
import warnings
import zipfile

from atomicwrites import atomic_write
import numpy as np
import pandas as pd

from snps import SNPs
from snps.resources import Resources, ReferenceSequence
from snps.utils import gzip_file
from tests import BaseSNPsTestCase


class TestResources(BaseSNPsTestCase):
    def setUp(self):
        self.resource = Resources(resources_dir="resources")
        self.del_output_dir_helper()

    def test_get_assembly_mapping_data(self):
        assembly_mapping_data = self.resource.get_assembly_mapping_data(
            "NCBI36", "GRCh37"
        )
        self.assertEqual(len(assembly_mapping_data), 25)

    def test_get_gsa_resources(self):
        gsa_resources = self.resource.get_gsa_resources()
        self.assertEqual(len(gsa_resources["rsid_map"]), 618541)
        self.assertEqual(len(gsa_resources["chrpos_map"]), 665609)

    def test_get_all_resources(self):
        resources = self.resource.get_all_resources()

        for k, v in resources.items():
            self.assertGreater(len(v), 0)

    def test_download_example_datasets(self):
        paths = self.resource.download_example_datasets()

        for path in paths:
            if not path or not os.path.exists(path):
                warnings.warn("Example dataset(s) not currently available")
                return

    def test_get_paths_reference_sequences_invalid_assembly(self):
        assembly, chroms, urls, paths = self.resource._get_paths_reference_sequences(
            assembly="36"
        )
        self.assertFalse(assembly)
        self.assertFalse(chroms)
        self.assertFalse(urls)
        self.assertFalse(paths)

    def test_create_reference_sequences_NCBI36(self):
        assembly, chroms, urls, paths = self.resource._get_paths_reference_sequences(
            assembly="NCBI36", chroms=["MT"]
        )
        seqs = self.resource._create_reference_sequences(assembly, chroms, urls, paths)
        self.assertEqual(len(seqs), 1)
        self.assertEqual(
            seqs["MT"].__repr__(), "ReferenceSequence(assembly='NCBI36', ID='MT')"
        )
        self.assertEqual(seqs["MT"].ID, "MT")
        self.assertEqual(seqs["MT"].chrom, "MT")
        self.assertEqual(
            seqs["MT"].url,
            "ftp://ftp.ensembl.org/pub/release-54/fasta/homo_sapiens/dna/Homo_sapiens.NCBI36.54.dna.chromosome.MT.fa.gz",
        )
        self.assertEqual(
            seqs["MT"].path,
            "resources/fasta/NCBI36/Homo_sapiens.NCBI36.54.dna.chromosome.MT.fa.gz",
        )
        self.assertTrue(os.path.exists(seqs["MT"].path))
        self.assertEqual(seqs["MT"].assembly, "NCBI36")
        self.assertEqual(seqs["MT"].build, "B36")
        self.assertEqual(seqs["MT"].species, "Homo sapiens")
        self.assertEqual(seqs["MT"].taxonomy, "x")

    def test_create_reference_sequences_GRCh37(self):
        assembly, chroms, urls, paths = self.resource._get_paths_reference_sequences(
            assembly="GRCh37", chroms=["MT"]
        )
        seqs = self.resource._create_reference_sequences(assembly, chroms, urls, paths)
        self.assertEqual(len(seqs), 1)
        self.assertEqual(
            seqs["MT"].__repr__(), "ReferenceSequence(assembly='GRCh37', ID='MT')"
        )
        self.assertEqual(seqs["MT"].ID, "MT")
        self.assertEqual(seqs["MT"].chrom, "MT")
        self.assertEqual(
            seqs["MT"].url,
            "ftp://ftp.ensembl.org/pub/grch37/release-96/fasta/homo_sapiens/dna/Homo_sapiens.GRCh37.dna.chromosome.MT.fa.gz",
        )
        self.assertEqual(
            seqs["MT"].path,
            "resources/fasta/GRCh37/Homo_sapiens.GRCh37.dna.chromosome.MT.fa.gz",
        )
        self.assertTrue(os.path.exists(seqs["MT"].path))
        self.assertEqual(seqs["MT"].assembly, "GRCh37")
        self.assertEqual(seqs["MT"].build, "B37")
        self.assertEqual(seqs["MT"].species, "Homo sapiens")
        self.assertEqual(seqs["MT"].taxonomy, "x")

    def test_create_reference_sequences_GRCh38(self):
        assembly, chroms, urls, paths = self.resource._get_paths_reference_sequences(
            assembly="GRCh38", chroms=["MT"]
        )
        seqs = self.resource._create_reference_sequences(assembly, chroms, urls, paths)
        self.assertEqual(len(seqs), 1)
        self.assertEqual(
            seqs["MT"].__repr__(), "ReferenceSequence(assembly='GRCh38', ID='MT')"
        )
        self.assertEqual(seqs["MT"].ID, "MT")
        self.assertEqual(seqs["MT"].chrom, "MT")
        self.assertEqual(
            seqs["MT"].url,
            "ftp://ftp.ensembl.org/pub/release-96/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.chromosome.MT.fa.gz",
        )
        self.assertEqual(
            seqs["MT"].path,
            "resources/fasta/GRCh38/Homo_sapiens.GRCh38.dna.chromosome.MT.fa.gz",
        )
        self.assertTrue(os.path.exists(seqs["MT"].path))
        self.assertEqual(seqs["MT"].assembly, "GRCh38")
        self.assertEqual(seqs["MT"].build, "B38")
        self.assertEqual(seqs["MT"].species, "Homo sapiens")
        self.assertEqual(seqs["MT"].taxonomy, "x")

    def test_create_reference_sequences_invalid_path(self):
        assembly, chroms, urls, paths = self.resource._get_paths_reference_sequences(
            assembly="GRCh38", chroms=["MT"]
        )
        paths[0] = ""
        seqs = self.resource._create_reference_sequences(assembly, chroms, urls, paths)
        self.assertEqual(len(seqs), 0)

    def test_download_file_socket_timeout(self):
        mock = Mock(side_effect=socket.timeout)
        with patch("urllib.request.urlopen", mock):
            path = self.resource._download_file("http://url", "test.txt")
        self.assertEqual(path, "")

    def test_get_reference_sequences(self):
        seqs = self.resource.get_reference_sequences(chroms=["MT"])
        self.assertEqual(len(seqs), 1)
        self.assertEqual(
            seqs["MT"].__repr__(), "ReferenceSequence(assembly='GRCh37', ID='MT')"
        )
        self.assertEqual(seqs["MT"].ID, "MT")
        self.assertEqual(seqs["MT"].chrom, "MT")
        self.assertEqual(
            seqs["MT"].url,
            "ftp://ftp.ensembl.org/pub/grch37/release-96/fasta/homo_sapiens/dna/Homo_sapiens.GRCh37.dna.chromosome.MT.fa.gz",
        )
        self.assertEqual(
            seqs["MT"].path,
            "resources/fasta/GRCh37/Homo_sapiens.GRCh37.dna.chromosome.MT.fa.gz",
        )
        self.assertTrue(os.path.exists(seqs["MT"].path))
        self.assertEqual(seqs["MT"].assembly, "GRCh37")
        self.assertEqual(seqs["MT"].build, "B37")
        self.assertEqual(seqs["MT"].species, "Homo sapiens")
        self.assertEqual(seqs["MT"].taxonomy, "x")

    def test_get_all_reference_sequences(self):
        seqs = self.resource.get_all_reference_sequences(chroms=["MT"])
        self.assertEqual(len(seqs), 3)
        self.assertEqual(len(seqs["NCBI36"]), 1)
        self.assertEqual(
            seqs["NCBI36"]["MT"].path,
            "resources/fasta/NCBI36/Homo_sapiens.NCBI36.54.dna.chromosome.MT.fa.gz",
        )
        self.assertEqual(len(seqs["GRCh37"]), 1)
        self.assertEqual(
            seqs["GRCh37"]["MT"].path,
            "resources/fasta/GRCh37/Homo_sapiens.GRCh37.dna.chromosome.MT.fa.gz",
        )
        self.assertEqual(len(seqs["GRCh38"]), 1)
        self.assertEqual(
            seqs["GRCh38"]["MT"].path,
            "resources/fasta/GRCh38/Homo_sapiens.GRCh38.dna.chromosome.MT.fa.gz",
        )

    def test_get_reference_sequences_invalid_assembly(self):
        seqs = self.resource.get_reference_sequences(assembly="36")
        self.assertEqual(len(seqs), 0)

    def test_get_reference_sequences_chrom_not_available(self):
        self.resource.get_reference_sequences(chroms=["MT"])
        del self.resource._reference_sequences["GRCh37"]["MT"]
        seqs = self.resource.get_reference_sequences(chroms=["MT"])
        self.assertEqual(len(seqs), 1)
        self.assertEqual(
            seqs["MT"].__repr__(), "ReferenceSequence(assembly='GRCh37', ID='MT')"
        )
        self.assertEqual(seqs["MT"].ID, "MT")
        self.assertEqual(seqs["MT"].chrom, "MT")
        self.assertEqual(
            seqs["MT"].url,
            "ftp://ftp.ensembl.org/pub/grch37/release-96/fasta/homo_sapiens/dna/Homo_sapiens.GRCh37.dna.chromosome.MT.fa.gz",
        )
        self.assertEqual(
            seqs["MT"].path,
            "resources/fasta/GRCh37/Homo_sapiens.GRCh37.dna.chromosome.MT.fa.gz",
        )
        self.assertTrue(os.path.exists(seqs["MT"].path))
        self.assertEqual(seqs["MT"].assembly, "GRCh37")
        self.assertEqual(seqs["MT"].build, "B37")
        self.assertEqual(seqs["MT"].species, "Homo sapiens")
        self.assertEqual(seqs["MT"].taxonomy, "x")

    def test_reference_sequence_load_sequence(self):
        seqs = self.resource.get_reference_sequences(chroms=["MT"])
        self.assertEqual(len(seqs["MT"].sequence), 16569)
        self.assertEqual(seqs["MT"].md5, "c68f52674c9fb33aef52dcf399755519")
        self.assertEqual(seqs["MT"].start, 1)
        self.assertEqual(seqs["MT"].end, 16569)
        self.assertEqual(seqs["MT"].length, 16569)

        seqs["MT"].clear()
        self.assertEqual(seqs["MT"]._sequence.size, 0)
        self.assertEqual(seqs["MT"]._md5, "")
        self.assertEqual(seqs["MT"]._start, 0)
        self.assertEqual(seqs["MT"]._end, 0)
        self.assertEqual(seqs["MT"]._length, 0)

        self.assertEqual(len(seqs["MT"].sequence), 16569)
        self.assertEqual(seqs["MT"].md5, "c68f52674c9fb33aef52dcf399755519")
        self.assertEqual(seqs["MT"].start, 1)
        self.assertEqual(seqs["MT"].end, 16569)
        self.assertEqual(seqs["MT"].length, 16569)

    def test_reference_sequence_generic_load_sequence(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            dest = os.path.join(tmpdir, "generic.fa.gz")
            gzip_file("tests/input/generic.fa", dest)

            seq = ReferenceSequence(ID="1", path=dest)
            self.assertEqual(seq.ID, "1")
            self.assertEqual(seq.chrom, "1")
            self.assertEqual(seq.path, dest)
            np.testing.assert_array_equal(
                seq.sequence,
                np.array(
                    bytearray(
                        "NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNAGGCCGGACNNNNNNNN",
                        encoding="utf-8",
                        errors="strict",
                    ),
                    dtype=np.uint8,
                ),
            )
            self.assertListEqual(
                list("AGGCCGGAC"), list(map(chr, seq.sequence[100:109]))
            )
            self.assertEqual(seq.md5, "6ac6176535ad0e38aba2d05d786c39b6")
            self.assertEqual(seq.start, 1)
            self.assertEqual(seq.end, 117)
            self.assertEqual(seq.length, 117)

    def test_get_opensnp_datadump_filenames(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            # temporarily set resources dir to tests
            self.resource._resources_dir = tmpdir

            # write test openSNP datadump zip
            with atomic_write(
                os.path.join(tmpdir, "opensnp_datadump.current.zip"),
                mode="wb",
                overwrite=True,
            ) as f:
                with zipfile.ZipFile(f, "w") as f_zip:
                    f_zip.write("tests/input/generic.csv", arcname="generic1.csv")
                    f_zip.write("tests/input/generic.csv", arcname="generic2.csv")

            filenames = self.resource.get_opensnp_datadump_filenames()

            self.assertListEqual(filenames, ["generic1.csv", "generic2.csv"])

            self.resource._resources_dir = "resources"

    def test_load_opensnp_datadump_file(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            # temporarily set resources dir to tests
            self.resource._resources_dir = tmpdir

            # write test openSNP datadump zip
            with atomic_write(
                os.path.join(tmpdir, "opensnp_datadump.current.zip"),
                mode="wb",
                overwrite=True,
            ) as f:
                with zipfile.ZipFile(f, "w") as f_zip:
                    f_zip.write("tests/input/generic.csv", arcname="generic1.csv")
                    f_zip.write("tests/input/generic.csv", arcname="generic2.csv")

            snps1 = SNPs(self.resource.load_opensnp_datadump_file("generic1.csv"))
            snps2 = SNPs(self.resource.load_opensnp_datadump_file("generic2.csv"))

            pd.testing.assert_frame_equal(snps1.snps, self.generic_snps())
            pd.testing.assert_frame_equal(snps2.snps, self.generic_snps())

            self.resource._resources_dir = "resources"
