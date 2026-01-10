import gzip
import os
import socket
import tempfile
import urllib.error
from unittest.mock import Mock, mock_open, patch

import numpy as np

from snps import SNPs
from snps.resources import ReferenceSequence, Resources
from snps.utils import gzip_file
from tests import BaseSNPsTestCase


class TestResources(BaseSNPsTestCase):
    def _reset_resource(self):
        self.resource._init_resource_attributes()

    def run(self, result=None):
        # set resources directory based on if downloads are being performed
        # https://stackoverflow.com/a/11180583

        self.resource = Resources()
        self._reset_resource()
        if self.downloads_enabled:
            self.resource._resources_dir = "resources"
            super().run(result)
        else:
            # use a temporary directory for test resource data
            with tempfile.TemporaryDirectory() as tmpdir:
                self.resource._resources_dir = tmpdir
                super().run(result)
                self.resource._resources_dir = "resources"

    def test_get_assembly_mapping_data(self):
        def f():
            effects = [{"mappings": []} for _ in range(1, 26)]
            for k, v in self.NCBI36_GRCh37().items():
                effects[int(k) - 1] = v
            mock = Mock(side_effect=effects)
            with patch("snps.ensembl.EnsemblRestClient.perform_rest_action", mock):
                return self.resource.get_assembly_mapping_data("NCBI36", "GRCh37")

        assembly_mapping_data = (
            self.resource.get_assembly_mapping_data("NCBI36", "GRCh37")
            if self.downloads_enabled
            else f()
        )

        self.assertEqual(len(assembly_mapping_data), 25)

    def test_get_gsa_resources(self):
        def f():
            # mock download of test data for each resource
            self._generate_test_gsa_resources()
            # load test resources saved to `tmpdir`
            return self.resource.get_gsa_resources()

        gsa_resources = (
            self.resource.get_gsa_resources() if self.downloads_enabled else f()
        )

        self.assertEqual(len(gsa_resources["rsid_map"]), 618540)
        self.assertEqual(len(gsa_resources["chrpos_map"]), 665608)
        self.assertEqual(len(gsa_resources["dbsnp_151_37_reverse"]), 2393418)

    def _generate_test_gsa_resources(self):
        lines = ["Name\tRsID"]

        for i in range(1, 618541):
            lines.append(f"rs{i}\trs{i}")

        s = "\n".join(lines)
        mock = mock_open(read_data=gzip.compress(s.encode()))
        with patch("urllib.request.urlopen", mock):
            self.resource.get_gsa_rsid()

        lines = ["Name\tChr\tMapInfo\tdeCODE(cM)"]

        for i in range(1, 665609):
            lines.append(f"rs{i}\t1\t{i}\t0.0000")

        s = "\n".join(lines)

        mock = mock_open(read_data=gzip.compress(s.encode()))
        with patch("urllib.request.urlopen", mock):
            self.resource.get_gsa_chrpos()

        lines = ["# comment", "rs1 0.0 0.0 0.0 0.0"]

        for i in range(2, 2393419):
            lines.append(f"rs{i}")

        s = "\n".join(lines)

        mock = mock_open(read_data=gzip.compress(s.encode()))
        with patch("urllib.request.urlopen", mock):
            self.resource.get_dbsnp_151_37_reverse()

    def test_get_all_resources(self):
        def f():
            # mock download of test data for each resource
            self._generate_test_gsa_resources()
            self._generate_test_chip_clusters()
            self._generate_test_low_quality_snps()

            # generate test data for permutations of remapping data
            effects = [{"mappings": []} for _ in range(1, 26)]
            for k, v in self.NCBI36_GRCh37().items():
                effects[int(k) - 1] = v
            mock = Mock(side_effect=effects * 6)
            with patch("snps.ensembl.EnsemblRestClient.perform_rest_action", mock):
                return self.resource.get_all_resources()

        resources = self.resource.get_all_resources() if self.downloads_enabled else f()

        for k, v in resources.items():
            self.assertGreater(len(v), 0)

    def test_get_paths_reference_sequences_invalid_assembly(self):
        assembly, chroms, urls, paths = self.resource._get_paths_reference_sequences(
            assembly="36"
        )
        self.assertFalse(assembly)
        self.assertFalse(chroms)
        self.assertFalse(urls)
        self.assertFalse(paths)

    def run_reference_sequences_test(self, f, assembly="GRCh37"):
        if self.downloads_enabled:
            f()
        else:
            s = f">MT dna:chromosome chromosome:{assembly}:MT:1:16569:1 REF\n"
            for i in range(276):
                s += "A" * 60
                s += "\n"
            s += "A" * 9
            s += "\n"
            with patch(
                "urllib.request.urlopen", mock_open(read_data=gzip.compress(s.encode()))
            ):
                f()

    def run_create_reference_sequences_test(self, assembly_expect, url_expect):
        def f():
            (
                assembly,
                chroms,
                urls,
                paths,
            ) = self.resource._get_paths_reference_sequences(
                assembly=assembly_expect, chroms=["MT"]
            )
            seqs = self.resource._create_reference_sequences(
                assembly, chroms, urls, paths
            )
            self.assertEqual(len(seqs), 1)
            self.assertEqual(
                seqs["MT"].__repr__(),
                f"ReferenceSequence(assembly='{assembly_expect}', ID='MT')",
            )
            self.assertEqual(seqs["MT"].ID, "MT")
            self.assertEqual(seqs["MT"].chrom, "MT")
            self.assertEqual(seqs["MT"].url, f"{url_expect}")
            self.assertEqual(
                seqs["MT"].path,
                os.path.relpath(
                    f"{os.path.join(self.resource._resources_dir, 'fasta', assembly_expect, os.path.basename(url_expect))}"
                ),
            )
            self.assertTrue(os.path.exists(seqs["MT"].path))
            self.assertEqual(seqs["MT"].assembly, assembly_expect)
            self.assertEqual(seqs["MT"].build, f"B{assembly_expect[-2:]}")
            self.assertEqual(seqs["MT"].species, "Homo sapiens")
            self.assertEqual(seqs["MT"].taxonomy, "x")

        self.run_reference_sequences_test(f, assembly_expect)

    def test_create_reference_sequences_NCBI36(self):
        self.run_create_reference_sequences_test(
            "NCBI36",
            "ftp://ftp.ensembl.org/pub/release-54/fasta/homo_sapiens/dna/Homo_sapiens.NCBI36.54.dna.chromosome.MT.fa.gz",
        )

    def test_create_reference_sequences_GRCh37(self):
        self.run_create_reference_sequences_test(
            "GRCh37",
            "ftp://ftp.ensembl.org/pub/grch37/release-96/fasta/homo_sapiens/dna/Homo_sapiens.GRCh37.dna.chromosome.MT.fa.gz",
        )

    def test_create_reference_sequences_GRCh38(self):
        self.run_create_reference_sequences_test(
            "GRCh38",
            "ftp://ftp.ensembl.org/pub/release-96/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.chromosome.MT.fa.gz",
        )

    def test_create_reference_sequences_invalid_path(self):
        def f():
            (
                assembly,
                chroms,
                urls,
                paths,
            ) = self.resource._get_paths_reference_sequences(
                assembly="GRCh37", chroms=["MT"]
            )
            paths[0] = ""
            seqs = self.resource._create_reference_sequences(
                assembly, chroms, urls, paths
            )
            self.assertEqual(len(seqs), 0)

        self.run_reference_sequences_test(f)

    def test_download_file_socket_timeout(self):
        mock = Mock(side_effect=socket.timeout)
        with patch("urllib.request.urlopen", mock):
            path = self.resource._download_file("http://url", "test.txt")
        self.assertEqual(path, "")

    def test_download_file_URL_error(self):
        mock = Mock(side_effect=urllib.error.URLError("test error"))
        with patch("urllib.request.urlopen", mock):
            path1 = self.resource._download_file("http://url", "test.txt")
            path2 = self.resource._download_file("ftp://url", "test.txt")
        self.assertEqual(path1, "")
        self.assertEqual(path2, "")

    def test_get_reference_sequences(self):
        def f():
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
                os.path.relpath(
                    f"{os.path.join(self.resource._resources_dir, 'fasta', 'GRCh37', 'Homo_sapiens.GRCh37.dna.chromosome.MT.fa.gz')}"
                ),
            )
            self.assertTrue(os.path.exists(seqs["MT"].path))
            self.assertEqual(seqs["MT"].assembly, "GRCh37")
            self.assertEqual(seqs["MT"].build, "B37")
            self.assertEqual(seqs["MT"].species, "Homo sapiens")
            self.assertEqual(seqs["MT"].taxonomy, "x")

        self.run_reference_sequences_test(f)

    def test_get_all_reference_sequences(self):
        def f():
            seqs = self.resource.get_all_reference_sequences(chroms=["MT"])
            self.assertEqual(len(seqs), 3)
            self.assertEqual(len(seqs["NCBI36"]), 1)
            self.assertEqual(
                seqs["NCBI36"]["MT"].path,
                os.path.relpath(
                    os.path.join(
                        self.resource._resources_dir,
                        "fasta",
                        "NCBI36",
                        "Homo_sapiens.NCBI36.54.dna.chromosome.MT.fa.gz",
                    )
                ),
            )
            self.assertEqual(len(seqs["GRCh37"]), 1)
            self.assertEqual(
                seqs["GRCh37"]["MT"].path,
                os.path.relpath(
                    os.path.join(
                        self.resource._resources_dir,
                        "fasta",
                        "GRCh37",
                        "Homo_sapiens.GRCh37.dna.chromosome.MT.fa.gz",
                    )
                ),
            )
            self.assertEqual(len(seqs["GRCh38"]), 1)
            self.assertEqual(
                seqs["GRCh38"]["MT"].path,
                os.path.relpath(
                    os.path.join(
                        self.resource._resources_dir,
                        "fasta",
                        "GRCh38",
                        "Homo_sapiens.GRCh38.dna.chromosome.MT.fa.gz",
                    )
                ),
            )

        self.run_reference_sequences_test(f)

    def test_get_reference_sequences_invalid_assembly(self):
        seqs = self.resource.get_reference_sequences(assembly="36")
        self.assertEqual(len(seqs), 0)

    def test_get_reference_sequences_chrom_not_available(self):
        def f():
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
                os.path.relpath(
                    os.path.join(
                        self.resource._resources_dir,
                        "fasta",
                        "GRCh37",
                        "Homo_sapiens.GRCh37.dna.chromosome.MT.fa.gz",
                    )
                ),
            )
            self.assertTrue(os.path.exists(seqs["MT"].path))
            self.assertEqual(seqs["MT"].assembly, "GRCh37")
            self.assertEqual(seqs["MT"].build, "B37")
            self.assertEqual(seqs["MT"].species, "Homo sapiens")
            self.assertEqual(seqs["MT"].taxonomy, "x")

        self.run_reference_sequences_test(f)

    def run_reference_sequence_load_sequence_test(self, hash):
        def f():
            seqs = self.resource.get_reference_sequences(chroms=["MT"])
            self.assertEqual(len(seqs["MT"].sequence), 16569)
            self.assertEqual(seqs["MT"].md5, hash)
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
            self.assertEqual(seqs["MT"].md5, hash)
            self.assertEqual(seqs["MT"].start, 1)
            self.assertEqual(seqs["MT"].end, 16569)
            self.assertEqual(seqs["MT"].length, 16569)

        self.run_reference_sequences_test(f)

    def test_reference_sequence_load_sequence(self):
        if self.downloads_enabled:
            self.run_reference_sequence_load_sequence_test(
                "c68f52674c9fb33aef52dcf399755519"
            )
        else:
            self.run_reference_sequence_load_sequence_test(
                "d432324413a21aa9247321c56c300ad3"
            )

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

    def test_create_example_datasets(self):
        """Test creating synthetic example datasets."""
        with tempfile.TemporaryDirectory() as tmpdir:
            paths = self.resource.create_example_datasets(tmpdir)

            # Verify two files were created
            self.assertEqual(len(paths), 2)
            self.assertTrue(os.path.exists(paths[0]))
            self.assertTrue(os.path.exists(paths[1]))

            # Verify filenames
            self.assertTrue(paths[0].endswith("sample1.23andme.txt.gz"))
            self.assertTrue(paths[1].endswith("sample2.ftdna.csv.gz"))

            # Verify files can be loaded
            s1 = SNPs(paths[0])
            s2 = SNPs(paths[1])

            # Verify file sources are detected correctly
            self.assertEqual(s1.source, "23andMe")
            self.assertEqual(s2.source, "FTDNA")

            # Build detection should work with injected marker SNPs
            self.assertEqual(s1.build, 37)
            self.assertTrue(s1.build_detected)
            self.assertEqual(s2.build, 37)
            self.assertTrue(s2.build_detected)

            # Verify SNP counts are approximately correct
            self.assertGreater(s1.count, 900000)
            # FTDNA file only contains autosomal chromosomes (1-22), so count is lower
            self.assertGreater(s2.count, 650000)

            # Verify 23andMe file includes all chromosomes
            s1_chroms = set(s1.snps["chrom"].unique())
            self.assertIn("X", s1_chroms)
            self.assertIn("Y", s1_chroms)
            self.assertIn("MT", s1_chroms)

            # Verify FTDNA file only includes autosomal chromosomes (1-22)
            s2_chroms = set(s2.snps["chrom"].unique())
            self.assertNotIn("X", s2_chroms)
            self.assertNotIn("Y", s2_chroms)
            self.assertNotIn("MT", s2_chroms)

    def _generate_test_chip_clusters(self):
        s = "1:1\tc1\n" * 2135214
        mock = mock_open(read_data=gzip.compress(s.encode()))
        with patch("urllib.request.urlopen", mock):
            self.resource.get_chip_clusters()

    def test_get_chip_clusters(self):
        def f():
            # mock download of test data for chip clusters
            self._generate_test_chip_clusters()
            # load test resource
            return self.resource.get_chip_clusters()

        chip_clusters = (
            self.resource.get_chip_clusters() if self.downloads_enabled else f()
        )

        self.assertEqual(len(chip_clusters), 2135214)

    def _generate_test_low_quality_snps(self):
        s = "c1\t" + "1:1," * 56024 + "1:1\n"
        mock = mock_open(read_data=gzip.compress(s.encode()))
        with patch("urllib.request.urlopen", mock):
            self.resource.get_low_quality_snps()

    def test_get_low_quality_snps(self):
        def f():
            # mock download of test data for low quality SNPs
            self._generate_test_low_quality_snps()
            # load test resource
            return self.resource.get_low_quality_snps()

        low_quality_snps = (
            self.resource.get_low_quality_snps() if self.downloads_enabled else f()
        )

        self.assertEqual(len(low_quality_snps), 56025)
