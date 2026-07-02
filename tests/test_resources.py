import os
import tempfile
import unittest

import numpy as np

from snps.resources import ReferenceSequence, Resources
from snps.utils import gzip_file
from tests import BaseSNPsTestCase
from tests.support import FakeResources

# Live integration tests contact the real resource servers and run only when downloads
# are explicitly enabled (the dedicated CI live-integration job).
DOWNLOADS_ENABLED = os.getenv("DOWNLOADS_ENABLED") == "true"


class TestResources(BaseSNPsTestCase):
    # --- offline behavior, served by the fixture-backed provider --------------

    def test_get_assembly_mapping_data(self):
        data = FakeResources().get_assembly_mapping_data("NCBI36", "GRCh37")
        self.assertEqual(sorted(data.keys()), ["1", "3"])

    def test_get_gsa_resources(self):
        gsa_resources = FakeResources().get_gsa_resources()
        self.assertEqual(len(gsa_resources["rsid_map"]), 8)
        self.assertEqual(len(gsa_resources["chrpos_map"]), 8)
        self.assertEqual(len(gsa_resources["dbsnp_151_37_reverse"]), 2)

    def test_get_chip_clusters(self):
        self.assertEqual(len(FakeResources().get_chip_clusters()), 8)

    def test_get_low_quality_snps(self):
        self.assertEqual(len(FakeResources().get_low_quality_snps()), 3)

    def test_get_reference_sequences(self):
        seqs = FakeResources().get_reference_sequences(chroms=["1"])
        self.assertEqual(len(seqs), 1)
        self.assertEqual(seqs["1"].ID, "1")
        self.assertEqual(seqs["1"].chrom, "1")
        self.assertEqual(seqs["1"].assembly, "GRCh37")
        self.assertEqual(seqs["1"].build, "B37")
        self.assertEqual(len(seqs["1"].sequence), 117)
        self.assertEqual(seqs["1"].md5, "6ac6176535ad0e38aba2d05d786c39b6")

    def test_get_par_lookup(self):
        fr = FakeResources()
        self.assertEqual(fr.get_par_lookup("rs758419898")["refsnp_id"], "758419898")
        # a merged RefSNP resolves to the snapshot it was merged into
        self.assertEqual(fr.get_par_lookup("rs113378274")["refsnp_id"], "72608386")

    def test_get_reference_sequences_invalid_assembly(self):
        self.assertEqual(FakeResources().get_reference_sequences(assembly="36"), {})

    def test_get_paths_reference_sequences_invalid_assembly(self):
        assembly, chroms, urls, paths = Resources(
            resources_dir="resources"
        )._get_paths_reference_sequences(assembly="36")
        self.assertFalse(assembly)
        self.assertFalse(chroms)
        self.assertFalse(urls)
        self.assertFalse(paths)

    def test_create_example_datasets(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            paths = Resources(resources_dir=tmpdir).create_example_datasets(tmpdir)

            self.assertEqual(len(paths), 2)
            self.assertTrue(paths[0].endswith("sample1.23andme.txt.gz"))
            self.assertTrue(paths[1].endswith("sample2.ftdna.csv.gz"))
            self.assertTrue(os.path.exists(paths[0]))
            self.assertTrue(os.path.exists(paths[1]))

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

    # --- cached files load without re-downloading (offline) -------------------

    def test_cache_layout_resolves_without_download(self):
        # An already-cached resource must load without downloading; this pins the cache
        # paths/filenames so existing user caches keep working.
        with tempfile.TemporaryDirectory() as tmpdir:
            r = Resources(resources_dir=tmpdir)

            gzip_file(
                "tests/resources/chip_clusters.tsv",
                os.path.join(tmpdir, "chip_clusters.tsv.gz"),
            )
            gzip_file(
                "tests/resources/gsa_rsid_map.txt",
                os.path.join(tmpdir, "gsa_rsid_map.txt.gz"),
            )
            os.makedirs(os.path.join(tmpdir, "fasta", "GRCh37"))
            gzip_file(
                "tests/input/generic.fa",
                os.path.join(
                    tmpdir,
                    "fasta",
                    "GRCh37",
                    "Homo_sapiens.GRCh37.dna.chromosome.1.fa.gz",
                ),
            )

            self.assertEqual(len(r.get_chip_clusters()), 8)
            self.assertEqual(len(r.get_gsa_rsid()), 8)
            seqs = r.get_reference_sequences(chroms=["1"])
            self.assertEqual(seqs["1"].ID, "1")
            self.assertEqual(len(seqs["1"].sequence), 117)

    # --- live integration (real servers; gated to the live-integration job) ---

    @unittest.skipUnless(DOWNLOADS_ENABLED, "live downloads disabled")
    def test_get_assembly_mapping_data_live(self):
        data = Resources().get_assembly_mapping_data("NCBI36", "GRCh37")
        self.assertEqual(len(data), 25)

    @unittest.skipUnless(DOWNLOADS_ENABLED, "live downloads disabled")
    def test_get_gsa_resources_live(self):
        gsa_resources = Resources().get_gsa_resources()
        self.assertEqual(len(gsa_resources["rsid_map"]), 618540)
        self.assertEqual(len(gsa_resources["chrpos_map"]), 665608)
        self.assertEqual(len(gsa_resources["dbsnp_151_37_reverse"]), 2393418)

    @unittest.skipUnless(DOWNLOADS_ENABLED, "live downloads disabled")
    def test_get_chip_clusters_live(self):
        self.assertEqual(len(Resources().get_chip_clusters()), 2135214)

    @unittest.skipUnless(DOWNLOADS_ENABLED, "live downloads disabled")
    def test_get_low_quality_snps_live(self):
        self.assertEqual(len(Resources().get_low_quality_snps()), 56025)

    @unittest.skipUnless(DOWNLOADS_ENABLED, "live downloads disabled")
    def test_get_par_lookup_live(self):
        response = Resources().get_par_lookup("rs28736870")
        self.assertEqual(response["refsnp_id"], "28736870")
        self.assertIn("primary_snapshot_data", response)

    @unittest.skipUnless(DOWNLOADS_ENABLED, "live downloads disabled")
    def test_get_reference_sequences_live(self):
        # each assembly uses a distinct Ensembl FASTA base URL, so verify all three
        expected_urls = {
            "NCBI36": "ftp://ftp.ensembl.org/pub/release-54/fasta/homo_sapiens/dna/Homo_sapiens.NCBI36.54.dna.chromosome.MT.fa.gz",
            "GRCh37": "ftp://ftp.ensembl.org/pub/grch37/release-96/fasta/homo_sapiens/dna/Homo_sapiens.GRCh37.dna.chromosome.MT.fa.gz",
            "GRCh38": "ftp://ftp.ensembl.org/pub/release-96/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.chromosome.MT.fa.gz",
        }
        r = Resources()
        for assembly, url in expected_urls.items():
            seqs = r.get_reference_sequences(assembly=assembly, chroms=["MT"])
            self.assertEqual(seqs["MT"].url, url)
            self.assertEqual(seqs["MT"].assembly, assembly)
            self.assertEqual(seqs["MT"].build, f"B{assembly[-2:]}")
            self.assertGreater(len(seqs["MT"].sequence), 0)

        # GRCh37 MT is the rCRS; pin its exact size and hash
        grch37 = r.get_reference_sequences(assembly="GRCh37", chroms=["MT"])
        self.assertEqual(len(grch37["MT"].sequence), 16569)
        self.assertEqual(grch37["MT"].md5, "c68f52674c9fb33aef52dcf399755519")
