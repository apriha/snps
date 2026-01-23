"""Test synthetic SNP data generator."""

from __future__ import annotations

import gzip
import os
import tempfile

import numpy as np
import pandas as pd

from snps import SNPs
from snps.build_constants import BUILD_MARKER_SNPS
from snps.io.generator import SyntheticSNPGenerator
from tests import BaseSNPsTestCase


class TestSyntheticSNPGenerator(BaseSNPsTestCase):
    def test_init_build_37(self) -> None:
        """Test generator initialization with Build 37."""
        gen = SyntheticSNPGenerator(build=37, seed=123)
        self.assertEqual(gen.build, 37)
        self.assertEqual(gen.seed, 123)
        self.assertIsNotNone(gen.chrom_sizes)

    def test_init_build_36(self) -> None:
        """Test generator initialization with Build 36."""
        gen = SyntheticSNPGenerator(build=36, seed=456)
        self.assertEqual(gen.build, 36)
        self.assertEqual(gen.seed, 456)

    def test_init_build_38(self) -> None:
        """Test generator initialization with Build 38."""
        gen = SyntheticSNPGenerator(build=38, seed=789)
        self.assertEqual(gen.build, 38)
        self.assertEqual(gen.seed, 789)

    def test_init_invalid_build(self) -> None:
        """Test generator initialization with invalid build."""
        with self.assertRaises(ValueError):
            SyntheticSNPGenerator(build=99)

    def test_generate_snps_empty_chromosomes(self) -> None:
        """Test SNP generation with empty chromosomes list."""
        gen = SyntheticSNPGenerator(build=37, seed=123)
        df = gen.generate_snps(num_snps=1000, chromosomes=[])

        # Verify empty DataFrame with correct structure
        self.assertIsInstance(df, pd.DataFrame)
        self.assertEqual(len(df), 0)
        self.assertEqual(df.index.name, "rsid")
        self.assertListEqual(list(df.columns), ["chrom", "pos", "genotype"])

        # Verify dtypes
        self.assertEqual(df["chrom"].dtype, object)
        self.assertEqual(df["pos"].dtype, np.uint32)
        self.assertEqual(df["genotype"].dtype, object)

    def test_generate_snps_basic(self) -> None:
        """Test basic SNP generation."""
        gen = SyntheticSNPGenerator(build=37, seed=123)
        df = gen.generate_snps(num_snps=1000)

        # Verify DataFrame structure
        self.assertIsInstance(df, pd.DataFrame)
        self.assertEqual(df.index.name, "rsid")
        self.assertListEqual(list(df.columns), ["chrom", "pos", "genotype"])

        # Verify dtypes
        self.assertEqual(df["chrom"].dtype, object)
        self.assertEqual(df["pos"].dtype, np.uint32)
        self.assertEqual(df["genotype"].dtype, object)

        # Verify approximate SNP count (may vary due to distribution)
        self.assertGreater(len(df), 900)
        self.assertLess(len(df), 1100)

    def test_generate_snps_chromosomes(self) -> None:
        """Test SNP generation with specific chromosomes."""
        gen = SyntheticSNPGenerator(build=37, seed=123)
        df = gen.generate_snps(num_snps=500, chromosomes=["1", "2", "X"])

        # Verify only specified chromosomes are present
        unique_chroms = set(df["chrom"].unique())
        self.assertTrue(unique_chroms.issubset({"1", "2", "X"}))

    def test_generate_snps_missing_rate(self) -> None:
        """Test SNP generation with missing genotypes."""
        gen = SyntheticSNPGenerator(build=37, seed=123)
        df = gen.generate_snps(num_snps=10000, missing_rate=0.1)

        # Count missing genotypes
        missing_count = (df["genotype"] == "--").sum()
        total_count = len(df)

        # Verify approximately 10% are missing (with some tolerance)
        missing_rate = missing_count / total_count
        self.assertGreater(missing_rate, 0.05)  # At least 5%
        self.assertLess(missing_rate, 0.15)  # At most 15%

    def test_generate_snps_reproducibility(self) -> None:
        """Test that same seed produces same results."""
        gen1 = SyntheticSNPGenerator(build=37, seed=42)
        gen2 = SyntheticSNPGenerator(build=37, seed=42)

        df1 = gen1.generate_snps(num_snps=1000)
        df2 = gen2.generate_snps(num_snps=1000)

        # Verify identical output
        self.assert_frame_equal_with_string_index(df1, df2)

    def test_generate_snps_rsid_format(self) -> None:
        """Test that rsIDs are correctly formatted."""
        gen = SyntheticSNPGenerator(build=37, seed=123)
        df = gen.generate_snps(num_snps=100)

        # Verify all rsIDs start with "rs" followed by digits
        for rsid in df.index:
            self.assertTrue(rsid.startswith("rs"))
            self.assertTrue(rsid[2:].isdigit())

    def test_generate_snps_positions_sorted(self) -> None:
        """Test that positions within each chromosome are sorted."""
        gen = SyntheticSNPGenerator(build=37, seed=123)
        df = gen.generate_snps(num_snps=5000)

        # Check each chromosome separately
        for chrom in df["chrom"].unique():
            chrom_df = df[df["chrom"] == chrom]
            positions = chrom_df["pos"].values
            # Verify positions are sorted
            self.assertTrue(np.all(positions[:-1] <= positions[1:]))

    def test_generate_snps_genotype_format(self) -> None:
        """Test that genotypes are correctly formatted."""
        gen = SyntheticSNPGenerator(build=37, seed=123)
        df = gen.generate_snps(num_snps=1000, missing_rate=0.05)

        # Valid bases
        valid_bases = {"A", "C", "G", "T", "-"}

        # Verify all genotypes are valid
        for genotype in df["genotype"].unique():
            if pd.notna(genotype):
                # Should be either 2 characters or "--"
                self.assertIn(len(genotype), [2])
                for base in genotype:
                    self.assertIn(base, valid_bases)

    def test_save_as_23andme(self) -> None:
        """Test saving SNPs in 23andMe format."""
        gen = SyntheticSNPGenerator(build=37, seed=123)

        with tempfile.TemporaryDirectory() as tmpdir:
            output_path = os.path.join(tmpdir, "test_23andme.txt.gz")
            result_path = gen.save_as_23andme(output_path, num_snps=10000)

            # Verify file was created
            self.assertEqual(result_path, output_path)
            self.assertTrue(os.path.exists(output_path))

            # Verify file can be read
            with gzip.open(output_path, "rt") as f:
                lines = f.readlines()

            # Verify header
            self.assertTrue(lines[0].startswith("# 23andMe"))

            # Verify format (tab-separated)
            data_lines = [line for line in lines if not line.startswith("#")]
            self.assertGreater(len(data_lines), 9000)  # Should have ~10000 SNPs
            first_data = data_lines[0].strip().split("\t")
            self.assertEqual(len(first_data), 4)  # rsid, chrom, pos, genotype

            # Verify SNPs can load the file
            s = SNPs(output_path)
            self.assertEqual(s.source, "23andMe")
            # Build detection should work with injected marker SNPs
            self.assertEqual(s.build, 37)
            self.assertTrue(s.build_detected)

    def test_save_as_ancestry(self) -> None:
        """Test saving SNPs in AncestryDNA format."""
        gen = SyntheticSNPGenerator(build=37, seed=456)

        with tempfile.TemporaryDirectory() as tmpdir:
            output_path = os.path.join(tmpdir, "test_ancestry.txt.gz")
            result_path = gen.save_as_ancestry(output_path, num_snps=5000)

            # Verify file was created
            self.assertEqual(result_path, output_path)
            self.assertTrue(os.path.exists(output_path))

            # Verify file can be read
            with gzip.open(output_path, "rt") as f:
                lines = f.readlines()

            # Verify header
            self.assertTrue(lines[0].startswith("#AncestryDNA"))

            # Verify format (allele1 and allele2 columns)
            header_line = [line for line in lines if line.startswith("rsid")][0]
            self.assertIn("allele1", header_line)
            self.assertIn("allele2", header_line)

            # Verify SNPs can load the file
            s = SNPs(output_path)
            self.assertEqual(s.source, "AncestryDNA")

    def test_save_as_ftdna(self) -> None:
        """Test saving SNPs in FTDNA format."""
        gen = SyntheticSNPGenerator(build=36, seed=789)

        with tempfile.TemporaryDirectory() as tmpdir:
            output_path = os.path.join(tmpdir, "test_ftdna.csv.gz")
            result_path = gen.save_as_ftdna(output_path, num_snps=7000)

            # Verify file was created
            self.assertEqual(result_path, output_path)
            self.assertTrue(os.path.exists(output_path))

            # Verify file can be read
            with gzip.open(output_path, "rt") as f:
                lines = f.readlines()

            # Verify CSV format with quotes
            header = lines[0].strip()
            self.assertEqual(header, "RSID,CHROMOSOME,POSITION,RESULT")

            # Verify data line format
            self.assertTrue(lines[1].startswith('"rs'))

            # Verify SNPs can load the file
            s = SNPs(output_path)
            self.assertEqual(s.source, "FTDNA")
            # Build detection should work with injected marker SNPs
            self.assertEqual(s.build, 36)
            self.assertTrue(s.build_detected)

    def test_save_as_generic_csv(self) -> None:
        """Test saving SNPs in generic CSV format."""
        gen = SyntheticSNPGenerator(build=37, seed=202)

        with tempfile.TemporaryDirectory() as tmpdir:
            output_path = os.path.join(tmpdir, "test.csv.gz")
            result_path = gen.save_as_generic(output_path, format="csv", num_snps=2000)

            # Verify file was created
            self.assertEqual(result_path, output_path)
            self.assertTrue(os.path.exists(output_path))

            # Verify file can be read
            with gzip.open(output_path, "rt") as f:
                lines = f.readlines()

            # Verify CSV header
            header = lines[0].strip()
            self.assertEqual(header, "rsid,chrom,pos,genotype")

            # Verify data format (comma-separated with 4 fields)
            data_lines = [line.strip() for line in lines[1:] if line.strip()]
            self.assertGreater(len(data_lines), 0, "No data lines found")
            self.assertEqual(len(data_lines[0].split(",")), 4)

            # Verify SNPs can load the file
            s = SNPs(output_path)
            self.assertGreater(s.count, 1900)

    def test_save_as_generic_tsv(self) -> None:
        """Test saving SNPs in generic TSV format."""
        gen = SyntheticSNPGenerator(build=37, seed=303)

        with tempfile.TemporaryDirectory() as tmpdir:
            output_path = os.path.join(tmpdir, "test.tsv")
            result_path = gen.save_as_generic(output_path, format="tsv", num_snps=1500)

            # Verify file was created (not gzipped)
            self.assertEqual(result_path, output_path)
            self.assertTrue(os.path.exists(output_path))

            # Verify file can be read
            with open(output_path, "r") as f:
                lines = f.readlines()

            # Verify TSV header
            header = lines[0].strip()
            self.assertEqual(header, "rsid\tchrom\tpos\tgenotype")

            # Verify data format (tab-separated with 4 fields)
            data_lines = [line.strip() for line in lines[1:] if line.strip()]
            self.assertGreater(len(data_lines), 0, "No data lines found")
            self.assertEqual(len(data_lines[0].split("\t")), 4)

    def test_save_creates_directory(self) -> None:
        """Test that save methods create output directory if needed."""
        gen = SyntheticSNPGenerator(build=37, seed=404)

        with tempfile.TemporaryDirectory() as tmpdir:
            # Create a nested path that doesn't exist
            output_path = os.path.join(tmpdir, "subdir1", "subdir2", "test.txt.gz")
            gen.save_as_23andme(output_path, num_snps=1000)

            # Verify file was created
            self.assertTrue(os.path.exists(output_path))

    def test_large_dataset(self) -> None:
        """Test generating a large dataset."""
        gen = SyntheticSNPGenerator(build=37, seed=505)

        with tempfile.TemporaryDirectory() as tmpdir:
            output_path = os.path.join(tmpdir, "large.txt.gz")
            gen.save_as_23andme(output_path, num_snps=1000000)

            # Verify file was created
            self.assertTrue(os.path.exists(output_path))

            # Verify SNPs can load it
            s = SNPs(output_path)
            self.assertGreater(s.count, 950000)  # Should be close to 1M
            self.assertLess(s.count, 1050000)

    def test_different_builds_have_different_chromosome_sizes(self) -> None:
        """Test that different builds use different chromosome size ranges."""
        gen36 = SyntheticSNPGenerator(build=36)
        gen37 = SyntheticSNPGenerator(build=37)
        gen38 = SyntheticSNPGenerator(build=38)

        # Verify different chromosome sizes are used
        self.assertNotEqual(gen36.chrom_sizes["1"], gen37.chrom_sizes["1"])
        self.assertNotEqual(gen37.chrom_sizes["1"], gen38.chrom_sizes["1"])
        self.assertNotEqual(gen36.chrom_sizes["1"], gen38.chrom_sizes["1"])

    def test_build_marker_injection(self) -> None:
        """Test that build marker SNPs are injected correctly."""
        # Test Build 37
        gen37 = SyntheticSNPGenerator(build=37, seed=123)
        df37 = gen37.generate_snps(num_snps=10000, inject_build_markers=True)

        # Check that at least one marker SNP was injected
        found_markers = [rsid for rsid in BUILD_MARKER_SNPS if rsid in df37.index]
        self.assertGreater(len(found_markers), 0)

        # Verify marker positions are correct for Build 37
        if "rs3094315" in df37.index:
            self.assertEqual(df37.loc["rs3094315"]["pos"], 752566)
            self.assertEqual(df37.loc["rs3094315"]["chrom"], "1")

        # Test Build 36
        gen36 = SyntheticSNPGenerator(build=36, seed=456)
        df36 = gen36.generate_snps(num_snps=10000, inject_build_markers=True)

        # Verify marker position is different for Build 36
        if "rs3094315" in df36.index:
            self.assertEqual(df36.loc["rs3094315"]["pos"], 742429)
            self.assertEqual(df36.loc["rs3094315"]["chrom"], "1")

    def test_build_detection_with_markers(self) -> None:
        """Test that SNPs with injected markers can be correctly build-detected."""
        gen37 = SyntheticSNPGenerator(build=37, seed=789)

        with tempfile.TemporaryDirectory() as tmpdir:
            output_path = os.path.join(tmpdir, "test_build37.txt.gz")
            gen37.save_as_23andme(output_path, num_snps=50000)

            # Load and verify build detection works
            s = SNPs(output_path)
            self.assertEqual(s.build, 37)
            self.assertTrue(s.build_detected)

    def test_generate_without_markers(self) -> None:
        """Test generating SNPs without build markers."""
        gen = SyntheticSNPGenerator(build=37, seed=999)
        df = gen.generate_snps(num_snps=1000, inject_build_markers=False)

        # Verify no marker SNPs are present
        found_markers = [rsid for rsid in BUILD_MARKER_SNPS if rsid in df.index]
        self.assertEqual(len(found_markers), 0)

    def test_inject_build_markers_no_matching_chromosomes(self) -> None:
        """Test _inject_build_markers when chromosomes have no marker SNPs.

        Marker SNPs are on chromosomes 1, 2, and 20. Using chromosomes
        without markers should skip marker injection (empty marker_rows branch).
        """
        gen = SyntheticSNPGenerator(build=37, seed=123)

        # Generate SNPs only on chromosomes that don't have marker SNPs
        # Marker SNPs are on chroms: 1, 2, 20
        df = gen.generate_snps(
            num_snps=500, chromosomes=["MT", "X", "Y"], inject_build_markers=True
        )

        # Verify DataFrame was created with correct structure
        self.assertIsInstance(df, pd.DataFrame)
        self.assertGreater(len(df), 0)
        self.assertEqual(df.index.name, "rsid")

        # Verify only the specified chromosomes are present
        unique_chroms = set(df["chrom"].unique())
        self.assertTrue(unique_chroms.issubset({"MT", "X", "Y"}))

        # Verify no marker SNPs were injected (since they're on other chromosomes)
        found_markers = [rsid for rsid in BUILD_MARKER_SNPS if rsid in df.index]
        self.assertEqual(len(found_markers), 0)

    def test_create_example_dataset_pair(self) -> None:
        """Test creating a pair of example datasets for merge demonstration."""
        gen = SyntheticSNPGenerator(build=37, seed=47)

        with tempfile.TemporaryDirectory() as tmpdir:
            path1, path2 = gen.create_example_dataset_pair(tmpdir)

            # Verify files were created
            self.assertTrue(os.path.exists(path1))
            self.assertTrue(os.path.exists(path2))
            self.assertTrue(path1.endswith("sample1.23andme.txt.gz"))
            self.assertTrue(path2.endswith("sample2.ftdna.csv.gz"))

            # Load both files
            s1 = SNPs(path1)
            s2 = SNPs(path2)

            # Verify sources
            self.assertEqual(s1.source, "23andMe")
            self.assertEqual(s2.source, "FTDNA")

            # Verify builds
            self.assertEqual(s1.build, 37)
            self.assertEqual(s2.build, 37)

            # Verify file 1 includes all chromosomes (1-22, X, Y, MT)
            s1_chroms = set(s1.snps["chrom"].unique())
            self.assertIn("X", s1_chroms)
            self.assertIn("Y", s1_chroms)
            self.assertIn("MT", s1_chroms)

            # Verify file 2 only has autosomal chromosomes (1-22)
            s2_chroms = set(s2.snps["chrom"].unique())
            self.assertNotIn("X", s2_chroms)
            self.assertNotIn("Y", s2_chroms)
            self.assertNotIn("MT", s2_chroms)

            # Verify approximate SNP counts
            # File 1: ~700K base + ~292K unique = ~992K
            self.assertGreater(s1.count, 950000)
            self.assertLess(s1.count, 1050000)

            # File 2: autosomal base SNPs + ~15K unique (no X, Y, MT)
            self.assertGreater(s2.count, 600000)
            self.assertLess(s2.count, 750000)

            # Verify shared RSIDs exist (base SNPs rs1-rs700000)
            common_rsids = set(s1.snps.index) & set(s2.snps.index)
            self.assertGreater(len(common_rsids), 600000)

    def test_unique_positions_small_chromosome(self) -> None:
        """Test that positions are unique even on small chromosomes like MT."""
        gen = SyntheticSNPGenerator(build=37, seed=42)

        # MT chromosome is only ~16,569 bp, generate many SNPs
        df = gen.generate_snps(num_snps=10000, chromosomes=["MT"])

        # Check for duplicate positions within chromosome
        duplicates = df.duplicated(subset=["chrom", "pos"], keep=False)
        self.assertEqual(duplicates.sum(), 0, "Found duplicate positions on MT")

        # Verify all positions are within valid range
        self.assertTrue((df["pos"] >= 1).all())
        self.assertTrue((df["pos"] <= gen.chrom_sizes["MT"]).all())

    def test_unique_positions_large_dataset(self) -> None:
        """Test that positions are unique in large datasets."""
        gen = SyntheticSNPGenerator(build=37, seed=123)

        # Generate 100K SNPs across all chromosomes
        df = gen.generate_snps(num_snps=100000)

        # Check for duplicate positions within each chromosome
        duplicates = df.duplicated(subset=["chrom", "pos"], keep=False)
        self.assertEqual(duplicates.sum(), 0, "Found duplicate positions")
