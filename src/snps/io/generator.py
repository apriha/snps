"""Generate synthetic genotype data for testing and examples."""

from __future__ import annotations

import gzip
import logging
import os
from collections.abc import Callable
from io import StringIO
from typing import Any

import numpy as np
import pandas as pd
from atomicwrites import atomic_write
from numpy.typing import NDArray
from pandas.api.types import CategoricalDtype

from snps.build_constants import BUILD_MARKER_SNPS, CHROM_SIZES, VALID_BUILDS
from snps.constants import REFERENCE_SEQUENCE_CHROMS
from snps.io.reader import get_empty_snps_dataframe

logger = logging.getLogger(__name__)


class SyntheticSNPGenerator:
    """Generate realistic synthetic genotype data.

    This class generates synthetic SNP data that mimics real genotype files
    from various DNA testing companies. The generated data is suitable for
    testing, examples, and documentation.

    Parameters
    ----------
    build : int
        Genome build (36, 37, or 38), default is 37
    seed : int, optional
        Random seed for reproducibility

    Examples
    --------
    >>> gen = SyntheticSNPGenerator(build=37, seed=123)
    >>> gen.save_as_23andme("output.txt", num_snps=10000)
    'output.txt'
    """

    def __init__(self, build: int = 37, seed: int | None = None) -> None:
        if build not in VALID_BUILDS:
            raise ValueError(f"Unsupported build: {build}. Use {VALID_BUILDS}.")
        self.build: int = build
        self.seed: int | None = seed
        self.rng: np.random.Generator = np.random.default_rng(seed)
        self.chrom_sizes: dict[str, int] = CHROM_SIZES[build]

    def generate_snps(
        self,
        num_snps: int = 10000,
        chromosomes: list[str] | None = None,
        missing_rate: float = 0.01,
        inject_build_markers: bool = True,
    ) -> pd.DataFrame:
        """Generate a DataFrame of synthetic SNPs.

        Parameters
        ----------
        num_snps : int
            Approximate number of SNPs to generate
        chromosomes : list of str, optional
            Chromosomes to include (default: all autosomes plus X, Y, MT)
        missing_rate : float
            Proportion of SNPs with missing genotypes (default: 0.01)
        inject_build_markers : bool
            Inject known marker SNPs for build detection (default: True)

        Returns
        -------
        pd.DataFrame
            DataFrame with columns: rsid (index), chrom, pos, genotype
        """
        if chromosomes is None:
            chromosomes = list(REFERENCE_SEQUENCE_CHROMS)

        snps_per_chrom = self._calculate_snps_per_chromosome(num_snps, chromosomes)
        all_snps = []
        rsid_counter = 1

        for chrom in chromosomes:
            n = snps_per_chrom[chrom]
            chrom_size = self.chrom_sizes[chrom]

            # Cap SNP count to available positions (chromosome size - 1)
            max_positions = chrom_size - 1
            n = min(n, max_positions)

            # Generate unique positions efficiently
            positions = self._generate_unique_positions(n, chrom_size)
            genotypes = self._generate_genotypes(n, missing_rate)
            rsids = [f"rs{rsid_counter + i}" for i in range(n)]
            rsid_counter += n

            all_snps.append(
                pd.DataFrame(
                    {
                        "rsid": rsids,
                        "chrom": chrom,
                        "pos": positions,
                        "genotype": genotypes,
                    }
                )
            )

        if not all_snps:
            return get_empty_snps_dataframe()

        df = pd.concat(all_snps, ignore_index=True).set_index("rsid")
        df["chrom"] = df["chrom"].astype(object)
        df["pos"] = df["pos"].astype(np.uint32)
        df["genotype"] = df["genotype"].astype(object)

        if inject_build_markers:
            df = self._inject_build_markers(df, chromosomes)

        # Renumber RSIDs sequentially while preserving build markers
        df = self._renumber_rsids_sequentially(df)

        return df

    def create_example_dataset_pair(self, output_dir: str = ".") -> tuple[str, str]:
        """Create a pair of realistic example datasets suitable for merging.

        Generates two correlated genotype files that share a large number of
        common SNPs, with some discrepancies to demonstrate merge functionality.

        Parameters
        ----------
        output_dir : str
            Directory for output files

        Returns
        -------
        tuple of (str, str)
            Paths to (file1_23andme, file2_ftdna)
        """
        # Generate base dataset (~700K SNPs shared between files)
        base_snps = self._generate_base_snps(700000)

        # File 1: 23andMe format with ~991K SNPs
        file1_df = self._create_file1_dataframe(base_snps, 291786)

        # File 2: FTDNA format with ~715K SNPs, including discrepancies
        file2_df = self._create_file2_dataframe(base_snps, 15194)

        # Save files
        path1 = os.path.join(output_dir, "sample1.23andme.txt.gz")
        path2 = os.path.join(output_dir, "sample2.ftdna.csv.gz")

        logger.info(f"Creating {os.path.relpath(path1)}")
        logger.info(f"Creating {os.path.relpath(path2)}")

        self._write_snps_as_23andme(file1_df, path1)
        self._write_snps_as_ftdna(file2_df, path2)

        return path1, path2

    def _generate_base_snps(self, num_snps: int) -> pd.DataFrame:
        """Generate base SNPs with sequential rsIDs."""
        base = self.generate_snps(num_snps=num_snps, inject_build_markers=False)
        base = self._sort_snps_dataframe(base).reset_index(drop=True)
        base.index = pd.Index([f"rs{i + 1}" for i in range(len(base))], name="rsid")
        return base

    def _create_file1_dataframe(
        self, base_snps: pd.DataFrame, unique_count: int
    ) -> pd.DataFrame:
        """Create file 1 DataFrame (23andMe format).

        Note: RSIDs are NOT renumbered here to preserve RSID correspondence
        with file2 for merge functionality. Unique SNPs use offset 10_000_001
        to clearly distinguish them from base SNPs (rs1-rs700000).
        """
        # Use large offset to clearly separate unique SNPs from base SNPs
        unique = self._generate_unique_snps(unique_count, 10_000_000)
        df = self._sort_snps_dataframe(pd.concat([base_snps, unique]))
        return self._inject_build_markers(df, list(REFERENCE_SEQUENCE_CHROMS))

    def _create_file2_dataframe(
        self, base_snps: pd.DataFrame, unique_count: int
    ) -> pd.DataFrame:
        """Create file 2 DataFrame (FTDNA format) with discrepancies.

        Note: RSIDs are NOT renumbered here to preserve RSID correspondence
        with file1 for merge functionality. Unique SNPs use offset 20_000_001
        to clearly distinguish them from base SNPs and file1's unique SNPs.

        This file only includes autosomal chromosomes (1-22), no X, Y, or MT.
        """
        # Only include autosomal chromosomes (1-22) for FTDNA format
        autosomal_chroms = [str(i) for i in range(1, 23)]

        # Filter base SNPs to autosomal chromosomes only
        file2_base = base_snps[base_snps["chrom"].isin(autosomal_chroms)].copy()

        # Create position discrepancies (27 SNPs)
        pos_disc_indices = self._introduce_position_discrepancies(file2_base, 27)

        # Create genotype discrepancies (151 SNPs, excluding position-discrepant ones)
        self._introduce_genotype_discrepancies(
            file2_base, 151, exclude=pos_disc_indices
        )

        # Use large offset to clearly separate unique SNPs from base SNPs
        # Generate unique SNPs only for autosomal chromosomes
        unique = self._generate_unique_snps(
            unique_count, 20_000_000, chromosomes=autosomal_chroms
        )
        df = self._sort_snps_dataframe(pd.concat([file2_base, unique]))
        return self._inject_build_markers(df, autosomal_chroms)

    def _generate_unique_snps(
        self, count: int, rsid_offset: int, chromosomes: list[str] | None = None
    ) -> pd.DataFrame:
        """Generate unique SNPs with rsIDs starting after offset."""
        unique = self.generate_snps(
            num_snps=count, inject_build_markers=False, chromosomes=chromosomes
        )
        unique = unique.reset_index(drop=True)
        unique.index = pd.Index(
            [f"rs{rsid_offset + i + 1}" for i in range(len(unique))], name="rsid"
        )
        return unique

    def _introduce_position_discrepancies(
        self, df: pd.DataFrame, count: int
    ) -> set[Any]:
        """Introduce position discrepancies in a DataFrame. Returns affected indices."""
        indices = self.rng.choice(df.index, size=count, replace=False)
        for idx in indices:
            shift = int(self.rng.integers(-1000, 1000))
            chrom = str(df.loc[idx, "chrom"])
            max_pos = self.chrom_sizes[chrom]
            new_pos = max(1, min(int(df.loc[idx, "pos"]) + shift, max_pos))  # type: ignore[arg-type]
            df.loc[idx, "pos"] = np.uint32(new_pos)
        return set(indices)

    def _introduce_genotype_discrepancies(
        self, df: pd.DataFrame, count: int, exclude: set[Any] | None = None
    ) -> None:
        """Introduce genotype discrepancies in a DataFrame, optionally excluding indices."""
        valid_indices = [
            idx for idx in df.index if exclude is None or idx not in exclude
        ]
        indices = self.rng.choice(valid_indices, size=count, replace=False)
        bases = ["A", "C", "G", "T"]
        for idx in indices:
            current = str(df.loc[idx, "genotype"])
            if current != "--":
                new_allele = self.rng.choice([b for b in bases if b not in current])
                df.loc[idx, "genotype"] = new_allele + new_allele

    def _calculate_snps_per_chromosome(
        self, num_snps: int, chromosomes: list[str]
    ) -> dict[str, int]:
        """Calculate proportional SNP distribution across chromosomes."""
        total_size = sum(self.chrom_sizes[c] for c in chromosomes)
        return {
            chrom: max(1, int(num_snps * self.chrom_sizes[chrom] / total_size))
            for chrom in chromosomes
        }

    def _generate_unique_positions(self, n: int, chrom_size: int) -> NDArray[np.uint32]:
        """Generate n unique positions within chromosome efficiently.

        Uses fast randint with deduplication for large chromosomes (where
        collisions are rare) and np.random.choice for small chromosomes
        (where collisions would be frequent).
        """
        max_pos = chrom_size - 1

        # For small chromosomes or high density, use choice (guaranteed unique)
        # Threshold: if requesting more than 1% of positions, use choice
        if n > max_pos * 0.01:
            positions = self.rng.choice(max_pos, size=n, replace=False) + 1
        else:
            # For large chromosomes, use randint with deduplication (faster)
            # Generate extra to account for potential duplicates
            oversample = int(n * 1.1) + 100
            positions = self.rng.integers(1, chrom_size, size=oversample)
            positions = np.unique(positions)

            # If we don't have enough unique positions, sample more
            while len(positions) < n:
                extra = self.rng.integers(1, chrom_size, size=n - len(positions) + 100)
                positions = np.unique(np.concatenate([positions, extra]))

            # Take exactly n positions
            positions = self.rng.choice(positions, size=n, replace=False)

        return np.sort(positions).astype(np.uint32)

    def _generate_genotypes(self, n: int, missing_rate: float) -> NDArray[np.object_]:
        """Generate n genotypes using vectorized operations."""
        bases = np.array(["A", "C", "G", "T"])

        missing_mask = self.rng.random(n) < missing_rate
        homo_mask = self.rng.random(n) < 0.7

        homo_bases = self.rng.integers(0, 4, size=n)
        het_base1 = self.rng.integers(0, 4, size=n)
        het_base2 = (het_base1 + self.rng.integers(1, 4, size=n)) % 4

        genotypes = np.empty(n, dtype=object)
        genotypes[missing_mask] = "--"

        homo_idx = ~missing_mask & homo_mask
        homo_alleles = bases[homo_bases[homo_idx]]
        genotypes[homo_idx] = np.char.add(homo_alleles, homo_alleles)

        het_idx = ~missing_mask & ~homo_mask
        allele1, allele2 = bases[het_base1[het_idx]], bases[het_base2[het_idx]]
        genotypes[het_idx] = np.array(
            ["".join(sorted([a1, a2])) for a1, a2 in zip(allele1, allele2)]
        )

        return genotypes

    def _inject_build_markers(
        self, df: pd.DataFrame, chromosomes: list[str]
    ) -> pd.DataFrame:
        """Inject known marker SNPs at their correct positions for build detection."""
        marker_rows = []
        for rsid, marker_data in BUILD_MARKER_SNPS.items():
            chrom = marker_data["chrom"]
            if chrom in chromosomes and self.build in marker_data:
                marker_rows.append(
                    {
                        "rsid": rsid,
                        "chrom": chrom,
                        "pos": np.uint32(marker_data[self.build]),
                        "genotype": self._random_genotype(),
                    }
                )

        if marker_rows:
            marker_df = pd.DataFrame(marker_rows).set_index("rsid")
            df = df[~df.index.isin(marker_df.index)]
            df = self._sort_snps_dataframe(pd.concat([df, marker_df]))

        return df

    def _renumber_rsids_sequentially(self, df: pd.DataFrame) -> pd.DataFrame:
        """Renumber RSIDs sequentially while preserving build marker RSIDs.

        After sorting by chromosome and position, RSIDs are renumbered to be
        sequential (rs1, rs2, rs3, ...) for a cleaner output. Build marker
        RSIDs are preserved to ensure build detection continues to work.
        """
        # Get the set of build marker rsids that must be preserved
        marker_rsids = set(BUILD_MARKER_SNPS.keys())

        # Create new index with sequential rsids, preserving markers
        new_index = []
        counter = 1
        for rsid in df.index:
            if rsid in marker_rsids:
                new_index.append(rsid)
            else:
                new_index.append(f"rs{counter}")
                counter += 1

        df.index = pd.Index(new_index, name="rsid")
        return df

    def _sort_snps_dataframe(self, df: pd.DataFrame) -> pd.DataFrame:
        """Sort SNPs dataframe by chromosome and position.

        Uses the same chromosome ordering as SNPs.sort(): numeric chromosomes
        in natural order (1, 2, ..., 22), then X, Y, and MT at the end.

        Parameters
        ----------
        df : pd.DataFrame
            DataFrame with 'chrom' and 'pos' columns

        Returns
        -------
        pd.DataFrame
            Sorted DataFrame with chrom as object dtype
        """
        # Get unique chromosomes and determine sort order
        unique_chroms = df["chrom"].unique()

        # Use REFERENCE_SEQUENCE_CHROMS order, filtering to only include present chroms
        chrom_order = [c for c in REFERENCE_SEQUENCE_CHROMS if c in unique_chroms]

        # Add any chromosomes not in standard list (e.g., PAR)
        for c in unique_chroms:
            if c not in chrom_order:
                chrom_order.append(c)

        # Convert to categorical for proper sorting
        df["chrom"] = df["chrom"].astype(
            CategoricalDtype(categories=chrom_order, ordered=True)
        )

        # Sort by chromosome and position
        df = df.sort_values(["chrom", "pos"])

        # Convert back to object dtype
        df["chrom"] = df["chrom"].astype(object)

        return df

    def _random_genotype(self) -> str:
        """Generate a single random genotype."""
        bases = ["A", "C", "G", "T"]
        if self.rng.random() < 0.7:
            base = self.rng.choice(bases)
            return base + base
        alleles = self.rng.choice(bases, size=2, replace=False)
        return "".join(sorted(alleles))

    # -------------------------------------------------------------------------
    # File Format Writers
    # -------------------------------------------------------------------------

    def _save_as(
        self,
        output_path: str,
        writer_method: Callable[[pd.DataFrame, str], None],
        num_snps: int,
        **kwargs: Any,
    ) -> str:
        """Generate SNPs and save using the specified writer method."""
        df = self.generate_snps(num_snps=num_snps, **kwargs)
        logger.info(f"Creating {os.path.relpath(output_path)}")
        writer_method(df, output_path)
        return output_path

    def save_as_23andme(
        self, output_path: str, num_snps: int = 991786, **kwargs: Any
    ) -> str:
        """Save SNPs in 23andMe format."""
        return self._save_as(
            output_path, self._write_snps_as_23andme, num_snps, **kwargs
        )

    def save_as_ancestry(
        self, output_path: str, num_snps: int = 700000, **kwargs: Any
    ) -> str:
        """Save SNPs in AncestryDNA format."""
        return self._save_as(
            output_path, self._write_snps_as_ancestry, num_snps, **kwargs
        )

    def save_as_ftdna(
        self, output_path: str, num_snps: int = 715194, **kwargs: Any
    ) -> str:
        """Save SNPs in Family Tree DNA (FTDNA) format."""
        return self._save_as(output_path, self._write_snps_as_ftdna, num_snps, **kwargs)

    def save_as_generic(
        self,
        output_path: str,
        format: str = "csv",
        num_snps: int = 10000,
        **kwargs: Any,
    ) -> str:
        """Save SNPs in generic CSV or TSV format."""
        df = self.generate_snps(num_snps=num_snps, **kwargs)
        logger.info(f"Creating {os.path.relpath(output_path)}")
        os.makedirs(os.path.dirname(output_path) or ".", exist_ok=True)
        sep = "," if format == "csv" else "\t"
        df_out = df.reset_index()
        if output_path.endswith(".gz"):
            with gzip.open(output_path, "wt") as f:
                df_out.to_csv(f, sep=sep, index=False)
        else:
            df_out.to_csv(output_path, sep=sep, index=False)
        return output_path

    def _write_snps_as_23andme(self, df: pd.DataFrame, output_path: str) -> None:
        """Write a DataFrame of SNPs in 23andMe format."""
        header = (
            "# 23andMe\n#\n"
            "# This is synthetic genotype data generated for demonstration purposes.\n#\n"
            "# rsid\tchromosome\tposition\tgenotype\n"
        )
        buffer = StringIO()
        buffer.write(header)
        df.reset_index().to_csv(buffer, sep="\t", header=False, index=False)
        self._write_file(output_path, buffer.getvalue())

    def _write_snps_as_ftdna(self, df: pd.DataFrame, output_path: str) -> None:
        """Write a DataFrame of SNPs in FTDNA format."""
        df_out = df.reset_index()
        df_out.columns = ["RSID", "CHROMOSOME", "POSITION", "RESULT"]
        buffer = StringIO()
        buffer.write("RSID,CHROMOSOME,POSITION,RESULT\n")
        df_out.to_csv(buffer, index=False, header=False, quoting=1)
        self._write_file(output_path, buffer.getvalue())

    def _write_snps_as_ancestry(self, df: pd.DataFrame, output_path: str) -> None:
        """Write a DataFrame of SNPs in AncestryDNA format."""
        header = (
            "#AncestryDNA\n#\n"
            "# This is synthetic genotype data generated for demonstration purposes.\n#\n"
            "rsid\tchromosome\tposition\tallele1\tallele2\n"
        )
        df_out = df.reset_index()
        genotypes = df_out["genotype"].values

        allele1 = np.where(
            genotypes == "--",
            "0",
            np.array([g[0] if len(g) >= 1 else "0" for g in genotypes]),
        )
        allele2 = np.where(
            genotypes == "--",
            "0",
            np.array([g[1] if len(g) >= 2 else "0" for g in genotypes]),
        )

        result = pd.DataFrame(
            {
                "rsid": df_out["rsid"],
                "chrom": df_out["chrom"],
                "pos": df_out["pos"],
                "allele1": allele1,
                "allele2": allele2,
            }
        )
        buffer = StringIO()
        buffer.write(header)
        result.to_csv(buffer, sep="\t", header=False, index=False)
        self._write_file(output_path, buffer.getvalue())

    def _write_file(self, path: str, content: str) -> None:
        """Write content to file atomically, handling gzip compression if needed."""
        os.makedirs(os.path.dirname(path) or ".", exist_ok=True)
        if path.endswith(".gz"):
            with atomic_write(path, mode="wb", overwrite=True) as f:
                with gzip.open(f, "wt", compresslevel=1) as f_gzip:
                    f_gzip.write(content)
        else:
            with atomic_write(path, mode="w", overwrite=True) as f:
                f.write(content)
