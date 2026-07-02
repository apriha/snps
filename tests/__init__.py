import os
import shutil
import tempfile
from unittest import TestCase

import numpy as np
import pandas as pd
from pandas.api.types import is_object_dtype, is_string_dtype, is_unsigned_integer_dtype

from snps import SNPs
from snps.resources import set_default_provider
from snps.testing import SNPsTestMixin, create_simulated_snp_df
from snps.utils import gzip_file, zip_file
from tests.support import FakeResources, chip_clusters_df, low_quality_snps_df
from tests.support import data as _data


class BaseSNPsTestCase(SNPsTestMixin, TestCase):
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
        s._snps = create_simulated_snp_df(
            chrom=chrom,
            pos_start=pos_start,
            pos_max=pos_max,
            pos_step=pos_step,
            pos_dtype=np.uint32,
            genotype=genotype,
            insert_nulls=insert_nulls,
            null_snp_step=null_snp_step,
            complement_genotype_one_allele=complement_genotype_one_chrom,
            complement_genotype_two_alleles=complement_genotype_two_chroms,
            complement_snp_step=complement_snp_step,
        )
        return s

    def load_assign_PAR_SNPs(self, path):
        """Load and assign PAR SNPs offline via the fixture-backed resource provider.

        Parameters
        ----------
        path : str

        Returns
        -------
        SNPs

        References
        ----------
        1. National Center for Biotechnology Information, Variation Services, RefSNP,
           https://api.ncbi.nlm.nih.gov/variation/v0/
        4. Sherry ST, Ward MH, Kholodov M, Baker J, Phan L, Smigielski EM, Sirotkin K.
           dbSNP: the NCBI database of genetic variation. Nucleic Acids Res. 2001 Jan 1;
           29(1):308-11.
        5. Database of Single Nucleotide Polymorphisms (dbSNP). Bethesda (MD): National Center
           for Biotechnology Information, National Library of Medicine. dbSNP accession:
           rs28736870, rs113313554, rs758419898, and rs113378274 (dbSNP Build ID: 151).
           Available from: http://www.ncbi.nlm.nih.gov/SNP/
        """
        # PAR lookups are served by the active (fixture-backed) default provider; do not
        # replace it here so callers can inject specific data (e.g., remap mappings).
        return SNPs(path, assign_par_snps=True, deduplicate_XY_chrom=False)

    def _get_test_assembly_mapping_data(self, source, target, strands, mappings):
        return _data.get_test_assembly_mapping_data(source, target, strands, mappings)

    def NCBI36_GRCh37(self):
        return _data.NCBI36_GRCh37()

    def GRCh37_NCBI36(self):
        return _data.GRCh37_NCBI36()

    def GRCh37_GRCh38(self):
        return _data.GRCh37_GRCh38()

    def GRCh37_GRCh38_PAR(self):
        return _data.GRCh37_GRCh38_PAR()

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

    def snps_GRCh38(self):
        return self.create_snp_df(
            rsid=["rs3094315", "rsIndelTest", "rs2500347", "rs11928389"],
            chrom=["1", "1", "1", "3"],
            pos=[817186, 148946168, 148946169, 50889578],
            genotype=["AA", "ID", np.nan, "TC"],
        )

    def snps_GRCh37_PAR(self):
        return self.create_snp_df(
            rsid=["rs28736870", "rs113378274", "rs113313554", "rs758419898"],
            chrom=["X", "X", "Y", "PAR"],
            pos=[220770, 91941056, 535258, 1],
            genotype=["AA", "AA", "AA", "AA"],
        )

    def snps_GRCh38_PAR(self):
        return self.create_snp_df(
            rsid=["rs28736870", "rs113378274", "rs113313554"],
            chrom=["X", "X", "Y"],
            pos=[304103, 92686057, 624523],
            genotype=["AA", "AA", "AA"],
        )

    def generic_snps_vcf(self):
        df = self.generic_snps()
        return pd.concat(
            [
                df,
                self.create_snp_df(
                    rsid=["rs" + str(i) for i in range(12, 18)],
                    chrom=["1"] * 6,
                    pos=list(range(112, 118)),
                    genotype=[np.nan] * 6,
                ),
            ]
        )

    def run_parsing_tests(
        self, file, source, phased=False, build=37, build_detected=False, snps_df=None
    ):
        self.make_parsing_assertions(
            self.parse_file(file), source, phased, build, build_detected, snps_df
        )
        self.make_parsing_assertions(
            self.parse_bytes(file), source, phased, build, build_detected, snps_df
        )

        with tempfile.TemporaryDirectory() as tmpdir:
            base = os.path.basename(file)
            dest = os.path.join(tmpdir, f"{base}.gz")
            gzip_file(file, dest)
            self.make_parsing_assertions(
                self.parse_file(dest), source, phased, build, build_detected, snps_df
            )
            self.make_parsing_assertions(
                self.parse_bytes(dest), source, phased, build, build_detected, snps_df
            )
            # remove .gz extension
            shutil.move(dest, dest[:-3])
            self.make_parsing_assertions(
                self.parse_file(dest[:-3]),
                source,
                phased,
                build,
                build_detected,
                snps_df,
            )

            dest = os.path.join(tmpdir, f"{base}.zip")
            zip_file(file, dest, base)
            self.make_parsing_assertions(
                self.parse_file(dest), source, phased, build, build_detected, snps_df
            )
            self.make_parsing_assertions(
                self.parse_bytes(dest), source, phased, build, build_detected, snps_df
            )
            # remove .zip extension
            shutil.move(dest, dest[:-4])
            self.make_parsing_assertions(
                self.parse_file(dest[:-4]),
                source,
                phased,
                build,
                build_detected,
                snps_df,
            )

    def run_parsing_tests_vcf(
        self,
        file,
        source="vcf",
        phased=False,
        unannotated=False,
        rsids=(),
        build=37,
        build_detected=False,
        snps_df=None,
    ):
        # https://samtools.github.io/hts-specs/VCFv4.3.pdf
        # this tests for homozygous snps, heterozygous snps, multiallelic snps,
        # phased snps, and snps with missing rsID
        self.make_parsing_assertions_vcf(
            self.parse_file(file, rsids),
            source,
            phased,
            unannotated,
            rsids,
            build,
            build_detected,
            snps_df,
        )
        self.make_parsing_assertions_vcf(
            self.parse_bytes(file, rsids),
            source,
            phased,
            unannotated,
            rsids,
            build,
            build_detected,
            snps_df,
        )

        with tempfile.TemporaryDirectory() as tmpdir:
            base = os.path.basename(file)
            dest = os.path.join(tmpdir, f"{base}.gz")
            gzip_file(file, dest)
            self.make_parsing_assertions_vcf(
                self.parse_file(dest, rsids),
                source,
                phased,
                unannotated,
                rsids,
                build,
                build_detected,
                snps_df,
            )
            self.make_parsing_assertions_vcf(
                self.parse_bytes(dest, rsids),
                source,
                phased,
                unannotated,
                rsids,
                build,
                build_detected,
                snps_df,
            )
            # remove .gz extension
            shutil.move(dest, dest[:-3])
            self.make_parsing_assertions_vcf(
                self.parse_file(dest[:-3], rsids),
                source,
                phased,
                unannotated,
                rsids,
                build,
                build_detected,
                snps_df,
            )

    def make_normalized_dataframe_assertions(self, df):
        self.assertEqual(df.index.name, "rsid")
        # Accept both object dtype and StringDtype (used in Python 3.14+)
        self.assertTrue(is_string_dtype(df.index.dtype))
        self.assertTrue(
            is_string_dtype(df.chrom.dtype) or is_object_dtype(df.chrom.dtype)
        )
        self.assertTrue(is_unsigned_integer_dtype(df.pos.dtype))
        self.assertTrue(
            is_string_dtype(df.genotype.dtype) or is_object_dtype(df.genotype.dtype)
        )

    def parse_file(self, file, rsids=()):
        return SNPs(file, rsids=rsids)

    def parse_bytes(self, file, rsids=()):
        with open(file, "rb") as f:
            return SNPs(f.read(), rsids=rsids)

    def make_parsing_assertions(
        self, snps, source, phased, build, build_detected, snps_df
    ):
        if snps_df is None:
            snps_df = self.generic_snps()

        # these are useful for debugging if there is a problem
        print("Observed:")
        print(snps.snps)
        print(snps.snps.info())
        print("Expected:")
        print(snps_df)
        print(snps_df.info())

        self.assertEqual(snps.source, source)
        self.assert_frame_equal_with_string_index(snps.snps, snps_df, check_exact=True)
        self.assertTrue(snps.phased) if phased else self.assertFalse(snps.phased)
        self.assertEqual(snps.build, build)
        (
            self.assertTrue(snps.build_detected)
            if build_detected
            else self.assertFalse(snps.build_detected)
        )
        self.make_normalized_dataframe_assertions(snps.snps)

    def make_parsing_assertions_vcf(
        self, snps, source, phased, unannotated, rsids, build, build_detected, snps_df
    ):
        if snps_df is None:
            snps_df = self.generic_snps_vcf()

        self.assertEqual(snps.source, source)

        if unannotated:
            self.assertTrue(snps.unannotated_vcf)
            self.assertEqual(0, snps.count)
        else:
            self.assertFalse(snps.unannotated_vcf)
            (
                self.assert_frame_equal_with_string_index(
                    snps.snps, snps_df.loc[rsids], check_exact=True
                )
                if rsids
                else self.assert_frame_equal_with_string_index(
                    snps.snps, snps_df, check_exact=True
                )
            )

        self.assertTrue(snps.phased) if phased else self.assertFalse(snps.phased)
        self.assertEqual(snps.build, build)
        (
            self.assertTrue(snps.build_detected)
            if build_detected
            else self.assertFalse(snps.build_detected)
        )
        self.make_normalized_dataframe_assertions(snps.snps)

    def get_low_quality_snps(self, pos=(104, 106, 1001), cluster="c1"):
        return low_quality_snps_df(pos=pos, cluster=cluster)

    def run_low_quality_snps_test(self, f, low_quality_snps, cluster="c1"):
        # Cluster detection runs the real overlap computation (no mocks) against these
        # chip clusters. The generic test SNPs are rs1-rs8 at positions 101-108.
        if cluster:
            # clusters overlap the SNPs, so `cluster` is detected
            chip_clusters = chip_clusters_df(tuple(range(101, 109)), cluster, 8)
        else:
            # clusters do not overlap the SNPs, so no cluster is detected
            chip_clusters = chip_clusters_df(tuple(range(1001, 1009)), "c1", 8)
        set_default_provider(
            FakeResources(
                chip_clusters=chip_clusters, low_quality_snps=low_quality_snps
            )
        )
        f()
