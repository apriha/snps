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
import tempfile
from unittest import TestCase
from unittest.mock import Mock, patch

import numpy as np
import pandas as pd

from snps import SNPs
from snps.utils import gzip_file, zip_file


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

    @property
    def downloads_enabled(self):
        """ Property indicating if downloads are enabled.

        Only download from external resources when an environment variable named
        "DOWNLOADS_ENABLED" is set to "true".

        Returns
        -------
        bool
        """
        return True if os.getenv("DOWNLOADS_ENABLED") == "true" else False

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

    def load_assign_PAR_SNPs(self, path):
        """ Load and assign PAR SNPs.

        If downloads are not enabled, use a minimal subset of the real responses.

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
        2. Yates et. al. (doi:10.1093/bioinformatics/btu613),
           `<http://europepmc.org/search/?query=DOI:10.1093/bioinformatics/btu613>`_
        3. Zerbino et. al. (doi.org/10.1093/nar/gkx1098), https://doi.org/10.1093/nar/gkx1098
        4. Sherry ST, Ward MH, Kholodov M, Baker J, Phan L, Smigielski EM, Sirotkin K.
           dbSNP: the NCBI database of genetic variation. Nucleic Acids Res. 2001 Jan 1;
           29(1):308-11.
        5. Database of Single Nucleotide Polymorphisms (dbSNP). Bethesda (MD): National Center
           for Biotechnology Information, National Library of Medicine. dbSNP accession:
           rs28736870, rs113313554, rs758419898, and rs113378274 (dbSNP Build ID: 151).
           Available from: http://www.ncbi.nlm.nih.gov/SNP/
        """
        effects = [
            {
                "refsnp_id": "758419898",
                "create_date": "2015-04-1T22:25Z",
                "last_update_date": "2019-07-14T04:19Z",
                "last_update_build_id": "153",
                "primary_snapshot_data": {
                    "placements_with_allele": [
                        {
                            "seq_id": "NC_000024.9",
                            "placement_annot": {
                                "seq_id_traits_by_assembly": [
                                    {"assembly_name": "GRCh37.p13"}
                                ]
                            },
                            "alleles": [
                                {
                                    "allele": {
                                        "spdi": {
                                            "seq_id": "NC_000024.9",
                                            "position": 7364103,
                                        }
                                    }
                                }
                            ],
                        }
                    ]
                },
            },
            {
                "refsnp_id": "28736870",
                "create_date": "2005-05-24T14:43Z",
                "last_update_date": "2019-07-14T04:18Z",
                "last_update_build_id": "153",
                "primary_snapshot_data": {
                    "placements_with_allele": [
                        {
                            "seq_id": "NC_000023.10",
                            "placement_annot": {
                                "seq_id_traits_by_assembly": [
                                    {"assembly_name": "GRCh37.p13"}
                                ]
                            },
                            "alleles": [
                                {
                                    "allele": {
                                        "spdi": {
                                            "seq_id": "NC_000023.10",
                                            "position": 220769,
                                        }
                                    }
                                }
                            ],
                        }
                    ]
                },
            },
            {
                "refsnp_id": "113313554",
                "create_date": "2010-07-4T18:13Z",
                "last_update_date": "2019-07-14T04:18Z",
                "last_update_build_id": "153",
                "primary_snapshot_data": {
                    "placements_with_allele": [
                        {
                            "seq_id": "NC_000024.9",
                            "placement_annot": {
                                "seq_id_traits_by_assembly": [
                                    {"assembly_name": "GRCh37.p13"}
                                ]
                            },
                            "alleles": [
                                {
                                    "allele": {
                                        "spdi": {
                                            "seq_id": "NC_000024.9",
                                            "position": 535257,
                                        }
                                    }
                                }
                            ],
                        }
                    ]
                },
            },
            {
                "refsnp_id": "113378274",
                "create_date": "2010-07-4T18:14Z",
                "last_update_date": "2016-03-3T10:51Z",
                "last_update_build_id": "147",
                "merged_snapshot_data": {"merged_into": ["72608386"]},
            },
            {
                "refsnp_id": "72608386",
                "create_date": "2009-02-14T01:08Z",
                "last_update_date": "2019-07-14T04:05Z",
                "last_update_build_id": "153",
                "primary_snapshot_data": {
                    "placements_with_allele": [
                        {
                            "seq_id": "NC_000023.10",
                            "placement_annot": {
                                "seq_id_traits_by_assembly": [
                                    {"assembly_name": "GRCh37.p13"}
                                ]
                            },
                            "alleles": [
                                {
                                    "allele": {
                                        "spdi": {
                                            "seq_id": "NC_000023.10",
                                            "position": 91941055,
                                        }
                                    }
                                }
                            ],
                        }
                    ]
                },
            },
        ]

        if self.downloads_enabled:
            return SNPs(path, assign_par_snps=True, deduplicate_XY_chrom=False)
        else:
            mock = Mock(side_effect=effects)
            with patch("snps.ensembl.EnsemblRestClient.perform_rest_action", mock):
                return SNPs(path, assign_par_snps=True, deduplicate_XY_chrom=False)

    def _get_test_assembly_mapping_data(self, source, target, strands, mappings):
        return {
            "1": {
                "mappings": [
                    {
                        "original": {
                            "seq_region_name": "1",
                            "strand": strands[0],
                            "start": mappings[0],
                            "end": mappings[0],
                            "assembly": "{}".format(source),
                        },
                        "mapped": {
                            "seq_region_name": "1",
                            "strand": strands[1],
                            "start": mappings[1],
                            "end": mappings[1],
                            "assembly": "{}".format(target),
                        },
                    },
                    {
                        "original": {
                            "seq_region_name": "1",
                            "strand": strands[2],
                            "start": mappings[2],
                            "end": mappings[2],
                            "assembly": "{}".format(source),
                        },
                        "mapped": {
                            "seq_region_name": "1",
                            "strand": strands[3],
                            "start": mappings[3],
                            "end": mappings[3],
                            "assembly": "{}".format(target),
                        },
                    },
                    {
                        "original": {
                            "seq_region_name": "1",
                            "strand": strands[4],
                            "start": mappings[4],
                            "end": mappings[4],
                            "assembly": "{}".format(source),
                        },
                        "mapped": {
                            "seq_region_name": "1",
                            "strand": strands[5],
                            "start": mappings[5],
                            "end": mappings[5],
                            "assembly": "{}".format(target),
                        },
                    },
                ]
            },
            "3": {
                "mappings": [
                    {
                        "original": {
                            "seq_region_name": "3",
                            "strand": strands[6],
                            "start": mappings[6],
                            "end": mappings[6],
                            "assembly": "{}".format(source),
                        },
                        "mapped": {
                            "seq_region_name": "3",
                            "strand": strands[7],
                            "start": mappings[7],
                            "end": mappings[7],
                            "assembly": "{}".format(target),
                        },
                    }
                ]
            },
        }

    def NCBI36_GRCh37(self):
        return self._get_test_assembly_mapping_data(
            "NCBI36",
            "GRCh37",
            [1, 1, 1, 1, 1, 1, 1, -1],
            [
                742429,
                752566,
                143649677,
                144938320,
                143649678,
                144938321,
                50908372,
                50927009,
            ],
        )

    def GRCh37_NCBI36(self):
        return self._get_test_assembly_mapping_data(
            "GRCh37",
            "NCBI36",
            [1, 1, 1, 1, 1, 1, 1, -1],
            [
                752566,
                742429,
                144938320,
                143649677,
                144938321,
                143649678,
                50927009,
                50908372,
            ],
        )

    def GRCh37_GRCh38(self):
        return self._get_test_assembly_mapping_data(
            "GRCh37",
            "GRCh38",
            [1, 1, 1, -1, 1, -1, 1, 1],
            [
                752566,
                817186,
                144938320,
                148946169,
                144938321,
                148946168,
                50927009,
                50889578,
            ],
        )

    def GRCh37_GRCh38_PAR(self):
        return {
            "X": {
                "mappings": [
                    {
                        "original": {
                            "seq_region_name": "X",
                            "strand": 1,
                            "start": 220770,
                            "end": 220770,
                            "assembly": "GRCh37",
                        },
                        "mapped": {
                            "seq_region_name": "X",
                            "strand": 1,
                            "start": 304103,
                            "end": 304103,
                            "assembly": "GRCh38",
                        },
                    },
                    {
                        "original": {
                            "seq_region_name": "X",
                            "strand": 1,
                            "start": 91941056,
                            "end": 91941056,
                            "assembly": "GRCh37",
                        },
                        "mapped": {
                            "seq_region_name": "X",
                            "strand": 1,
                            "start": 92686057,
                            "end": 92686057,
                            "assembly": "GRCh38",
                        },
                    },
                ]
            },
            "Y": {
                "mappings": [
                    {
                        "original": {
                            "seq_region_name": "Y",
                            "strand": 1,
                            "start": 535258,
                            "end": 535258,
                            "assembly": "GRCh37",
                        },
                        "mapped": {
                            "seq_region_name": "Y",
                            "strand": 1,
                            "start": 624523,
                            "end": 624523,
                            "assembly": "GRCh38",
                        },
                    }
                ]
            },
        }

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
            dest = os.path.join(tmpdir, "{}.gz".format(base))
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

            dest = os.path.join(tmpdir, "{}.zip".format(base))
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
    ):
        # https://samtools.github.io/hts-specs/VCFv4.2.pdf
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
        )
        self.make_parsing_assertions_vcf(
            self.parse_bytes(file, rsids),
            source,
            phased,
            unannotated,
            rsids,
            build,
            build_detected,
        )

        with tempfile.TemporaryDirectory() as tmpdir:
            base = os.path.basename(file)
            dest = os.path.join(tmpdir, "{}.gz".format(base))
            gzip_file(file, dest)
            self.make_parsing_assertions_vcf(
                self.parse_file(dest, rsids),
                source,
                phased,
                unannotated,
                rsids,
                build,
                build_detected,
            )
            self.make_parsing_assertions_vcf(
                self.parse_bytes(dest, rsids),
                source,
                phased,
                unannotated,
                rsids,
                build,
                build_detected,
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

        self.assertEqual(snps.source, source)
        pd.testing.assert_frame_equal(snps.snps, snps_df, check_exact=True)
        self.assertTrue(snps.phased) if phased else self.assertFalse(snps.phased)
        self.assertEqual(snps.build, build)
        self.assertTrue(snps.build_detected) if build_detected else self.assertFalse(
            snps.build_detected
        )

    def make_parsing_assertions_vcf(
        self, snps, source, phased, unannotated, rsids, build, build_detected
    ):
        self.assertEqual(snps.source, source)

        if unannotated:
            self.assertTrue(snps.unannotated_vcf)
            self.assertEqual(0, snps.snp_count)
        else:
            self.assertFalse(snps.unannotated_vcf)
            pd.testing.assert_frame_equal(
                snps.snps, self.generic_snps_vcf().loc[rsids], check_exact=True
            ) if rsids else pd.testing.assert_frame_equal(
                snps.snps, self.generic_snps_vcf(), check_exact=True
            )

        self.assertTrue(snps.phased) if phased else self.assertFalse(snps.phased)
        self.assertEqual(snps.build, build)
        self.assertTrue(snps.build_detected) if build_detected else self.assertFalse(
            snps.build_detected
        )
