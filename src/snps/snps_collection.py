""" `SNPsCollection` extends `SNPs` to merge genotype / raw data files.

"""

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

import logging
import os

import numpy as np
import pandas as pd

from snps.snps import SNPs
from snps.utils import save_df_as_csv, clean_str

logger = logging.getLogger(__name__)


class SNPsCollection(SNPs):
    def __init__(self, raw_data=None, output_dir="output", name="", **kwargs):
        """ Object used to merge genotype / raw data files.

        Parameters
        ----------
        raw_data : list or str
            path(s) to file(s) with raw genotype data
        output_dir : str
            path to output directory
        name : str
            name for this ``SNPsCollection``
        """
        super().__init__(file="", output_dir=output_dir, **kwargs)

        self._source = []
        self._discrepant_positions_file_count = 0
        self._discrepant_genotypes_file_count = 0
        self._discrepant_positions = pd.DataFrame()
        self._discrepant_genotypes = pd.DataFrame()
        self._name = name

        if raw_data is not None:
            self.load_snps(raw_data)

    def __repr__(self):
        return "SNPsCollection(name={!r})".format(self._name)

    @property
    def source(self):
        """ Summary of the SNP data source for ``SNPs``.

        Returns
        -------
        str
        """
        return ", ".join(self._source)

    @property
    def discrepant_positions(self):
        """ SNPs with discrepant positions discovered while loading SNPs.

        Returns
        -------
        pandas.DataFrame
        """
        return self._discrepant_positions

    @property
    def discrepant_genotypes(self):
        """ SNPs with discrepant genotypes discovered while loading SNPs.

        Returns
        -------
        pandas.DataFrame
        """
        return self._discrepant_genotypes

    @property
    def discrepant_snps(self):
        """ SNPs with discrepant positions and / or genotypes discovered while loading SNPs.

        Returns
        -------
        pandas.DataFrame
        """
        df = self._discrepant_positions.append(self._discrepant_genotypes)
        if len(df) > 1:
            df = df.drop_duplicates()
        return df

    def load_snps(
        self,
        raw_data,
        discrepant_snp_positions_threshold=100,
        discrepant_genotypes_threshold=500,
        save_output=False,
    ):
        """ Load raw genotype data.

        Parameters
        ----------
        raw_data : list or str
            path(s) to file(s) with raw genotype data
        discrepant_snp_positions_threshold : int
            threshold for discrepant SNP positions between existing data and data to be loaded,
            a large value could indicate mismatched genome assemblies
        discrepant_genotypes_threshold : int
            threshold for discrepant genotype data between existing data and data to be loaded,
            a large value could indicated mismatched individuals
        save_output : bool
            specifies whether to save discrepant SNP output to CSV files in the output directory
        """
        if type(raw_data) is list:
            for file in raw_data:
                self._load_snps_helper(
                    file,
                    discrepant_snp_positions_threshold,
                    discrepant_genotypes_threshold,
                    save_output,
                )
        elif type(raw_data) is str:
            self._load_snps_helper(
                raw_data,
                discrepant_snp_positions_threshold,
                discrepant_genotypes_threshold,
                save_output,
            )
        else:
            raise TypeError("invalid filetype")

    def _load_snps_helper(
        self,
        file,
        discrepant_snp_positions_threshold,
        discrepant_genotypes_threshold,
        save_output,
    ):
        logger.info("Loading {}".format(os.path.relpath(file)))
        discrepant_positions, discrepant_genotypes = self._add_snps(
            SNPs(file),
            discrepant_snp_positions_threshold,
            discrepant_genotypes_threshold,
            save_output,
        )

        self._discrepant_positions = self._discrepant_positions.append(
            discrepant_positions, sort=True
        )
        self._discrepant_genotypes = self._discrepant_genotypes.append(
            discrepant_genotypes, sort=True
        )

    def save_snps(self, filename="", vcf=False, atomic=True, **kwargs):
        """ Save SNPs to file.

        Parameters
        ----------
        filename : str
            filename for file to save
        vcf : bool
            flag to save file as VCF
        atomic : bool
            atomically write output to a file on local filesystem
        **kwargs
            additional parameters to `pandas.DataFrame.to_csv`

        Returns
        -------
        str
            path to file in output directory if SNPs were saved, else empty str
        """
        if not self._name:
            prefix = ""
        else:
            prefix = "{}_".format(clean_str(self._name))

        if not filename:
            if vcf:
                ext = ".vcf"
            else:
                ext = ".csv"

            filename = "{}{}{}".format(prefix, self.assembly, ext)
        return super().save_snps(filename=filename, vcf=vcf, atomic=atomic, **kwargs)

    def save_discrepant_positions(self, filename=""):
        """ Save SNPs with discrepant positions to file.

        Parameters
        ----------
        filename : str
            filename for file to save

        Returns
        -------
        str
            path to file in output directory if SNPs were saved, else empty str
        """
        return self._save_discrepant_snps_file(
            self.discrepant_positions, "discrepant_positions", filename
        )

    def save_discrepant_genotypes(self, filename=""):
        """ Save SNPs with discrepant genotypes to file.

        Parameters
        ----------
        filename : str
            filename for file to save

        Returns
        -------
        str
            path to file in output directory if SNPs were saved, else empty str
        """
        return self._save_discrepant_snps_file(
            self.discrepant_genotypes, "discrepant_genotypes", filename
        )

    def save_discrepant_snps(self, filename=""):
        """ Save SNPs with discrepant positions and / or genotypes to file.

        Parameters
        ----------
        filename : str
            filename for file to save

        Returns
        -------
        str
            path to file in output directory if SNPs were saved, else empty str
        """
        return self._save_discrepant_snps_file(
            self.discrepant_snps, "discrepant_snps", filename
        )

    def _save_discrepant_snps_file(self, df, discrepant_snps_type, filename):
        if not filename:
            if not self._name:
                filename = "{}.csv".format(discrepant_snps_type)
            else:
                filename = "{}_{}.csv".format(
                    clean_str(self._name), discrepant_snps_type
                )

        return save_df_as_csv(
            df,
            self._output_dir,
            filename,
            comment="# Source(s): {}\n".format(self.source),
        )

    def _add_snps(
        self,
        snps,
        discrepant_snp_positions_threshold,
        discrepant_genotypes_threshold,
        save_output,
    ):
        """ Add SNPs to this ``SNPsCollection``.

        Parameters
        ----------
        snps : SNPs
            SNPs to add
        discrepant_snp_positions_threshold : int
            see above
        discrepant_genotypes_threshold : int
            see above
        save_output
            see above

        Returns
        -------
        discrepant_positions : pandas.DataFrame
        discrepant_genotypes : pandas.DataFrame
        """
        discrepant_positions = pd.DataFrame()
        discrepant_genotypes = pd.DataFrame()

        if snps._snps.empty:
            return discrepant_positions, discrepant_genotypes

        build = snps._build
        source = [s.strip() for s in snps._source.split(",")]

        if not snps._build_detected:
            logger.warning("build not detected, assuming build {}".format(snps._build))

        if not self._build:
            self._build = build
        elif self._build != build:
            logger.warning(
                "build / assembly mismatch between current build of SNPs and SNPs being loaded"
            )

        # ensure there are always two X alleles
        snps = self._double_single_alleles(snps._snps, "X")

        if self._snps.empty:
            self._source.extend(source)
            self._snps = snps
        else:
            common_snps = self._snps.join(snps, how="inner", rsuffix="_added")

            discrepant_positions = common_snps.loc[
                (common_snps["chrom"] != common_snps["chrom_added"])
                | (common_snps["pos"] != common_snps["pos_added"])
            ]

            if not self._name:
                prefix = ""
            else:
                prefix = "{}_".format(clean_str(self._name))

            if 0 < len(discrepant_positions) < discrepant_snp_positions_threshold:
                logger.warning(
                    "{} SNP positions were discrepant; keeping original positions".format(
                        str(len(discrepant_positions))
                    )
                )

                if save_output:
                    self._discrepant_positions_file_count += 1
                    save_df_as_csv(
                        discrepant_positions,
                        self._output_dir,
                        "{}discrepant_positions_{}{}".format(
                            prefix, str(self._discrepant_positions_file_count), ".csv"
                        ),
                    )
            elif len(discrepant_positions) >= discrepant_snp_positions_threshold:
                logger.warning(
                    "too many SNPs differ in position; ensure same genome build is being used"
                )
                return discrepant_positions, discrepant_genotypes

            # remove null genotypes
            common_snps = common_snps.loc[
                ~common_snps["genotype"].isnull()
                & ~common_snps["genotype_added"].isnull()
            ]

            # discrepant genotypes are where alleles are not equivalent (i.e., alleles are not the
            # same and not swapped)
            discrepant_genotypes = common_snps.loc[
                (
                    (common_snps["genotype"].str.len() == 1)
                    & (common_snps["genotype_added"].str.len() == 1)
                    & ~(
                        common_snps["genotype"].str[0]
                        == common_snps["genotype_added"].str[0]
                    )
                )
                | (
                    (common_snps["genotype"].str.len() == 2)
                    & (common_snps["genotype_added"].str.len() == 2)
                    & ~(
                        (
                            common_snps["genotype"].str[0]
                            == common_snps["genotype_added"].str[0]
                        )
                        & (
                            common_snps["genotype"].str[1]
                            == common_snps["genotype_added"].str[1]
                        )
                    )
                    & ~(
                        (
                            common_snps["genotype"].str[0]
                            == common_snps["genotype_added"].str[1]
                        )
                        & (
                            common_snps["genotype"].str[1]
                            == common_snps["genotype_added"].str[0]
                        )
                    )
                )
            ]

            if 0 < len(discrepant_genotypes) < discrepant_genotypes_threshold:
                logger.warning(
                    "{} SNP genotypes were discrepant; marking those as null".format(
                        str(len(discrepant_genotypes))
                    )
                )

                if save_output:
                    self._discrepant_genotypes_file_count += 1
                    save_df_as_csv(
                        discrepant_genotypes,
                        self._output_dir,
                        "{}discrepant_genotypes_{}{}".format(
                            prefix, str(self._discrepant_genotypes_file_count), ".csv"
                        ),
                    )
            elif len(discrepant_genotypes) >= discrepant_genotypes_threshold:
                logger.warning(
                    "too many SNPs differ in their genotype; ensure file is for same "
                    "individual"
                )
                return discrepant_positions, discrepant_genotypes

            # add new SNPs
            self._source.extend(source)
            self._snps = self._snps.combine_first(snps)
            self._snps.loc[discrepant_genotypes.index, "genotype"] = np.nan

            # combine_first converts position to float64, so convert it back to int64
            self._snps["pos"] = self._snps["pos"].astype(np.int64)

        self.sort_snps()

        return discrepant_positions, discrepant_genotypes

    @staticmethod
    def _double_single_alleles(df, chrom):
        """ Double any single alleles in the specified chromosome.

        Parameters
        ----------
        df : pandas.DataFrame
            SNPs
        chrom : str
            chromosome of alleles to double

        Returns
        -------
        df : pandas.DataFrame
            SNPs with specified chromosome's single alleles doubled
        """
        # find all single alleles of the specified chromosome
        single_alleles = np.where(
            (df["chrom"] == chrom) & (df["genotype"].str.len() == 1)
        )[0]

        # double those alleles
        df.iloc[single_alleles, 2] = df.iloc[single_alleles, 2] * 2

        return df
