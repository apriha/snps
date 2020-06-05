""" Class for writing SNPs.

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

import datetime
import logging

import numpy as np
import pandas as pd

import snps
from snps.utils import save_df_as_csv, clean_str

logger = logging.getLogger(__name__)


class Writer:
    """ Class for writing SNPs to files. """

    def __init__(self, snps=None, filename="", vcf=False, atomic=True, **kwargs):
        """ Initialize a `Writer`.

        Parameters
        ----------
        snps : SNPs
            SNPs to save to file or write to buffer
        filename : str or buffer
            filename for file to save or buffer to write to
        vcf : bool
            flag to save file as VCF
        atomic : bool
            atomically write output to a file on local filesystem
        **kwargs
            additional parameters to `pandas.DataFrame.to_csv`
        """
        self._snps = snps
        self._filename = filename
        self._vcf = vcf
        self._atomic = atomic
        self._kwargs = kwargs

    def __call__(self):
        if self._vcf:
            return self._write_vcf()
        else:
            return self._write_csv()

    @classmethod
    def write_file(cls, snps=None, filename="", vcf=False, atomic=True, **kwargs):
        """ Save SNPs to file.

        Parameters
        ----------
        snps : SNPs
            SNPs to save to file or write to buffer
        filename : str or buffer
            filename for file to save or buffer to write to
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
        w = cls(snps=snps, filename=filename, vcf=vcf, atomic=atomic, **kwargs)
        return w()

    def _write_csv(self):
        """ Write SNPs to a CSV file.

        Returns
        -------
        str
            path to file in output directory if SNPs were saved, else empty str
        """
        filename = self._filename
        if not filename:
            filename = "{}_{}{}".format(
                clean_str(self._snps._source), self._snps.assembly, ".csv"
            )

        comment = (
            "# Source(s): {}\n"
            "# Assembly: {}\n"
            "# Phased: {}\n"
            "# SNPs: {}\n"
            "# Chromosomes: {}\n".format(
                self._snps.source,
                self._snps.assembly,
                self._snps.phased,
                self._snps.snp_count,
                self._snps.chromosomes_summary,
            )
        )
        if "header" in self._kwargs:
            if isinstance(self._kwargs["header"], bool):
                if self._kwargs["header"]:
                    self._kwargs["header"] = ["chromosome", "position", "genotype"]
        else:
            self._kwargs["header"] = ["chromosome", "position", "genotype"]

        return save_df_as_csv(
            self._snps._snps,
            self._snps._output_dir,
            filename,
            comment=comment,
            atomic=self._atomic,
            **self._kwargs
        )

    def _write_vcf(self):
        """ Write SNPs to a VCF file.

        References
        ----------
        1. The Variant Call Format (VCF) Version 4.2 Specification, 8 Mar 2019,
           https://samtools.github.io/hts-specs/VCFv4.2.pdf

        Returns
        -------
        str
            path to file in output directory if SNPs were saved, else empty str
        """
        filename = self._filename
        if not filename:
            filename = "{}_{}{}".format(
                clean_str(self._snps._source), self._snps.assembly, ".vcf"
            )

        comment = (
            "##fileformat=VCFv4.2\n"
            "##fileDate={}\n"
            '##source="{}; snps v{}; https://pypi.org/project/snps/"\n'.format(
                datetime.datetime.utcnow().strftime("%Y%m%d"),
                self._snps._source,
                snps.__version__,
            )
        )

        reference_sequence_chroms = (
            "1",
            "2",
            "3",
            "4",
            "5",
            "6",
            "7",
            "8",
            "9",
            "10",
            "11",
            "12",
            "13",
            "14",
            "15",
            "16",
            "17",
            "18",
            "19",
            "20",
            "21",
            "22",
            "X",
            "Y",
            "MT",
        )

        df = self._snps.snps

        tasks = []

        # skip insertions and deletions
        df = df.drop(
            df.loc[
                df["genotype"].notnull()
                & (
                    (df["genotype"].str[0] == "I")
                    | (df["genotype"].str[0] == "D")
                    | (df["genotype"].str[1] == "I")
                    | (df["genotype"].str[1] == "D")
                )
            ].index
        )

        chroms_to_drop = []
        for chrom in df["chrom"].unique():
            if chrom not in reference_sequence_chroms:
                chroms_to_drop.append(chrom)
                continue

            tasks.append(
                {
                    "resources": self._snps._resources,
                    "assembly": self._snps.assembly,
                    "chrom": chrom,
                    "snps": pd.DataFrame(df.loc[(df["chrom"] == chrom)]),
                }
            )

        # drop chromosomes without reference sequence data (e.g., unassigned PAR)
        for chrom in chroms_to_drop:
            df = df.drop(df.loc[df["chrom"] == chrom].index)

        # create the VCF representation for SNPs
        results = map(self._create_vcf_representation, tasks)

        contigs = []
        vcf = [pd.DataFrame()]
        for result in list(results):
            contigs.append(result["contig"])
            vcf.append(result["vcf"])

        vcf = pd.concat(vcf)

        comment += "".join(contigs)
        comment += '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n'
        comment += "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE\n"

        return save_df_as_csv(
            vcf,
            self._snps._output_dir,
            filename,
            comment=comment,
            prepend_info=False,
            header=False,
            index=False,
            na_rep=".",
            sep="\t",
        )

    def _create_vcf_representation(self, task):
        resources = task["resources"]
        assembly = task["assembly"]
        chrom = task["chrom"]
        snps = task["snps"]

        if len(snps.loc[snps["genotype"].notnull()]) == 0:
            return {"contig": "", "vcf": pd.DataFrame()}

        seqs = resources.get_reference_sequences(assembly, [chrom])
        seq = seqs[chrom]

        contig = '##contig=<ID={},URL={},length={},assembly={},md5={},species="{}">\n'.format(
            seq.ID, seq.url, seq.length, seq.build, seq.md5, seq.species
        )

        snps = snps.reset_index()

        df = pd.DataFrame(
            columns=[
                "CHROM",
                "POS",
                "ID",
                "REF",
                "ALT",
                "QUAL",
                "FILTER",
                "INFO",
                "FORMAT",
                "SAMPLE",
            ]
        )
        df = df.astype(
            {
                "CHROM": object,
                "POS": np.int64,
                "ID": object,
                "REF": object,
                "ALT": object,
                "QUAL": np.int64,
                "FILTER": object,
                "INFO": object,
                "FORMAT": object,
                "SAMPLE": object,
            }
        )

        df["CHROM"] = snps["chrom"]
        df["POS"] = snps["pos"]
        df["ID"] = snps["rsid"]

        # https://stackoverflow.com/a/24838429
        df["REF"] = list(map(chr, seq.sequence[snps.pos - seq.start]))

        df["FORMAT"] = "GT"

        seq.clear()

        df["genotype"] = snps["genotype"]

        temp = df.loc[df["genotype"].notnull()]

        # https://stackoverflow.com/a/19976286
        df.loc[df["genotype"].notnull(), "ALT"] = np.vectorize(self._compute_alt)(
            temp["REF"], temp["genotype"]
        )

        temp = df.loc[df["genotype"].notnull()]

        df.loc[df["genotype"].notnull(), "SAMPLE"] = np.vectorize(
            self._compute_genotype
        )(temp["REF"], temp["ALT"], temp["genotype"])

        df.loc[df["SAMPLE"].isnull(), "SAMPLE"] = "./."

        del df["genotype"]

        return {"contig": contig, "vcf": df}

    def _compute_alt(self, ref, genotype):
        genotype_alleles = list(set(genotype))

        if ref in genotype_alleles:
            if len(genotype_alleles) == 1:
                return "N"
            else:
                genotype_alleles.remove(ref)
                return genotype_alleles.pop(0)
        else:
            return ",".join(genotype_alleles)

    def _compute_genotype(self, ref, alt, genotype):
        alleles = [ref]

        if self._snps.phased:
            separator = "|"
        else:
            separator = "/"

        if pd.notna(alt):
            alleles.extend(alt.split(","))

        if len(genotype) == 2:
            return "{}{}{}".format(
                alleles.index(genotype[0]), separator, alleles.index(genotype[1])
            )
        else:
            return "{}".format(alleles.index(genotype[0]))
