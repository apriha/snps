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
        discrepant_vcf_position : pd.DataFrame
            SNPs with discrepant positions discovered while saving VCF
        """

        if vcf:
            return cls._write_vcf(snps, filename, atomic)
        else:
            return (cls._write_csv(snps, filename, atomic, kwargs),)

    @classmethod
    def _write_csv(cls, s, filename, atomic, kwargs):
        """ Write SNPs to a CSV file.

        Returns
        -------
        str
            path to file in output directory if SNPs were saved, else empty str
        """
        if not filename:
            ext = ".txt"

            if "sep" in kwargs and kwargs["sep"] == ",":
                ext = ".csv"

            filename = f"{clean_str(s.source)}_{s.assembly}{ext}"

        comment = (
            f"# Source(s): {s.source}\n"
            f"# Build: {s.build}\n"
            f"# Build Detected: {s.build_detected}\n"
            f"# Phased: {s.phased}\n"
            f"# SNPs: {s.count}\n"
            f"# Chromosomes: {s.chromosomes_summary}\n"
        )
        if "header" in kwargs:
            if isinstance(kwargs["header"], bool):
                if kwargs["header"]:
                    kwargs["header"] = ["chromosome", "position", "genotype"]
        else:
            kwargs["header"] = ["chromosome", "position", "genotype"]

        return save_df_as_csv(
            s.snps, s._output_dir, filename, comment=comment, atomic=atomic, **kwargs,
        )

    @classmethod
    def _write_vcf(cls, s, filename, atomic):
        """ Write SNPs to a VCF file.

        References
        ----------
        1. The Variant Call Format (VCF) Version 4.2 Specification, 8 Mar 2019,
           https://samtools.github.io/hts-specs/VCFv4.2.pdf

        Returns
        -------
        str
            path to file in output directory if SNPs were saved, else empty str
        discrepant_vcf_position : pd.DataFrame
            SNPs with discrepant positions discovered while saving VCF
        """
        if not filename:
            filename = f"{clean_str(s.source)}_{s.assembly}{'.vcf'}"

        comment = (
            f"##fileformat=VCFv4.2\n"
            f'##fileDate={datetime.datetime.utcnow().strftime("%Y%m%d")}\n'
            f'##source="{s.source}; snps v{snps.__version__}; '
            f'https://pypi.org/project/snps/"\n'
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

        df = s.snps

        p = s._parallelizer
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
                    "resources": s._resources,
                    "assembly": s.assembly,
                    "chrom": chrom,
                    "snps_df": pd.DataFrame(df.loc[(df["chrom"] == chrom)]),
                    "phased": s.phased,
                }
            )

        # drop chromosomes without reference sequence data (e.g., unassigned PAR)
        for chrom in chroms_to_drop:
            df = df.drop(df.loc[df["chrom"] == chrom].index)

        # create the VCF representation for SNPs
        results = p(cls._create_vcf_representation, tasks)

        contigs = []
        vcf = [pd.DataFrame()]
        discrepant_vcf_position = [pd.DataFrame()]
        for result in list(results):
            contigs.append(result["contig"])
            vcf.append(result["vcf"])
            discrepant_vcf_position.append(result["discrepant_vcf_position"])

        vcf = pd.concat(vcf)
        discrepant_vcf_position = pd.concat(discrepant_vcf_position)

        comment += "".join(contigs)
        comment += '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n'
        comment += "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE\n"

        return (
            save_df_as_csv(
                vcf,
                s._output_dir,
                filename,
                comment=comment,
                atomic=atomic,
                prepend_info=False,
                header=False,
                index=False,
                na_rep=".",
                sep="\t",
            ),
            discrepant_vcf_position,
        )

    @classmethod
    def _create_vcf_representation(cls, task):
        resources = task["resources"]
        assembly = task["assembly"]
        chrom = task["chrom"]
        snps_df = task["snps_df"]
        phased = task["phased"]

        if len(snps_df.loc[snps_df["genotype"].notnull()]) == 0:
            return {
                "contig": "",
                "vcf": pd.DataFrame(),
                "discrepant_vcf_position": pd.DataFrame(),
            }

        seqs = resources.get_reference_sequences(assembly, [chrom])
        seq = seqs[chrom]

        contig = f'##contig=<ID={seq.ID},URL={seq.url},length={seq.length},assembly={seq.build},md5={seq.md5},species="{seq.species}">\n'

        snps_df = snps_df.reset_index()

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
                "POS": np.uint32,
                "ID": object,
                "REF": object,
                "ALT": object,
                "QUAL": np.float32,
                "FILTER": object,
                "INFO": object,
                "FORMAT": object,
                "SAMPLE": object,
            }
        )

        df["CHROM"] = snps_df["chrom"]
        df["POS"] = snps_df["pos"]
        df["ID"] = snps_df["rsid"]

        # drop SNPs with discrepant positions (outside reference sequence)
        discrepant_vcf_position = snps_df.loc[
            (snps_df.pos - seq.start < 0) | (snps_df.pos - seq.start > seq.length - 1)
        ]
        df.drop(discrepant_vcf_position.index, inplace=True)

        # https://stackoverflow.com/a/24838429
        df["REF"] = list(map(chr, seq.sequence[df.POS - seq.start]))

        df["FORMAT"] = "GT"

        seq.clear()

        df["genotype"] = snps_df["genotype"]

        temp = df.loc[df["genotype"].notnull()]

        # https://stackoverflow.com/a/19976286
        df.loc[df["genotype"].notnull(), "ALT"] = np.vectorize(cls._compute_alt)(
            temp["REF"], temp["genotype"]
        )

        temp = df.loc[df["genotype"].notnull()]

        df.loc[df["genotype"].notnull(), "SAMPLE"] = np.vectorize(
            cls._compute_genotype
        )(temp["REF"], temp["ALT"], temp["genotype"], phased)

        df.loc[df["SAMPLE"].isnull(), "SAMPLE"] = "./."

        del df["genotype"]

        return {
            "contig": contig,
            "vcf": df,
            "discrepant_vcf_position": discrepant_vcf_position,
        }

    @staticmethod
    def _compute_alt(ref, genotype):
        genotype_alleles = list(set(genotype))

        if ref in genotype_alleles:
            if len(genotype_alleles) == 1:
                return "N"
            else:
                genotype_alleles.remove(ref)
                return genotype_alleles.pop(0)
        else:
            genotype_alleles.sort()
            return ",".join(genotype_alleles)

    @staticmethod
    def _compute_genotype(ref, alt, genotype, phased):
        alleles = [ref]

        if phased:
            separator = "|"
        else:
            separator = "/"

        if pd.notna(alt):
            alleles.extend(alt.split(","))

        if len(genotype) == 2:
            return (
                f"{alleles.index(genotype[0])}{separator}{alleles.index(genotype[1])}"
            )
        else:
            return f"{alleles.index(genotype[0])}"
