"""Class for writing SNPs."""

import logging
import warnings

import numpy as np
import pandas as pd

import snps
from snps.constants import REFERENCE_SEQUENCE_CHROMS
from snps.io import get_empty_snps_dataframe
from snps.utils import clean_str, get_utc_now, save_df_as_csv

logger = logging.getLogger(__name__)


class Writer:
    """Class for writing SNPs to files."""

    def __init__(
        self,
        snps=None,
        filename="",
        vcf=False,
        atomic=True,
        vcf_alt_unavailable=".",
        vcf_chrom_prefix="",
        vcf_qc_only=False,
        vcf_qc_filter=False,
        **kwargs,
    ):
        """Initialize a `Writer`.

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
        vcf_alt_unavailable : str
            representation of VCF ALT allele when ALT is not able to be determined
        vcf_chrom_prefix : str
            prefix for chromosomes in VCF CHROM column
        vcf_qc_only : bool
            for VCF, output only SNPs that pass quality control
        vcf_qc_filter : bool
            for VCF, populate VCF FILTER column based on quality control results
        **kwargs
            additional parameters to `pandas.DataFrame.to_csv`
        """
        self._snps = snps
        self._filename = filename
        self._vcf = vcf
        self._atomic = atomic
        self._vcf_alt_unavailable = vcf_alt_unavailable
        self._vcf_chrom_prefix = vcf_chrom_prefix
        self._vcf_qc_only = vcf_qc_only
        self._vcf_qc_filter = vcf_qc_filter
        self._kwargs = kwargs

    def write(self):
        """Write SNPs to file or buffer.

        Returns
        -------
        str
            path to file in output directory if SNPs were saved, else empty str
        discrepant_vcf_position : pd.DataFrame
            SNPs with discrepant positions discovered while saving VCF
        """
        if self._vcf:
            return self._write_vcf()
        else:
            return (self._write_csv(),)

    @classmethod
    def write_file(
        cls,
        snps=None,
        filename="",
        vcf=False,
        atomic=True,
        vcf_alt_unavailable=".",
        vcf_qc_only=False,
        vcf_qc_filter=False,
        **kwargs,
    ):
        warnings.warn(
            "This method will be removed in a future release.", DeprecationWarning
        )
        w = cls(
            snps=snps,
            filename=filename,
            vcf=vcf,
            atomic=atomic,
            vcf_alt_unavailable=vcf_alt_unavailable,
            vcf_qc_only=vcf_qc_only,
            vcf_qc_filter=vcf_qc_filter,
            **kwargs,
        )
        return w.write()

    def _write_csv(self):
        """Write SNPs to a CSV file.

        Returns
        -------
        str
            path to file in output directory if SNPs were saved, else empty str
        """
        filename = self._filename
        if not filename:
            ext = ".txt"

            if "sep" in self._kwargs and self._kwargs["sep"] == ",":
                ext = ".csv"

            filename = f"{clean_str(self._snps.source)}_{self._snps.assembly}{ext}"

        comment = [
            f"# Source(s): {self._snps.source}",
            f"# Build: {self._snps.build}",
            f"# Build Detected: { self._snps.build_detected}",
            f"# Phased: {self._snps.phased}",
            f"# SNPs: {self._snps.count}",
            f"# Chromosomes: {self._snps.chromosomes_summary}",
        ]
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
            comment="\n".join(comment) + "\n",
            atomic=self._atomic,
            **self._kwargs,
        )

    def _write_vcf(self):
        """Write SNPs to a VCF file.

        References
        ----------
        1. The Variant Call Format (VCF) Version 4.3 Specification, 27 Nov 2022,
           https://samtools.github.io/hts-specs/VCFv4.3.pdf

        Returns
        -------
        str
            path to file in output directory if SNPs were saved, else empty str
        discrepant_vcf_position : pd.DataFrame
            SNPs with discrepant positions discovered while saving VCF
        """
        filename = self._filename
        if not filename:
            filename = f"{clean_str(self._snps.source)}_{self._snps.assembly}{'.vcf'}"

        comment = [
            "##fileformat=VCFv4.3",
            f'##fileDate={get_utc_now().strftime("%Y%m%d")}',
            f'##source="snps v{snps.__version__}; https://pypi.org/project/snps/"',
            f'##detectedCompany="{self._snps.source}"',
        ]

        if self._snps.build_original:
            comment.append(f"##detectedOriginalBuild={self._snps.build_original}")

        if self._snps.determine_sex():
            comment.append(f"##detectedSex={self._snps.determine_sex()}")

        if self._vcf_qc_only or self._vcf_qc_filter:
            chip_version = ""
            if self._snps.chip_version:
                chip_version = f" {self._snps.chip_version}"

            if self._snps.chip:
                comment.append(
                    f'##detectedChip="{self._snps.chip}{chip_version} per Lu et al.: https://doi.org/10.1016/j.csbj.2021.06.040"'
                )

        df = self._snps.snps

        p = self._snps._parallelizer
        tasks = []

        chroms_to_drop = []
        for chrom in df["chrom"].unique():
            if chrom not in REFERENCE_SEQUENCE_CHROMS:
                chroms_to_drop.append(chrom)
                continue

            tasks.append(
                {
                    "resources": self._snps._resources,
                    "assembly": self._snps.assembly,
                    "chrom": chrom,
                    "snps": pd.DataFrame(df.loc[(df["chrom"] == chrom)]),
                    "cluster": (
                        self._snps.cluster
                        if self._vcf_qc_only or self._vcf_qc_filter
                        else ""
                    ),
                    "low_quality_snps": (
                        self._snps.low_quality
                        if self._vcf_qc_only or self._vcf_qc_filter
                        else get_empty_snps_dataframe()
                    ),
                    "sex": self._snps.determine_sex(),
                }
            )

        # drop chromosomes without reference sequence data (e.g., unassigned PAR)
        for chrom in chroms_to_drop:
            df = df.drop(df.loc[df["chrom"] == chrom].index)

        # Check for the presence of insertions or deletions
        has_ins = df["genotype"].str.contains("I", na=False).any()
        has_del = df["genotype"].str.contains("D", na=False).any()

        # create the VCF representation for SNPs
        results = p(self._create_vcf_representation, tasks)

        contigs = []
        vcf = [pd.DataFrame()]
        discrepant_vcf_position = [pd.DataFrame()]
        for result in list(results):
            if result["contig"]:
                contigs.append(result["contig"])
            vcf.append(result["vcf"])
            discrepant_vcf_position.append(result["discrepant_vcf_position"])

        vcf = pd.concat(vcf)
        discrepant_vcf_position = pd.concat(discrepant_vcf_position)

        comment.extend(contigs)

        if has_del:
            comment.append(
                '##ALT=<ID=DEL,Description="Deletion relative to the reference">'
            )
        if has_ins:
            comment.append(
                '##ALT=<ID=INS,Description="Insertion of novel sequence relative to the reference">'
            )

        if has_ins or has_del:
            comment.append(
                '##INFO=<ID=SVTYPE,Number=.,Type=String,Description="Type of structural variant: INS (Insertion), DEL (Deletion)">'
            )
            comment.append(
                '##INFO=<ID=IMPRECISE,Number=0,Type=Flag,Description="Imprecise structural variation">'
            )

        if self._vcf_qc_filter and self._snps.cluster:
            comment.append(
                '##FILTER=<ID=lq,Description="Low quality SNP per Lu et al.: https://doi.org/10.1016/j.csbj.2021.06.040">'
            )

        comment.append('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">')
        comment.append("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE")

        return (
            save_df_as_csv(
                vcf,
                self._snps._output_dir,
                filename,
                comment="\n".join(comment) + "\n",
                prepend_info=False,
                header=False,
                index=False,
                na_rep=".",
                sep="\t",
            ),
            discrepant_vcf_position,
        )

    def _create_vcf_representation(self, task):
        resources = task["resources"]
        assembly = task["assembly"]
        chrom = task["chrom"]
        snps = task["snps"]
        cluster = task["cluster"]
        low_quality_snps = task["low_quality_snps"]
        sex = task["sex"]

        if len(snps.loc[snps["genotype"].notnull()]) == 0:
            return {
                "contig": "",
                "vcf": pd.DataFrame(),
                "discrepant_vcf_position": pd.DataFrame(),
            }

        seqs = resources.get_reference_sequences(assembly, [chrom])
        seq = seqs[chrom]

        contig = f'##contig=<ID={self._vcf_chrom_prefix}{seq.ID},URL={seq.url},length={seq.length},assembly={seq.build},md5={seq.md5},species="{seq.species}">'

        if self._vcf_qc_only and cluster:
            # drop low quality SNPs if SNPs object maps to a cluster
            snps = snps.drop(snps.index.intersection(low_quality_snps.index))

        if self._vcf_qc_filter and cluster:
            # initialize filter for  all SNPs if SNPs object maps to a cluster,
            snps["filter"] = "PASS"
            # then indicate SNPs that were identified as low quality
            snps.loc[snps.index.intersection(low_quality_snps.index), "filter"] = "lq"

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

        df["CHROM"] = self._vcf_chrom_prefix + snps["chrom"]
        df["POS"] = snps["pos"]
        df["ID"] = snps["rsid"]

        if self._vcf_qc_filter and cluster:
            df["FILTER"] = snps["filter"]

        # drop SNPs with discrepant positions (outside reference sequence)
        discrepant_vcf_position = snps.loc[
            (snps.pos - seq.start < 0) | (snps.pos - seq.start > seq.length - 1)
        ]
        df.drop(discrepant_vcf_position.index, inplace=True)

        # https://stackoverflow.com/a/24838429
        df["REF"] = list(map(chr, seq.sequence[df.POS - seq.start]))

        df["FORMAT"] = "GT"

        seq.clear()

        df["genotype"] = snps["genotype"]

        temp = df.loc[df["genotype"].notnull()]

        # https://stackoverflow.com/a/19976286
        df.loc[df["genotype"].notnull(), "ALT"] = np.vectorize(self._compute_alt)(
            temp["REF"], temp["genotype"]
        )

        # Populate INFO field
        df["INFO"] = df["ALT"].apply(self._compute_info)

        temp = df.loc[df["genotype"].notnull()]

        df.loc[df["genotype"].notnull(), "SAMPLE"] = np.vectorize(
            self._compute_genotype
        )(temp["REF"], temp["ALT"], temp["genotype"])

        if sex == "Female":
            haploid_chroms = ["Y", "MT"]
        else:
            haploid_chroms = ["X", "Y", "MT"]

        # populate null values for haploid chromosomes
        df.loc[
            (df["SAMPLE"].isnull())
            & (df["CHROM"].str.contains("|".join(haploid_chroms))),
            "SAMPLE",
        ] = "."

        df.loc[df["SAMPLE"].isnull(), "SAMPLE"] = "./."

        del df["genotype"]

        return {
            "contig": contig,
            "vcf": df,
            "discrepant_vcf_position": discrepant_vcf_position,
        }

    def _replace_genotype_indels(self, genotype):
        # Replace 'I' and 'D' with '<INS>' and '<DEL>'
        return [
            "<INS>" if allele == "I" else "<DEL>" if allele == "D" else allele
            for allele in genotype
        ]

    def _compute_alt(self, ref, genotype):
        genotype_alleles = list(set(genotype))

        genotype_alleles = self._replace_genotype_indels(genotype_alleles)

        if ref in genotype_alleles:
            if len(genotype_alleles) == 1:
                return self._vcf_alt_unavailable
            else:
                genotype_alleles.remove(ref)
                return genotype_alleles.pop(0)
        else:
            genotype_alleles.sort()
            return ",".join(genotype_alleles)

    def _compute_genotype(self, ref, alt, genotype):
        genotype = list(genotype)

        genotype = self._replace_genotype_indels(genotype)

        alleles = [ref]

        if self._snps.phased:
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

    def _compute_info(self, alt):
        """Generate the INFO field based on ALT values."""
        if pd.isna(alt):
            return "."

        alt_values = alt.split(",")
        svtypes = []
        for alt_value in alt_values:
            if alt_value == "<INS>":
                svtypes.append("INS")
            elif alt_value == "<DEL>":
                svtypes.append("DEL")

        if not svtypes:
            return "."

        svtype_str = ",".join(svtypes)
        return f"SVTYPE={svtype_str};IMPRECISE" if svtype_str else "."
