""" Classes for reading and writing SNPs.

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
import os
import io
import gzip
import zipfile
import binascii
from copy import deepcopy

import numpy as np
import pandas as pd

import snps
from snps.utils import save_df_as_csv, clean_str

import logging

logger = logging.getLogger(__name__)


class Reader:
    """ Class for reading and parsing raw data / genotype files. """

    def __init__(self, file="", only_detect_source=False, resources=None, rsids=()):
        """ Initialize a `Reader`.

        Parameters
        ----------
        file : str or bytes
            path to file to load or bytes to load
        only_detect_source : bool
            only detect the source of the data
        resources : Resources
            instance of Resources
        rsids : tuple, optional
            rsids to extract if loading a VCF file

        """
        self._file = file
        self._only_detect_source = only_detect_source
        self._resources = resources
        self._rsids = rsids

    def __call__(self):
        """ Read and parse a raw data / genotype file.

        Returns
        -------
        dict
            dict with the following items:

            snps (pandas.DataFrame)
                dataframe of parsed SNPs
            source (str)
                detected source of SNPs
            phased (bool)
                flag indicating if SNPs are phased
        """
        file = self._file
        compression = "infer"
        d = {"snps": pd.DataFrame(), "source": "", "phased": False}

        # peek into files to determine the data format
        if isinstance(file, str) and os.path.exists(file):

            if ".zip" in file:
                with zipfile.ZipFile(file) as z:
                    with z.open(z.namelist()[0], "r") as f:
                        first_line, comments, data = self._extract_comments(
                            f, decode=True
                        )
            elif ".gz" in file:
                with gzip.open(file, "rt") as f:
                    first_line, comments, data = self._extract_comments(f)
            else:
                with open(file, "r") as f:
                    first_line, comments, data = self._extract_comments(f)

        elif isinstance(file, bytes):

            first_line, comments, data, compression = self._handle_bytes_data(file)
            file = io.BytesIO(file)

        else:
            return d

        if "23andMe" in first_line:
            return self.read_23andme(file, compression)
        elif "Ancestry" in first_line:
            return self.read_ancestry(file, compression)
        elif first_line.startswith("RSID"):
            return self.read_ftdna(file, compression)
        elif "famfinder" in first_line:
            return self.read_ftdna_famfinder(file, compression)
        elif "MyHeritage" in first_line:
            return self.read_myheritage(file, compression)
        elif "Living DNA" in first_line:
            return self.read_livingdna(file, compression)
        elif "SNP Name	rsID	Sample.ID	Allele1...Top" in first_line:
            return self.read_mapmygenome(file, compression)
        elif "lineage" in first_line or "snps" in first_line:
            return self.read_snps_csv(file, comments, compression)
        elif first_line.startswith("rsid"):
            return self.read_generic_csv(file, compression)
        elif "vcf" in comments.lower() or "##contig" in comments.lower():
            return self.read_vcf(file, compression, self._rsids)
        elif ("Genes for Good" in comments) | ("PLINK" in comments):
            return self.read_genes_for_good(file, compression)
        elif "CODIGO46" in comments:
            return self.read_codigo46(file)
        elif "SANO" in comments:
            return self.read_sano(file)
        else:
            return d

    @classmethod
    def read_file(cls, file, only_detect_source, resources, rsids):
        """ Read `file`.

        Parameters
        ----------
        file : str or bytes
            path to file to load or bytes to load
        only_detect_source : bool
            only detect the source of the data
        resources : Resources
            instance of Resources
        rsids : tuple
            rsids to extract if loading a VCF file

        Returns
        -------
        dict
            dict with the following items:

            snps (pandas.DataFrame)
                dataframe of parsed SNPs
            source (str)
                detected source of SNPs
            phased (bool)
                flag indicating if SNPs are phased
        """
        r = cls(file, only_detect_source, resources, rsids)
        return r()

    def _extract_comments(self, f, decode=False, include_data=False):
        line = self._read_line(f, decode)

        first_line = line
        comments = ""
        data = ""

        if first_line.startswith("#"):
            while line.startswith("#"):
                comments += line
                line = self._read_line(f, decode)
            if include_data:
                while line:
                    data += line
                    line = self._read_line(f, decode)

        elif first_line.startswith("[Header]"):
            while not line.startswith("[Data]"):
                comments += line
                line = self._read_line(f, decode)
            # Ignore the [Data] row
            line = self._read_line(f, decode)
            if include_data:
                while line:
                    data += line
                    line = self._read_line(f, decode)
        if not isinstance(f, zipfile.ZipExtFile):
            f.seek(0)
        return first_line, comments, data

    def _handle_bytes_data(self, file, include_data=False):
        compression = "infer"
        if self.is_zip(file):
            compression = "zip"
            with zipfile.ZipFile(io.BytesIO(file)) as z:
                namelist = z.namelist()
                key = "GFG_filtered_unphased_genotypes_23andMe.txt"
                key_search = [key in name for name in namelist]

                if any(key_search):
                    filename = namelist[key_search.index(True)]
                else:
                    filename = namelist[0]

                with z.open(filename, "r") as f:
                    first_line, comments, data = self._extract_comments(
                        f, decode=True, include_data=include_data
                    )

        elif self.is_gzip(file):
            compression = "gzip"

            with gzip.open(io.BytesIO(file), "rb") as f:
                first_line, comments, data = self._extract_comments(
                    f, decode=True, include_data=include_data
                )

        else:
            file = io.BytesIO(file)
            first_line, comments, data = self._extract_comments(
                deepcopy(file), decode=True, include_data=include_data
            )
            file.seek(0)
        return first_line, comments, data, compression

    @staticmethod
    def is_zip(bytes_data):
        """Check whether or not a bytes_data file is a valid Zip file."""
        return zipfile.is_zipfile(io.BytesIO(bytes_data))

    @staticmethod
    def is_gzip(bytes_data):
        """Check whether or not a bytes_data file is a valid gzip file."""
        return binascii.hexlify(bytes_data[:2]) == b"1f8b"

    @staticmethod
    def _read_line(f, decode):
        if decode:
            # https://stackoverflow.com/a/606199
            return f.readline().decode("utf-8")
        else:
            return f.readline()

    def read_helper(self, source, parser):
        """ Generic method to help read files.

        Parameters
        ----------
        source : str
            name of data source
        parser : func
            parsing function, which returns a tuple with the following items:

            0 (pandas.DataFrame)
                dataframe of parsed SNPs (empty if only detecting source)
            1 (bool), optional
                flag indicating if SNPs are phased

        Returns
        -------
        dict
            dict with the following items:

            snps (pandas.DataFrame)
                dataframe of parsed SNPs
            source (str)
                detected source of SNPs
            phased (bool)
                flag indicating if SNPs are phased

        References
        ----------
        1. Fluent Python by Luciano Ramalho (O'Reilly). Copyright 2015 Luciano Ramalho,
           978-1-491-94600-8.
        """
        phased = False

        if self._only_detect_source:
            df = pd.DataFrame()
        else:
            df, *extra = parser()

            if len(extra) == 1:
                phased = extra[0]

        return {"snps": df, "source": source, "phased": phased}

    def read_23andme(self, file, compression):
        """ Read and parse 23andMe file.

        https://www.23andme.com

        Parameters
        ----------
        file : str
            path to file

        Returns
        -------
        dict
            result of `read_helper`
        """

        def parser():
            return (
                pd.read_csv(
                    file,
                    comment="#",
                    sep="\t",
                    na_values="--",
                    names=["rsid", "chrom", "pos", "genotype"],
                    index_col=0,
                    dtype={"chrom": object},
                    compression=compression,
                ),
            )

        return self.read_helper("23andMe", parser)

    def read_ftdna(self, file, compression):
        """Read and parse Family Tree DNA (FTDNA) file.

        https://www.familytreedna.com

        Parameters
        ----------
        file : str
            path to file

        Returns
        -------
        dict
            result of `read_helper`
        """

        def parser():
            df = pd.read_csv(
                file,
                skiprows=1,
                na_values="--",
                names=["rsid", "chrom", "pos", "genotype"],
                index_col=0,
                dtype={"chrom": object},
                compression=compression,
            )

            # remove incongruous data
            df = df.drop(
                df.loc[df.index == "RSID"].index
            )  # second header for concatenated data

            # if second header existed, pos dtype will be object (should be np.int64)
            df["pos"] = df["pos"].astype(np.int64)

            return (df,)

        return self.read_helper("FTDNA", parser)

    def read_ftdna_famfinder(self, file, compression):
        """ Read and parse Family Tree DNA (FTDNA) "famfinder" file.

        https://www.familytreedna.com

        Parameters
        ----------
        file : str
            path to file

        Returns
        -------
        dict
            result of `read_helper`
        """

        def parser():
            df = pd.read_csv(
                file,
                comment="#",
                na_values="-",
                names=["rsid", "chrom", "pos", "allele1", "allele2"],
                index_col=0,
                dtype={"chrom": object},
                compression=compression,
            )

            # create genotype column from allele columns
            df["genotype"] = df["allele1"] + df["allele2"]

            # delete allele columns
            # http://stackoverflow.com/a/13485766
            del df["allele1"]
            del df["allele2"]

            return (df,)

        return self.read_helper("FTDNA", parser)

    def read_ancestry(self, file, compression):
        """ Read and parse Ancestry.com file.

        http://www.ancestry.com

        Parameters
        ----------
        file : str
            path to file

        Returns
        -------
        dict
            result of `read_helper`
        """

        def parser():
            df = pd.read_csv(
                file,
                comment="#",
                header=0,
                sep="\t",
                na_values=0,
                names=["rsid", "chrom", "pos", "allele1", "allele2"],
                index_col=0,
                dtype={"chrom": object},
                compression=compression,
            )

            # create genotype column from allele columns
            df["genotype"] = df["allele1"] + df["allele2"]

            # delete allele columns
            # http://stackoverflow.com/a/13485766
            del df["allele1"]
            del df["allele2"]

            # https://redd.it/5y90un
            df.iloc[np.where(df["chrom"] == "23")[0], 0] = "X"
            df.iloc[np.where(df["chrom"] == "24")[0], 0] = "Y"
            df.iloc[np.where(df["chrom"] == "25")[0], 0] = "PAR"
            df.iloc[np.where(df["chrom"] == "26")[0], 0] = "MT"

            return (df,)

        return self.read_helper("AncestryDNA", parser)

    def read_myheritage(self, file, compression):
        """ Read and parse MyHeritage file.

        https://www.myheritage.com

        Parameters
        ----------
        file : str
            path to file

        Returns
        -------
        dict
            result of `read_helper`
        """

        def parser():
            return (
                pd.read_csv(
                    file,
                    comment="#",
                    header=0,
                    na_values="--",
                    names=["rsid", "chrom", "pos", "genotype"],
                    index_col=0,
                    dtype={"chrom": object, "pos": np.int64},
                    compression=compression,
                ),
            )

        return self.read_helper("MyHeritage", parser)

    def read_livingdna(self, file, compression):
        """ Read and parse LivingDNA file.

        https://livingdna.com/

        Parameters
        ----------
        file : str
            path to file

        Returns
        -------
        dict
            result of `read_helper`
        """

        def parser():
            return (
                pd.read_csv(
                    file,
                    comment="#",
                    sep="\t",
                    na_values="--",
                    names=["rsid", "chrom", "pos", "genotype"],
                    index_col=0,
                    dtype={"chrom": object},
                    compression=compression,
                ),
            )

        return self.read_helper("LivingDNA", parser)

    def read_mapmygenome(self, file, compression):
        """ Read and parse Mapmygenome file.

        https://mapmygenome.in

        Parameters
        ----------
        file : str
            path to file

        Returns
        -------
        dict
            result of `read_helper`
        """

        def parser():
            df = pd.read_csv(
                file,
                comment="#",
                sep="\t",
                na_values="--",
                header=0,
                index_col=1,
                dtype={"Chr": object},
                compression=compression,
            )

            df["genotype"] = df["Allele1...Top"] + df["Allele2...Top"]
            df.rename(columns={"Chr": "chrom", "Position": "pos"}, inplace=True)
            df.index.name = "rsid"
            df = df[["chrom", "pos", "genotype"]]

            return (df,)

        return self.read_helper("Mapmygenome", parser)

    def read_genes_for_good(self, file, compression):
        """ Read and parse Genes For Good file.

        https://genesforgood.sph.umich.edu/readme/readme1.2.txt

        Parameters
        ----------
        file : str
            path to file

        Returns
        -------
        dict
            result of `read_helper`
        """

        def parser():
            return (
                pd.read_csv(
                    file,
                    comment="#",
                    sep="\t",
                    na_values="--",
                    names=["rsid", "chrom", "pos", "genotype"],
                    index_col=0,
                    dtype={"chrom": object},
                    compression=compression,
                ),
            )

        return self.read_helper("GenesForGood", parser)

    def _read_gsa_helper(self, file, source, strand, na_values="--"):
        def parser():
            gsa_resources = self._resources.get_gsa_resources()

            if isinstance(file, str):
                with open(file, "rb") as f:
                    first_line, comments, data = self._extract_comments(
                        f, decode=True, include_data=True
                    )
            else:
                first_line, comments, data, compression = self._handle_bytes_data(
                    file.read(), include_data=True
                )

            df = pd.read_csv(io.StringIO(data), sep="\t", na_values=na_values)

            def map_rsids(x):
                return gsa_resources["rsid_map"].get(x)

            def map_chr(x):
                chrpos = gsa_resources["chrpos_map"].get(x)
                return chrpos.split(":")[0] if chrpos else None

            def map_pos(x):
                chrpos = gsa_resources["chrpos_map"].get(x)
                return chrpos.split(":")[1] if chrpos else None

            df["rsid"] = df["SNP Name"].apply(map_rsids)
            df["chrom"] = df["SNP Name"].apply(map_chr)
            df["pos"] = df["SNP Name"].apply(map_pos)
            df["genotype"] = (
                df["Allele1 - {}".format(strand)] + df["Allele2 - {}".format(strand)]
            )
            df.dropna(subset=["rsid", "chrom", "pos"], inplace=True)

            df = df.astype({"chrom": object, "pos": np.int64})
            df = df[["rsid", "chrom", "pos", "genotype"]]
            df.set_index(["rsid"], inplace=True)

            return (df,)

        return self.read_helper(source, parser)

    def read_codigo46(self, file):
        """ Read and parse Codigo46 files.

        https://codigo46.com.mx

        Parameters
        ----------
        data : str
            data string

        Returns
        -------
        dict
            result of `read_helper`
        """
        return self._read_gsa_helper(file, "Codigo46", "Plus")

    def read_sano(self, file):
        """ Read and parse Sano Genetics files.

        https://sanogenetics.com

        Parameters
        ----------
        data : str
            data string

        Returns
        -------
        dict
            result of `read_helper`
        """
        return self._read_gsa_helper(file, "Sano", "Forward", na_values="-")

    def read_snps_csv(self, file, comments, compression):
        """ Read and parse CSV file generated by `snps`.

        Parameters
        ----------
        file : str
            path to file
        comments : str
            comments at beginning of file

        Returns
        -------
        pandas.DataFrame
            genetic data normalized for use with `snps`
        str
            name of data source(s)
        """
        source = ""
        phased = False

        for comment in comments.split("\n"):
            if "Source(s):" in comment:
                source = comment.split("Source(s):")[1].strip()
                break

        for comment in comments.split("\n"):
            if "Phased:" in comment:
                phased_str = comment.split("Phased:")[1].strip()
                if phased_str == "True":
                    phased = True
                break

        def parser():
            return (
                pd.read_csv(
                    file,
                    comment="#",
                    header=0,
                    na_values="--",
                    names=["rsid", "chrom", "pos", "genotype"],
                    index_col=0,
                    dtype={"chrom": object, "pos": np.int64},
                    compression=compression,
                ),
                phased,
            )

        return self.read_helper(source, parser)

    def read_generic_csv(self, file, compression):
        """ Read and parse generic CSV file.

        Notes
        -----
        Assumes columns are 'rsid', 'chrom' / 'chromosome', 'pos' / 'position', and 'genotype';
        values are comma separated; unreported genotypes are indicated by '--'; and one header row
        precedes data. For example:

            rsid,chromosome,position,genotype
            rs1,1,1,AA
            rs2,1,2,CC
            rs3,1,3,--

        Parameters
        ----------
        file : str
            path to file

        Returns
        -------
        dict
            result of `read_helper`
        """

        def parser():
            return (
                pd.read_csv(
                    file,
                    skiprows=1,
                    na_values="--",
                    names=["rsid", "chrom", "pos", "genotype"],
                    index_col=0,
                    dtype={"chrom": object, "pos": np.int64},
                    compression=compression,
                ),
            )

        return self.read_helper("generic", parser)

    def read_vcf(self, file, compression, rsids=()):
        """ Read and parse VCF file.

        Notes
        -----
        This method attempts to read and parse a VCF file or buffer, optionally
        compressed with gzip. Some assumptions are made throughout this process:

            * SNPs that are not annotated with an RSID are skipped
            * If the VCF contains multiple samples, only the first sample is used to
              lookup the genotype
            * Insertions and deletions are skipped
            * If a sample allele is not specified, the genotype is reported as NaN
            * If a sample allele refers to a REF or ALT allele that is not specified,
              the genotype is reported as NaN

        Parameters
        ----------
        file : str or bytes
            path to file or bytes to load
        rsids : tuple, optional
            rsids to extract if loading a VCF file

        Returns
        -------
        dict
            result of `read_helper`
        """

        def parser():
            if not isinstance(file, io.BytesIO):
                with open(file, "rb") as f:
                    df, phased = self._parse_vcf(f, rsids)
            else:
                df, phased = self._parse_vcf(file, rsids)

            return (df, phased)

        return self.read_helper("vcf", parser)

    def _parse_vcf(self, buffer, rsids):
        rows = []
        phased = True
        first_four_bytes = buffer.read(4)
        buffer.seek(0)

        if self.is_gzip(first_four_bytes):
            f = gzip.open(buffer)
        else:
            f = buffer

        with io.TextIOWrapper(io.BufferedReader(f)) as file:

            for line in file:

                line_strip = line.strip("\n")
                if line_strip.startswith("#"):
                    continue
                rsid = line_strip.split("\t")[2]
                # skip SNPs with missing rsIDs.
                if rsid == ".":
                    continue
                if rsids:
                    if rsid not in rsids:
                        continue

                line_split = line_strip.split("\t")

                # snps does not yet support multi-sample vcf.
                if len(line_split) > 10:
                    logger.info("Multiple samples detected in the vcf file")

                ref = line_split[3]
                alt = line_split[4]
                zygote = line_split[9]
                zygote = zygote.split(":")[0]

                ref_alt = [ref] + alt.split(",")

                # skip insertions and deletions
                if sum(map(len, ref_alt)) > len(ref_alt):
                    continue

                zygote1, zygote2 = zygote.replace("|", " ").replace("/", " ").split(" ")
                if zygote1 == "." or zygote2 == ".":
                    # assign null genotypes if either allele is None
                    genotype = np.nan
                elif (zygote1 == "0" or zygote2 == "0") and ref == ".":
                    # sample allele specifies REF allele, which is None
                    genotype = np.nan
                elif (zygote1 == "1" or zygote2 == "1") and alt == ".":
                    # sample allele specifies ALT allele, which is None
                    genotype = np.nan
                else:
                    # Could capture full genotype, if REF is None, but genotype is 1/1 or
                    # if ALT is None, but genotype is 0/0
                    genotype = ref_alt[int(zygote1)] + ref_alt[int(zygote2)]

                if "/" in zygote and pd.notna(genotype):
                    phased = False

                record_array = [
                    rsid,
                    "{}".format(line_split[0]).strip("chr"),
                    line_split[1],
                    genotype,
                ]
                rows.append(record_array)

            df = pd.DataFrame(rows, columns=["rsid", "chrom", "pos", "genotype"])
            df = df.astype(
                {"rsid": object, "chrom": object, "pos": np.int64, "genotype": object}
            )

            df.set_index("rsid", inplace=True, drop=True)

        return (df, phased)


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
        vcf = []
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
