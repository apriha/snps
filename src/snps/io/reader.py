""" Class for reading SNPs.

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

import binascii
from copy import deepcopy
import gzip
import io
import logging
import os
import re
import zipfile
import zlib

import numpy as np
import pandas as pd


logger = logging.getLogger(__name__)

NORMALIZED_DTYPES = {
    "rsid": object,
    "chrom": object,
    "pos": np.uint32,
    "genotype": object,
}

TWO_ALLELE_DTYPES = {
    "rsid": object,
    "chrom": object,
    "pos": np.uint32,
    "allele1": object,
    "allele2": object,
}


def get_empty_snps_dataframe():
    """ Get empty dataframe normalized for usage with ``snps``.

    Returns
    -------
    pd.DataFrame
    """
    df = pd.DataFrame(columns=["rsid", "chrom", "pos", "genotype"])
    df = df.astype(NORMALIZED_DTYPES)
    df.set_index("rsid", inplace=True)
    return df


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
        self._rsids = frozenset(rsids)

    def read(self):
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
        d = {
            "snps": get_empty_snps_dataframe(),
            "source": "",
            "phased": False,
            "build": 0,
        }

        # peek into files to determine the data format
        if isinstance(file, str) and os.path.exists(file):

            if ".zip" in file:
                with zipfile.ZipFile(file) as z:
                    with z.open(z.namelist()[0], "r") as f:
                        first_line, comments, data = self._extract_comments(
                            f, decode=True
                        )
                    compression = "zip"
            elif ".gz" in file:
                with gzip.open(file, "rt") as f:
                    first_line, comments, data = self._extract_comments(f)
                compression = "gzip"
            else:
                with open(file, "rb") as f:
                    first_line, comments, data, compression = self._handle_bytes_data(
                        f.read()
                    )

        elif isinstance(file, bytes):

            first_line, comments, data, compression = self._handle_bytes_data(file)
            file = io.BytesIO(file)

        else:
            return d

        if "23andMe" in first_line:
            d = self.read_23andme(file, compression)
        elif "Ancestry" in first_line:
            d = self.read_ancestry(file, compression)
        elif first_line.startswith("RSID"):
            d = self.read_ftdna(file, compression)
        elif "famfinder" in first_line:
            d = self.read_ftdna_famfinder(file, compression)
        elif "MyHeritage" in first_line:
            d = self.read_myheritage(file, compression)
        elif "Living DNA" in first_line:
            d = self.read_livingdna(file, compression)
        elif "SNP Name\trsID" in first_line or "SNP.Name\tSample.ID" in first_line:
            d = self.read_mapmygenome(file, compression, first_line)
        elif "lineage" in first_line or "snps" in first_line:
            d = self.read_snps_csv(file, comments, compression)
        elif "Chromosome" in first_line:
            d = self.read_tellmegen(file, compression)
        elif re.match("^#*[ \t]*rsid[, \t]*chr", first_line):
            d = self.read_generic(file, compression)
        elif re.match("^rs[0-9]*[, \t]{1}[1]", first_line):
            d = self.read_generic(file, compression, skip=0)
        elif "vcf" in comments.lower() or "##contig" in comments.lower():
            d = self.read_vcf(file, compression, "vcf", self._rsids)
        elif ("Genes for Good" in comments) | ("PLINK" in comments):
            d = self.read_genes_for_good(file, compression)
        elif "DNA.Land" in comments:
            d = self.read_dnaland(file, compression)
        elif "CODIGO46" in comments:
            d = self.read_codigo46(file)
        elif "SANO" in comments:
            d = self.read_sano(file)

        # detect build from comments if build was not already detected from `read` method
        if not d["build"]:
            d.update({"build": self._detect_build_from_comments(comments, d["source"])})

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
        return r.read()

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
        if not data and include_data:
            data = f.read()
            if decode:
                data = data.decode()
        if not isinstance(f, zipfile.ZipExtFile):
            f.seek(0)
        return first_line, comments, data

    def _detect_build_from_comments(self, comments, source):
        comments = comments.lower()
        if "build 37" in comments:
            return 37
        elif "build 36" in comments:
            return 36

        # allow more variations for VCF
        if source == "vcf":
            if "https://pypi.org/project/snps/" in comments:  # remove `snps` version
                comments = f"{comments[: comments.find('snps v')]}{comments[comments.find('https://pypi.org/project/snps/'):]}"
            if "hg19" in comments:
                return 37
            elif "ncbi36" in comments:
                return 36
            elif "grch38" in comments:
                return 38
            elif "build 38" in comments:
                return 38
            elif "b37" in comments:
                return 37
            elif "hg38" in comments:
                return 38
            elif "b38" in comments:
                return 38
        return 0

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
            compression = None
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
            2 (int), optional
                detected build of SNPs

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
            build (int)
                detected build of SNPs

        References
        ----------
        1. Fluent Python by Luciano Ramalho (O'Reilly). Copyright 2015 Luciano Ramalho,
           978-1-491-94600-8.
        """
        phased = False
        build = 0

        if self._only_detect_source:
            df = get_empty_snps_dataframe()
        else:
            df, *extra = parser()

            if len(extra) == 1:
                phased = extra[0]
            elif len(extra) == 2:
                phased = extra[0]
                build = extra[1]

        return {"snps": df, "source": source, "phased": phased, "build": build}

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
                    dtype=NORMALIZED_DTYPES,
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
            try:
                df = pd.read_csv(
                    file,
                    skiprows=1,
                    na_values="--",
                    names=["rsid", "chrom", "pos", "genotype"],
                    index_col=0,
                    dtype=NORMALIZED_DTYPES,
                    compression=compression,
                )
            except ValueError:
                # read files with second header for concatenated data
                if isinstance(file, io.BytesIO):
                    file.seek(0)
                    (*data,) = self._handle_bytes_data(file.read(), include_data=True)
                    file.seek(0)
                else:
                    with open(file, "rb") as f:
                        (*data,) = self._handle_bytes_data(f.read(), include_data=True)
                # reconstruct file content from `_handle_bytes_data` results
                lines = data[0] + data[2]
                lines = [line.strip() for line in lines.split("\n")]
                # find index of second header
                second_header_idx = lines.index("RSID,CHROMOSOME,POSITION,RESULT", 1)

                df = pd.read_csv(
                    file,
                    skiprows=[0, second_header_idx],
                    na_values="--",
                    names=["rsid", "chrom", "pos", "genotype"],
                    index_col=0,
                    dtype=NORMALIZED_DTYPES,
                    compression=compression,
                )
            except OSError:
                # read concatenated gzip files with extra data
                if isinstance(file, io.BytesIO):
                    file.seek(0)
                    data = file.getbuffer()
                else:
                    with open(file, "rb") as f:
                        data = f.read()

                # https://stackoverflow.com/q/4928560
                # https://stackoverflow.com/a/37042747
                decompressor = zlib.decompressobj(31)

                # decompress data from first concatenated gzip file
                data = decompressor.decompress(data)

                # decompress data from second concatenated gzip file
                additional_data = zlib.decompress(decompressor.unused_data, 31)
                data += additional_data[33:]  # skip over second header

                new_file = io.BytesIO(data)

                df = pd.read_csv(
                    new_file,
                    skiprows=1,
                    na_values="--",
                    names=["rsid", "chrom", "pos", "genotype"],
                    index_col=0,
                    dtype=NORMALIZED_DTYPES,
                    compression=None,  # already decompressed
                )

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
                dtype=TWO_ALLELE_DTYPES,
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
                engine="c",
                sep="\s+",
                # delim_whitespace=True,  # https://stackoverflow.com/a/15026839
                na_values=0,
                names=["rsid", "chrom", "pos", "allele1", "allele2"],
                index_col=0,
                dtype=TWO_ALLELE_DTYPES,
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

            if isinstance(file, str):
                with open(file, "rb") as f:
                    first_line, comments, data, comrpession = self._handle_bytes_data(
                        f.read(), include_data=True
                    )
            else:
                first_line, comments, data, compression = self._handle_bytes_data(
                    file.read(), include_data=True
                )

            file_string_in = io.StringIO(data)
            file_string_out = io.StringIO()
            for line in file_string_in:
                # user the number of quotes in a line to tell old from new
                if line.count('"') == 14:
                    # extra-quoted new variant file
                    # can all be stripped so pandas can read it normally
                    line = line.replace('"', "")
                    # take it apart and put it back together so it looks
                    # like the older MyHeritage files
                    line = '"' + '","'.join(line.strip().split(",")) + '"\n'
                file_string_out.write(line)

            return (
                pd.read_csv(
                    io.StringIO(file_string_out.getvalue()),
                    comment="#",
                    header=0,
                    na_values="--",
                    names=["rsid", "chrom", "pos", "genotype"],
                    index_col=0,
                    dtype=NORMALIZED_DTYPES,
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
                    dtype=NORMALIZED_DTYPES,
                    compression=compression,
                ),
            )

        return self.read_helper("LivingDNA", parser)

    def read_mapmygenome(self, file, compression, header):
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
            def parse(rsid_col_name, rsid_col_idx):
                return pd.read_csv(
                    file,
                    comment="#",
                    sep="\t",
                    na_values="--",
                    header=0,
                    index_col=rsid_col_idx,
                    dtype={
                        rsid_col_name: object,
                        "Chr": object,
                        "Position": np.uint32,
                        "Allele1...Top": object,
                        "Allele2...Top": object,
                    },
                    compression=compression,
                )

            if "rsID" in header:
                df = parse("rsID", 1)
            else:
                df = parse("SNP.Name", 0)

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
                    dtype=NORMALIZED_DTYPES,
                    compression=compression,
                ),
            )

        return self.read_helper("GenesForGood", parser)

    def _read_gsa_helper(self, file, source, strand, dtypes, na_values="--"):
        def parser():
            gsa_resources = self._resources.get_gsa_resources()

            if isinstance(file, str):
                try:
                    with open(file, "rb") as f:
                        first_line, comments, data = self._extract_comments(
                            f, decode=True, include_data=True
                        )
                except UnicodeDecodeError:
                    # compressed file on filesystem
                    with open(file, "rb") as f:
                        (
                            first_line,
                            comments,
                            data,
                            compression,
                        ) = self._handle_bytes_data(f.read(), include_data=True)
            else:
                first_line, comments, data, compression = self._handle_bytes_data(
                    file.read(), include_data=True
                )

            data_io = io.StringIO(data)
            df = pd.read_csv(data_io, sep="\t", dtype=dtypes, na_values=na_values)

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
            df["genotype"] = df[f"Allele1 - {strand}"] + df[f"Allele2 - {strand}"]
            df.dropna(subset=["rsid", "chrom", "pos"], inplace=True)

            df = df.astype(NORMALIZED_DTYPES)
            df = df[["rsid", "chrom", "pos", "genotype"]]
            df.set_index(["rsid"], inplace=True)

            return (df,)

        return self.read_helper(source, parser)

    def read_tellmegen(self, file, compression):
        """ Read and parse tellmeGen files.

        https://www.tellmegen.com/

        Parameters
        ----------
        data : str
            data string

        Returns
        -------
        dict
            result of `read_helper`
        """

        def parser():
            gsa_resources = self._resources.get_gsa_resources()

            df = pd.read_csv(
                file,
                sep="\t",
                skiprows=1,
                na_values="--",
                names=["rsid", "chrom", "pos", "genotype"],
                index_col=0,
                dtype=NORMALIZED_DTYPES,
                compression=compression,
            )
            df.rename(index=gsa_resources["rsid_map"], inplace=True)
            return (df,)

        return self.read_helper("tellmeGen", parser)

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
        return self._read_gsa_helper(file, "Codigo46", "Plus", {})

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
        dtype = {"Chr": object, "Position": np.uint32}
        return self._read_gsa_helper(file, "Sano", "Forward", dtype, na_values="-")

    def read_dnaland(self, file, compression):
        """ Read and parse DNA.land files.

        https://dna.land/

        Parameters
        ----------
        data : str
            data string

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
                    dtype=NORMALIZED_DTYPES,
                    compression=compression,
                ),
            )

        return self.read_helper("DNA.Land", parser)

    def read_snps_csv(self, file, comments, compression):
        """ Read and parse CSV file generated by ``snps``.

        https://pypi.org/project/snps/

        Parameters
        ----------
        file : str or buffer
            path to file or buffer to read
        comments : str
            comments at beginning of file

        Returns
        -------
        dict
            result of `read_helper`
        """
        source = ""
        phased = False
        build = 0

        comment_lines = comments.split("\n")
        for comment1 in comment_lines:
            if "Source(s):" in comment1:
                source = comment1.split("Source(s):")[1].strip()
            elif "Phased:" in comment1:
                if comment1.split("Phased:")[1].strip() == "True":
                    phased = True
            elif "Build:" in comment1:
                temp = int(comment1.split("Build:")[1].strip())
                for comment2 in comment_lines:
                    if "Build Detected:" in comment2:
                        # only assign build if it was detected
                        if comment2.split("Build Detected:")[1].strip() == "True":
                            build = temp
                        break

        def parser():
            def parse_csv(sep):
                return pd.read_csv(
                    file,
                    sep=sep,
                    comment="#",
                    header=0,
                    na_values="--",
                    names=["rsid", "chrom", "pos", "genotype"],
                    index_col=0,
                    dtype=NORMALIZED_DTYPES,
                    compression=compression,
                )

            try:
                return (parse_csv(","), phased)
            except pd.errors.ParserError:
                if isinstance(file, io.BufferedIOBase):
                    file.seek(0)

                return (parse_csv("\t"), phased, build)

        return self.read_helper(source, parser)

    def read_generic(self, file, compression, skip=1):
        """ Read and parse generic CSV or TSV file.

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
            def parse(sep):
                return pd.read_csv(
                    file,
                    sep=sep,
                    skiprows=skip,
                    na_values="--",
                    names=["rsid", "chrom", "pos", "genotype"],
                    index_col=0,
                    dtype=NORMALIZED_DTYPES,
                    compression=compression,
                )

            try:
                df = parse(",")
            except ValueError:
                try:
                    if isinstance(file, io.BufferedIOBase):
                        file.seek(0)

                    df = parse("\t")
                except ValueError:
                    if isinstance(file, io.BufferedIOBase):
                        file.seek(0)

                    df = pd.read_csv(
                        file,
                        sep=None,
                        na_values="--",
                        skiprows=skip,
                        engine="python",
                        names=["rsid", "chrom", "pos", "genotype"],
                        usecols=[0, 1, 2, 3],
                        index_col=0,
                        dtype=NORMALIZED_DTYPES,
                        compression=compression,
                    )
            return (df,)

        return self.read_helper("generic", parser)

    def read_vcf(self, file, compression, provider, rsids=()):
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

        return self.read_helper(provider, parser)

    def _parse_vcf(self, buffer, rsids):
        rows = []
        phased = True
        first_four_bytes = buffer.read(4)
        buffer.seek(0)

        if self.is_gzip(first_four_bytes):
            f = gzip.open(buffer)
        else:
            f = buffer

        logged_multi_sample = False

        with io.TextIOWrapper(io.BufferedReader(f)) as file:

            for line in file:
                line_strip = line.strip("\n")

                # skip blank lines
                if not line_strip:
                    continue

                # skip comment lines
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
                if not logged_multi_sample and len(line_split) > 10:
                    logger.info("Multiple samples detected in the vcf file")
                    logged_multi_sample = True

                ref = line_split[3]
                alt = line_split[4]
                if len(alt.split(",")) > 1 and alt.split(",")[1] == "<NON_REF>":
                    alt = alt.split(",")[0]

                ref_alt = [ref] + alt.split(",")

                # skip insertions and deletions
                if sum(map(len, ref_alt)) > len(ref_alt):
                    continue

                # GT (genotype) is usually the first sample-specific field
                # | = diploid phased
                # / = diploid unphased
                # or haploid e.g. male sex chromosome
                genotype = ""
                zygote = line_split[9]
                zygote = zygote.split(":")[0]
                for z in zygote.replace("|", "/").split("/"):
                    if z == ".":
                        # missing genotype
                        genotype = np.nan
                        break
                    z = int(z)
                    if z >= len(ref_alt):
                        # invalid genotype number
                        genotype = np.nan
                        break
                    elif ref_alt[z] == ".":
                        # missing genotype in ref or alt
                        genotype = np.nan
                        break
                    else:
                        genotype = genotype + ref_alt[z]

                if "/" in zygote and pd.notna(genotype):
                    phased = False

                record_array = [
                    rsid,
                    f"{line_split[0]}".strip("chr"),
                    line_split[1],
                    genotype,
                ]
                rows.append(record_array)

            if len(rows) == 0:
                phased = False

            df = pd.DataFrame(rows, columns=["rsid", "chrom", "pos", "genotype"])
            df = df.astype(NORMALIZED_DTYPES)

            df.set_index("rsid", inplace=True, drop=True)

        return (df, phased)
