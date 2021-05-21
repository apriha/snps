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

    @classmethod
    def read_file(cls, file, only_detect_source=False, resources=None, rsids=()):
        """ Read and parse a raw data / genotype file.

        Notes
        -----
        This method is called when instantiating a ``SNPs`` object, and calling it directly is
        not normally required.

        If `file` is on the file system, it can optionally include an extension. Also, `file`
        can optionally be compressed with zip or gzip.

        If `file` is a VCF file or buffer, some assumptions are made:

            * VCF is uncompressed, or compressed with gzip
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
            build (int)
                detected build of SNPs from parser or comments
        """
        rsids = frozenset(rsids)

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
                        first_line, comments, data = cls._extract_comments(
                            f, decode=True
                        )
                    compression = "zip"
            elif ".gz" in file:
                with gzip.open(file, "rt") as f:
                    first_line, comments, data = cls._extract_comments(f)
                compression = "gzip"
            else:
                with open(file, "rb") as f:
                    first_line, comments, data, compression = cls._handle_bytes_data(
                        f.read()
                    )

        elif isinstance(file, bytes):

            first_line, comments, data, compression = cls._handle_bytes_data(file)
            file = io.BytesIO(file)

        else:
            return d

        parse_kwargs = {
            "file": file,
            "compression": compression,
            "only_detect_source": only_detect_source,
        }

        if "23andMe" in first_line:
            d.update(_23andMe.parse(**parse_kwargs))
        elif "Ancestry" in first_line:
            d.update(_AncestryDNA.parse(**parse_kwargs))
        elif first_line.startswith("RSID"):
            d.update(_FTDNA.parse(**parse_kwargs))
        elif "famfinder" in first_line:
            d.update(_FTDNAFamFinder.parse(**parse_kwargs))
        elif "MyHeritage" in first_line:
            d.update(_MyHeritage.parse(**parse_kwargs))
        elif "Living DNA" in first_line:
            d.update(_LivingDNA.parse(**parse_kwargs))
        elif "SNP Name\trsID" in first_line or "SNP.Name\tSample.ID" in first_line:
            d.update(_Mapmygenome.parse(header=first_line, **parse_kwargs))
        elif "lineage" in first_line or "snps" in first_line:
            d.update(_snps.parse(comments=comments, **parse_kwargs))
        elif "Chromosome" in first_line:
            d.update(_tellmeGen.parse(resources=resources, **parse_kwargs))
        elif re.match("^#*[ \t]*rsid[, \t]*chr", first_line):
            d.update(_Generic.parse(**parse_kwargs))
        elif re.match("^rs[0-9]*[, \t]{1}[1]", first_line):
            d.update(_Generic.parse(skip=0, **parse_kwargs))
        elif "vcf" in comments.lower() or "##contig" in comments.lower():
            d.update(_VCF.parse(provider="vcf", rsids=rsids, **parse_kwargs))
        elif ("Genes for Good" in comments) | ("PLINK" in comments):
            d.update(_GenesForGood.parse(**parse_kwargs))
        elif "DNA.Land" in comments:
            d.update(_DNALand.parse(**parse_kwargs))
        elif "CODIGO46" in comments:
            d.update(_Codigo46.parse(resources=resources, **parse_kwargs))
        elif "SANO" in comments:
            d.update(_Sano.parse(resources=resources, **parse_kwargs))

        # detect build from comments if build was not already detected by parser
        if not d["build"]:
            d.update({"build": cls._detect_build_from_comments(comments, d["source"])})

        return d

    @classmethod
    def _extract_comments(cls, f, decode=False, include_data=False):
        line = cls._read_line(f, decode)

        first_line = line
        comments = ""
        data = ""

        if first_line.startswith("#"):
            while line.startswith("#"):
                comments += line
                line = cls._read_line(f, decode)
            if include_data:
                while line:
                    data += line
                    line = cls._read_line(f, decode)

        elif first_line.startswith("[Header]"):
            while not line.startswith("[Data]"):
                comments += line
                line = cls._read_line(f, decode)
            # Ignore the [Data] row
            line = cls._read_line(f, decode)
            if include_data:
                while line:
                    data += line
                    line = cls._read_line(f, decode)
        if not data and include_data:
            data = f.read()
            if decode:
                data = data.decode()
        if not isinstance(f, zipfile.ZipExtFile):
            f.seek(0)
        return first_line, comments, data

    @staticmethod
    def _detect_build_from_comments(comments, source):
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

    @classmethod
    def _handle_bytes_data(cls, file, include_data=False):
        if cls.is_zip(file):
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
                    first_line, comments, data = cls._extract_comments(
                        f, decode=True, include_data=include_data
                    )

        elif cls.is_gzip(file):
            compression = "gzip"

            with gzip.open(io.BytesIO(file), "rb") as f:
                first_line, comments, data = cls._extract_comments(
                    f, decode=True, include_data=include_data
                )

        else:
            compression = None
            file = io.BytesIO(file)
            first_line, comments, data = cls._extract_comments(
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

    @classmethod
    def _read_helper(cls, source, parser, only_detect_source):
        """ Generic method to help read files.

        Parameters
        ----------
        source : str
            name of data source
        parser : func
            parsing function, which returns a dict with any of the the following items:

            snps (pandas.DataFrame)
                dataframe of parsed SNPs
            source (str)
                detected source of SNPs
            phased (bool)
                flag indicating if SNPs are phased
            build (int)
                detected build of SNPs from parser

        Returns
        -------
        dict
            dict with any of the following items:

            snps (pandas.DataFrame)
                dataframe of parsed SNPs
            source (str)
                detected source of SNPs
            phased (bool)
                flag indicating if SNPs are phased
            build (int)
                detected build of SNPs from parser

        """
        d = {"source": source}

        if not only_detect_source:
            d.update(parser())

        return d

    @staticmethod
    def _parse_4_cols(
        file, compression, **kwargs,
    ):
        return pd.read_csv(
            file,
            names=["rsid", "chrom", "pos", "genotype"],
            index_col=0,
            dtype=NORMALIZED_DTYPES,
            na_values="--",
            compression=compression,
            comment="#",
            **kwargs,
        )

    @staticmethod
    def _parse_5_cols(file, compression, **kwargs):
        return pd.read_csv(
            file,
            names=["rsid", "chrom", "pos", "allele1", "allele2"],
            index_col=0,
            dtype=TWO_ALLELE_DTYPES,
            engine="c",
            compression=compression,
            comment="#",
            **kwargs,
        )


class _23andMe(Reader):
    @classmethod
    def parse(cls, file, compression, only_detect_source):
        """ Read and parse 23andMe file.

        https://www.23andme.com
        """

        def f():
            return {"snps": cls._parse_4_cols(file, compression, sep="\t")}

        return cls._read_helper("23andMe", f, only_detect_source)


class _FTDNA(Reader):
    @classmethod
    def parse(cls, file, compression, only_detect_source):
        """Read and parse Family Tree DNA (FTDNA) file.

        https://www.familytreedna.com
        """

        def f():
            try:
                df = cls._parse_4_cols(file, compression, sep=",", skiprows=1)
            except ValueError:
                # read files with second header for concatenated data
                if isinstance(file, io.BytesIO):
                    file.seek(0)
                    (*data,) = cls._handle_bytes_data(file.read(), include_data=True)
                    file.seek(0)
                else:
                    with open(file, "rb") as f:
                        (*data,) = cls._handle_bytes_data(f.read(), include_data=True)
                # reconstruct file content from `_handle_bytes_data` results
                lines = data[0] + data[2]
                lines = [line.strip() for line in lines.split("\n")]
                # find index of second header
                second_header_idx = lines.index("RSID,CHROMOSOME,POSITION,RESULT", 1)

                df = cls._parse_4_cols(
                    file, compression, sep=",", skiprows=[0, second_header_idx]
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

                df = cls._parse_4_cols(new_file, compression=None, sep=",", skiprows=1)

            return {"snps": df}

        return cls._read_helper("FTDNA", f, only_detect_source)


class _FTDNAFamFinder(Reader):
    @classmethod
    def parse(cls, file, compression, only_detect_source):
        """ Read and parse Family Tree DNA (FTDNA) "famfinder" file.

        https://www.familytreedna.com
        """

        def f():
            df = cls._parse_5_cols(
                file, compression, sep=",", na_values="-", header=None
            )

            # create genotype column from allele columns
            df["genotype"] = df["allele1"] + df["allele2"]

            # delete allele columns
            # http://stackoverflow.com/a/13485766
            del df["allele1"]
            del df["allele2"]

            return {"snps": df}

        return cls._read_helper("FTDNA", f, only_detect_source)


class _AncestryDNA(Reader):
    @classmethod
    def parse(cls, file, compression, only_detect_source):
        """ Read and parse Ancestry.com file.

        http://www.ancestry.com
        """

        def f():
            df = cls._parse_5_cols(file, compression, sep=r"\s+", na_values=0, header=0)

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

            return {"snps": df}

        return cls._read_helper("AncestryDNA", f, only_detect_source)


class _MyHeritage(Reader):
    @classmethod
    def parse(cls, file, compression, only_detect_source):
        """ Read and parse MyHeritage file.

        https://www.myheritage.com
        """

        def f():
            if isinstance(file, str):
                with open(file, "rb") as f:
                    first_line, comments, data, compression = cls._handle_bytes_data(
                        f.read(), include_data=True
                    )
            else:
                first_line, comments, data, compression = cls._handle_bytes_data(
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

            return {
                "snps": cls._parse_4_cols(
                    io.StringIO(file_string_out.getvalue()),
                    compression=None,
                    header=0,
                    sep=",",
                )
            }

        return cls._read_helper("MyHeritage", f, only_detect_source)


class _LivingDNA(Reader):
    @classmethod
    def parse(cls, file, compression, only_detect_source):
        """ Read and parse LivingDNA file.

        https://livingdna.com/
        """

        def f():
            return {"snps": cls._parse_4_cols(file, compression, sep="\t")}

        return cls._read_helper("LivingDNA", f, only_detect_source)


class _Mapmygenome(Reader):
    @classmethod
    def parse(cls, file, compression, only_detect_source, header):
        """ Read and parse Mapmygenome file.

        https://mapmygenome.in
        """

        def f():
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

            return {"snps": df}

        return cls._read_helper("Mapmygenome", f, only_detect_source)


class _GenesForGood(Reader):
    @classmethod
    def parse(cls, file, compression, only_detect_source):
        """ Read and parse Genes For Good file.

        https://genesforgood.sph.umich.edu/readme/readme1.2.txt
        """

        def f():
            return {"snps": cls._parse_4_cols(file, compression, sep="\t")}

        return cls._read_helper("GenesForGood", f, only_detect_source)


class _GSA(Reader):
    @classmethod
    def _read_gsa_helper(
        cls, file, source, strand, dtypes, resources, only_detect_source, na_values="--"
    ):
        def f():
            gsa_resources = resources.get_gsa_resources()

            if isinstance(file, str):
                try:
                    with open(file, "rb") as f:
                        first_line, comments, data = cls._extract_comments(
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
                        ) = cls._handle_bytes_data(f.read(), include_data=True)
            else:
                first_line, comments, data, compression = cls._handle_bytes_data(
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

            return {"snps": df}

        return cls._read_helper(source, f, only_detect_source)


class _tellmeGen(Reader):
    @classmethod
    def parse(cls, file, compression, only_detect_source, resources):
        """ Read and parse tellmeGen files.

        https://www.tellmegen.com/
        """

        def f():
            gsa_resources = resources.get_gsa_resources()

            df = cls._parse_4_cols(file, compression, sep="\t", skiprows=1)
            df.rename(index=gsa_resources["rsid_map"], inplace=True)
            return {"snps": df}

        return cls._read_helper("tellmeGen", f, only_detect_source)


class _Codigo46(_GSA):
    @classmethod
    def parse(cls, file, compression, only_detect_source, resources):
        """ Read and parse Codigo46 files.

        https://codigo46.com.mx
        """
        return cls._read_gsa_helper(
            file, "Codigo46", "Plus", {}, resources, only_detect_source
        )


class _Sano(_GSA):
    @classmethod
    def parse(cls, file, compression, only_detect_source, resources):
        """ Read and parse Sano Genetics files.

        https://sanogenetics.com
        """
        dtype = {"Chr": object, "Position": np.uint32}
        return cls._read_gsa_helper(
            file,
            "Sano",
            "Forward",
            dtype,
            resources,
            only_detect_source,
            na_values="-",
        )


class _DNALand(Reader):
    @classmethod
    def parse(cls, file, compression, only_detect_source):
        """ Read and parse DNA.land files.

        https://dna.land/
        """

        def f():
            return {"snps": cls._parse_4_cols(file, compression, sep="\t")}

        return cls._read_helper("DNA.Land", f, only_detect_source)


class _snps(Reader):
    @classmethod
    def parse(cls, file, compression, only_detect_source, comments):
        """ Read and parse CSV file generated by ``snps``.

        https://pypi.org/project/snps/
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

        def f():
            def parse_csv(sep):
                return cls._parse_4_cols(file, compression, sep=sep, header=0)

            try:
                return {"snps": parse_csv(","), "phased": phased, "build": build}
            except pd.errors.ParserError:
                if isinstance(file, io.BufferedIOBase):
                    file.seek(0)

                return {"snps": parse_csv("\t"), "phased": phased, "build": build}

        return cls._read_helper(source, f, only_detect_source)


class _Generic(Reader):
    @classmethod
    def parse(cls, file, compression, only_detect_source, skip=1):
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
        """

        def f():
            def parse(sep):
                return cls._parse_4_cols(file, compression, sep=sep, skiprows=skip)

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

                    df = cls._parse_4_cols(
                        file,
                        compression,
                        sep=None,
                        skiprows=skip,
                        engine="python",
                        usecols=[0, 1, 2, 3],
                    )
            return {"snps": df}

        return cls._read_helper("generic", f, only_detect_source)


class _VCF(Reader):
    @classmethod
    def parse(cls, file, compression, provider, only_detect_source, rsids=()):
        """ Read and parse VCF file.

        See ``Reader`` notes for more details on VCF parsing capability.
        """

        def f():
            if not isinstance(file, io.BytesIO):
                with open(file, "rb") as f:
                    df, phased = cls._parse_vcf(f, rsids)
            else:
                df, phased = cls._parse_vcf(file, rsids)

            return {"snps": df, "phased": phased}

        return cls._read_helper(provider, f, only_detect_source)

    @classmethod
    def _parse_vcf(cls, buffer, rsids):
        rows = []
        phased = True
        first_four_bytes = buffer.read(4)
        buffer.seek(0)

        if cls.is_gzip(first_four_bytes):
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
