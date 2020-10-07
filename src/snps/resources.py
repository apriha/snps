""" Class for downloading and loading required external resources.

References
----------
1. International Human Genome Sequencing Consortium. Initial sequencing and
   analysis of the human genome. Nature. 2001 Feb 15;409(6822):860-921.
   http://dx.doi.org/10.1038/35057062
2. hg19 (GRCh37): Hiram Clawson, Brooke Rhead, Pauline Fujita, Ann Zweig, Katrina
   Learned, Donna Karolchik and Robert Kuhn, https://genome.ucsc.edu/cgi-bin/hgGateway?db=hg19
3. Yates et. al. (doi:10.1093/bioinformatics/btu613),
   `<http://europepmc.org/search/?query=DOI:10.1093/bioinformatics/btu613>`_
4. Zerbino et. al. (doi.org/10.1093/nar/gkx1098), https://doi.org/10.1093/nar/gkx1098

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

import gzip
import hashlib
import itertools
import json
import logging
import os
import socket
import tarfile
import tempfile
import urllib.error
import urllib.request
import zipfile

from atomicwrites import atomic_write
import numpy as np

from snps.ensembl import EnsemblRestClient
from snps.utils import create_dir, Singleton

logger = logging.getLogger(__name__)


class Resources(metaclass=Singleton):
    """ Object used to manage resources required by `snps`. """

    def __init__(self, resources_dir="resources"):
        """ Initialize a ``Resources`` object.

        Parameters
        ----------
        resources_dir : str
            name / path of resources directory
        """
        self._resources_dir = os.path.abspath(resources_dir)
        self._ensembl_rest_client = EnsemblRestClient()
        self._reference_sequences = {}
        self._gsa_resources = {}
        self._opensnp_datadump_filenames = []

    def get_reference_sequences(
        self,
        assembly="GRCh37",
        chroms=(
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
        ),
    ):
        """ Get Homo sapiens reference sequences for `chroms` of `assembly`.

        Notes
        -----
        This function can download over 800MB of data for each assembly.

        Parameters
        ----------
        assembly : {'NCBI36', 'GRCh37', 'GRCh38'}
            reference sequence assembly
        chroms : list of str
            reference sequence chromosomes

        Returns
        -------
        dict
            dict of ReferenceSequence, else {}
        """
        valid_assemblies = ["NCBI36", "GRCh37", "GRCh38"]

        if assembly not in valid_assemblies:
            logger.warning("Invalid assembly")
            return {}

        if not self._reference_chroms_available(assembly, chroms):
            self._reference_sequences[assembly] = self._create_reference_sequences(
                *self._get_paths_reference_sequences(assembly=assembly, chroms=chroms)
            )

        return self._reference_sequences[assembly]

    def _reference_chroms_available(self, assembly, chroms):
        if assembly in self._reference_sequences:
            for chrom in chroms:
                if chrom not in self._reference_sequences[assembly]:
                    return False
            return True
        else:
            return False

    def get_assembly_mapping_data(self, source_assembly, target_assembly):
        """ Get assembly mapping data.

        Parameters
        ----------
        source_assembly : {'NCBI36', 'GRCh37', 'GRCh38'}
            assembly to remap from
        target_assembly : {'NCBI36', 'GRCh37', 'GRCh38'}
            assembly to remap to

        Returns
        -------
        dict
            dict of json assembly mapping data if loading was successful, else {}
        """
        return self._load_assembly_mapping_data(
            self._get_path_assembly_mapping_data(source_assembly, target_assembly)
        )

    def download_example_datasets(self):
        """ Download example datasets from `openSNP <https://opensnp.org>`_.

        Per openSNP, "the data is donated into the public domain using `CC0 1.0
        <http://creativecommons.org/publicdomain/zero/1.0/>`_."

        Returns
        -------
        paths : list of str or empty str
            paths to example datasets

        References
        ----------
        1. Greshake B, Bayer PE, Rausch H, Reda J (2014), "openSNP-A Crowdsourced Web Resource
           for Personal Genomics," PLOS ONE, 9(3): e89204,
           https://doi.org/10.1371/journal.pone.0089204
        """
        paths = []
        paths.append(
            self._download_file(
                "https://opensnp.org/data/662.23andme.340",
                "662.23andme.340.txt.gz",
                compress=True,
            )
        )
        paths.append(
            self._download_file(
                "https://opensnp.org/data/662.ftdna-illumina.341",
                "662.ftdna-illumina.341.csv.gz",
                compress=True,
            )
        )

        return paths

    def get_all_resources(self):
        """ Get / download all resources used throughout `snps`.

        Notes
        -----
        This function does not download reference sequences and the openSNP datadump,
        due to their large sizes.

        Returns
        -------
        dict
            dict of resources
        """
        resources = {}
        for source, target in itertools.permutations(["NCBI36", "GRCh37", "GRCh38"], 2):
            resources[source + "_" + target] = self.get_assembly_mapping_data(
                source, target
            )
        resources["gsa_resources"] = self.get_gsa_resources()
        return resources

    def get_all_reference_sequences(self, **kwargs):
        """ Get Homo sapiens reference sequences for Builds 36, 37, and 38 from Ensembl.

        Notes
        -----
        This function can download over 2.5GB of data.

        Returns
        -------
        dict
            dict of ReferenceSequence, else {}
        """
        for assembly in ("NCBI36", "GRCh37", "GRCh38"):
            self.get_reference_sequences(assembly=assembly, **kwargs)
        return self._reference_sequences

    def get_gsa_resources(self):
        """ Get resources for reading Global Screening Array files.

        https://support.illumina.com/downloads/infinium-global-screening-array-v2-0-product-files.html

        Returns
        -------
        dict
        """
        if not self._gsa_resources:
            self._gsa_resources = self._load_gsa_resources(
                self._get_path_gsa_rsid_map(), self._get_path_gsa_chrpos_map()
            )
        return self._gsa_resources

    def get_opensnp_datadump_filenames(self):
        """ Get filenames internal to the `openSNP <https://opensnp.org>`_ datadump zip.

        Per openSNP, "the data is donated into the public domain using `CC0 1.0
        <http://creativecommons.org/publicdomain/zero/1.0/>`_."

        Notes
        -----
        This function can download over 27GB of data. If the download is not successful,
        try using a different tool like `wget` or `curl` to download the file and move it
        to the resources directory (see `_get_path_opensnp_datadump`).

        Returns
        -------
        filenames : list of str
            filenames internal to the openSNP datadump

        References
        ----------
        1. Greshake B, Bayer PE, Rausch H, Reda J (2014), "openSNP-A Crowdsourced Web Resource
           for Personal Genomics," PLOS ONE, 9(3): e89204,
           https://doi.org/10.1371/journal.pone.0089204
        """
        if not self._opensnp_datadump_filenames:
            self._opensnp_datadump_filenames = self._get_opensnp_datadump_filenames(
                self._get_path_opensnp_datadump()
            )
        return self._opensnp_datadump_filenames

    def load_opensnp_datadump_file(self, filename):
        """ Load the specified file from the openSNP datadump.

        Per openSNP, "the data is donated into the public domain using `CC0 1.0
        <http://creativecommons.org/publicdomain/zero/1.0/>`_."

        Parameters
        ----------
        filename : str
            filename internal to the openSNP datadump

        Returns
        -------
        bytes
            content of specified file internal to the openSNP datadump

        References
        ----------
        1. Greshake B, Bayer PE, Rausch H, Reda J (2014), "openSNP-A Crowdsourced Web Resource
           for Personal Genomics," PLOS ONE, 9(3): e89204,
           https://doi.org/10.1371/journal.pone.0089204
        """
        if self._get_path_opensnp_datadump():
            with zipfile.ZipFile(self._get_path_opensnp_datadump()) as z:
                with z.open(filename, "r") as f:
                    return f.read()
        else:
            return bytes()

    @staticmethod
    def _get_opensnp_datadump_filenames(filename):
        """ Get list of filenames internal to the openSNP datadump zip.

        Parameters
        ----------
        filename : str
            path to openSNP datadump

        Returns
        -------
        filenames : list of str
            filenames internal to the openSNP datadump
        """
        if filename:
            with zipfile.ZipFile(filename) as z:
                return z.namelist()
        else:
            return []

    @staticmethod
    def _write_data_to_gzip(f, data):
        """ Write `data` to `f` in `gzip` format.

        Parameters
        ----------
        f : file object opened with `mode="wb"`
        data : `bytes` object
        """
        with gzip.open(f, "wb") as f_gzip:
            f_gzip.write(data)

    @staticmethod
    def _load_assembly_mapping_data(filename):
        """ Load assembly mapping data.

        Parameters
        ----------
        filename : str
            path to compressed archive with assembly mapping data

        Returns
        -------
        assembly_mapping_data : dict
            dict of assembly maps

        Notes
        -----
        Keys of returned dict are chromosomes and values are the corresponding assembly map.
        """
        assembly_mapping_data = {}

        with tarfile.open(filename, "r") as tar:
            # http://stackoverflow.com/a/2018576
            for member in tar.getmembers():
                if ".json" in member.name:
                    with tar.extractfile(member) as tar_file:
                        tar_bytes = tar_file.read()
                    # https://stackoverflow.com/a/42683509/4727627
                    assembly_mapping_data[member.name.split(".")[0]] = json.loads(
                        tar_bytes.decode("utf-8")
                    )

        return assembly_mapping_data

    def _get_paths_reference_sequences(
        self, sub_dir="fasta", assembly="GRCh37", chroms=()
    ):
        """ Get local paths to Homo sapiens reference sequences from Ensembl.

        Notes
        -----
        This function can download over 800MB of data for each assembly.

        Parameters
        ----------
        sub_dir : str
            directory under resources to store reference sequence data
        assembly : {'NCBI36', 'GRCh37', 'GRCh38'}
            reference sequence assembly
        chroms : list of str
            reference sequence chromosomes

        Returns
        -------
        assembly : str
            reference sequence assembly
        chroms : list of str
            reference sequence chromosomes
        urls : list of str
            urls to Ensembl reference sequences
        paths : list of str
            paths to local reference sequences

        References
        ----------
        1. Daniel R. Zerbino, Premanand Achuthan, Wasiu Akanni, M. Ridwan Amode,
           Daniel Barrell, Jyothish Bhai, Konstantinos Billis, Carla Cummins, Astrid Gall,
           Carlos García Giro´n, Laurent Gil, Leo Gordon, Leanne Haggerty, Erin Haskell,
           Thibaut Hourlier, Osagie G. Izuogu, Sophie H. Janacek, Thomas Juettemann,
           Jimmy Kiang To, Matthew R. Laird, Ilias Lavidas, Zhicheng Liu, Jane E. Loveland,
           Thomas Maurel, William McLaren, Benjamin Moore, Jonathan Mudge, Daniel N. Murphy,
           Victoria Newman, Michael Nuhn, Denye Ogeh, Chuang Kee Ong, Anne Parker,
           Mateus Patricio, Harpreet Singh Riat, Helen Schuilenburg, Dan Sheppard,
           Helen Sparrow, Kieron Taylor, Anja Thormann, Alessandro Vullo, Brandon Walts,
           Amonida Zadissa, Adam Frankish, Sarah E. Hunt, Myrto Kostadima, Nicholas Langridge,
           Fergal J. Martin, Matthieu Muffato, Emily Perry, Magali Ruffier, Dan M. Staines,
           Stephen J. Trevanion, Bronwen L. Aken, Fiona Cunningham, Andrew Yates, Paul Flicek
           Ensembl 2018.
           PubMed PMID: 29155950.
           doi:10.1093/nar/gkx1098
        2. NCBI 36, Oct 2005, Ensembl release 54, Database version: 54.36p
        3. GRCh37.p13 (Genome Reference Consortium Human Reference 37),
           INSDC Assembly GCA_000001405.14, Feb 2009, Ensembl GRCh37 release 96, Database
           version: 96.37
        4. GRCh38.p12 (Genome Reference Consortium Human Build 38),
           INSDC Assembly GCA_000001405.27, Dec 2013, Ensembl release 96, Database
           version: 96.38
        """
        release = ""

        # https://www.biostars.org/p/374149/#374219
        if assembly == "GRCh37":
            base = "ftp://ftp.ensembl.org/pub/grch37/release-96/fasta/homo_sapiens/dna/"
        elif assembly == "NCBI36":
            base = "ftp://ftp.ensembl.org/pub/release-54/fasta/homo_sapiens/dna/"
            release = "54."
        elif assembly == "GRCh38":
            base = "ftp://ftp.ensembl.org/pub/release-96/fasta/homo_sapiens/dna/"
        else:
            return ("", [], [], [])

        filenames = [
            f"Homo_sapiens.{assembly}.{release}dna.chromosome.{chrom}.fa.gz"
            for chrom in chroms
        ]

        urls = [f"{base}{filename}" for filename in filenames]

        local_filenames = [
            f"{sub_dir}{os.sep}{assembly}{os.sep}{filename}" for filename in filenames
        ]

        return (
            assembly,
            chroms,
            urls,
            list(map(self._download_file, urls, local_filenames)),
        )

    def _create_reference_sequences(self, assembly, chroms, urls, paths):
        # https://samtools.github.io/hts-specs/VCFv4.2.pdf
        seqs = {}

        for i, path in enumerate(paths):
            if not path:
                continue

            d = {}
            d["ID"] = chroms[i]
            d["url"] = urls[i]
            d["path"] = os.path.relpath(path)
            d["assembly"] = assembly
            d["species"] = "Homo sapiens"
            d["taxonomy"] = "x"
            seqs[chroms[i]] = ReferenceSequence(**d)

        return seqs

    def _get_path_assembly_mapping_data(
        self, source_assembly, target_assembly, retries=10
    ):
        """ Get local path to assembly mapping data, downloading if necessary.

        Parameters
        ----------
        source_assembly : {'NCBI36', 'GRCh37', 'GRCh38'}
            assembly to remap from
        target_assembly : {'NCBI36', 'GRCh37', 'GRCh38'}
            assembly to remap to
        retries : int
            number of retries per chromosome to download assembly mapping data

        Returns
        -------
        str
            path to <source_assembly>_<target_assembly>.tar.gz

        References
        ----------
        1. Ensembl, Assembly Information Endpoint,
           https://rest.ensembl.org/documentation/info/assembly_info
        2. Ensembl, Assembly Map Endpoint,
           http://rest.ensembl.org/documentation/info/assembly_map

        """

        if not create_dir(self._resources_dir):
            return ""

        chroms = [str(i) for i in range(1, 23)]
        chroms.extend(["X", "Y", "MT"])

        assembly_mapping_data = source_assembly + "_" + target_assembly
        destination = os.path.join(
            self._resources_dir, assembly_mapping_data + ".tar.gz"
        )

        if not os.path.exists(destination):
            logger.info(f"Downloading {os.path.relpath(destination)}")

            self._download_assembly_mapping_data(
                destination, chroms, source_assembly, target_assembly, retries
            )

        return destination

    def _download_assembly_mapping_data(
        self, destination, chroms, source_assembly, target_assembly, retries
    ):
        with atomic_write(destination, mode="wb", overwrite=True) as f:
            with tarfile.open(fileobj=f, mode="w:gz") as out_tar:
                for chrom in chroms:
                    file = chrom + ".json"

                    map_endpoint = (
                        f"/map/human/{source_assembly}/{chrom}/{target_assembly}?"
                    )

                    # get assembly mapping data
                    response = None
                    retry = 0
                    while response is None and retry < retries:
                        response = self._ensembl_rest_client.perform_rest_action(
                            map_endpoint
                        )
                        retry += 1

                    if response is not None:
                        # open temp file, save json response to file, close temp file
                        with tempfile.NamedTemporaryFile(
                            delete=False, mode="w"
                        ) as f_tmp:
                            json.dump(response, f_tmp)

                        # add temp file to archive
                        out_tar.add(f_tmp.name, arcname=file)

                        # remove temp file
                        os.remove(f_tmp.name)

    def _load_gsa_resources(self, rsid_map, chrpos_map):
        d = {}

        with gzip.open(rsid_map, "rb") as f:
            gsa_rsid_map = f.read().decode("utf-8")

        d["rsid_map"] = dict(
            (x.split("\t")[0], x.split("\t")[1]) for x in gsa_rsid_map.split("\n")[:-1]
        )

        with gzip.open(chrpos_map, "rb") as f:
            gsa_chrpos_map = f.read().decode("utf-8")

        d["chrpos_map"] = dict(
            (x.split("\t")[0], x.split("\t")[1] + ":" + x.split("\t")[2])
            for x in gsa_chrpos_map.split("\n")[:-1]
        )

        return d

    def _get_path_gsa_rsid_map(self):
        return self._download_file(
            "https://sano-public.s3.eu-west-2.amazonaws.com/gsa_rsid_map.txt.gz",
            "gsa_rsid_map.txt.gz",
        )

    def _get_path_gsa_chrpos_map(self):
        return self._download_file(
            "https://sano-public.s3.eu-west-2.amazonaws.com/gsa_chrpos_map.txt.gz",
            "gsa_chrpos_map.txt.gz",
        )

    def _get_path_opensnp_datadump(self):
        return self._download_file(
            "https://opensnp.org/data/zip/opensnp_datadump.current.zip",
            "opensnp_datadump.current.zip",
        )

    def _download_file(self, url, filename, compress=False, timeout=30):
        """ Download a file to the resources folder.

        Download data from `url`, save as `filename`, and optionally compress with gzip.

        Parameters
        ----------
        url : str
            URL to download data from
        filename : str
            name of file to save; if compress, ensure '.gz' is appended
        compress : bool
            compress with gzip
        timeout : int
            seconds for timeout of download request

        Returns
        -------
        str
            path to downloaded file, empty str if error
        """
        if compress and filename[-3:] != ".gz":
            filename += ".gz"

        destination = os.path.join(self._resources_dir, filename)

        if not create_dir(os.path.relpath(os.path.dirname(destination))):
            return ""

        if not os.path.exists(destination):
            try:
                # get file if it hasn't already been downloaded
                # http://stackoverflow.com/a/7244263
                with urllib.request.urlopen(
                    url, timeout=timeout
                ) as response, atomic_write(destination, mode="wb") as f:
                    self._print_download_msg(destination)
                    data = response.read()  # a `bytes` object

                    if compress:
                        self._write_data_to_gzip(f, data)
                    else:
                        f.write(data)
            except urllib.error.URLError as err:
                logger.warning(err)
                destination = ""
                # try HTTP if an FTP error occurred
                if "ftp://" in url:
                    destination = self._download_file(
                        url.replace("ftp://", "http://"),
                        filename,
                        compress=compress,
                        timeout=timeout,
                    )
            except socket.timeout:
                logger.warning(f"Timeout downloading {url}")
                destination = ""

        return destination

    @staticmethod
    def _print_download_msg(path):
        """ Print download message.

        Parameters
        ----------
        path : str
            path to file being downloaded
        """
        logger.info(f"Downloading {os.path.relpath(path)}")


class ReferenceSequence:
    """ Object used to represent and interact with a reference sequence. """

    def __init__(self, ID="", url="", path="", assembly="", species="", taxonomy=""):
        """ Initialize a ``ReferenceSequence`` object.

        Parameters
        ----------
        ID : str
            reference sequence chromosome
        url : str
            url to Ensembl reference sequence
        path : str
            path to local reference sequence
        assembly : str
            reference sequence assembly (e.g., "GRCh37")
        species : str
            reference sequence species
        taxonomy : str
            reference sequence taxonomy

        References
        ----------
        1. The Variant Call Format (VCF) Version 4.2 Specification, 8 Mar 2019,
           https://samtools.github.io/hts-specs/VCFv4.2.pdf
        """
        self._ID = ID
        self._url = url
        self._path = path
        self._assembly = assembly
        self._species = species
        self._taxonomy = taxonomy
        self._sequence = np.array([], dtype=np.uint8)
        self._md5 = ""
        self._start = 0
        self._end = 0
        self._length = 0

    def __repr__(self):
        return f"ReferenceSequence(assembly={self._assembly!r}, ID={self._ID!r})"

    @property
    def ID(self):
        """ Get reference sequence chromosome.

        Returns
        -------
        str
        """
        return self._ID

    @property
    def chrom(self):
        """ Get reference sequence chromosome.

        Returns
        -------
        str
        """
        return self._ID

    @property
    def url(self):
        """ Get URL to Ensembl reference sequence.

        Returns
        -------
        str
        """
        return self._url

    @property
    def path(self):
        """ Get path to local reference sequence.

        Returns
        -------
        str
        """
        return self._path

    @property
    def assembly(self):
        """ Get reference sequence assembly.

        Returns
        -------
        str
        """
        return self._assembly

    @property
    def build(self):
        """ Get reference sequence build.

        Returns
        -------
        str
            e.g., "B37"
        """
        return f"B{self._assembly[-2:]}"

    @property
    def species(self):
        """ Get reference sequence species.

        Returns
        -------
        str
        """
        return self._species

    @property
    def taxonomy(self):
        """ Get reference sequence taxonomy.

        Returns
        -------
        str
        """
        return self._taxonomy

    @property
    def sequence(self):
        """ Get reference sequence.

        Returns
        -------
        np.array(dtype=np.uint8)
        """
        self._load_sequence()
        return self._sequence

    @property
    def md5(self):
        """ Get reference sequence MD5 hash.

        Returns
        -------
        str
        """
        self._load_sequence()
        return self._md5

    @property
    def start(self):
        """ Get reference sequence start position (1-based).

        Returns
        -------
        int
        """
        self._load_sequence()
        return self._start

    @property
    def end(self):
        """ Get reference sequence end position (1-based).

        Returns
        -------
        int
        """
        self._load_sequence()
        return self._end

    @property
    def length(self):
        """ Get reference sequence length.

        Returns
        -------
        int
        """
        self._load_sequence()
        return self._sequence.size

    def clear(self):
        """ Clear reference sequence. """
        self._sequence = np.array([], dtype=np.uint8)
        self._md5 = ""
        self._start = 0
        self._end = 0
        self._length = 0

    def _load_sequence(self):
        if not self._sequence.size:
            # decompress and read file
            with gzip.open(self._path, "rb") as f:
                data = f.read()

            # convert bytes to str
            data = str(data, encoding="utf-8", errors="strict")

            data = data.split("\n")

            self._start, self._end = self._parse_first_line(data[0])

            # convert str (FASTA sequence) to bytes
            data = bytearray("".join(data[1:]), encoding="utf-8", errors="strict")

            # get MD5 of FASTA sequence
            self._md5 = hashlib.md5(data).hexdigest()

            # store FASTA sequence as `np.uint8` array
            self._sequence = np.array(data, dtype=np.uint8)

    def _parse_first_line(self, first_line):
        items = first_line.split(":")
        return (
            int(items[items.index(self._ID) + 1]),
            int(items[items.index(self._ID) + 2]),
        )
