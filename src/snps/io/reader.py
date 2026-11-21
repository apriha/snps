"""Class for reading SNPs."""

import gzip
import io
import os
import warnings
import zipfile

from snps.data import get_empty_snps_dataframe
from snps.io.readers import formats
from snps.io.readers._common import (
    detect_build_from_comments,
    extract_comments,
    handle_bytes_data,
)
from snps.io.readers._registry import Header, ParseInput, detect


class Reader:
    """Class for reading and parsing raw data / genotype files."""

    def __init__(self, file="", only_detect_source=False, resources=None, rsids=()):
        """Initialize a `Reader`.

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
        """Read and parse a raw data / genotype file.

        Returns
        -------
        dict
            dict with the following items:

            snps (pandas.DataFrame)
                dataframe of parsed SNPs; may carry a boolean ``phased`` column
                saying whether each SNP is phased
            source (str)
                detected source of SNPs
            phased (bool)
                flag indicating if the whole file is phased
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
                        first_line, comments, data = extract_comments(f, decode=True)
                    compression = "zip"
            elif ".gz" in file:
                with gzip.open(file, "rt") as f:
                    first_line, comments, data = extract_comments(f)
                compression = "gzip"
            else:
                with open(file, "rb") as f:
                    first_line, comments, data, compression = handle_bytes_data(
                        f.read()
                    )

        elif isinstance(file, bytes):
            first_line, comments, data, compression = handle_bytes_data(file)
            file = io.BytesIO(file)

        else:
            return d

        header = Header(first_line, comments)
        fmt = detect(header, formats())

        if fmt is not None:
            source = fmt.source(header) if fmt.source is not None else fmt.name
            if source is None:
                return d

            if self._only_detect_source:
                d["source"] = source
            else:
                parsed = fmt.parse(
                    ParseInput(file, compression, header, self._resources, self._rsids)
                )
                if parsed is None:
                    return d

                df = parsed.frame
                # a parser that reports the whole file as phased yields phased calls
                if parsed.phased and "phased" not in df.columns:
                    df = df.assign(phased=df["genotype"].notna())

                d = {
                    "snps": df,
                    "source": source,
                    "phased": parsed.phased,
                    "build": parsed.build,
                }

        # otherwise use a build named in the comments
        if not d["build"]:
            d.update({"build": detect_build_from_comments(comments, d["source"])})

        return d

    @classmethod
    def read_file(cls, file, only_detect_source, resources, rsids):
        warnings.warn(
            "This method will be removed in a future release.", DeprecationWarning
        )
        r = cls(file, only_detect_source, resources, rsids)
        return r.read()
