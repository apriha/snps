"""Class for reading SNPs."""

import warnings

from snps.data import get_empty_snps_dataframe
from snps.exceptions import UnsupportedFormatError
from snps.io.readers import formats
from snps.io.readers._common import detect_build_from_comments
from snps.io.readers._registry import ParseInput, detect
from snps.io.readers._source import InputSource, can_read
from snps.resources import get_default_provider


class Reader:
    """Class for reading and parsing raw data / genotype files."""

    def __init__(self, file="", only_detect_source=False, resources=None, rsids=()):
        """Initialize a `Reader`.

        Parameters
        ----------
        file : str, os.PathLike, bytes or binary file object
            file to read, as for ``SNPs``
        only_detect_source : bool
            only detect the source of the data
        resources : Resources
            instance of Resources
        rsids : tuple, optional
            rsids to extract if loading a VCF file

        """
        self._file = file
        self._only_detect_source = only_detect_source
        self._resources = resources if resources is not None else get_default_provider()
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
        d = {
            "snps": get_empty_snps_dataframe(),
            "source": "",
            "phased": False,
            "build": 0,
        }

        # peek into files to determine the data format
        if not can_read(self._file):
            return d
        try:
            input_source = InputSource(self._file)
        except UnsupportedFormatError:
            return d

        with input_source:
            header = input_source.header()
            comments = header.comments
            fmt = detect(header, formats())

            if fmt is not None:
                source = fmt.source(header) if fmt.source is not None else fmt.name
                if source is None:
                    return d

                if self._only_detect_source:
                    d["source"] = source
                else:
                    parsed = fmt.parse(
                        ParseInput(input_source, header, self._resources, self._rsids)
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
