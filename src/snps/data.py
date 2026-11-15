"""Immutable value type for normalized genotype data."""

from __future__ import annotations

import re
from collections import Counter
from dataclasses import dataclass, field, replace

import pandas as pd
from pandas.api.types import CategoricalDtype

from snps.build_constants import BUILD_MARKER_SNPS
from snps.constants import NORMALIZED_DTYPES

_NUMBER_OR_TEXT = re.compile("([0-9]+)|([^0-9]+)")


def _natural_sort_key(name: str) -> list[tuple[int, int | str]]:
    """Key for sorting names in natural order: 1, 2, 10, X (not 1, 10, 2, X).

    A name is split into numbers and text; numbers sort by value and before text.
    """
    return [
        (0, int(number)) if number else (1, text.lower())
        for number, text in _NUMBER_OR_TEXT.findall(name)
    ]


def _consecutive_runs(numbers: list[int]) -> list[str]:
    """Format each run of integers increasing by one as ``first-last``, or ``first``."""
    runs: list[list[int]] = []
    for n in numbers:
        if runs and n == runs[-1][-1] + 1:
            runs[-1].append(n)
        else:
            runs.append([n])
    return [f"{run[0]}-{run[-1]}" if len(run) > 1 else f"{run[0]}" for run in runs]


def get_empty_snps_dataframe() -> pd.DataFrame:
    """Get empty dataframe normalized for usage with ``snps``.

    Returns
    -------
    pandas.DataFrame
        empty normalized ``snps`` dataframe (``rsid`` index; ``chrom``, ``pos``,
        ``genotype`` columns)
    """
    df = pd.DataFrame(columns=["rsid", "chrom", "pos", "genotype"])
    df = df.astype(NORMALIZED_DTYPES)
    df = df.set_index("rsid")
    return df


@dataclass(frozen=True, eq=False)
class SNPData:
    """Immutable container for a normalized SNP frame and its core metadata.

    ``SNPData`` is a frozen dataclass: its fields cannot be rebound after
    construction. Transformations (:meth:`sort`, :meth:`with_frame`) return a *new*
    ``SNPData`` rather than modifying this one.

    The SNPs are held in one DataFrame, by reference and not copied, with an internal
    boolean ``phased`` column beside the columns of :attr:`frame`. :attr:`frame` and
    :attr:`phased` share its data; under copy-on-write, writing to them leaves the held
    data unchanged.

    Instances compare equal only if they are the same object, and the frame is
    omitted from ``repr``.

    Parameters
    ----------
    _frame : pandas.DataFrame
        normalized ``snps`` dataframe (``rsid`` index; ``chrom``, ``pos``, ``genotype``
        columns), optionally with the internal boolean ``phased`` column (read it
        through :attr:`phased`); without that column, every SNP is unphased
    build : int
        build / assembly of the SNPs (36, 37, 38), or 0 if unknown
    build_detected : bool
        whether the build was detected (as opposed to assumed)
    build_original : int
        build of the SNPs before any remapping
    source : tuple of str
        data source(s) the SNPs were parsed from
    """

    _frame: pd.DataFrame = field(default_factory=get_empty_snps_dataframe, repr=False)
    build: int = 0
    build_detected: bool = False
    build_original: int = 0
    source: tuple[str, ...] = ()

    def __post_init__(self) -> None:
        # a frame without a phased column holds unphased SNPs
        if "phased" not in self._frame.columns:
            object.__setattr__(self, "_frame", self._frame.assign(phased=False))

    # -- frame and phased -----------------------------------------------------

    @property
    def frame(self) -> pd.DataFrame:
        """Normalized SNPs.

        Returns
        -------
        pandas.DataFrame
            normalized ``snps`` dataframe, per :attr:`snps.snps.SNPs.snps`
        """
        return self._frame[["chrom", "pos", "genotype"]]

    @property
    def phased(self) -> pd.Series[bool]:
        """Whether each SNP is phased.

        A no-call is never phased.

        Returns
        -------
        pandas.Series
            boolean, indexed like :attr:`frame`
        """
        return self._frame["phased"]

    # -- transformations (return a new SNPData) -------------------------------

    def with_frame(self, frame: pd.DataFrame) -> SNPData:
        """Return a copy of this ``SNPData`` with a replaced frame.

        Parameters
        ----------
        frame : pandas.DataFrame
            normalized ``snps`` dataframe, optionally with the internal boolean
            ``phased`` column; without that column, every SNP is unphased

        Returns
        -------
        SNPData
        """
        return replace(self, _frame=frame)

    def sort(self) -> SNPData:
        """Return a copy sorted by ordered chromosome list and position.

        Chromosomes are sorted in natural order, with ``PAR`` and ``MT`` moved to the
        end. The input frame is not mutated.

        Returns
        -------
        SNPData
        """
        sorted_list = sorted(
            (str(x) for x in self._frame["chrom"].unique()), key=_natural_sort_key
        )

        # move PAR and MT to the end of the dataframe
        if "PAR" in sorted_list:
            sorted_list.remove("PAR")
            sorted_list.append("PAR")

        if "MT" in sorted_list:
            sorted_list.remove("MT")
            sorted_list.append("MT")

        # convert chrom column to category for sorting; the shallow copy shares
        # column data with self._frame, and replacing its chrom column leaves
        # self._frame unchanged
        snps = self._frame.copy(deep=False)
        snps["chrom"] = snps["chrom"].astype(
            CategoricalDtype(categories=sorted_list, ordered=True)
        )

        # sort based on ordered chromosome list and position
        snps = snps.sort_values(["chrom", "pos"])

        # convert chromosome back to object
        snps["chrom"] = snps["chrom"].astype(object)

        return replace(self, _frame=snps)

    # -- derived scalars ------------------------------------------------------

    @property
    def valid(self) -> bool:
        """Whether this ``SNPData`` holds any SNPs.

        Returns
        -------
        bool
        """
        return not self._frame.empty

    @property
    def count(self) -> int:
        """Count of SNPs.

        Returns
        -------
        int
        """
        return self.get_count()

    def get_count(self, chrom: str = "") -> int:
        """Count of SNPs, optionally for a single chromosome.

        Parameters
        ----------
        chrom : str, optional
            chromosome (e.g., "1", "X", "MT")

        Returns
        -------
        int
        """
        return len(self._filter(chrom))

    @property
    def assembly(self) -> str:
        """Assembly of SNPs.

        Returns
        -------
        str
        """
        if self.build == 37:
            return "GRCh37"
        elif self.build == 36:
            return "NCBI36"
        elif self.build == 38:
            return "GRCh38"
        else:
            return ""

    @property
    def chromosomes(self) -> list[str]:
        """Chromosomes of SNPs.

        Returns
        -------
        list of str
            e.g., ['1', '2', '3', 'MT']; empty list if no chromosomes
        """
        if not self._frame.empty:
            return list(pd.unique(self._frame["chrom"]))
        else:
            return []

    @property
    def chromosomes_summary(self) -> str:
        """Summary of the chromosomes of SNPs.

        Returns
        -------
        str
            human-readable listing (e.g., '1-3, MT'); empty str if no chromosomes
        """
        if not self._frame.empty:
            chroms = list(pd.unique(self._frame["chrom"]))

            int_chroms: list[int] = []
            str_chroms: list[str] = []
            for chrom in chroms:
                if chrom.isascii() and chrom.isdigit():
                    int_chroms.append(int(chrom))
                else:
                    str_chroms.append(chrom)

            int_summary = ", ".join(_consecutive_runs(int_chroms))
            str_summary = ", ".join(str_chroms)

            if int_summary and str_summary:
                return f"{int_summary}, {str_summary}"
            return int_summary + str_summary
        else:
            return ""

    @property
    def sex(self) -> str:
        """Sex derived from SNPs.

        Returns
        -------
        str
            'Male' or 'Female' if detected, else empty str
        """
        sex = self.determine_sex(chrom="X")
        if not sex:
            sex = self.determine_sex(chrom="Y")
        return sex

    # -- frame views ----------------------------------------------------------

    def _filter(self, chrom: str = "") -> pd.DataFrame:
        frame = self.frame
        return frame.loc[frame.chrom == chrom] if chrom else frame

    def heterozygous(self, chrom: str = "") -> pd.DataFrame:
        """Get heterozygous SNPs.

        Parameters
        ----------
        chrom : str, optional
            chromosome (e.g., "1", "X", "MT")

        Returns
        -------
        pandas.DataFrame
            normalized ``snps`` dataframe
        """
        df = self._filter(chrom)

        return df.loc[
            (df.genotype.notnull())
            & (df.genotype.str.len() == 2)
            & (df.genotype.str[0] != df.genotype.str[1])
        ]

    def homozygous(self, chrom: str = "") -> pd.DataFrame:
        """Get homozygous SNPs.

        Parameters
        ----------
        chrom : str, optional
            chromosome (e.g., "1", "X", "MT")

        Returns
        -------
        pandas.DataFrame
            normalized ``snps`` dataframe
        """
        df = self._filter(chrom)

        return df.loc[
            (df.genotype.notnull())
            & (df.genotype.str.len() == 2)
            & (df.genotype.str[0] == df.genotype.str[1])
        ]

    def notnull(self, chrom: str = "") -> pd.DataFrame:
        """Get not null genotype SNPs.

        Parameters
        ----------
        chrom : str, optional
            chromosome (e.g., "1", "X", "MT")

        Returns
        -------
        pandas.DataFrame
            normalized ``snps`` dataframe
        """
        df = self._filter(chrom)

        return df.loc[df.genotype.notnull()]

    # -- sex / build detection ------------------------------------------------

    def determine_sex(
        self,
        heterozygous_x_snps_threshold: float = 0.03,
        y_snps_not_null_threshold: float = 0.3,
        chrom: str = "X",
    ) -> str:
        """Determine sex from SNPs using thresholds.

        Parameters
        ----------
        heterozygous_x_snps_threshold : float
            percentage heterozygous X SNPs; above this threshold, Female is determined
        y_snps_not_null_threshold : float
            percentage Y SNPs that are not null; above this threshold, Male is determined
        chrom : {"X", "Y"}
            use X or Y chromosome SNPs to determine sex

        Returns
        -------
        str
            'Male' or 'Female' if detected, else empty str
        """
        if not self._frame.empty:
            if chrom == "X":
                return self._determine_sex_X(heterozygous_x_snps_threshold)
            elif chrom == "Y":
                return self._determine_sex_Y(y_snps_not_null_threshold)
        return ""

    def _determine_sex_X(self, threshold: float) -> str:
        x_snps = self.get_count("X")

        if x_snps > 0:
            if len(self.heterozygous("X")) / x_snps > threshold:
                return "Female"
            else:
                return "Male"
        else:
            return ""

    def _determine_sex_Y(self, threshold: float) -> str:
        y_snps = self.get_count("Y")

        if y_snps > 0:
            if len(self.notnull("Y")) / y_snps > threshold:
                return "Male"
            else:
                return "Female"
        else:
            return ""

    def _build_votes(self) -> Counter[int]:
        # Tally the build matched by each marker SNP present, in marker order; see
        # SNPs.detect_build for the marker SNPs and references
        votes: Counter[int] = Counter()
        for rsid, marker_data in BUILD_MARKER_SNPS.items():
            if rsid in self._frame.index:
                # rsids repeat if not deduplicated; use the first row, which
                # deduplication keeps
                pos = self._frame.loc[[rsid], "pos"].iloc[0]
                for build in (37, 38, 36):  # Check most common builds first
                    if marker_data[build] == pos:
                        votes[build] += 1
                        break
        return votes
