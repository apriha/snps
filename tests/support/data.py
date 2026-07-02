"""Committed test data served by :class:`FakeResources` and shared with test cases.

This is the single source for the canned resource data the suite uses offline:
assembly mapping data, PAR (RefSNP) lookups, and small chip-cluster / low-quality
DataFrames. Keeping it here avoids duplicating the data across test modules.
"""

import numpy as np
import pandas as pd


def get_test_assembly_mapping_data(source, target, strands, mappings):
    """Build assembly mapping data (chroms "1" and "3") for `source` -> `target`."""
    return {
        "1": {
            "mappings": [
                {
                    "original": {
                        "seq_region_name": "1",
                        "strand": strands[0],
                        "start": mappings[0],
                        "end": mappings[0],
                        "assembly": f"{source}",
                    },
                    "mapped": {
                        "seq_region_name": "1",
                        "strand": strands[1],
                        "start": mappings[1],
                        "end": mappings[1],
                        "assembly": f"{target}",
                    },
                },
                {
                    "original": {
                        "seq_region_name": "1",
                        "strand": strands[2],
                        "start": mappings[2],
                        "end": mappings[2],
                        "assembly": f"{source}",
                    },
                    "mapped": {
                        "seq_region_name": "1",
                        "strand": strands[3],
                        "start": mappings[3],
                        "end": mappings[3],
                        "assembly": f"{target}",
                    },
                },
                {
                    "original": {
                        "seq_region_name": "1",
                        "strand": strands[4],
                        "start": mappings[4],
                        "end": mappings[4],
                        "assembly": f"{source}",
                    },
                    "mapped": {
                        "seq_region_name": "1",
                        "strand": strands[5],
                        "start": mappings[5],
                        "end": mappings[5],
                        "assembly": f"{target}",
                    },
                },
            ]
        },
        "3": {
            "mappings": [
                {
                    "original": {
                        "seq_region_name": "3",
                        "strand": strands[6],
                        "start": mappings[6],
                        "end": mappings[6],
                        "assembly": f"{source}",
                    },
                    "mapped": {
                        "seq_region_name": "3",
                        "strand": strands[7],
                        "start": mappings[7],
                        "end": mappings[7],
                        "assembly": f"{target}",
                    },
                }
            ]
        },
    }


def NCBI36_GRCh37():
    return get_test_assembly_mapping_data(
        "NCBI36",
        "GRCh37",
        [1, 1, 1, 1, 1, 1, 1, -1],
        [
            742429,
            752566,
            143649677,
            144938320,
            143649678,
            144938321,
            50908372,
            50927009,
        ],
    )


def GRCh37_NCBI36():
    return get_test_assembly_mapping_data(
        "GRCh37",
        "NCBI36",
        [1, 1, 1, 1, 1, 1, 1, -1],
        [
            752566,
            742429,
            144938320,
            143649677,
            144938321,
            143649678,
            50927009,
            50908372,
        ],
    )


def GRCh37_GRCh38():
    return get_test_assembly_mapping_data(
        "GRCh37",
        "GRCh38",
        [1, 1, 1, -1, 1, -1, 1, 1],
        [
            752566,
            817186,
            144938320,
            148946169,
            144938321,
            148946168,
            50927009,
            50889578,
        ],
    )


def GRCh37_GRCh38_PAR():
    return {
        "X": {
            "mappings": [
                {
                    "original": {
                        "seq_region_name": "X",
                        "strand": 1,
                        "start": 220770,
                        "end": 220770,
                        "assembly": "GRCh37",
                    },
                    "mapped": {
                        "seq_region_name": "X",
                        "strand": 1,
                        "start": 304103,
                        "end": 304103,
                        "assembly": "GRCh38",
                    },
                },
                {
                    "original": {
                        "seq_region_name": "X",
                        "strand": 1,
                        "start": 91941056,
                        "end": 91941056,
                        "assembly": "GRCh37",
                    },
                    "mapped": {
                        "seq_region_name": "X",
                        "strand": 1,
                        "start": 92686057,
                        "end": 92686057,
                        "assembly": "GRCh38",
                    },
                },
            ]
        },
        "Y": {
            "mappings": [
                {
                    "original": {
                        "seq_region_name": "Y",
                        "strand": 1,
                        "start": 535258,
                        "end": 535258,
                        "assembly": "GRCh37",
                    },
                    "mapped": {
                        "seq_region_name": "Y",
                        "strand": 1,
                        "start": 624523,
                        "end": 624523,
                        "assembly": "GRCh38",
                    },
                }
            ]
        },
    }


def standard_assembly_mappings():
    """Default ``(source, target) -> mapping data`` map served by ``FakeResources``."""
    return {
        ("NCBI36", "GRCh37"): NCBI36_GRCh37(),
        ("GRCh37", "NCBI36"): GRCh37_NCBI36(),
        ("GRCh37", "GRCh38"): GRCh37_GRCh38(),
    }


# RefSNP snapshots for PAR SNPs, keyed by RefSNP id (incl. one merged snapshot).
PAR_EFFECTS = [
    {
        "refsnp_id": "758419898",
        "create_date": "2015-04-1T22:25Z",
        "last_update_date": "2019-07-14T04:19Z",
        "last_update_build_id": "153",
        "primary_snapshot_data": {
            "placements_with_allele": [
                {
                    "seq_id": "NC_000024.9",
                    "placement_annot": {
                        "seq_id_traits_by_assembly": [{"assembly_name": "GRCh37.p13"}]
                    },
                    "alleles": [
                        {
                            "allele": {
                                "spdi": {"seq_id": "NC_000024.9", "position": 7364103}
                            }
                        }
                    ],
                }
            ]
        },
    },
    {
        "refsnp_id": "28736870",
        "create_date": "2005-05-24T14:43Z",
        "last_update_date": "2019-07-14T04:18Z",
        "last_update_build_id": "153",
        "primary_snapshot_data": {
            "placements_with_allele": [
                {
                    "seq_id": "NC_000023.10",
                    "placement_annot": {
                        "seq_id_traits_by_assembly": [{"assembly_name": "GRCh37.p13"}]
                    },
                    "alleles": [
                        {
                            "allele": {
                                "spdi": {"seq_id": "NC_000023.10", "position": 220769}
                            }
                        }
                    ],
                }
            ]
        },
    },
    {
        "refsnp_id": "113313554",
        "create_date": "2010-07-4T18:13Z",
        "last_update_date": "2019-07-14T04:18Z",
        "last_update_build_id": "153",
        "primary_snapshot_data": {
            "placements_with_allele": [
                {
                    "seq_id": "NC_000024.9",
                    "placement_annot": {
                        "seq_id_traits_by_assembly": [{"assembly_name": "GRCh37.p13"}]
                    },
                    "alleles": [
                        {
                            "allele": {
                                "spdi": {"seq_id": "NC_000024.9", "position": 535257}
                            }
                        }
                    ],
                }
            ]
        },
    },
    {
        "refsnp_id": "113378274",
        "create_date": "2010-07-4T18:14Z",
        "last_update_date": "2016-03-3T10:51Z",
        "last_update_build_id": "147",
        "merged_snapshot_data": {"merged_into": ["72608386"]},
    },
    {
        "refsnp_id": "72608386",
        "create_date": "2009-02-14T01:08Z",
        "last_update_date": "2019-07-14T04:05Z",
        "last_update_build_id": "153",
        "primary_snapshot_data": {
            "placements_with_allele": [
                {
                    "seq_id": "NC_000023.10",
                    "placement_annot": {
                        "seq_id_traits_by_assembly": [{"assembly_name": "GRCh37.p13"}]
                    },
                    "alleles": [
                        {
                            "allele": {
                                "spdi": {"seq_id": "NC_000023.10", "position": 91941055}
                            }
                        }
                    ],
                }
            ]
        },
    },
]


def chip_clusters_df(pos=tuple(range(101, 109)), cluster="c1", length=8):
    """Build a parsed chip-clusters DataFrame (as returned by ``get_chip_clusters``)."""
    df = pd.DataFrame(
        {"chrom": ["1"] * length, "pos": pos, "clusters": [cluster] * length},
        columns=["chrom", "pos", "clusters"],
    )
    df.chrom = df.chrom.astype(pd.CategoricalDtype(ordered=False))
    df.pos = df.pos.astype(np.uint32)
    df.clusters = df.clusters.astype(pd.CategoricalDtype(ordered=False))
    return df


def low_quality_snps_df(pos=(104, 106, 1001), cluster="c1"):
    """Build a parsed low-quality-SNPs DataFrame (as returned by ``get_low_quality_snps``)."""
    df = pd.DataFrame(
        {"chrom": ["1"] * len(pos), "pos": pos, "cluster": [cluster] * len(pos)},
        columns=["chrom", "pos", "cluster"],
    )
    df.chrom = df.chrom.astype(pd.CategoricalDtype(ordered=False))
    df.pos = df.pos.astype(np.uint32)
    df.cluster = df.cluster.astype(pd.CategoricalDtype(ordered=False))
    return df
