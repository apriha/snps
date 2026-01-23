"""Centralized constants for genome builds and marker SNPs.

This module provides single source of truth for:
- Valid genome builds
- Marker SNPs for build detection with positions
- Chromosome sizes for each build

References
----------
1. Genome Reference Consortium, https://www.ncbi.nlm.nih.gov/grc
   - GRCh38: https://www.ncbi.nlm.nih.gov/grc/human/data?asm=GRCh38
   - GRCh37: https://www.ncbi.nlm.nih.gov/grc/human/data?asm=GRCh37
2. Database of Single Nucleotide Polymorphisms (dbSNP). Bethesda (MD): National Center
   for Biotechnology Information, National Library of Medicine. dbSNP accession: rs3094315,
   rs11928389, rs2500347, rs964481, rs2341354, rs3850290, and rs1329546
   (dbSNP Build ID: 151). Available from: http://www.ncbi.nlm.nih.gov/SNP/
3. UCSC Genome Browser chromosome size data:
   - NCBI36/hg18: https://hgdownload.soe.ucsc.edu/goldenPath/hg18/bigZips/hg18.chrom.sizes
   - GRCh37/hg19: https://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/hg19.chrom.sizes
   - GRCh38/hg38: https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.chrom.sizes
"""

from typing import Dict, Final, Tuple, Union

# Valid genome builds supported by this library
VALID_BUILDS: Final[Tuple[int, int, int]] = (36, 37, 38)

# Build name mappings
BUILD_ASSEMBLY_MAP: Final[Dict[int, str]] = {
    36: "NCBI36",
    37: "GRCh37",
    38: "GRCh38",
}

# Marker SNPs for build detection
# These SNPs have known positions across different builds and are used
# to automatically detect which build a genotype file uses
# Format: {rsid: {"chrom": str, 36: pos, 37: pos, 38: pos}}
BUILD_MARKER_SNPS: Final[Dict[str, Dict[Union[str, int], Union[str, int]]]] = {
    "rs3094315": {"chrom": "1", 36: 742429, 37: 752566, 38: 817186},
    "rs11928389": {"chrom": "1", 36: 50908372, 37: 50927009, 38: 50889578},
    "rs2500347": {"chrom": "1", 36: 143649677, 37: 144938320, 38: 148946169},
    "rs964481": {"chrom": "20", 36: 27566744, 37: 27656823, 38: 27638706},
    "rs2341354": {"chrom": "1", 36: 908436, 37: 918573, 38: 983193},
    "rs3850290": {"chrom": "2", 36: 22315141, 37: 23245301, 38: 22776092},
    "rs1329546": {"chrom": "1", 36: 135302086, 37: 135474420, 38: 136392261},
}

# Chromosome sizes indexed by build, then chromosome
# Values are approximate sizes in base pairs
CHROM_SIZES: Final[Dict[int, Dict[str, int]]] = {
    36: {
        "1": 247249719,
        "2": 242951149,
        "3": 199501827,
        "4": 191273063,
        "5": 180857866,
        "6": 170899992,
        "7": 158821424,
        "8": 146274826,
        "9": 140273252,
        "10": 135374737,
        "11": 134452384,
        "12": 132349534,
        "13": 114142980,
        "14": 106368585,
        "15": 100338915,
        "16": 88827254,
        "17": 78774742,
        "18": 76117153,
        "19": 63811651,
        "20": 62435964,
        "21": 46944323,
        "22": 49691432,
        "X": 154913754,
        "Y": 57772954,
        "MT": 16571,
    },
    37: {
        "1": 249250621,
        "2": 243199373,
        "3": 198022430,
        "4": 191154276,
        "5": 180915260,
        "6": 171115067,
        "7": 159138663,
        "8": 146364022,
        "9": 141213431,
        "10": 135534747,
        "11": 135006516,
        "12": 133851895,
        "13": 115169878,
        "14": 107349540,
        "15": 102531392,
        "16": 90354753,
        "17": 81195210,
        "18": 78077248,
        "19": 59128983,
        "20": 63025520,
        "21": 48129895,
        "22": 51304566,
        "X": 155270560,
        "Y": 59373566,
        "MT": 16571,
    },
    38: {
        "1": 248956422,
        "2": 242193529,
        "3": 198295559,
        "4": 190214555,
        "5": 181538259,
        "6": 170805979,
        "7": 159345973,
        "8": 145138636,
        "9": 138394717,
        "10": 133797422,
        "11": 135086622,
        "12": 133275309,
        "13": 114364328,
        "14": 107043718,
        "15": 101991189,
        "16": 90338345,
        "17": 83257441,
        "18": 80373285,
        "19": 58617616,
        "20": 64444167,
        "21": 46709983,
        "22": 50818468,
        "X": 156040895,
        "Y": 57227415,
        "MT": 16569,
    },
}
