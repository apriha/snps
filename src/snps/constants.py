import numpy as np

REFERENCE_SEQUENCE_CHROMS = tuple(str(i) for i in range(1, 23)) + ("X", "Y", "MT")

# dtypes of the columns of a normalized ``snps`` dataframe (``rsid`` is the index)
NORMALIZED_DTYPES: dict[str, type] = {
    "rsid": object,
    "chrom": object,
    "pos": np.uint32,
    "genotype": object,
}
