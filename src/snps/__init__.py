"""`snps`

tools for reading, writing, generating, merging, and remapping SNPs

"""

from importlib.metadata import version

from snps.snps import SNPs as SNPs

__version__ = version("snps")

__all__ = ["SNPs", "__version__"]
