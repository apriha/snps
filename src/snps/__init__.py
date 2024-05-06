"""`snps`

tools for reading, writing, merging, and remapping SNPs

"""

from snps.snps import SNPs as SNPs

# set version string with Versioneer
from . import _version

__version__ = _version.get_versions()["version"]
