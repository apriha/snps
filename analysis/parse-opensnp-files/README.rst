parse-opensnp-files
===================
scripts to load and debug parsing of openSNP datadump files

Method
------
Attempt to parse each file in the `openSNP <https://opensnp.org>`_ datadump by creating a
``SNPs`` object. For files where SNPs were loaded, save summary statistics to a dataframe and
output as a CSV. For files where no SNPs were loaded, save a message for each file indicating
the issue and optionally extract these files from the datadump for debugging.

Results
-------
As of May 2020, ``snps`` can parse ~96.6% of the genotype files in the datadump. Additionally,
``snps`` can detect the build in ~99.9% of those files.
