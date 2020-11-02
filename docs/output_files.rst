Output Files
============
The various output files produced by ``snps`` are detailed below. Output files are saved in the
output directory, which is defined at the instantiation of a :class:`~snps.snps.SNPs` object.

Save SNPs
---------
SNPs can be saved with :meth:`SNPs.save <snps.snps.SNPs.save>`. By default, one tab-separated
``.txt`` or ``.vcf`` file (``vcf=True``) is output when SNPs are saved. If comma is specified as
the separator (``sep=","``), the default extension is ``.csv``.

The content of non-VCF files (after comment lines, which start with ``#``) is as follows:

==========  ===========
Column      Description
==========  ===========
rsid        SNP ID
chromosome  Chromosome of SNP
position    Position of SNP
genotype    Genotype of SNP
==========  ===========

When ``filename`` is not specified, default filenames are used as described below.

:meth:`SNPs.save <snps.snps.SNPs.save>`
```````````````````````````````````````

<source>_<assembly>.txt / <source>_<assembly>.csv
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Where ``source`` is the detected source(s) of SNPs data and ``assembly`` is the assembly of the
SNPs being saved.
