xy-chrom-snp-ratios
===================
an analysis of heterozygous X SNP and not-null Y SNP ratios

Method
------
All files in the `openSNP <https://opensnp.org>`_ data dump that ``snps`` could read were
analyzed for heterozygous X SNP and not-null Y SNP ratios.

Specifically, for each sample, the total number of X snps (`x_snps`) and heterozygous X snps
(`heterozygous_x_snps`) were calculated. The `heterozygous_x_snps` to `x_snps` ratio was
computed to determine if an optimum ratio could be used to separate Male and Female genotypes.

Similarly, for each sample, the total number of Y snps (`y_snps`) and not-null Y snps
(`y_snps_not_null`) were calculated. The `y_snps_not_null` to `y_snps` ratio was
computed to determine if an optimum ratio could be used to separate Male and Female genotypes.

Results
-------
Over nearly 4000 samples, heterozygous X chromosome SNPs appear to separate Male and Female
genotypes at a ratio of 0.03 (i.e., there are two separate distributions). Similarly, not-null Y
SNPs can be used to separate Male and Female genotypes, with a ratio of 0.3.

.. image:: https://raw.githubusercontent.com/apriha/snps/master/analysis/xy-chrom-snp-ratios/xy-chrom-snp-ratios.png
