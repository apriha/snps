Output Files
============
The various output files produced by ``snps`` are detailed below. Output files are saved in the
output directory, which is defined at the instantiation of a :class:`~snps.snps.SNPs` or
:class:`~snps.snps_collection.SNPsCollection` object.

Load SNPs
---------
Multiple raw data files can be loaded when a :class:`~snps.snps_collection.SNPsCollection` is
instantiated.

When loading multiple raw data files, if there are any discrepancies between the existing data
and the new data, those are noted. Specifically, discrepancies in SNP positions and genotypes
are output as CSV files. Output of these files is enabled via the ``save_output=True`` argument
to :meth:`~snps.snps_collection.SNPsCollection.load_snps`.

discrepant_positions_<num>.csv
``````````````````````````````
Where ``num`` indicates the file count for discrepant positions files.

==============  ===========
Column          Description
==============  ===========
rsid            SNP ID
chrom           Chromosome of existing SNP
pos             Position of existing SNP
genotype        Genotype of existing SNP
chrom_added     Chromosome of added SNP
pos_added       Position of added SNP (discrepant with pos)
genotype_added  Genotype of added SNP
==============  ===========

A large number of discrepant positions could indicate that the files contain SNPs from different
builds / assemblies.

discrepant_genotypes_<num>.csv
``````````````````````````````
Where ``num`` indicates the file count for discrepant genotypes files.

===============  ===========
Column           Description
===============  ===========
rsid             SNP ID
chrom            Chromosome of existing SNP
pos              Position of existing SNP
genotype         Genotype of existing SNP
chrom_added      Chromosome of added SNP
pos_added        Position of added SNP
genotype_added   Genotype of added SNP (discrepant with genotype)
===============  ===========

A large number of discrepant genotypes could indicate that the files contain SNPs for different
individuals.

Discrepant SNPs
---------------
Summaries can be saved of the discrepant SNPs found while loading files.

discrepant_positions.csv
````````````````````````
SNPs with discrepant positions can be saved with
:meth:`~snps.snps_collection.SNPsCollection.save_discrepant_positions`.

==============  ===========
Column          Description
==============  ===========
rsid            SNP ID
chrom           Chromosome of existing SNP
pos             Position of existing SNP
genotype        Genotype of existing SNP
chrom_added     Chromosome of added SNP
pos_added       Position of added SNP (discrepant with pos)
genotype_added  Genotype of added SNP
==============  ===========

discrepant_genotypes.csv
````````````````````````
SNPs with discrepant genotypes can be saved with
:meth:`~snps.snps_collection.SNPsCollection.save_discrepant_genotypes`.

===============  ===========
Column           Description
===============  ===========
rsid             SNP ID
chrom            Chromosome of existing SNP
pos              Position of existing SNP
genotype         Genotype of existing SNP
chrom_added      Chromosome of added SNP
pos_added        Position of added SNP
genotype_added   Genotype of added SNP (discrepant with genotype)
===============  ===========

discrepant_snps.csv
```````````````````
SNPs with discrepant positions and / or genotypes can be saved with
:meth:`~snps.snps_collection.SNPsCollection.save_discrepant_snps`.

===============  ===========
Column           Description
===============  ===========
rsid             SNP ID
chrom            Chromosome of existing SNP
pos              Position of existing SNP
genotype         Genotype of existing SNP
chrom_added      Chromosome of added SNP
pos_added        Position of added SNP (possibly discrepant with pos)
genotype_added   Genotype of added SNP (possibly discrepant with genotype)
===============  ===========

Save SNPs
---------
SNPs can be saved with :meth:`SNPs.save_snps <snps.snps.SNPs.save_snps>` or
:meth:`SNPsCollection.save_snps <snps.snps_collection.SNPsCollection.save_snps>`. By default, one
tab-separated ``.txt`` or ``.vcf`` file (``vcf=True``) is output when SNPs are saved. If comma
is specified as the separator (``sep=","``), the default extension is ``.csv``.

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

:meth:`SNPs.save_snps <snps.snps.SNPs.save_snps>`
`````````````````````````````````````````````````

<source>_<assembly>.txt / <source>_<assembly>.csv
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Where ``source`` is the detected source of SNPs data and ``assembly`` is the assembly of the
SNPs being saved.


:meth:`SNPsCollection.save_snps <snps.snps_collection.SNPsCollection.save_snps>`
````````````````````````````````````````````````````````````````````````````````

<name>_<assembly>.txt / <name>_<assembly>.csv
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Where ``name`` is the name (if any) for the :class:`~snps.snps_collection.SNPsCollection` and
``assembly`` is the assembly of the SNPs being saved. If name is not specified, ``<name>_`` is
not included in the filename.
