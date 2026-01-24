API Reference
=============

This section documents the complete API for the ``snps`` library.

Core Classes
------------

SNPs
~~~~

The main class for reading, writing, and analyzing genotype data.

.. automodule:: snps.snps
   :members:
   :undoc-members:
   :show-inheritance:
   :exclude-members: sort_snps, remap_snps, save_snps, snp_count, get_snp_count, not_null_snps,
      get_summary, get_assembly, get_chromosomes, get_chromosomes_summary, duplicate_snps,
      discrepant_XY_snps, heterozygous_MT_snps, heterozygous_snps, homozygous_snps,
      discrepant_positions, discrepant_genotypes, discrepant_snps, is_valid, save

I/O Operations
--------------

Modules for reading, writing, and generating SNP data files.

snps.io
~~~~~~~

.. automodule:: snps.io
   :members:
   :undoc-members:
   :show-inheritance:

snps.io.reader
~~~~~~~~~~~~~~

File format readers for various genotype data sources.

.. automodule:: snps.io.reader
   :members:
   :undoc-members:
   :show-inheritance:
   :exclude-members: read_file

snps.io.writer
~~~~~~~~~~~~~~

File format writers for exporting genotype data.

.. automodule:: snps.io.writer
   :members:
   :undoc-members:
   :show-inheritance:
   :exclude-members: write_file

snps.io.generator
~~~~~~~~~~~~~~~~~

Synthetic SNP data generation utilities.

.. automodule:: snps.io.generator
   :members:
   :undoc-members:
   :show-inheritance:

Data Resources
--------------

snps.ensembl
~~~~~~~~~~~~

Interface to Ensembl REST API for genomic data.

.. automodule:: snps.ensembl
   :members:
   :undoc-members:
   :show-inheritance:

snps.resources
~~~~~~~~~~~~~~

Resource management for reference data and assembly mappings.

.. automodule:: snps.resources
   :members:
   :undoc-members:
   :show-inheritance:

Utilities
---------

snps.utils
~~~~~~~~~~

Helper functions and utilities.

.. automodule:: snps.utils
   :members:
   :undoc-members:
   :show-inheritance:
