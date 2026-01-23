.. image:: https://raw.githubusercontent.com/apriha/snps/main/docs/images/snps_banner.png

|ci| |codecov| |docs| |pypi| |python| |downloads| |ruff|

snps
====
tools for reading, writing, generating, merging, and remapping SNPs 🧬

``snps`` *strives to be an easy-to-use and accessible open-source library for working with
genotype data*

Features
--------
Input / Output
``````````````
- Read raw data (genotype) files from a variety of direct-to-consumer (DTC) DNA testing
  sources with a `SNPs <https://snps.readthedocs.io/en/stable/snps.html#snps.snps.SNPs>`_
  object
- Read and write VCF files (e.g., convert `23andMe <https://www.23andme.com>`_ to VCF)
- Merge raw data files from different DNA tests, identifying discrepant SNPs in the process
- Read data in a variety of formats (e.g., files, bytes, compressed with `gzip` or `zip`)
- Handle several variations of file types, historically validated using
  data from `openSNP <https://opensnp.org>`_
- Generate synthetic genotype data for testing and examples

Build / Assembly Detection and Remapping
````````````````````````````````````````
- Detect the build / assembly of SNPs (supports builds 36, 37, and 38)
- Remap SNPs between builds / assemblies

Data Cleaning
`````````````
- Perform quality control (QC) / filter low quality SNPs based on `chip clusters <https://doi.org/10.1016/j.csbj.2021.06.040>`_
- Fix several common issues when loading SNPs
- Sort SNPs based on chromosome and position
- Deduplicate RSIDs
- Deduplicate alleles in the non-PAR regions of the X and Y chromosomes for males
- Deduplicate alleles on MT
- Assign PAR SNPs to the X or Y chromosome

Analysis
````````
- Derive sex from SNPs
- Detect deduced genotype / chip array and chip version based on `chip clusters <https://doi.org/10.1016/j.csbj.2021.06.040>`_
- Predict ancestry from SNPs (when installed with `ezancestry <https://github.com/arvkevi/ezancestry>`_)

Supported Genotype Files
------------------------
``snps`` supports `VCF <https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3137218/>`_ files and
genotype files from the following DNA testing sources:

- `23andMe <https://www.23andme.com>`_
- `23Mofang <https://www.23mofang.com>`_
- `Ancestry <https://www.ancestry.com>`_
- `CircleDNA <https://circledna.com/>`_
- `Código 46 <https://codigo46.com.mx>`_
- `DNA.Land <https://dna.land>`_
- `Family Tree DNA <https://www.familytreedna.com>`_
- `Genes for Good <https://genesforgood.sph.umich.edu>`_
- `LivingDNA <https://livingdna.com>`_
- `Mapmygenome <https://mapmygenome.in>`_
- `MyHeritage <https://www.myheritage.com>`_
- `PLINK <https://www.cog-genomics.org/plink/>`_
- `Sano Genetics <https://sanogenetics.com>`_
- `SelfDecode <https://selfdecode.com>`_
- `tellmeGen <https://www.tellmegen.com>`_

Additionally, ``snps`` can read a variety of "generic" CSV and TSV files.

Dependencies
------------
``snps`` requires `Python <https://www.python.org>`_ 3.9+ and the following Python
packages:

- `numpy <http://www.numpy.org>`_
- `pandas <http://pandas.pydata.org>`_
- `atomicwrites <https://github.com/untitaker/python-atomicwrites>`_

Installation
------------
``snps`` is `available <https://pypi.org/project/snps/>`_ on the
`Python Package Index <https://pypi.org>`_. Install ``snps`` (and its required
Python dependencies) via ``pip``::

    $ pip install snps

For `ancestry prediction <https://snps.readthedocs.io/en/stable/snps.html#snps.snps.SNPs.predict_ancestry>`_
capability, ``snps`` can be installed with `ezancestry <https://github.com/arvkevi/ezancestry>`_::

    $ pip install snps[ezancestry]

Examples
--------

To try these examples, first generate some sample data:

>>> from snps.resources import Resources
>>> paths = Resources().create_example_datasets()

Load a Raw Data File
````````````````````
Load a raw data file exported from a DNA testing source (e.g.,
`23andMe <https://www.23andme.com>`_, `AncestryDNA <https://www.ancestry.com>`_,
`Family Tree DNA <https://www.familytreedna.com>`_):

>>> from snps import SNPs
>>> s = SNPs("resources/sample1.23andme.txt.gz")

``snps`` automatically detects the source format and `normalizes
<https://snps.readthedocs.io/en/stable/snps.html#snps.snps.SNPs.snps>`_ the data:

>>> s.source
'23andMe'
>>> s.count
991767
>>> s.build
37
>>> s.assembly
'GRCh37'

The SNPs are available as a ``pandas.DataFrame``:

>>> df = s.snps
>>> df.columns.values
array(['chrom', 'pos', 'genotype'], dtype=object)
>>> len(df)
991767

Merge Raw Data Files
````````````````````
Combine SNPs from multiple files (e.g., combine data from different testing companies):

>>> results = s.merge([SNPs("resources/sample2.ftdna.csv.gz")])
>>> s.count
1006949

SNPs are compared during the merge. Position and genotype discrepancies are identified and
can be inspected via properties of the ``SNPs`` object:

>>> len(s.discrepant_merge_positions)
27
>>> len(s.discrepant_merge_genotypes)
156

Remap SNPs
``````````
Convert SNPs between genome assemblies (Build 36/NCBI36, Build 37/GRCh37, Build 38/GRCh38):

>>> chromosomes_remapped, chromosomes_not_remapped = s.remap(38)
>>> s.assembly
'GRCh38'

Save SNPs
`````````
Save SNPs to common file formats:

>>> _ = s.to_tsv("output.txt")
>>> _ = s.to_csv("output.csv")

To save as VCF, ``snps`` automatically downloads the required reference sequences for the
assembly. This ensures the REF alleles in the VCF are accurate:

>>> _ = s.to_vcf("output.vcf")  # doctest: +SKIP

All output files are saved to the `output directory
<https://snps.readthedocs.io/en/stable/output_files.html>`_.

Generate Synthetic Data
```````````````````````
Generate synthetic genotype data for testing, examples, or demonstrations:

>>> from snps.io import SyntheticSNPGenerator
>>> gen = SyntheticSNPGenerator(build=37, seed=123)
>>> gen.save_as_23andme("synthetic_23andme.txt.gz", num_snps=10000)
'synthetic_23andme.txt.gz'

The generator supports multiple output formats (23andMe, AncestryDNA, FTDNA) and
automatically injects build-specific marker SNPs to ensure accurate build detection.

Documentation
-------------
Documentation is available `here <https://snps.readthedocs.io/>`_.

Acknowledgements
----------------
Thanks to Mike Agostino, Padma Reddy, Kevin Arvai, `Open Humans <https://www.openhumans.org>`_,
and `Sano Genetics <https://sanogenetics.com>`_. This project was historically validated using
data from `openSNP <https://opensnp.org>`_.

``snps`` incorporates code and concepts generated with the assistance of various
generative AI tools (including but not limited to `ChatGPT <https://chatgpt.com>`_,
`Grok <https://grok.com>`_, and `Claude <https://claude.ai>`_). ✨

License
-------
``snps`` is licensed under the `BSD 3-Clause License <https://github.com/apriha/snps/blob/main/LICENSE.txt>`_.

.. https://github.com/rtfd/readthedocs.org/blob/master/docs/badges.rst
.. |ci| image:: https://github.com/apriha/snps/actions/workflows/ci.yml/badge.svg?branch=main
   :target: https://github.com/apriha/snps/actions/workflows/ci.yml
.. |codecov| image:: https://codecov.io/gh/apriha/snps/branch/main/graph/badge.svg
   :target: https://codecov.io/gh/apriha/snps
.. |docs| image:: https://readthedocs.org/projects/snps/badge/?version=stable
   :target: https://snps.readthedocs.io/
.. |pypi| image:: https://img.shields.io/pypi/v/snps.svg
   :target: https://pypi.python.org/pypi/snps
.. |python| image:: https://img.shields.io/pypi/pyversions/snps.svg
   :target: https://www.python.org
.. |downloads| image:: https://pepy.tech/badge/snps
   :target: https://pepy.tech/project/snps
.. |ruff| image:: https://img.shields.io/endpoint?url=https://raw.githubusercontent.com/astral-sh/ruff/main/assets/badge/v2.json
   :target: https://github.com/astral-sh/ruff
   :alt: Ruff
