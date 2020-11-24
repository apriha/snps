.. image:: https://raw.githubusercontent.com/apriha/snps/master/docs/images/snps_banner.png

|build| |codecov| |docs| |pypi| |python| |downloads| |black|

snps
====
tools for reading, writing, merging, and remapping SNPs 🧬

``snps`` *strives to be an easy-to-use and accessible open-source library for working with
genotype data*

Features
--------
Input / Output
``````````````
- Read raw data (genotype) files from a variety of direct-to-consumer (DTC) DNA testing
  sources with a `SNPs <https://snps.readthedocs.io/en/latest/snps.html#snps.snps.SNPs>`_
  object
- Read and write VCF files (e.g., convert `23andMe <https://www.23andme.com>`_ to VCF)
- Merge raw data files from different DNA tests, identifying discrepant SNPs in the process
- Read data in a variety of formats (e.g., files, bytes, compressed with `gzip` or `zip`)
- Handle several variations of file types, validated via
  `openSNP parsing analysis <https://github.com/apriha/snps/tree/master/analysis/parse-opensnp-files>`_

Build / Assembly Detection and Remapping
````````````````````````````````````````
- Detect the build / assembly of SNPs (supports builds 36, 37, and 38)
- Remap SNPs between builds / assemblies

Data Cleaning
`````````````
- Fix several common issues when loading SNPs
- Sort SNPs based on chromosome and position
- Deduplicate RSIDs
- Deduplicate alleles in the non-PAR regions of the X and Y chromosomes for males
- Deduplicate alleles on MT
- Assign PAR SNPs to the X or Y chromosome

Supported Genotype Files
------------------------
``snps`` supports `VCF <https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3137218/>`_ files and
genotype files from the following DNA testing sources:

- `23andMe <https://www.23andme.com>`_
- `Ancestry <https://www.ancestry.com>`_
- `Código 46 <https://codigo46.com.mx>`_
- `DNA.Land <https://dna.land>`_
- `Family Tree DNA <https://www.familytreedna.com>`_
- `Genes for Good <https://genesforgood.sph.umich.edu>`_
- `LivingDNA <https://livingdna.com>`_
- `Mapmygenome <https://mapmygenome.in>`_
- `MyHeritage <https://www.myheritage.com>`_
- `Sano Genetics <https://sanogenetics.com>`_
- `tellmeGen <https://www.tellmegen.com>`_

Additionally, ``snps`` can read a variety of "generic" CSV and TSV files.

Dependencies
------------
``snps`` requires `Python <https://www.python.org>`_ 3.6.1+ and the following Python
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

Examples
--------
Download Example Data
`````````````````````
First, let's setup logging to get some helpful output:

>>> import logging, sys
>>> logger = logging.getLogger()
>>> logger.setLevel(logging.INFO)
>>> logger.addHandler(logging.StreamHandler(sys.stdout))

Now we're ready to download some example data from `openSNP <https://opensnp.org>`_:

>>> from snps.resources import Resources
>>> r = Resources()
>>> paths = r.download_example_datasets()
Downloading resources/662.23andme.340.txt.gz
Downloading resources/662.ftdna-illumina.341.csv.gz

Load Raw Data
`````````````
Load a `23andMe <https://www.23andme.com>`_ raw data file:

>>> from snps import SNPs
>>> s = SNPs("resources/662.23andme.340.txt.gz")
>>> s.source
'23andMe'
>>> s.count
991786

The ``SNPs`` class accepts a path to a file or a bytes object. A ``Reader`` class attempts to
infer the data source and load the SNPs. The loaded SNPs are
`normalized <https://snps.readthedocs.io/en/latest/snps.html#snps.snps.SNPs.snps>`_ and
available via a ``pandas.DataFrame``:

>>> df = s.snps
>>> df.columns.values
array(['chrom', 'pos', 'genotype'], dtype=object)
>>> df.index.name
'rsid'
>>> df.chrom.dtype.name
'object'
>>> df.pos.dtype.name
'uint32'
>>> df.genotype.dtype.name
'object'
>>> len(df)
991786

``snps`` also attempts to detect the build / assembly of the data:

>>> s.build
37
>>> s.build_detected
True
>>> s.assembly
'GRCh37'

Merge Raw Data Files
````````````````````
The dataset consists of raw data files from two different DNA testing sources - let's combine
these files. Specifically, we'll update the ``SNPs`` object with SNPs from a
`Family Tree DNA <https://www.familytreedna.com>`_ file.

>>> merge_results = s.merge([SNPs("resources/662.ftdna-illumina.341.csv.gz")])
Merging SNPs('resources/662.ftdna-illumina.341.csv.gz')
SNPs('resources/662.ftdna-illumina.341.csv.gz') has Build 36; remapping to Build 37
Downloading resources/NCBI36_GRCh37.tar.gz
27 SNP positions were discrepant; keeping original positions
151 SNP genotypes were discrepant; marking those as null
>>> s.source
'23andMe, FTDNA'
>>> s.count
1006960
>>> s.build
37
>>> s.build_detected
True

If the SNPs being merged have a build that differs from the destination build, the SNPs to merge
will be remapped automatically. After this example merge, the build is still detected, since the
build was detected for all ``SNPs`` objects that were merged.

As the data gets added, it's compared to the existing data, and SNP position and genotype
discrepancies are identified. (The discrepancy thresholds can be tuned via parameters.) These
discrepant SNPs are available for inspection after the merge via properties of the ``SNPs`` object.

>>> len(s.discrepant_merge_genotypes)
151

Additionally, any non-called / null genotypes will be updated during the merge, if the file
being merged has a called genotype for the SNP.

Moreover, ``merge`` takes a ``chrom`` parameter - this enables merging of only SNPs associated
with the specified chromosome (e.g., "Y" or "MT").

Finally, ``merge`` returns a list of ``dict``, where each ``dict`` has information corresponding
to the results of each merge (e.g., SNPs in common).

>>> sorted(list(merge_results[0].keys()))
['common_rsids', 'discrepant_genotype_rsids', 'discrepant_position_rsids', 'merged']
>>> merge_results[0]["merged"]
True
>>> len(merge_results[0]["common_rsids"])
692918

Remap SNPs
``````````
Now, let's remap the merged SNPs to change the assembly / build:

>>> s.snps.loc["rs3094315"].pos
752566
>>> chromosomes_remapped, chromosomes_not_remapped = s.remap(38)
Downloading resources/GRCh37_GRCh38.tar.gz
>>> s.build
38
>>> s.assembly
'GRCh38'
>>> s.snps.loc["rs3094315"].pos
817186

SNPs can be remapped between Build 36 (``NCBI36``), Build 37 (``GRCh37``), and Build 38
(``GRCh38``).

Save SNPs
`````````
Ok, so far we've merged the SNPs from two files (ensuring the same build in the process and
identifying discrepancies along the way). Then, we remapped the SNPs to Build 38. Now, let's save
the merged and remapped dataset consisting of 1M+ SNPs to a tab-separated values (TSV) file:

>>> saved_snps = s.save("out.txt")
Saving output/out.txt
>>> print(saved_snps)
output/out.txt

Moreover, let's get the reference sequences for this assembly and save the SNPs as a VCF file:

>>> saved_snps = s.save("out.vcf", vcf=True)
Downloading resources/fasta/GRCh38/Homo_sapiens.GRCh38.dna.chromosome.1.fa.gz
Downloading resources/fasta/GRCh38/Homo_sapiens.GRCh38.dna.chromosome.2.fa.gz
Downloading resources/fasta/GRCh38/Homo_sapiens.GRCh38.dna.chromosome.3.fa.gz
Downloading resources/fasta/GRCh38/Homo_sapiens.GRCh38.dna.chromosome.4.fa.gz
Downloading resources/fasta/GRCh38/Homo_sapiens.GRCh38.dna.chromosome.5.fa.gz
Downloading resources/fasta/GRCh38/Homo_sapiens.GRCh38.dna.chromosome.6.fa.gz
Downloading resources/fasta/GRCh38/Homo_sapiens.GRCh38.dna.chromosome.7.fa.gz
Downloading resources/fasta/GRCh38/Homo_sapiens.GRCh38.dna.chromosome.8.fa.gz
Downloading resources/fasta/GRCh38/Homo_sapiens.GRCh38.dna.chromosome.9.fa.gz
Downloading resources/fasta/GRCh38/Homo_sapiens.GRCh38.dna.chromosome.10.fa.gz
Downloading resources/fasta/GRCh38/Homo_sapiens.GRCh38.dna.chromosome.11.fa.gz
Downloading resources/fasta/GRCh38/Homo_sapiens.GRCh38.dna.chromosome.12.fa.gz
Downloading resources/fasta/GRCh38/Homo_sapiens.GRCh38.dna.chromosome.13.fa.gz
Downloading resources/fasta/GRCh38/Homo_sapiens.GRCh38.dna.chromosome.14.fa.gz
Downloading resources/fasta/GRCh38/Homo_sapiens.GRCh38.dna.chromosome.15.fa.gz
Downloading resources/fasta/GRCh38/Homo_sapiens.GRCh38.dna.chromosome.16.fa.gz
Downloading resources/fasta/GRCh38/Homo_sapiens.GRCh38.dna.chromosome.17.fa.gz
Downloading resources/fasta/GRCh38/Homo_sapiens.GRCh38.dna.chromosome.18.fa.gz
Downloading resources/fasta/GRCh38/Homo_sapiens.GRCh38.dna.chromosome.19.fa.gz
Downloading resources/fasta/GRCh38/Homo_sapiens.GRCh38.dna.chromosome.20.fa.gz
Downloading resources/fasta/GRCh38/Homo_sapiens.GRCh38.dna.chromosome.21.fa.gz
Downloading resources/fasta/GRCh38/Homo_sapiens.GRCh38.dna.chromosome.22.fa.gz
Downloading resources/fasta/GRCh38/Homo_sapiens.GRCh38.dna.chromosome.X.fa.gz
Downloading resources/fasta/GRCh38/Homo_sapiens.GRCh38.dna.chromosome.Y.fa.gz
Downloading resources/fasta/GRCh38/Homo_sapiens.GRCh38.dna.chromosome.MT.fa.gz
Saving output/out.vcf
1 SNP positions were found to be discrepant when saving VCF

When saving a VCF, if any SNPs have positions outside of the reference sequence, they are marked
as discrepant and are available via a property of the ``SNPs`` object.

All `output files <https://snps.readthedocs.io/en/latest/output_files.html>`_ are saved to the
output directory.

Documentation
-------------
Documentation is available `here <https://snps.readthedocs.io/>`_.

Acknowledgements
----------------
Thanks to Mike Agostino, Padma Reddy, Kevin Arvai, `openSNP <https://opensnp.org>`_,
`Open Humans <https://www.openhumans.org>`_, and `Sano Genetics <https://sanogenetics.com>`_.

.. https://github.com/rtfd/readthedocs.org/blob/master/docs/badges.rst
.. |build| image:: https://travis-ci.com/apriha/snps.svg?branch=master
   :target: https://travis-ci.com/apriha/snps
.. |codecov| image:: https://codecov.io/gh/apriha/snps/branch/master/graph/badge.svg
   :target: https://codecov.io/gh/apriha/snps
.. |docs| image:: https://readthedocs.org/projects/snps/badge/?version=latest
   :target: https://snps.readthedocs.io/
.. |pypi| image:: https://img.shields.io/pypi/v/snps.svg
   :target: https://pypi.python.org/pypi/snps
.. |python| image:: https://img.shields.io/pypi/pyversions/snps.svg
   :target: https://www.python.org
.. |downloads| image:: https://pepy.tech/badge/snps
   :target: https://pepy.tech/project/snps
.. |black| image:: https://img.shields.io/badge/code%20style-black-000000.svg
   :target: https://github.com/psf/black
