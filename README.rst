::

 GATCACAGGTCTATCAC     CCTATTAACCACTCAC     GGGAGCTCTCCATGCAT     TTGGTATTTTCGTCTGG
 GGGGTATGCACGCGATA     GCATTGCGAGACGCTG     GAGCCGGAGCACCCTAT     GTCGCAGTATCTGTCTT
 TGA                   TTC          CTG     CCT           CAT     CCT
 ATTATTTATCGCACCTA     CGT          TCA     ATATTACAGGCGAACAT     ACTTACTAAAGTGTGTT
 AATTAATTAATGCTTGT     AGG          ACA     TAATAATAACAATTGAA     TGTCTGCACAGCCACTT
               TCC     ACA          CAG     ACA                                 TCA
 TAACAAAAAATTTCCAC     CAA          ACC     CCC                   CCTCCCCCGCTTCTGGC
 CACAGCACTTAAACACA     TCT          CTG     CCA                   AACCCCAAAAACAAAGA

|build| |codecov| |docs| |pypi| |python|

snps
====
tools for reading, writing, merging, and remapping SNPs 🧬

Capabilities
------------
- Read raw data (genotype) files from a variety of direct-to-consumer (DTC) DNA testing companies
- Read and write VCF files for Builds 36, 37, and 38 (e.g., convert `23andMe <https://www.23andme.com>`_ to VCF)
- Merge raw data files from different DNA testing companies, identifying discrepant SNPs in the process
- Remap SNPs between assemblies / builds (e.g., convert SNPs from Build 36 to Build 37, etc.)

Supported Genotype Files
------------------------
``snps`` supports `VCF <https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3137218/>`_ files and
genotype files from the following DTC DNA testing companies:

- `23andMe <https://www.23andme.com>`_
- `Ancestry <https://www.ancestry.com>`_
- `Family Tree DNA <https://www.familytreedna.com>`_
- `MyHeritage <https://www.myheritage.com>`_

Dependencies
------------
``snps`` requires `Python <https://www.python.org>`_ 3.5+ and the following Python packages:

- `numpy <http://www.numpy.org>`_
- `pandas <http://pandas.pydata.org>`_
- `atomicwrites <https://github.com/untitaker/python-atomicwrites>`_
- `PyVCF <https://github.com/jamescasbon/PyVCF>`_

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
Let's download some example data from `openSNP <https://opensnp.org>`_:

>>> from snps.resources import Resources
>>> r = Resources()
>>> paths = r.download_example_datasets()
Downloading resources/662.23andme.340.txt.gz
Downloading resources/662.ftdna-illumina.341.csv.gz

Load Raw Data
`````````````
Load a `23andMe <https://www.23andme.com>`_ raw data file:

>>> from snps import SNPs
>>> s = SNPs('resources/662.23andme.340.txt.gz')

The loaded SNPs are available via a ``pandas.DataFrame``:

>>> df = s.snps
>>> df.columns.values
array(['chrom', 'pos', 'genotype'], dtype=object)
>>> df.index.name
'rsid'
>>> len(df)
991786

``snps`` also attempts to detect the build / assembly of the data:

>>> s.build
37
>>> s.build_detected
True
>>> s.assembly
'GRCh37'

Remap SNPs
``````````
Let's remap the SNPs to change the assembly / build:

>>> s.snps.loc["rs3094315"].pos
752566
>>> chromosomes_remapped, chromosomes_not_remapped = s.remap_snps(38)
Downloading resources/GRCh37_GRCh38.tar.gz
>>> s.build
38
>>> s.assembly
'GRCh38'
>>> s.snps.loc["rs3094315"].pos
817186

SNPs can be remapped between Build 36 (``NCBI36``), Build 37 (``GRCh37``), and Build 38
(``GRCh38``).

Merge Raw Data Files
````````````````````
The dataset consists of raw data files from two different DNA testing companies. Let's combine
these files using a ``SNPsCollection``.

>>> from snps import SNPsCollection
>>> sc = SNPsCollection("resources/662.ftdna-illumina.341.csv.gz", name="User662")
Loading resources/662.ftdna-illumina.341.csv.gz
>>> sc.build
36
>>> chromosomes_remapped, chromosomes_not_remapped = sc.remap_snps(37)
Downloading resources/NCBI36_GRCh37.tar.gz
>>> sc.snp_count
708092

As the data gets added, it's compared to the existing data, and SNP position and genotype
discrepancies are identified. (The discrepancy thresholds can be tuned via parameters.)

>>> sc.load_snps(["resources/662.23andme.340.txt.gz"], discrepant_genotypes_threshold=300)
Loading resources/662.23andme.340.txt.gz
27 SNP positions were discrepant; keeping original positions
151 SNP genotypes were discrepant; marking those as null
>>> len(sc.discrepant_snps)  # SNPs with discrepant positions and genotypes, dropping dups
169
>>> sc.snp_count
1006960

Save SNPs
`````````
Ok, so far we've remapped the SNPs to the same build and merged the SNPs from two files,
identifying discrepancies along the way. Let's save the merged dataset consisting of over 1M+
SNPs to a CSV file:

>>> saved_snps = sc.save_snps()
Saving output/User662_GRCh37.csv

Moreover, let's get the reference sequences for this assembly and save the SNPs as a VCF file:

>>> saved_snps = sc.save_snps(vcf=True)
Downloading resources/fasta/GRCh37/Homo_sapiens.GRCh37.dna.chromosome.1.fa.gz
Downloading resources/fasta/GRCh37/Homo_sapiens.GRCh37.dna.chromosome.2.fa.gz
Downloading resources/fasta/GRCh37/Homo_sapiens.GRCh37.dna.chromosome.3.fa.gz
Downloading resources/fasta/GRCh37/Homo_sapiens.GRCh37.dna.chromosome.4.fa.gz
Downloading resources/fasta/GRCh37/Homo_sapiens.GRCh37.dna.chromosome.5.fa.gz
Downloading resources/fasta/GRCh37/Homo_sapiens.GRCh37.dna.chromosome.6.fa.gz
Downloading resources/fasta/GRCh37/Homo_sapiens.GRCh37.dna.chromosome.7.fa.gz
Downloading resources/fasta/GRCh37/Homo_sapiens.GRCh37.dna.chromosome.8.fa.gz
Downloading resources/fasta/GRCh37/Homo_sapiens.GRCh37.dna.chromosome.9.fa.gz
Downloading resources/fasta/GRCh37/Homo_sapiens.GRCh37.dna.chromosome.10.fa.gz
Downloading resources/fasta/GRCh37/Homo_sapiens.GRCh37.dna.chromosome.11.fa.gz
Downloading resources/fasta/GRCh37/Homo_sapiens.GRCh37.dna.chromosome.12.fa.gz
Downloading resources/fasta/GRCh37/Homo_sapiens.GRCh37.dna.chromosome.13.fa.gz
Downloading resources/fasta/GRCh37/Homo_sapiens.GRCh37.dna.chromosome.14.fa.gz
Downloading resources/fasta/GRCh37/Homo_sapiens.GRCh37.dna.chromosome.15.fa.gz
Downloading resources/fasta/GRCh37/Homo_sapiens.GRCh37.dna.chromosome.16.fa.gz
Downloading resources/fasta/GRCh37/Homo_sapiens.GRCh37.dna.chromosome.17.fa.gz
Downloading resources/fasta/GRCh37/Homo_sapiens.GRCh37.dna.chromosome.18.fa.gz
Downloading resources/fasta/GRCh37/Homo_sapiens.GRCh37.dna.chromosome.19.fa.gz
Downloading resources/fasta/GRCh37/Homo_sapiens.GRCh37.dna.chromosome.20.fa.gz
Downloading resources/fasta/GRCh37/Homo_sapiens.GRCh37.dna.chromosome.21.fa.gz
Downloading resources/fasta/GRCh37/Homo_sapiens.GRCh37.dna.chromosome.22.fa.gz
Downloading resources/fasta/GRCh37/Homo_sapiens.GRCh37.dna.chromosome.X.fa.gz
Downloading resources/fasta/GRCh37/Homo_sapiens.GRCh37.dna.chromosome.Y.fa.gz
Downloading resources/fasta/GRCh37/Homo_sapiens.GRCh37.dna.chromosome.MT.fa.gz
Saving output/User662_GRCh37.vcf

All `output files <https://snps.readthedocs.io/en/latest/output_files.html>`_ are saved to the
output directory.

Documentation
-------------
Documentation is available `here <https://snps.readthedocs.io/>`_.

Origins
-------
Initial ``snps`` capability was developed as part of `lineage <https://github.com/apriha/lineage>`_.

Acknowledgements
----------------
Thanks to Mike Agostino, Padma Reddy, and `openSNP <https://opensnp.org>`_. Logo composed of
nucleotides from `GRCh38 mitochondrial DNA <https://www.ncbi.nlm.nih.gov/nuccore/NC_012920.1>`_.

.. https://github.com/rtfd/readthedocs.org/blob/master/docs/badges.rst
.. |build| image:: https://travis-ci.org/apriha/snps.svg?branch=master
   :target: https://travis-ci.org/apriha/snps
.. |codecov| image:: https://codecov.io/gh/apriha/snps/branch/master/graph/badge.svg
   :target: https://codecov.io/gh/apriha/snps
.. |docs| image:: https://readthedocs.org/projects/snps/badge/?version=latest
   :target: https://snps.readthedocs.io/
.. |pypi| image:: https://img.shields.io/pypi/v/snps.svg
   :target: https://pypi.python.org/pypi/snps
.. |python| image:: https://img.shields.io/pypi/pyversions/lineage.svg
   :target: https://www.python.org
