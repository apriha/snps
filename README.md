![snps](https://raw.githubusercontent.com/apriha/snps/main/docs/images/snps_banner.png)

[![CI](https://github.com/apriha/snps/actions/workflows/ci.yml/badge.svg?branch=main)](https://github.com/apriha/snps/actions/workflows/ci.yml)
[![codecov](https://codecov.io/gh/apriha/snps/branch/main/graph/badge.svg)](https://codecov.io/gh/apriha/snps)
[![docs](https://readthedocs.org/projects/snps/badge/?version=stable)](https://snps.readthedocs.io/)
[![pypi](https://img.shields.io/pypi/v/snps.svg)](https://pypi.python.org/pypi/snps)
[![python](https://img.shields.io/pypi/pyversions/snps.svg)](https://www.python.org)
[![downloads](https://pepy.tech/badge/snps)](https://pepy.tech/project/snps)
[![Ruff](https://img.shields.io/endpoint?url=https://raw.githubusercontent.com/astral-sh/ruff/main/assets/badge/v2.json)](https://github.com/astral-sh/ruff)

# snps

tools for reading, writing, generating, merging, and remapping SNPs 🧬

*`snps` strives to be an easy-to-use and accessible open-source library for working with
genotype data*

## Features

### Input / Output

- Read raw data (genotype) files from a variety of direct-to-consumer (DTC) DNA testing
  sources with a [SNPs](https://snps.readthedocs.io/en/stable/snps.html#snps.snps.SNPs)
  object
- Read and write VCF files (e.g., convert [23andMe](https://www.23andme.com) to VCF)
- Merge raw data files from different DNA tests, identifying discrepant SNPs in the process
- Read data in a variety of formats (e.g., files, bytes, compressed with `gzip` or `zip`)
- Handle several variations of file types, historically validated using
  data from [openSNP](https://opensnp.org)
- Generate synthetic genotype data for testing and examples

### Build / Assembly Detection and Remapping

- Detect the build / assembly of SNPs (supports builds 36, 37, and 38)
- Remap SNPs between builds / assemblies

### Data Cleaning

- Perform quality control (QC) / filter low quality SNPs based on [chip clusters](https://doi.org/10.1016/j.csbj.2021.06.040)
- Fix several common issues when loading SNPs
- Sort SNPs based on chromosome and position
- Deduplicate RSIDs
- Deduplicate alleles in the non-PAR regions of the X and Y chromosomes for males
- Deduplicate alleles on MT
- Assign PAR SNPs to the X or Y chromosome

### Analysis

- Derive sex from SNPs
- Detect deduced genotype / chip array and chip version based on [chip clusters](https://doi.org/10.1016/j.csbj.2021.06.040)
- Predict ancestry from SNPs (when installed with [ezancestry](https://github.com/arvkevi/ezancestry))

## Supported Genotype Files

`snps` supports [VCF](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3137218/) files and
genotype files from the following DNA testing sources:

- [23andMe](https://www.23andme.com)
- [23Mofang](https://www.23mofang.com)
- [Ancestry](https://www.ancestry.com)
- [CircleDNA](https://circledna.com/)
- [Código 46](https://codigo46.com.mx)
- [DNA.Land](https://dna.land)
- [Family Tree DNA](https://www.familytreedna.com)
- [Genes for Good](https://genesforgood.sph.umich.edu)
- [LivingDNA](https://livingdna.com)
- [Mapmygenome](https://mapmygenome.in)
- [MyHeritage](https://www.myheritage.com)
- [PLINK](https://www.cog-genomics.org/plink/)
- [Sano Genetics](https://sanogenetics.com)
- [SelfDecode](https://selfdecode.com)
- [tellmeGen](https://www.tellmegen.com)

Additionally, `snps` can read a variety of "generic" CSV and TSV files.

## Dependencies

`snps` requires [Python](https://www.python.org) 3.9+ and the following Python
packages:

- [numpy](http://www.numpy.org)
- [pandas](http://pandas.pydata.org)
- [atomicwrites](https://github.com/untitaker/python-atomicwrites)

## Installation

`snps` is [available](https://pypi.org/project/snps/) on the
[Python Package Index](https://pypi.org). Install `snps` (and its required
Python dependencies) via `pip`:

```bash
$ pip install snps
```

For [ancestry prediction](https://snps.readthedocs.io/en/stable/snps.html#snps.snps.SNPs.predict_ancestry)
capability, `snps` can be installed with [ezancestry](https://github.com/arvkevi/ezancestry):

```bash
$ pip install snps[ezancestry]
```

## Examples

To try these examples, first generate some sample data:

```python
>>> from snps.resources import Resources
>>> paths = Resources().create_example_datasets()
```

### Load a Raw Data File

Load a raw data file exported from a DNA testing source (e.g.,
[23andMe](https://www.23andme.com), [AncestryDNA](https://www.ancestry.com),
[Family Tree DNA](https://www.familytreedna.com)):

```python
>>> from snps import SNPs
>>> s = SNPs("resources/sample1.23andme.txt.gz")
```

`snps` automatically detects the source format and [normalizes](https://snps.readthedocs.io/en/stable/snps.html#snps.snps.SNPs.snps) the data:

```python
>>> s.source
'23andMe'
>>> s.count
991767
>>> s.build
37
>>> s.assembly
'GRCh37'
```

The SNPs are available as a `pandas.DataFrame`:

```python
>>> df = s.snps
>>> df.columns.tolist()
['chrom', 'pos', 'genotype']
>>> len(df)
991767
```

### Merge Raw Data Files

Combine SNPs from multiple files (e.g., combine data from different testing companies):

```python
>>> results = s.merge([SNPs("resources/sample2.ftdna.csv.gz")])
>>> s.count
1006949
```

SNPs are compared during the merge. Position and genotype discrepancies are identified and
can be inspected via properties of the `SNPs` object:

```python
>>> len(s.discrepant_merge_positions)
27
>>> len(s.discrepant_merge_genotypes)
156
```

### Remap SNPs

Convert SNPs between genome assemblies (Build 36/NCBI36, Build 37/GRCh37, Build 38/GRCh38):

```python
>>> chromosomes_remapped, chromosomes_not_remapped = s.remap(38)
>>> s.assembly
'GRCh38'
```

### Save SNPs

Save SNPs to common file formats:

```python
>>> _ = s.to_tsv("output.txt")
>>> _ = s.to_csv("output.csv")
```

To save as VCF, `snps` automatically downloads the required reference sequences for the
assembly. This ensures the REF alleles in the VCF are accurate:

```python
>>> _ = s.to_vcf("output.vcf")  # doctest: +SKIP
```

All output files are saved to the [output directory](https://snps.readthedocs.io/en/stable/output_files.html).

### Generate Synthetic Data

Generate synthetic genotype data for testing, examples, or demonstrations:

```python
>>> from snps.io import SyntheticSNPGenerator
>>> gen = SyntheticSNPGenerator(build=37, seed=123)
>>> gen.save_as_23andme("synthetic_23andme.txt.gz", num_snps=10000)
'synthetic_23andme.txt.gz'
```

The generator supports multiple output formats (23andMe, AncestryDNA, FTDNA) and
automatically injects build-specific marker SNPs to ensure accurate build detection.

## Documentation

Documentation is available [here](https://snps.readthedocs.io/).

## Acknowledgements

Thanks to Mike Agostino, Padma Reddy, Kevin Arvai, [Open Humans](https://www.openhumans.org),
and [Sano Genetics](https://sanogenetics.com). This project was historically validated using
data from [openSNP](https://opensnp.org).

`snps` incorporates code and concepts generated with the assistance of various
generative AI tools (including but not limited to [ChatGPT](https://chatgpt.com),
[Grok](https://grok.com), and [Claude](https://claude.ai)). ✨

## License

`snps` is licensed under the [BSD 3-Clause License](https://github.com/apriha/snps/blob/main/LICENSE.txt).
