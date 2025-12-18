# Nanopore barcoding pipeline

This R-based pipeline infers *specimen barcodes* from Nanopore sequencing data
with multiplexed samples.

It uses a [DADA2](https://benjjneb.github.io/dada2)-based
[clustering strategy](analysis/workflow.md#details-on-the-clustering) with some refinements with the aim to resolve
haplotypes and other sequence polymorphisms as well as possible.

The DADA2 denoising method requires that each amplicon is supported by at least a
few error-free sequencing reads, which is usually the case with the latest R10.4 chemistry,
given that sequences are not overly long (tested with 500-1000 bp amplicons).
In case of low sequencing depths (with low replication), a fixed-threshold clustering
procedure is applied instead.

## Features

* Recognizes and removes primers with sample-specific tags for one or **multiple pooled amplicons**
  and does quality filtering
* Automatic **taxonomic assignments** based on user-defined reference databases
* Matching of morphological identifications with the sequence-based taxonomy
  and a procedure to automatically **recognize/remove contaminants**
* Detailed **output**:
    - *detailed processing report* (generated from `analysis.Rmd`)
      ([example](analysis-example.html))
    - *Excel table* with the *consensus sequences*, taxonomy and possible issues
      (see [curation](curation.md))
    - *BAM/FASTA alignment files* that can be viewed during
      [manual curation](curation.md)

[More on the clustering workflow and other processing steps](analysis/workflow.md)

## Documentation of full barcoding workflow

Topics:

1. [Primer design](primers.md)
   ([example report](primer-design-example.html))
2. [PCR, library preparation and sequencing](lab.md)
3. [Infer barcode sequences](analysis) using this pipeline
   ([example report](analysis-example.html))
4. [Manual curation](curation.md)

## Useful references

- [Srivathsan & Meier (2024)](https://doi.org/10.1007/978-1-0716-3581-0_14)
- [Hebert et al. (2024)](https://doi.org/10.1111/1755-0998.14028) barcoded up to 100k samples with a single MinION run.
- [ONT barcoding efforts by the MycoMap network](https://www.researchgate.net/publication/393048029_Approaching_Full-Scale_DNA_Barcoding_for_North_American_Macrofungi_Highlights_from_the_MycoMap_Network) (see also [methods description](https://dx.doi.org/10.17504/protocols.io.dm6gpbm88lzp/v4))
