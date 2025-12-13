# Nanopore barcoding pipeline

This R-based pipeline infers *specimen barcodes* Nanopore sequencing data
with multiplexed samples.

It uses a [DADA2](https://benjjneb.github.io/dada2)-based
[clustering strategy](analysis/workflow.md#details-on-the-clustering)
and reports consensus sequences, separately for haplotypes and other sequence
polymorphisms as well as possible.
This requires that each amplicon is supported by at least a few error-free sequencing reads,
which is usually the case with the latest R10.4 chemistry, given that sequences
are not overly long (tested with 500-1000 bp amplicons).
In case of low sequencing depths (with low replication), a fixed-threshold clustering
procedure is applied instead.

**Other features**

- Automatic *taxonomic assignments* (currently implemented for the ITS marker)
  for validating the morphological identifications and auto-filtering of contaminants
- Detailed output:
    * *detailed processing report* (generated from `analysis.Rmd`) ([example](analysis-example.html))
    * *Excel table* with detailed overview of sequences and possible issues (see [curation](curation.md))
    * *BAM/FASTA files* that assist with [manual curation](curation.md)

[More on the clustering workflow and other processing steps](analysis/workflow.md)

## Documentation of full barcoding workflow

1. [Primer design](primers.md) ([example report](primer-design-example.html))
1. [PCR, library preparation and sequencing](lab.md)
2. [Infer barcode sequences](analysis) using this pipeline ([example report](analysis-example.html))
3. [Manual curation](curation.md)

## Useful references

- [Srivathsan & Meier (2024)](https://doi.org/10.1007/978-1-0716-3581-0_14)
- [Hebert et al. (2024)](https://doi.org/10.1111/1755-0998.14028) barcoded up to 100k samples with a single MinION run.
- [ONT barcoding efforts by the MycoMap network](https://www.researchgate.net/publication/393048029_Approaching_Full-Scale_DNA_Barcoding_for_North_American_Macrofungi_Highlights_from_the_MycoMap_Network) (see also [methods description](https://dx.doi.org/10.17504/protocols.io.dm6gpbm88lzp/v4))
