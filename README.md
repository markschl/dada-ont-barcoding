# DadaNanoBC: Nanopore barcoding pipeline

This R package provides functions for:

- *sample demultiplexing* of base-called Nanopore sequencing data 
- inferring specimen *barcode sequences* using a [DADA2](https://benjjneb.github.io/dada2)-based
  strategy to accurately resolve haplotypes/polymorphic sequences
- automatic *taxonomic assignents* and recognition (down-ranking) of contaminant taxa
- detailed output files (comprehensive Excel table for curation, HTML report, BAM alignment files)

A fully-fledged [targets](https://books.ropensci.org/targets) pipeline connects
these features to a comparehensive workflow.

## What type of data does it work with?

- *Base-called* amplicon reads with sample-specific tags attached to both ends
  (dual-indexed)
- Sufficient reads of sufficient quality (R10.4 chemistry)
- Amplicon lengths of 0.5-1.5kb have been successfully tested
- Works with fixed- and variable-length amplicons from different taxonomic groups
  taht can be multiplexed and analyzed separately
- different taxonomic database types can be auto-downloaded
  (more may be implemented in the future)

## System requirements

- An UNIX environment (Linux, OS X or Windows with WSL2) with *Bash*
- An [R](https://cran.r-project.org) installation
- Some functions rely on third-party command-line programs:
  [seqtool](https://github.com/markschl/seqtool) for primer trimming/demultiplexing,
  [samtools](https://www.htslib.org) and [minimap2](https://github.com/lh3/minimap2)
  for dealing with alignments and 
  [VSEARCH](https://github.com/torognes/vsearch) for sequence-based taxonomic assignments
  (see [installation of programs](#installation-of-programs))

The package can always be installed; error messages may appear when calling
functions that rely on missing components.

If [basecalling](https://nanoporetech.com/document/data-analysis#basecalling-overview)
of the Nanopore data has not yet been done: an Nvidia GPU with Cuda support
(details in [basecalling tutorial](basecalling.md)).


## Quickstart

Install the package:

```r
devtools::install_github("markschl/DadaNanoBC")
```

Copy the pipeline templates to the current directory

```r
# optional: setwd("path/to/some/directory")
...
```


## Installation of programs


(...)
