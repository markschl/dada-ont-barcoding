# DadaNanoBC: Nanopore barcoding pipeline

This R package provides functions for:

- *sample demultiplexing* of base-called Nanopore sequencing data
- inferring specimen *barcode sequences* using a
  [DADA2](https://benjjneb.github.io/dada2)-based [clustering
  strategy](https://markschl.github.io/DadaNanoBC/articles/workflow.html#clustering)
  to accurately resolve haplotypes/polymorphic sequences
- automatic *taxonomic assignents* and recognition (down-ranking) of
  contaminant taxa
- detailed output files (comprehensive Excel table for curation, HTML
  report, BAM alignment files)

A fully-fledged [targets](https://books.ropensci.org/targets) pipeline
connects these features to a comparehensive workflow.

[See here for a comprehensive
explanation](https://markschl.github.io/DadaNanoBC/articles/workflow.html)

## What type of data does it work with?

- *Basecalled* amplicon reads with sample-specific tags attached to both
  ends (dual-indexed)
- Sufficient reads of sufficient quality (R10.4 chemistry)
- Amplicon lengths of 0.5-1.5kb have been successfully tested
- Works with fixed- and variable-length amplicons from different
  taxonomic groups that can be multiplexed and analyzed separately
- different taxonomic database types can be auto-downloaded (more may be
  implemented in the future)

## System requirements

- An UNIX environment (Linux, OS X or Windows with WSL2) with *Bash*
- An [R](https://cran.r-project.org) installation
- Some functions rely on third-party command-line programs:
  [seqtool](https://github.com/markschl/seqtool) for primer
  trimming/demultiplexing, [samtools](https://www.htslib.org) and
  [minimap2](https://github.com/lh3/minimap2) for dealing with
  alignments and [VSEARCH](https://github.com/torognes/vsearch) for
  sequence-based taxonomic assignments (see [installation of
  programs](#installation-of-programs))

The package can always be installed; error messages may appear when
calling functions that rely on missing components.

If
[basecalling](https://nanoporetech.com/document/data-analysis#basecalling-overview)
of the Nanopore data has not yet been done: an Nvidia GPU with Cuda
support (details in [basecalling
tutorial](https://markschl.github.io/DadaNanoBC/basecalling.md)).

## Preparation: lab work and basecalling

1.  [Primer
    design](https://markschl.github.io/DadaNanoBC/articles/primers.html)
    ([example
    report](https://markschl.github.io/DadaNanoBC/primer-design-example.html))
2.  [PCR, library preparation and
    sequencing](https://markschl.github.io/DadaNanoBC/articles/lab.html)
3.  [Basecalling](https://markschl.github.io/DadaNanoBC/articles/basecalling.html)

## Running the pipeline

In the `R` console, install the necessary packages and initialize the
pipeline within the desired directory.

``` r
devtools::install_github("markschl/DadaNanoBC")
install.packages(c("targets", "tarchetypes", "crew"))
setwd("path/to/pipeline")
DadaNanoBC::init_pipeline()
```

Some files are created:

    path/to/pipeline/
     ├─ _targets.R               Pipeline R code; execute with targets::tar_make()
     ├─ infer_barcodes           Script to execute the pipeline from Bash
     └─ analysis/
        ├─ meta-ITS5-ITS4.xlsx   Metadata template (to be modified and renamed to 'meta.xlsx')
        ├─ config.yaml           Pipeline configuration (to be modified)
        └─ reads.fastq.gz        (base-called nanopore sequences to be placed here)

There may also be warnings about missing software (see also
[installation of programs](#installation-of-programs)).

Next:

- Modify the example metadata file `meta-ITS5-ITS4.xlsx` to contain your
  primer sequences and sample metadata, and rename it to `meta.xlsx`
- Copy or move your
  [base-called](https://markschl.github.io/DadaNanoBC/articles/basecalling.html)
  sequencing reads to `reads.fastq.gz`
- Adapt `config.yaml` to your needs; at the minimum, the ‘taxonony’
  section needs to contain the correct taxonomic database.

> *Note*: Consider running the pipeline at reduced read depth to check
> the data maybe adjust some settings afterwards. Inspect `report.html`

> *Note 2* Some settings such as a custom analysis directory, the number
> of parallel workers or program paths are configured using
> `r Sys.setenv(DadaNanoBC_<setting> = ...)`. You may as well edit the
> `_targets.R` file and modify some options below the *Pipeline
> definition* header.

Run the pipeline:

``` r
targets::tar_make()
```

## Output

Some more files will appear in the analysis directory:

- *report.html* ([example
  report](https://markschl.github.io/DadaNanoBC/analysis-example.html)):
  detailed HTML report, useful for troubleshooting
- *report.xlsx*: Excel report with sequences and taxonomic assignments,
  which can be inspected for issues (see
  [curation](https://markschl.github.io/DadaNanoBC/articles/curation.html))
- *top_alignments* (and/or *alignments*, *separate_alignments* depending
  on configuration): aligned sequencing reads in BAM format, may be
  imported into sequence viewing/editing software such as Geneious/CLC
  Workbench/UGENE/IGV, etc.
- *tmp*: contains temporary (but sometimes useful) data such as the
  demultiplexed FASTQ sample files (sometimes *separate_alignments* with
  many small BAM files, depending on the configuration); if not needed
  you may delete the files

After inspecting the HTML report and going through the Excel report,
sequences can be exported as [detailed in the curation
tutorial](https://markschl.github.io/DadaNanoBC/articles/curation.html).

## Running in the bash console

Instead of `targets::tar_make()`, there is also a script to execute the
pipeline from the Bash console. It expects the path to the analysis
directory containing `meta.xlsx`, `reads.fastq.gz` and `config.yaml`,
and optionally also the number of workers:

``` sh
./infer_barcodes analysis_dir [workers]
```

This example runs the analysis in the default *analysis* directory with
8 parallel workers:

``` sh
cd path/to/pipeline
./infer_barcodes analysis 8
```

## Installation of programs

The above call to
[`init_pipeline()`](https://markschl.github.io/DadaNanoBC/reference/init_pipeline.md)
checks for available software and gives the URLs of missing software.
The following does the same:

``` r
DadaNanoBC::check_system_requirements()
```

The message points to the download locations of missing tools.
*Samtools*, *minimap2* and *VSEARCH* are available from different
sources:

*Debian/Ubuntu* (including Windows WSL2 with these OSes):

``` sh
sudo apt install -y samtools minimap2 vsearch
```

*OS X with [Homebrew](https://brew.sh):*

``` sh
brew install samtools minimap2 vsearch
```

*Conda package manager*
(e.g. [Conda-forge](https://conda-forge.org/download)) (run in the [bash
console](#running-in-the-bash-console)):

``` sh
conda install -y -n DadaNanoBC -c bioconda samtools minimap2 vsearch
conda activate DadaNanoBC
#./infer_barcodes analysis_dir [workers]
```

Instructions for downloading *Seqtool* are given on the [Github
site](https://github.com/markschl/seqtool/releases).

### Custom program location

R needs to know about the tool’s locations. If not placed at a standard
location, either the `PATH` environment variable needs adjustment, or
the paths can be communicated to DadaNanoBC.

In R (example with VSEARCH):

``` r
Sys.setenv(DadaNanoBC_vsearch = '~/Downloads/vsearch-2.30.2-linux-x86_64/bin/vsearch')
tar_make()
```

In the Bash console:

``` sh
cd path/to/pipelinea
export DadaNanoBC_vsearch=~/Downloads/vsearch-2.30.2-linux-x86_64/bin/vsearch
./infer_barcodes analysis 8
```

## Useful references

- [Srivathsan & Meier
  (2024)](https://doi.org/10.1007/978-1-0716-3581-0_14)
- [Hebert et al. (2024)](https://doi.org/10.1111/1755-0998.14028)
  barcoded up to 100k samples with a single MinION run.
- [ONT barcoding efforts by the MycoMap
  network](https://www.researchgate.net/publication/393048029_Approaching_Full-Scale_DNA_Barcoding_for_North_American_Macrofungi_Highlights_from_the_MycoMap_Network)
  (see also [methods
  description](https://dx.doi.org/10.17504/protocols.io.dm6gpbm88lzp/v4))
  (software used: [NGSpeciesID](https://github.com/ksahlin/NGSpeciesID))
