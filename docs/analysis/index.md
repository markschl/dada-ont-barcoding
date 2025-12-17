# Sequence analysis tutorial

[See here](workflow.md) for a description of the pipeline.

## System requirements

- UNIX environment (Linux, OS X or Windows with WSL2) with Bash console
- [R](https://cran.r-project.org)
  ([RStudio](https://posit.co/download/rstudio-desktop) for running the analysis interactively)
- A few external programs are needed (see [installation](#installation))
- If [basecalling](https://nanoporetech.com/document/data-analysis#basecalling-overview) of the Nanopore data
  has not yet been done: an Nvidia GPU with Cuda support (details in [basecalling tutorial](basecalling.md))


### Installation

#### Code download

On the GitHub page https://github.com/markschl/dada-ont-barcoding, go to *Code* → *Download ZIP* → extract to some folder. Or, run `git clone https://github.com/markschl/dada-ont-barcoding.git`

#### Dependencies

The following standard tools are needed:
[samtools](https://www.htslib.org) and [minimap2](https://github.com/lh3/minimap2) for dealing with alignments,
and [VSEARCH](https://github.com/torognes/vsearch) for sequence-based taxonomic assignments.
Furthermore, [seqtool](https://github.com/markschl/seqtool) is needed for for demultiplexing.

If you plan to run the analysis interactively in RStudio, the easiest may be a
*system-wide installation* for at least part of the tools (not possible for *seqtool*).

Debian/Ubuntu (including Windows WSL2 with these OSes)

```sh
sudo apt install -y samtools minimap2 vsearch
```

OS X with [Homebrew](https://brew.sh):

```sh
brew install samtools minimap2 vsearch
```

> *Alternatively*, the conda package manager (e.g. [Conda-forge](https://conda-forge.org/download))
> can be used if [running non-interactively in Bash](#in-the-bash-console) (RStudio may be difficult).
> 
> ```sh
> conda install -y -n dada-ont-barcoding -c bioconda samtools minimap2 vsearch
> conda activate dada-ont-barcoding
> ```

For all the missing tools (including *seqtool*) there is an *installation script*,
which places them in the `bin` subdirectory. *Samtools* and *minimap2* need to be compiled
from source (UNIX build environment required).
Run the following in a Bash console:

```sh
cd path/to/dada-ont-barcoding
scripts/install-deps.sh
```

### Basecalling

The pipeline requires [base-called](https://nanoporetech.com/platform/technology/basecalling)
reads in the FASTQ format (GZIP-compressed) as input, placed in `analysis/<run_name>/reads.fastq.gz`.
If starting from the raw squiggle data, follow [this tutorial](basecalling.md)

## Metadata and configuration

- Assemble a *meta.xlsx* file with correct primer sequences and plate layout,
  as well as a *config.yaml* file.
  (example/templates here: [ITS5-ITS4](https://github.com/markschl/dada-ont-barcoding/templates/analysis/ITS5-ITS4)).
- Add sample information for the current Nanopore run according to the instructions
  in the *information* worksheet
- Modify *config.yaml* to your needs [example with comments here](https://github.com/markschl/dada-ont-barcoding/templates/analysis/config.commented.yaml)

The directory structure:

```
ont-barcoding/
 ├─ analysis.Rmd               Main R-Markdown document (edit e.g. in RStudio)
 ├─ run_analysis.R             Run the analysis without RStudio
 ├─ analysis/
 │  └─ <run_name>            <- assign a name
 │    ├─ config.yaml         <- required configuration file
 │    ├─ meta.xlsx           <- required metadata file
 │    ├─ reads.fastq.gz      <- base called sequences
 │    ├─ output
 │    │  ├─ report.xlsx      -> results file with sequences
 │    │  ├─ alignments/...       -> detailed alignments (for curation)
 │    └─ analysis.pdf/html       -> analysis report (placed here by run_analysis.R)
 ├─ taxdb/                       -> taxonomy database (auto-downloaded)
```


## Run the analysis

### In RStudio

- Open [`analysis.Rmd`](https://github.com/markschl/dada-ont-barcoding/analysis.Rmd) in RStudio
- Edit the run name in the header (`params: <run_name>`). 
  It has to be identical to the prepared run directory (`analysis/<run_name>`).
- Run all chunks, check the results, maybe adjust settings in `config.yaml`

This may take from 15-20 minutes for ~3M reads up to a one hour for a larger MinION run
with ~15M reads (make sure to set multiple `cores` in `config.yaml`).
Re-running will be faster thanks to the intermediate files in `analysis/<run_name>/tmp`.

> Intermediate files in `tmp` (such as the `demux` folder or `cluster.rds`)
> need to be deleted manually if settings are changed.

> Consider generating a HTML/PDF report (Knit button) and moving it to `analysis/<run_name>`.

### In the Bash console

```sh
cd path/to/dada-ont-barcoding
./run_analysis.R run_name
# or
# Rscript run_analysis.R run_name
```

The analysis report will appear in `analysis/<run_name>/analysis.html`, along with the other output files.

> The disadvantage of this strategy: it is not possible to e.g. flexibly adjust quality filtering
> settings (or others) in `config.yaml` depending on the data.
> 
> Analysis.Rmd *may* be split into two scripts in future versions.

### Final steps

- [Manual curation](../curation.md)
- Delete `analysis/<run_name>/tmp` once the analysis is done (it may take a lot of space).
