# Attempts splitting sequence cluster into two haplotypes.

Attempts splitting sequence cluster into two haplotypes.

## Usage

``` r
try_split_haplotypes(
  reads,
  d,
  prefix,
  min_identical = 4,
  max_ratio = 3,
  fast = FALSE,
  cores = 1,
  ...,
  samtools = "samtools",
  minimap2 = "minimap2",
  verbose = FALSE
)
```

## Arguments

- reads:

  XStringQuality object

- d:

  a data frame as returned by
  [dada2_denoise](https://markschl.github.io/DadaNanoBC/reference/dada2_denoise.md)
  and `ambig_consensus` with these columns:

  - id

  - top_uniques

  - consensus_ambigs

- prefix:

  Path prefix where the .bam alignment and .fasta reference is found
  (may be overwritten!)

- min_identical:

  minimum number of *identical* sequences supporting *both* of the
  dominant unique sequences

- max_ratio:

  maximum abundance ratio (mapped read abundances) for splitting

- fast:

  use faster compression (larger BAM file)

## Details

Splitting is done if the sum of the consensus ambiguities of the two
resulting sequence variants is smaller than the number of consensus
ambiguities of the parent cluster.

Overwrites the existing BAM files and FASTA references in 'prefix'.
