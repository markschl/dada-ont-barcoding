# Run DADA2 denoising

Runs DADA2 on a FASTQ file, removes chimeras and and returns an
abundance-sorted data frame of ASVs with some additional information

## Usage

``` r
dada2_denoise(
  x,
  dada_err,
  cores = 1,
  max_members = 1e+06,
  singleton_threshold = 5,
  ...
)
```

## Arguments

- x:

  FASTQ file path or *derep-class* object

- dada_err:

  error information returned by
  [dada_learn_errors](https://markschl.github.io/DadaNanoBC/reference/dada_learn_errors.md)

- max_members:

  do random subsampling of overly large ASVs to the specified number of
  sequences

- singleton_threshold:

  by default, singletons cannot form a new cluster, but for very
  low-depth samples (\< N duplicates of any sequence), DETECT_SINGLETONS
  is turned on (see
  [dada2::setDadaOpt](https://rdrr.io/pkg/dada2/man/setDadaOpt.html)).
  This increases sensitivity with InDel-rich Nanopore data and prevents
  that a large proportion of sequences are discarded.

- ...:

  Arguments passed to
  [dada2::dada](https://rdrr.io/pkg/dada2/man/dada.html)

## Value

a data frame with the following columns:

- *sequence*

- *top_uniques*: abundance of top unique sequences
