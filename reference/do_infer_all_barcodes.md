# Infer barcode sequences for all samples

High-level function calling
[infer_barcodes](https://markschl.github.io/DadaNanoBC/reference/infer_barcodes.md)
on every sample, Parallel processing possible with `cores` \> 1. Custom
parallel processing frameworks can be plugged in by providing
`parallel_lapply_fn`.

## Usage

``` r
do_infer_all_barcodes(
  seq_tab,
  dada_err,
  aln_out = NULL,
  tmp_dir = NULL,
  parallel_lapply_fn = NULL,
  ...,
  cores = 1
)
```

## Arguments

- seq_tab:

  Sequence metadata table returned by
  [do_trim_demux](https://markschl.github.io/DadaNanoBC/reference/do_trim_demux.md)
  or
  [do_demux](https://markschl.github.io/DadaNanoBC/reference/do_demux.md)

- dada_err:

  Result of
  [dada_learn_errors](https://markschl.github.io/DadaNanoBC/reference/dada_learn_errors.md)

- aln_out:

  Optional output directory for BAM alignments (none saved if `NULL`)

- tmp_dir:

  Optional temporary directory (e.g. [RAM
  drive](https://en.wikipedia.org/wiki/RAM_drive) such as
  `/run/user/1000/DadaNanoBC` on Linux, for faster processing)

## Value

Input table (`seq_tab`) with a new column `clustering`, which is a list
of data frames as returned by
[infer_barcodes](https://markschl.github.io/DadaNanoBC/reference/infer_barcodes.md)
