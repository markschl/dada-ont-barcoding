# Infer barcode sequences for all samples

High-level function calling
[infer_barcode](https://markschl.github.io/DadaNanoBC/reference/infer_barcode.md)
and
[compare_seqs](https://markschl.github.io/DadaNanoBC/reference/compare_seqs.md)
on every sample. Parallel processing possible with `cores` \> 1. Custom
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

  Optional temporary directory (default: user-specific temporary
  directory, see
  [set_global_opts](https://markschl.github.io/DadaNanoBC/reference/set_global_opts.md))

## Value

Input table (`seq_tab`) with a new column `clustering`, which is a list
of data frames as returned by
[infer_barcode](https://markschl.github.io/DadaNanoBC/reference/infer_barcode.md)
