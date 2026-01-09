# Compare consensus/ASV and known sequences

Maps already known sequences or inconsistent ASV sequences against the
consensus (also useful for manual inspection)

## Usage

``` r
compare_seqs(d, bam_out = NULL, tmp_dir = NULL, known_seq = NULL)
```

## Arguments

- data:

  frame returned by
  [infer_barcodes](https://markschl.github.io/DadaNanoBC/reference/infer_barcodes.md)
  (needs the *consensus* column)

## Value

Adds a `consensus_diffs` column to `d` (NA if not compared, Inf if not
mapped due to too many mismatches)
