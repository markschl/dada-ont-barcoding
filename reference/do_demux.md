# Group sequences by sample index combination

Group sequences by sample index combination

## Usage

``` r
do_demux(
  primer_search_fq,
  out_dir,
  sample_tab,
  error_threshold = 2.5,
  min_barcode_length = 50
)
```

## Arguments

- primer_search_fq:

  One or multiple paths from
  [do_primer_search](https://markschl.github.io/DadaNanoBC/reference/do_primer_search.md)
  (`trimmed_fq` entry in the returned list).

## Value

`sample_tab` with a `reads_path` column and additional rows for samples
found in the reads but not in the sample table

## Details

See
[do_trim_demux](https://markschl.github.io/DadaNanoBC/reference/do_trim_demux.md)
for details on the other arguments

The read files in `out_dir` are named as follows:
`fprimer-fidx--rprimer-ridx.fastq.gz`

Requires [seqtool](https://github.com/markschl/seqtool) v0.4 or higher
