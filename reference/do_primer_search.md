# Search and remove primer sequences

Requires [seqtool](https://github.com/markschl/seqtool) v0.4 or higher

See
[do_trim_demux](https://markschl.github.io/DadaNanoBC/reference/do_trim_demux.md)
for more details on the arguments

## Usage

``` r
do_primer_search(
  fq_paths,
  amplicon_primers,
  out_dir,
  primer_max_err = 0.2,
  idx_max_diffs = 0,
  min_barcode_length = 50,
  cores = 1
)
```

## Value

Returns trimmed FASTQ paths (Zstandard-compressed) and statistics:
`list(trimmed_fq = c(...), stats = list(...))`. More FASTQ files are
generated for reads without primers/indexes or reads that are too short
or potential concatenated products.

## Details

The function generates different FASTQ files in `out_dir`, whereby
*trimmed.fastq.zst* contains "valid" reads with sample indexes in the
sequence headers: `<id> fi=<fwd-idx-name> ri=<rev-idx-name>`
