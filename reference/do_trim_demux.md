# Search and remove primers and sample indexes

do_trim_demux is a wrapper function running
[do_primer_search](https://markschl.github.io/DadaNanoBC/reference/do_primer_search.md)
and
[do_demux](https://markschl.github.io/DadaNanoBC/reference/do_demux.md)
sequentially.

## Usage

``` r
do_trim_demux(
  fq_paths,
  amplicon_primers,
  sample_tab,
  out_dir,
  primer_max_err = 0.2,
  idx_max_diffs = 0,
  min_barcode_length = 50,
  error_threshold = 2.5,
  cores = 1,
  keep_trimmed = FALSE
)
```

## Arguments

- fq_paths:

  Character vector of FASTQ file paths (compressed or uncompressed)

- amplicon_primers:

  Amplicon specification, as returned by
  [parse_primer_tab](https://markschl.github.io/DadaNanoBC/reference/parse_primer_tab.md)

- sample_tab:

  Sample table as returned by
  [parse_sample_tab](https://markschl.github.io/DadaNanoBC/reference/parse_sample_tab.md)

- out_dir:

  Output directory containing the demultiplexed files
  (`fprimer_fidx-rprimer_ridx.fastq.gz`)

- primer_max_err:

  Maximum allowed error rate in primers (edit distance =
  substitutions/InDels). Default: 0.2 (20%)

- idx_max_diffs:

  Maximum allowed differences in forward/reverse sample index sequences;
  don't set too high to avoid tag switching ("index hopping") resulting
  in a "background noise" of unspecific sequences. Default: 0
  differences

- min_barcode_length:

  Remove trimmed sequences shorter than `min_barcode_length`; make sure
  that the value is shorter than your expected minimal barcode length
  for *all* amplicons (default: 50 bp)

- error_threshold:

  Maximum sequencing errors allowed per sequence (as estimated from
  quality scores). This threshold impacts how many sequences pass the
  filter in (see
  [do_demux](https://markschl.github.io/DadaNanoBC/reference/do_demux.md)).
  Setting it too low may remove too many sequences, but setting it too
  high may lead to noisy/errorneous results. The default value of 2.5
  seems to work well with R10.4 data and ~0.5-1.5 kb amplicons, removing
  about half of the trimmed sequences.

## Value

A list of:

- `seq_tab`: `sample_tab` with a `reads_path` column and additional rows
  for samples found in the reads but not in the sample table

- `trim_stats`: a list containing primer/index search statistics (used
  in the HTML report)

## Details

Requires [seqtool](https://github.com/markschl/seqtool) v0.4 or higher
installed either system-wide (in `PATH`) or locally (provide path with
`set_program_path('seqtool', 'path/to/st')`).
