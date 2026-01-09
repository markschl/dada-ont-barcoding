# Create an Excel report

Create an Excel report

## Usage

``` r
create_excel_report(
  seq_tab,
  outfile,
  low_abund_threshold = 20,
  min_seqs_unknown = 10,
  n_curate = 4,
  n_show_depth = 6
)
```

## Arguments

- seq_tab:

  data frame returned by
  [do_assign_taxonomy](https://markschl.github.io/DadaNanoBC/reference/do_taxonomy.md)

- outfile:

  the Excel output file

- low_abund_threshold:

  add 'low-coverage' to the issues list for samples with less reads

- n_curate:

  max. number of sequences (polymorphism/haplotypes) from the top taxon
  to list in the *curation* section

- n_show_depth:

  max. number of sequences (polymorphism/haplotypes) from the top taxon
  for which to output read numbers

- min_seqs:

  ignore samples found in reads (index combination) but not in the
  samples table if they have less than the given number of reads
