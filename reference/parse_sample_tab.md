# Parse a data frame with sample metadata

Parse a data frame with sample metadata

## Usage

``` r
parse_sample_tab(sample_tab)
```

## Arguments

- sample_tab:

  The sample table (data frame-like, see details)

## Value

a valid sample metadata table (data frame) with sample names
de-duplicated (if necessary)

## Details

Required columns:

- Plate

- *plate*: A01-H12)

- *amplicon*: `forward-reverse`, name as in 'primer_index' column of
  primer table; see
  [parse_primer_tab](https://markschl.github.io/DadaNanoBC/reference/parse_primer_tab.md)

- *indexes*: `forward-reverse`, name as in 'primer_index' column of
  primer table; see
  [parse_primer_tab](https://markschl.github.io/DadaNanoBC/reference/parse_primer_tab.md)

- *sample*: sample name, should not be duplicated

- *sample_type*: Sample type such as 'negative control', will show up in
  reports; can be left empty

Optional columns:

- *taxon*: Known or suspected taxon (any rank), will be compared to the
  sequence-based identification and used for detecting contamination

- *known sequence*: Already known sequence (if any), e.g. from a
  previous Sanger sequencing; will be compared with the Nanopore-derived
  sequence

### Note on amplicon multiplexing

Primers are searched in the order that amplicons appear in the sample
table With *nested amplicons*, the shorter one should be placed at the
end. Also place amplicons with *little data* before other amplicons with
similar primers to make sure that they don't get "swallowed" in case of
unspecific primer or sample index matching.
