# Parse a data frame with primer information

Parse a data frame with primer information

## Usage

``` r
parse_primer_tab(primer_tab, amplicons = NULL)
```

## Arguments

- primer_tab:

  Primer table (data frame-like; see details)

- amplicons:

  Character vector of amplicons pooled in the library; their names
  should be in the form `forward-reverse` (primer names should not
  contain dashes)

## Value

a named list of amplicon metadata (named by amplicons =
`forward-reverse`). Each amplicon-specific item is another list with
*two entries* (for forward and reverse), which are again lists with
following entries:

- *primer*: named vector of primer sequences

- *index*: named vector of index sequences

- *index_len*: the length of the index sequences (no variable length
  currently allowed, but lengths may vary between forward and reverse
  sample indexes)

## Details

Columns of `primer_tab`:

- *primer-index*: Primer name and sample index name, delimited by a
  dash; the primer names **must not contain dashes** and both *primer*
  and *index* should correspond to the names in the *amplicon* column of
  the sample table (see
  [parse_sample_tab](https://markschl.github.io/DadaNanoBC/reference/parse_sample_tab.md)).

- *index_seq* sample index sequence

- *primer_seq* primer sequence (repeated)

If `amplicon` is not provided, the forward/reverse primers are assumed
to be in order. `amplicon` may be obtained from the *amplicon* column of
the sample table (see
[parse_sample_tab](https://markschl.github.io/DadaNanoBC/reference/parse_sample_tab.md)).
