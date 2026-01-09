# Match taxonomic names/lineages against the GBIF taxonomy

This function uses the GBIF web API (through
[rgbif](https://docs.ropensci.org/rgbif/reference/name_backbone_checklist.html))
to retrieve the currently accepted taxon name and lineage according to
the GBIF [backbone
taxonomy](https://www.gbif.org/dataset/d7dddbf4-2cf0-4f39-9b2a-bb099caae36c).

## Usage

``` r
get_gbif_taxa(
  taxa,
  cache_file,
  known_kingdom = NA,
  likely_kingdom = NA,
  verbose = FALSE
)
```

## Arguments

- taxa:

  A character vector of taxonomic names or a matrix/data frame with
  lineages; columns must have rank names and the first column should be
  *name*, containing taxonomic names such as "Russula sp." or
  "Bacteria". The lineages don't need to be complete (NA or only few
  ranks allowed), but their presence may help with the matching.'

- cache_file:

  Path where past searches are stored in a tab-separated text format.
  This avoids havoing to query the GBIF API repeatedly for thousands of
  taxa. Delete the file (or individual lines) to re-query the GBIF API
  for the taxa.

- known_kingdom:

  provide a known kingdom as alternative to providing the kingdom in a
  matrix of lineages. Supplying the kingdom may help to resolve
  ambiguous matches.

- likely_kingdom:

  If the kingdom is not known for sure (but very likely), it can be
  provided here. The GBIF API will be searched without and with kingdom,
  and the result yielding a more complete lineage is returned

## Value

A matrix with lineages; rows are in the same order as in the `taxa`
vector or matrix
