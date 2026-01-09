# Run the taxonomy assignment and compare with provided taxonomic names of specimens

This function runs `do_assign_taxonomy` and `do_compare_taxonomy`
sequentially (functions not documented).

## Usage

``` r
do_assign_compare_taxonomy(
  seq_tab,
  db_file,
  gbif_cache_file,
  tmp_dir = NULL,
  confidence_threshold = 0.8,
  summary_ranks = NULL,
  known_contaminants = NULL,
  likely_kingdom = NULL,
  contam_rank_delta = 3,
  cores = 1
)

do_assign_taxonomy(
  seq_tab,
  db_file,
  tmp_dir = NULL,
  confidence_threshold = 0.8,
  summary_ranks = NULL,
  cores = 1
)

do_compare_taxonomy(
  seq_tab,
  seq_lineages,
  gbif_cache_file,
  known_contaminants = NULL,
  likely_kingdom = NULL,
  contam_rank_delta = 3
)
```

## Arguments

- seq_tab:

  Sequence table as returned by
  [do_infer_all_barcodes](https://markschl.github.io/DadaNanoBC/reference/do_infer_all_barcodes.md)
  or downstream functions

- db_file:

  Sequence database with UTAX-formatted headers (optionally compressed);
  see
  [load_taxdb](https://markschl.github.io/DadaNanoBC/reference/load_taxdb.md)

- gbif_cache_file:

  Path to cache file for saving GBIF taxonomy lookup results (may be
  shared across datasets) (see
  [get_gbif_taxa](https://markschl.github.io/DadaNanoBC/reference/get_gbif_taxa.md))

- tmp_dir:

  optional temporary directory for saving intermediate sequence files

- confidence_threshold:

  SINTAX bootstrap threshold (default: 0.8 / 80%)

- summary_ranks:

  optional named character vector with ranks that should show up in the
  lineages (default: determine automatically)

- known_contaminants:

  nested list specifying known contaminant names, such as
  `list(family=c('Aspergillaceae', ...), genus=c(...))`

- likely_kingdom:

  kingdom name providing some guidance for the GBIF name search (see
  [get_gbif_taxa](https://markschl.github.io/DadaNanoBC/reference/get_gbif_taxa.md))

- contam_rank_delta:

  Require at least N additional ranks to be matching between the
  provided (e.g. morphology-based) and the auto-assigned sequence-based
  taxonomic lineages for a taxon to be determined as the "correct" taxon
  in the mix, and the top taxon being down-ranked. This setting affects
  the "sensitivity" of the contaminant detection (lower: more sensitive,
  but also more false classifications, higher: some contamination may
  not be recognized). The default of 3 seems to work well in most cases.

## Value

`do_assign_compare_taxonomy` returns `seq_tab`, with the nested data
frames in `clustering` **reordered** as such that dominant taxa flagged
as contamination are moved to the bottom. The following new columns are
added:

- *unique_id*: unique sequence id in the form: `u1`, `u2`, etc.

- *taxon*: auto-inferred taxon name (*may not be correct, confirm with
  BLAST or other identification methods*)

- *short_lineage*: taxonomic lineage of the auto-inferred taxon
  (provided or auto-inferred `summary_ranks` shown)

- *is_contaminant*: (logical) `TRUE` for taxa recognized as
  contamination

- *matching_ranks*: Number of matching taxonomic ranks in the comparison
  between the provided and auto-assigned sequence-based taxonomy
  (usually: kingdom/domain, phylum, class, order, family, genus,
  species). The comparison can only be done up to the highest defined
  rank in either lineage.

- *mismatching_ranks*: Number of non-matching ranks in the comparison
  between the provided and auto-assigned sequence-based taxonomy

- *unspecific*: `TRUE` for all taxa that are not at the top in the
  clusters table after the contamination ranking

The main `seq_tab` table further contains some columns with data from
the top sequence: *top_seq_matching_ranks*, *top_seq_mismatching_ranks*,
*top_seq_taxon*, *top_seq_short_lineage*

## Details

Sequences are flagged as contamination (*is_contaminant*) if:

- there exists a less abundant taxon with the correct kingdom (compared
  to *taxon*, GBIF name comparison guided by `likely_kingdom`)

- a taxon has has at least `contam_rank_delta` more consistent
  (matching) taxonomic ranks in the GBIF lineage than the top taxon
  (explanation for *matching_ranks* above); undefined ranks are excluded
  from the comparison
