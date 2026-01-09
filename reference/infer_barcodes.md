# Infer barcode sequences from raw reads

Applies a clustering workflow to the provided raw reads, based on DADA2
denoising and/or fixed-threshold clustering to infer the correct barcode
sequence(s). Further groups the resulting sequences by putative taxon.

## Usage

``` r
infer_barcodes(
  fq,
  dada_err,
  alignment_prefix,
  id_prefix = NULL,
  tmp_dir = NULL,
  dada_omega_a = c(1e-20, 1e-10, 0.01),
  omegaA_iter_threshold = 1000,
  dada_min_identical = 2,
  dada_min_n0 = 4,
  min_seq_abund = 3,
  max_sample_depth = 5000,
  consensus_max_depth = 3000,
  consensus_threshold = 0.65,
  consensus_by_qual = TRUE,
  homopoly_fix_min_ident = dada_min_identical,
  min_homopoly_len = 6,
  fixed_cluster_threshold = 0.97,
  taxa_cluster_threshold = fixed_cluster_threshold,
  cluster_single_linkage = TRUE,
  min_variant_freq = 0.2,
  split_min_identical = dada_min_identical,
  max_split_ratio = 3,
  cores = 1,
  verbose = FALSE
)
```

## Arguments

- fq:

  demultiplexed FASTQ file

- dada_err:

  output from
  [`dada2::learnErrors()`](https://rdrr.io/pkg/dada2/man/learnErrors.html)

- alignment_prefix:

  Optional output prefix for BAM alignment files and FASTA references
  for manual inspection and/or downstream analyses

- id_prefix:

  Name prefix for sequence IDs in the consensus FASTA and BAM alignment
  files

- tmp_dir:

  Optional path to a temporary directory (default is
  [`tempdir()`](https://rdrr.io/r/base/tempfile.html))

- dada_omega_a:

  OMEGA_A parameter value(s) to use for the denoising (see
  ?dada2::setDadaOpt): Can be a vector of increasing values (the
  default), which are tried sequentially until a grouping with
  well-separated haplotypes and non-ambiguous consensus sequence(s) is
  found for the *top taxon* (*note*: others currently ignored)

- omegaA_iter_threshold:

  Only evaluate multiple `dada_omega_a` values if the read depth is not
  larger than the given threshold. For large samples, it is better to
  avoid multiple denoising rounds, and it is anyway likely that DADA2
  will already find all sequence variants in the first round (with
  lowest OMEGA_A). Haplotype variation with InDels may still be found
  with 'try_split_haplotypes', but not with DADA2 denoising at high
  sensitivity.

- dada_min_identical:

  Minimum number of identical sequences required to do a DADA2
  denoising. Below this threshold, switch to simple fixed-threshold
  clustering instead (at `fixed_cluster_threshold`) and report the
  consensus without further attempting any haplotype splitting.

- dada_min_n0:

  Minimum number of error-free reads needed to retain the DADA2
  clustering (`n0` in dada-class \$clustering information). Switch to
  fixed-threshold clustering + consensus method below this threshold
  (see `fixed_cluster_threshold`). *Note*: *error-free* does not mean
  *identical*, as DADA2 only considers substitutions, and there can
  still be InDels, so `n0` is always equal or higher to the number of
  identical sequences.

- min_seq_abund:

  minimum number of supporting reads needed for any barcode sequence
  variant to be included in the results

- max_sample_depth:

  read a maximum of `max_sample_depth` demultiplexed reads from the
  input file (`fq`) for the denoising/clustering

- consensus_max_depth:

  maximum number sequences mapped against the clusters to infer the
  consensus sequence (if there are more, a random sample is taken)

- consensus_threshold:

  require at least the given proportion of bases to be identical at
  every alignment column for an unambiguous consensus call (values below
  the default 60% might be problematic)

- consensus_by_qual:

  Consider the read quality scores when building the consensus with
  [samtools
  consensus](https://www.htslib.org/doc/samtools-consensus.html). If
  `TRUE`, the relative base frequencies are weighted by the Phred
  quality scores.

- homopoly_fix_min_ident:

  Minimum number of identical sequences required to attempt adjusting
  ambiguous homopolymer sequences in the consensus. The homopolymer run
  length of the most frequent sequence is chosen.

- min_homopoly_len:

  Minimum length a homopolymer region needs to have in order to attempt
  "fixing" an ambiguous consensus in that region (see also
  `homopoly_fix_min_ident`)

- fixed_cluster_threshold:

  Similarity threshold for grouping sequence variants per taxon.
  Single-linkage clustering is applied, with the default threshold of
  0.97, the maximum divergence between any two sequences can be 3%.

- taxa_cluster_threshold:

  Similarity threshold and for post-clustering of DADA2 ASVs to group
  potential haplotypes of the same sequenced organism together.

- cluster_single_linkage:

  Whether to apply single-linkage fixed-threshold clustering for
  low-depth samples. Clusters grow as long as any two sequences have at
  least `fixed_cluster_threshold` similarity. If `FALSE`, the all
  cluster members are compared to one centroid sequence and remain more
  concise (=complete-linkage clustering). On (`TRUE`) by default, as
  this appears to work well for low-coverage samples.

- min_variant_freq:

  frequency threshold to consider separate sequences (ASVs) of the same
  taxon (as clustered with `taxa_cluster_threshold`) as "real"
  polymorphisms (haplotypes), not noise.

- split_min_identical:

  min. number of identical sequences that needs to support *both* of top
  two unique sequences in a sample in order to attempt haplotype
  splitting (with these sequences as new references).

- max_split_ratio:

  only split an ASV into 2 parental haplotypes if the ratio of
  larger:smaller is up to `max_split_ratio`. While the haplotypes might
  not always have an exact 1:1 ratio, this constraint still enforces a
  certain balancing of haplotype abundances. More unbalanced ratios are
  only possible if DADA2 is rerun with higher sensitivity (see
  `dada_omega_a`).

## Value

Returns a data frame with at least the following columns:

- *id*: Descriptive sequence ID, e.g.: 'taxon1_seq2'

- *full_id*: The ID used in the BAM and FASTA output files: *id*
  prefixed with `id_prefix`

- *taxon_num*: Nth (putative) taxon (integer; see also
  `taxa_cluster_threshold`)

- *sequence*: DADA2 ASV or most abundant unique split sequence (`NA` in
  case of fixed-threshold clustering)

- *consensus*: consensus sequence of haplotype member sequences

- *consensus_ambigs*: number of ambiguous bases in the consensus (see
  `consensus_threshold`)

- *consensus_diffs*: edit distance between consensus and the
  representative sequence (DADA2 ASV) or dominant sequence (after
  haplotype splitting/fixed-threshold clustering)

- *homopolymer_adjustments*: number of homopolymer locations that were
  "fixed" by inferring the repeat number from the most abundant sequence
  instead of relying on the consensus (containing Ns)

- *abundance*: number of reads clustering with the given sequence (up to
  `max_sample_depth`)

- *n_mapped*: number of reads mapped to the barcode sequences (usually
  the same or very similar to *abundance*)

- *n0* from DADA2: number of error-free reads (no substitutions, there
  may still be InDels)

- *max_identical*: max. number of identical reads

- *method*: method from which the haplotype emerged: one of 'dada',
  'dada_split' (if haplotype splitting was done), 'fixed_cluster'

- *is_rare*: `TRUE` if the relative frequency of a sequence within the
  same putative taxon is below `min_variant_freq`
