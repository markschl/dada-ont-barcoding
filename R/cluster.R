
#' Infer barcode sequence(s) from raw reads in a single sample
#'
#' Applies a clustering workflow to the provided raw reads, based on
#' DADA2 denoising and/or fixed-threshold clustering to infer the correct
#' barcode sequence(s). Further groups the resulting sequences by putative
#' taxon.
#'
#' @param fq demultiplexed FASTQ file
#' @param dada_err output from [dada_learn_errors]
#' @param alignment_prefix Optional output prefix for BAM alignment files and
#'    FASTA references for manual inspection and/or downstream analyses
#' @param id_prefix Name prefix for sequence IDs in the consensus FASTA and BAM alignment files
#' @param tmp_dir Optional path to a temporary directory (details in [set_global_opts])
#' @param dada_omega_a OMEGA_A parameter value(s) to use for the denoising
#'     (see [dada2::setDadaOpt]):
#'    Can be a vector of increasing values (the default), which are tried
#'    sequentially until a grouping with well-separated haplotypes and
#'    non-ambiguous consensus sequence(s) is found for the *top taxon*
#'    (*note*: others currently ignored)
#' @param omegaA_iter_threshold Only evaluate multiple `dada_omega_a` values if the
#'    read depth is not larger than the given threshold. For large samples, it is
#'    better to avoid multiple denoising rounds, and it is anyway likely that DADA2
#'    will already find all sequence variants in the first round (with lowest OMEGA_A).
#'    Haplotype variation with InDels may still be found with 'try_split_haplotypes',
#'    but not with DADA2 denoising at high sensitivity.
#' @param dada_min_identical Minimum number of identical sequences required to do a
#'    DADA2 denoising. Below this threshold, switch to simple fixed-threshold
#'    clustering instead (at `fixed_cluster_threshold`) and report the consensus
#'    without further attempting any haplotype splitting.
#' @param dada_min_n0 Minimum number of error-free reads needed to retain the DADA2 clustering
#'   (`n0` in dada-class $clustering information).
#'   Switch to fixed-threshold clustering + consensus method below this threshold
#'   (see `fixed_cluster_threshold`).
#'   *Note*: *error-free* does not mean *identical*, as DADA2 only considers
#'   substitutions, and there can still be InDels, so `n0` is always equal or higher
#'   to the number of identical sequences.
#' @param min_seq_abund minimum number of supporting reads needed for any barcode
#'   sequence variant to be included in the results
#' @param min_variant_freq frequency threshold to consider separate sequences (ASVs)
#'    of the same taxon (as clustered with `taxa_cluster_threshold`) as "real" polymorphisms
#'    (haplotypes), not noise.
#' @param max_sample_depth read a maximum of `max_sample_depth` demultiplexed reads
#'    from the input file (`fq`) for the denoising/clustering
#' @param consensus_max_depth maximum number sequences mapped against the
#'    clusters to infer the consensus sequence (if there are more, a random sample
#'    is taken)
#' @param consensus_threshold require at least the given proportion of bases
#'    to be identical at every alignment column for an unambiguous consensus call
#'    (values below the default 60% might be problematic)
#' @param consensus_by_qual Consider the read quality scores when building the
#'    consensus with [samtools consensus](https://www.htslib.org/doc/samtools-consensus.html).
#'    If `TRUE`, the relative base frequencies are weighted by the Phred quality scores.
#' @param homopoly_fix_min_ident Minimum number of identical sequences required to
#'    attempt adjusting ambiguous homopolymer sequences in the consensus.
#'    The homopolymer run length of the most frequent sequence is chosen.
#' @param min_homopoly_len Minimum length a homopolymer region needs to have in order
#'    to attempt "fixing" an ambiguous consensus in that region
#'    (see also `homopoly_fix_min_ident`)
#' @param fixed_cluster_threshold Similarity threshold for grouping sequence variants
#'    per taxon. Single-linkage clustering is applied, with the default threshold of
#'    0.97, the maximum divergence between any two sequences can be 3%.
#' @param cluster_single_linkage Whether to apply single-linkage fixed-threshold
#'    clustering for low-depth samples.
#'    Clusters grow as long as any two sequences have at least `fixed_cluster_threshold`
#'    similarity. If `FALSE`, the all cluster members are compared to one centroid
#'    sequence and remain more concise (=complete-linkage clustering).
#'    On (`TRUE`) by default, as this appears to work well for low-coverage samples.
#' @param taxa_cluster_threshold Similarity threshold
#'    and for post-clustering of DADA2 ASVs to group potential haplotypes of the
#'    same sequenced organism together.
#' @param split_min_identical min. number of identical sequences that needs to
#'    support *both* of top two unique sequences in a sample in order to attempt
#'    haplotype splitting (with these sequences as new references).
#' @param max_split_ratio only split an ASV into 2 parental haplotypes if the ratio
#'    of larger:smaller is up to `max_split_ratio`.
#'    While the haplotypes might not always have an exact 1:1 ratio, this
#'    constraint still enforces a certain balancing of haplotype abundances.
#'    More unbalanced ratios are only possible if DADA2 is rerun with higher sensitivity
#'    (see `dada_omega_a`).
#'
#' @returns Returns a data frame with at least the following columns:
#'    - *id*: Descriptive sequence ID, e.g.: 'taxon1_seq2'
#'    - *full_id*: The ID used in the BAM and FASTA output files: *id* prefixed with `id_prefix`
#'    - *taxon_num*: Nth (putative) taxon (integer; see also `taxa_cluster_threshold`)
#'    - *sequence*: DADA2 ASV or most abundant unique split sequence
#'       (`NA` in case of fixed-threshold clustering)
#'    - *consensus*: consensus sequence of haplotype member sequences
#'    - *consensus_ambigs*: number of ambiguous bases in the consensus (see `consensus_threshold`)
#'    - *consensus_diffs*: edit distance between consensus and the representative sequence
#'      (DADA2 ASV) or dominant sequence (after haplotype splitting/fixed-threshold clustering)
#'    - *homopolymer_adjustments*: number of homopolymer locations that were "fixed"
#'       by inferring the repeat number from the most abundant sequence instead of relying
#'       on the consensus (containing Ns)
#'    - *abundance*: number of reads clustering with the given sequence
#'       (up to `max_sample_depth`)
#'    - *n_mapped*: number of reads mapped to the barcode sequences
#'       (usually the same or very similar to *abundance*)
#'    - *n0* from DADA2: number of error-free reads (no substitutions, there may still be InDels)
#'    - *max_identical*: max. number of identical reads
#'    - *method*: method from which the haplotype emerged:
#'       one of 'dada', 'dada_split' (if haplotype splitting was done), 'fixed_cluster'
#'    - *is_rare*: `TRUE` if the relative frequency of a sequence within the same
#'      putative taxon is below `min_variant_freq`
#'
#' @export
infer_barcode <- function(fq,
                          dada_err,
                          alignment_prefix,
                          id_prefix = NULL,
                          tmp_dir = NULL,
                          dada_omega_a = c(1e-20, 1e-10, 1e-2),
                          omegaA_iter_threshold = 1000,
                          dada_min_identical = 2,
                          dada_min_n0 = 4,
                          min_seq_abund = 3,
                          max_sample_depth = 5000,
                          consensus_max_depth = 3000,
                          consensus_threshold = 0.65,
                          consensus_by_qual = TRUE,
                          homopoly_fix_min_ident = 4,
                          min_homopoly_len = 6,
                          fixed_cluster_threshold = 0.97,
                          taxa_cluster_threshold = fixed_cluster_threshold,
                          cluster_single_linkage = TRUE,
                          min_variant_freq = 0.2,
                          split_min_identical = 4,
                          max_split_ratio = 3,
                          cores = 1,
                          verbose = FALSE) {
  # paths/global options
  samtools <- get_program('samtools')
  minimap2 <- get_program('minimap2')
  tmp_dir <- tmp_dir %||% get_tmp_dir()

  # temporary directory
  tmp <- tempfile('infer_barcode_', tmpdir = tmp_dir)
  dir.create(tmp, FALSE, TRUE)

  # We will need the reads with quality scores below
  # cat('read/derep '); tictoc::tic()
  reads <- Biostrings::readQualityScaledDNAStringSet(fq, nrec = max_sample_depth)
  if (length(reads) == max_sample_depth) {
    # DADA2 derepFastq does not offer a way to only read the top N reads,
    # so we need to prepare a new FASTQ file for this
    fq <- file.path(tmp, '_fq_limited.fastq')
    Biostrings::writeQualityScaledXStringSet(reads, fq)
  }

  # DADA2: dereplicate once for multiple rounds of denoising
  dada_derep <- dada2::derepFastq(fq, qualityType = 'FastqQuality')
  # tictoc::toc()

  d.prev = NULL
  stopifnot(length(dada_omega_a) >= 1)
  for (dada_attempt in seq_along(dada_omega_a)) {
    round_prefix <- file.path(tmp, paste0('round_', dada_attempt))
    # Run DADA2 denoising
    # note: denoising is done in every case even if the dada_min_identical threshold
    # is not met; DADA2 is anyway very fast with low-coverage samples
    # cat('dada '); tictoc::tic()
    d <- dada2_denoise(
      dada_derep,
      dada_err,
      max_members = consensus_max_depth,
      OMEGA_A = dada_omega_a[dada_attempt]
    )
    d$method = 'dada'
    # tictoc::toc()
    dada_detail <- attributes(d)[c('n_reads', 'n_singletons')]

    # Output should never be empty with DADA2 if DETECT_SINGLETONS is on
    stopifnot(nrow(d) > 0)

    # low-coverage (and/or low-quality) sample?
    top_is_replicated <- d$max_identical[1] >= dada_min_identical &&
      d$n0[1] >= dada_min_n0
    # Switch to fixed-threshold clustering if there are very few identical
    # reads or or few reads without substitution errors (n0)
    # https://github.com/benjjneb/dada2/issues/943#issuecomment-585219967
    if (!top_is_replicated) {
      # cat('cluster fixed '); tictoc::tic()
      cl <- cluster_fixed(reads,
                          fixed_cluster_threshold,
                          single_linkage = cluster_single_linkage,
                          abund_order = TRUE)

      # For each cluster, get the most abundant unique read (and its count),
      # used as reference for read mapping + consensus.
      # If only singletons, we have to use the first-occurring sequence.
      reads_split <- split(as.character(reads), cl)
      top_uniq <- lapply(reads_split, function(r) sort(table(r), TRUE)[1])
      ref_seq <- sapply(top_uniq, names)
      stopifnot(names(reads_split) == names(ref_seq))
      d <- data.frame(
        row.names = names(ref_seq),
        id = as.integer(names(ref_seq)),
        sequence = ref_seq,
        abundance = sapply(reads_split, length),
        max_identical = sapply(top_uniq, max),
        method = 'cluster_fixed',
        check.names = FALSE
      )
      d$seq_indices = split(seq_along(reads), cl)[as.character(d$id)]
      # tictoc::toc()
    }

    # colum containing rarely-occurring messages
    d$message = ''

    # depth filter
    # TODO: initial cluster abundance may differ from mapped read abundance,
    #   but we use this for filtering
    d <- d[d$abundance >= min_seq_abund, ]

    if (nrow(d) > 0) {
      # Obtain consensus by mapping all ASV members against all ASV sequences
      # *Note*: reads not corrected by DADA2 are not mapped (their number depends on
      #  OMEGA_C, see ?setDadaOpt).
      # In some cases, members are mapped to another ASV than the one they belong to,
      # but this should rarely happen
      # cat('align/consensus all '); tictoc::tic()
      sel_reads <- reads[unlist(d$seq_indices)]  # mapped to any ASV
      ref_seq <- with(d, setNames(sequence, id))
      cons <- ambig_consensus(
        sel_reads,
        ref_seq,
        out_prefix = round_prefix,
        consensus_threshold = consensus_threshold,
        consensus_by_qual = consensus_by_qual,
        homopoly_fix = d$max_identical >= dada_min_identical,
        min_homopoly_len = min_homopoly_len,
        fast = TRUE,
        cores = cores,
        minimap2 = minimap2,
        samtools = samtools
      )
      stopifnot(row.names(cons) == d$id)
      d <- cbind(d, cons)
      # tictoc::toc()
      # Read numbers may (rarely) not be exactly the same as reads can map to the
      # "wrong" ASV (mapping not restricted) -> we use mapped reads
      # TODO: is this a good choice?
      stopifnot(names(d$consensus) == d$id)

      if (top_is_replicated && any(d$consensus_ambigs > 0)) {
        # In case of ambiguities, try further splitting
        # into 2 haplotypes
        # cat('attempt split '); tictoc::tic()
        d <- try_split_haplotypes(
          reads,
          d,
          prefix = round_prefix,
          min_identical = split_min_identical,
          max_ratio = max_split_ratio,
          consensus_threshold = consensus_threshold,
          consensus_by_qual = consensus_by_qual,
          min_homopoly_len = min_homopoly_len,
          fast = TRUE,
          cores = cores,
          minimap2 = minimap2,
          samtools = samtools,
          verbose = verbose
        )
        # tictoc::toc()
      }

      # post-cluster to group by (putative) taxon
      d$taxon_num = cluster_fixed(d$consensus,
                                  taxa_cluster_threshold,
                                  single_linkage = TRUE)

      # flag rare sequence variation
      d$is_rare <- as.logical(ave(d$n_mapped, d$taxon_num,
                                  FUN=function(x) x / sum(x) < min_variant_freq))

      # reorder taxon clusters by total abundance (keep members together)
      taxon_abund <- ave(ifelse(d$is_rare, 0, d$n_mapped), d$taxon_num, FUN=sum)
      d <- d[order(-taxon_abund, d$taxon_num, -d$n_mapped), ]
      d$taxon_num <- match(d$taxon_num, unique(d$taxon_num))

      # determine number of haplotypes
      is_cl1 <- d$taxon_num == d$taxon_num[1] & !d$is_rare
      n_seqs <- sum(is_cl1)
      if (verbose) {
        cat(
          sprintf(
            '%s\n    seq: %d | clust: %s | map: %s | n0: %s; ident: %s | diff: %s | ambig: %s | homop-adj: %s\n',
            if (dada_attempt == 1) {
              sprintf('%s (N = %d) | %s', basename(fq), dada_detail$n_reads, d$method[1])
            } else {
              paste("... omegaA =", dada_omega_a[dada_attempt])
            },
            n_seqs,
            paste(d$abundance[1:n_seqs], collapse = '/'),
            paste(d$n_mapped[1:n_seqs], collapse = '/'),
            paste(d$n0[1:n_seqs], collapse = '/'),
            paste(d$max_identical[1:n_seqs], collapse = '/'),
            paste(d$consensus_diffs[1:n_seqs], collapse = '/'),
            paste(d$consensus_ambigs[1:n_seqs], collapse = '/'),
            paste(d$homopolymer_adjustments[1:n_seqs], collapse = '/')
          ),
          file = stderr()
        )
      }
      stopifnot(length(unique(d$taxon_num[1:n_seqs])) == 1)
    }

    # go back to previous clustering if:
    # - all sequences have too little reads (below the abundance threshold)
    # - the total reads included in all clusters are >2x lower than last time
    #   (happens for very low-depth samples)
    # TODO: more conditions?
    if (!is.null(d.prev) &&
        (nrow(d) == 0 ||
         sum(d.prev$n_mapped) / sum(d$n_mapped) > 2)) {
      if (verbose)
        cat("...back\n", file = stderr())
      d <- d.prev
      dada_attempt <- dada_attempt - 1
      break
    }
    # stop with current clustering if there are no ambiguities in consensus of *top* taxon
    # OR the sample depth is large
    # (in which case it is unlikely that the haplotypes are not found with the initial omegaA)
    if (!top_is_replicated ||
        nrow(d) == 0 ||
        sum(d$consensus_ambigs[1:n_seqs]) == 0 ||
        dada_detail$n_reads >= omegaA_iter_threshold) {
      break
    }
    d.prev = d

  } # end of loop

  if (nrow(d) == 0) {
    unlink(tmp, TRUE)
    return(invisible(NULL))
  }

  # Delete entries of 'sequence' that are not well supported by duplicates
  # (if fixed clustering was applied, all 'sequence' entries will be NA,
  # otherwise only rare ones)
  # *note* sequence is not directly used later, leaving it may cause confusion
  d$sequence[d$max_identical < dada_min_identical] = NA

  # Assemble final BAM files:
  # - select clusters/haplotypes >= min_seq_abund (!d$is_rare)
  # - remap reads to consensus if it differs from the reference sequence
  #   (as it is the consensus that we will report)
  round_prefix <- file.path(tmp, paste0('round_', dada_attempt))
  d.sel <- d[!d$is_rare,]
  # We have to re-map to the consensus if there are ambiguous bases or
  # any differences to the dominant ASV/unique split haplotype
  remap_cons <- d.sel$consensus_diffs > 0 | d.sel$consensus_ambigs > 0
  remap_prefix <- file.path(tmp, 'remap_cons')
  # cat('assemble output '); tictoc::tic()
  if (all(!remap_cons)) {
    # nothing to remap -> just rename files or samtools view
    if (!any(d$is_rare)) {
      move_bam(round_prefix, remap_prefix)
    } else {
      subset_combine_bam(
        out_prefix = remap_prefix,
        list(list(round_prefix, d.sel$id)),
        do_index = FALSE,
        cores = cores,
        samtools = samtools
      )
    }
  } else {
    # some IDs require remapping
    sel_reads <- reads[unlist(d.sel$seq_indices[remap_cons])]
    # Max. 10 rounds of re-mapping to consensus:
    # Stop if 'samtools consensus' returns a consistent in two successive mappings.
    # Usually only 1 round necessary, but InDel-rich reads might require
    # multiple rounds.
    for (i in 1:10) {
      ref_seq <- setNames(d.sel$consensus[remap_cons], d.sel$id[remap_cons])
      cons <- ambig_consensus(
        sel_reads,
        ref_seq,
        out_prefix = remap_prefix,
        consensus_threshold = consensus_threshold,
        consensus_by_qual = consensus_by_qual,
        homopoly_fix = d.sel$max_identical[remap_cons] >= dada_min_identical,
        min_homopoly_len = min_homopoly_len,
        fast = !all(remap_cons),
        cores = cores,
        minimap2 = minimap2,
        samtools = samtools
      )
      stopifnot(row.names(cons) == d.sel$id[remap_cons])
      cons_same <- d.sel$consensus[remap_cons] == cons$consensus
      d.sel[remap_cons, names(cons)] <- cons
      if (all(cons_same)) {
        break
      }
      if (verbose)
        cat('Repeating consensus re-mapping\n')
    }
    if (any(!cons_same)) {
      # issue warning if no consistent consensus found after 10 rounds
      d.sel$message[remap_cons][!cons_same] <- paste0(
        d.sel$message[remap_cons][!cons_same],
        'Re-mapping to consensus gives another consensus! '
      )
    }
    if (!all(remap_cons)) {
      remap_prefix0 <- remap_prefix
      remap_prefix <- paste0(remap_prefix, '.2')
      subset_combine_bam(
        out_prefix = remap_prefix,
        list(
          list(remap_prefix0, d.sel$id[remap_cons]),
          list(round_prefix, d.sel$id[!remap_cons])
        ),
        do_index = FALSE,
        cores = cores,
        samtools = samtools
      )
    }
  }
  # tictoc::toc()

  # Finally: assign intuitively understandable names to sequences
  orig_id <- d$id
  d$id <- paste0('taxon', d$taxon_num)
  sel <- ave(d$id, d$taxon_num, FUN=length) > 1
  d$id[sel] <- sprintf(
    '%s_seq%s',
    d$id[sel],
    ave(d$id[sel], d$taxon_num[sel], FUN=seq_along)
  )
  id_prefix <- if (is.null(id_prefix)) {
    ''
  } else {
    stopifnot(!grepl(' ', id_prefix))
    paste0(id_prefix, '_')
  }
  d$full_id <- paste0(id_prefix, d$id)

  # Create BAM/FASTA output files (with IDs renamed in headers)
  if (!is.null(alignment_prefix)) {
    commented_id <- sprintf(
      '%s reads: %d | ambigs: %d | identical: %d | subst-free (n0): %d | %s',
      d$full_id,
      d$n_mapped,
      d$consensus_ambigs,
      d$max_identical,
      d$n0 %||% NA,
      d$message
    )
    rename_bam_refs(remap_prefix, alignment_prefix,
                    id_map = setNames(commented_id, orig_id),
                    ref_seq = setNames(d$consensus, orig_id),
                    do_index = TRUE,
                    samtools = samtools)
  }

  # clean up
  d$seq_indices = d$top_uniques = NULL
  unlink(tmp, TRUE)

  # calculate homopolymer stretch length
  l <- rle(utf8ToInt(d$consensus[1]))$lengths
  d$max_homopoly_len = max(l)

  attributes(d)[names(dada_detail)] = dada_detail
  attr(d, 'omega_a') = dada_omega_a[dada_attempt]
  d
}


#' Compare consensus/ASV and known sequences
#'
#' Maps already known sequences or inconsistent ASV sequences against the consensus
#' and creates very small BAM alignment files containing just inconsistent
#' consensus <-> cluster sequence comparisons and consensus <-> known sequence comparisons
#' (also useful for manual inspection)
#'
#' @param data frame returned by [infer_barcode] (needs the *consensus* column).
#'
#' @returns
#' Adds a `consensus_diffs` column to `d` (NA if not compared, Inf if not mapped
#' due to too many mismatches)
#' TODO: ambiguous bases unfortunately lead to mismatches
#'
#' @export
compare_seqs <- function(d,
                         bam_out = NULL,
                         tmp_dir = NULL,
                         known_seq = NULL) {
  samtools <- get_program('samtools')
  minimap2 <- get_program('minimap2')
  tmp_dir <- tmp_dir %||% get_tmp_dir()
  known_seq <- known_seq %||% NA
  stopifnot(length(known_seq) == 1)
  d$known_seq_diffs <- NA_integer_
  sel <- !is.na(d$sequence) & !d$is_rare & d$sequence != d$consensus
  map_seqs <- c()
  if (any(sel)) {
    map_seqs <- c(
      map_seqs,
      setNames(d$sequence[sel], paste0(d$full_id[sel], '_dominant_seq'))
    )
  }
  if (!is.na(known_seq)) {
    sel <- rep(TRUE, nrow(d))
    map_seqs <- c(map_seqs, setNames(known_seq, 'known_sequence'))
  }
  do_cmp <- attr(d, 'has_seq_comparison') <- length(map_seqs) > 0
  if (do_cmp) {
    stopifnot(any(sel))
    tmp_dir <- tmp_dir %||% tempdir()
    seq_file <- tempfile('cmp_seq', tmp_dir, fileext = '.fasta')
    ref_file <- tempfile('cmp_ref', tmp_dir, fileext = '.fasta')
    write_dna(map_seqs, seq_file)
    write_dna(setNames(d$consensus[sel], d$full_id[sel]), ref_file)
    if (no_bam <- is.null(bam_out)) {
      bam_out <- tempfile('cmp', tmp_dir, fileext = '.bam')
    }
    args <- c(
      seq_file,
      ref_file,
      bam_out,
      1,
      minimap2,
      samtools
    )
    sam <- run_bash_script('map_ref_simple.sh', args, stdout = TRUE)
    if (no_bam)
      file.remove(bam_out)
    sam <- if (length(sam) > 0) {
      sam <- read.delim(textConnection(sam), header = FALSE, colClasses = 'character')
      sam[sam[[1]] == 'known_sequence',]
    } else {
      data.frame()
    }
    if (nrow(sam) == 0) {
      if (!is.na(known_seq)) {
        # nothing mapped -> must be large number
        d$known_seq_diffs[1] <- Inf
      }
    } else {
      stopifnot(nrow(sam) == 1)
      dist_col <- which(grepl('NM:i:', unlist(sam), fixed=TRUE))[1]
      stopifnot(!is.na(dist_col))
      mapped_i <- gsub(' .*', '', d$full_id) == sam[[3]]
      stopifnot(sum(mapped_i) == 1)
      d$known_seq_diffs[mapped_i] <- as.integer(gsub('NM:i:', '', sam[[dist_col]]))
    }
    invisible(file.remove(seq_file, ref_file))
  }
  d
}

#' Learn error rates
#'
#' Calls [dada2::learnErrors] with settings adjusted for Nanopore data and the clustering
#' procedure ([infer_barcode])
#'
#' @export
dada_learn_errors <- function(fq_paths, omega_a = 1e-20, cores = 1, ...) {
  # To be sure we increase the NW alignment band size,
  # which is recommended for technologies with potentially some InDels (454, PacBio)
  # and also applies for Nanopore.
  dada2::learnErrors(fq_paths,
                     BAND_SIZE = 32,
                     qualityType = 'FastqQuality',
                     OMEGA_A = omega_a,
                     # this setting is currently hard-coded in denoise_dada2,
                     # so we also need it here
                     OMEGA_C = 1e-10,
                     multithread = cores,
                     ...,
                     verbose = 0)

}

#' Run DADA2 denoising
#'
#' Runs DADA2 on a FASTQ file, removes chimeras and and returns
#' an abundance-sorted data frame of ASVs with some additional information
#'
#' @param x FASTQ file path or *derep-class* object
#' @param dada_err error information returned by [dada_learn_errors]
#' @param max_members do random subsampling of overly large ASVs to
#' the specified number of sequences
#' @param singleton_threshold by default, singletons cannot form a new cluster,
#' but for very low-depth samples (< N duplicates of any sequence),
#' DETECT_SINGLETONS is turned on (see [dada2::setDadaOpt]).
#' This increases sensitivity with InDel-rich Nanopore data and prevents that
#' a large proportion of sequences are discarded.
#' @param ... Arguments passed to [dada2::dada]
#'
#' @returns a data frame with the following columns:
#' - *sequence*
#' - *top_uniques*: abundance of top unique sequences
dada2_denoise <- function(x,
                          dada_err,
                          cores = 1,
                          max_members = 1e6,
                          singleton_threshold = 5,
                          ...) {
  # DADA2 denoise
  derep <- if (inherits(x, 'character')) {
    dada2::derepFastq(x, qualityType = 'FastqQuality')
  } else {
    stopifnot(inherits(x, 'derep'))
    x
  }
  dd <- dada2::dada(
    derep,
    err = dada_err,
    BAND_SIZE = 32,
    DETECT_SINGLETONS = max(derep$uniques) < singleton_threshold,
    # controls how many seqs. end up in ASV
    # (higher than default 1e-40, meaning that more seqs. end up as ASV members;
    # appears to be a good value for low-coverage samples)
    # TODO: should it be configurable?
    OMEGA_C = 1e-10,
    multithread = cores,
    verbose = 0,
    ...
  )

  d <- dd$clustering[, c('sequence', 'abundance', 'n0', 'nunq'), drop = FALSE]
  # assign an ID, as we will reorder/subset later, but need the cluster number
  d$id = 1:nrow(d)

  # # Create Phred quality string with base 33 for ASV sequence
  # # TODO: should be fine like this?
  # # dada-class$quality has "average quality scores for each cluster (row) by position (col)"
  # d$quality = sapply(1:nrow(dd$clustering), function(i) {
  #   intToUtf8(33 + as.integer(round(dd$quality[i, 1:nchar(dd$clustering$sequence[i])])))
  # })

  # Chimera check (we do this independently for each barcode)
  if (nrow(d) > 2) {
    # create a mini count table for this sample (normally very few ASVs)
    st <- with(d, matrix(
      abundance,
      nrow = 1,
      dimnames = list(NULL, sequence)
    ))
    bim <- dada2::isBimeraDenovo(st, multithread = cores)
    d <- d[!bim, ]
  }

  # Sort by decreasing abundance
  d <- d[order(-d$abundance), ]

  # Remember the read indices for each cluster in the FASTQ file, used later in ambig_consensus()
  stopifnot(length(dd$map) == length(derep$uniques))
  # ASV cluster ID for *all* (non-dereplicated) input sequences (same order)
  cl <- dd$map[derep$map]
  cl <- split(seq_along(cl), cl)
  # random subsampling for overly large ASVs
  cl <- lapply(cl, function(i) {
    if (length(i) > max_members) {
      i <- sample(i, max_members)
    }
    i
  })
  d$seq_indices = cl[as.character(d$id)]

  # In addition, remember the two most abundant sequences, used later for haplotype splitting
  d$top_uniques = lapply(d$seq_indices, function(i) {
    uniq_freq <- tabulate(derep$map[i])
    i <- head(order(-uniq_freq), 2)
    setNames(uniq_freq[i], names(derep$uniques)[i])
  })

  # max. number of identical sequences (usually identical to inferred ASV sequence)
  d$max_identical = sapply(d$top_uniques, max)

  # attach more information
  attr(d, 'n_reads') = length(derep$map)
  attr(d, 'n_singletons') = sum(tabulate(derep$map) == 1)
  d
}


#' Attempts splitting sequence cluster into two haplotypes.
#'
#' @param reads XStringQuality object
#' @param d a data frame as returned by [dada2_denoise] and [ambig_consensus] with these columns:
#' - id
#' - top_uniques
#' - consensus_ambigs
#' @param prefix Path prefix where the .bam alignment and .fasta reference is found
#' (may be overwritten!)
#' @param min_identical minimum number of *identical* sequences supporting *both*
#' of the dominant unique sequences
#' @param max_ratio maximum abundance ratio (mapped read abundances) for splitting
#' @param fast use faster compression (larger BAM file)
#'
#' @details
#' Splitting is done if the sum of the consensus ambiguities of the two resulting
#' sequence variants is smaller than the number of consensus ambiguities of the parent
#' cluster.
#'
#' Overwrites the existing BAM files and FASTA references in 'prefix'.
try_split_haplotypes <- function(reads,
                                 d,
                                 prefix,
                                 min_identical = 4,
                                 max_ratio = 3,
                                 fast = FALSE,
                                 cores = 1,
                                 ...,
                                 samtools = 'samtools',
                                 minimap2 = 'minimap2',
                                 verbose = FALSE) {



  d$do_split = FALSE

  # for which ASVs should we *attempt* splitting?
  check_split <- d$consensus_ambigs > 0 &
    # at least min_identical uniques are required for both resulting haplotypes
    sapply(d$top_uniques, min) >= min_identical &
    # simple approximate pre-filtering:
    # enforce 2.5*max_ratio for uniques, based on the fact that unique abundances
    # correlate with final mapped read abundances
    # (might not be precise at low numbers, but there anyway we want to be cautious
    # not to over-split)
    sapply(d$top_uniques, function(x)
      max(x) / min(x)) <= 2.5 * max_ratio
  stopifnot(!is.na(check_split))

  if (any(check_split)) {
    # attempt splitting for those ambiguous ASVs that pass the thresholds
    # for each ASV, select the two top unique sequences as references for mapping + consensus
    prefix_split <- paste0(prefix, '_split')
    sel_reads <- reads[unlist(d$seq_indices[check_split])]
    ref_seq <- unlist(setNames(lapply(d$top_uniques[check_split], function(x) {
      setNames(names(x), seq_along(x))
    }), d$id[check_split]))
    # do mapping + consensus for all candidates at once
    # (rarely reads might switch to a different reference)
    cons <- ambig_consensus(
      sel_reads,
      ref_seq ,
      out_prefix = prefix_split,
      homopoly_fix = TRUE,
      fast = fast,
      cores = cores,
      samtools = samtools,
      minimap2 = minimap2,
      ...
    )
    # for each ambiguous cluster, check if splitting leads to a better result
    orig_id <- gsub('\\.[^.]+$', '', row.names(cons))
    grouped_cons <- split(cons, orig_id)
    stopifnot(sort(row.names(grouped_cons)) == sort(as.character(d$id[check_split])))
    split_d <- do.call(rbind, lapply(which(check_split), function(i) {
      d.s <- d[i,]
      cons <- grouped_cons[[as.character(d.s$id)]]
      # required for splitting:
      # 1) sum of split ambiguities < non-split ambiguities
      # 2) the abundance ratio is not too high
      # TODO: enforcing an 1:1 ratio with a binomial test seems too strict
      # ratio.p = binom.test(cons$n_mapped[1], sum(cons$n_mapped), p = 0.5)$p.value
      abund.ratio = max(cons$n_mapped) / min(cons$n_mapped)
      is_better <- !any(is.na(cons$consensus)) &&
        sum(cons$consensus_ambigs) < d.s$consensus_ambigs &&
        abund.ratio <= max_ratio
      if (verbose) {
        cat(
          sprintf(
            ' split %s %s: has cons: %s | ambigs %d -> %d/%d = %d | uniq %d:%d | hap %d:%d | ratio=%.2f\n',
            d.s$id,
            is_better,
            !any(is.na(cons$consensus)),
            d.s$consensus_ambigs,
            cons$consensus_ambigs[1],
            cons$consensus_ambigs[2],
            sum(cons$consensus_ambigs),
            d.s$top_uniques[[1]][1],
            d.s$top_uniques[[1]][2],
            cons$n_mapped[1],
            cons$n_mapped[2],
            abund.ratio
          ),
          file = stderr()
        )
      }
      if (is_better) {
        uniq <- d.s$top_uniques[[1]]
        # one of the most abundant sequences *should* be identical to the ASV sequence
        # (but unfortunately not always is, probably due to little duplication
        # of low-abundance taxa)
        # stopifnot(is.null(d.s$sequence) | any(names(uniq) == d.s$sequence))
        stopifnot(paste(d.s$id, 1:2, sep='.') == row.names(cons))
        d2 <- cbind(
          id = row.names(cons),
          sequence = names(uniq),
          cons,
          max_identical = unname(uniq),
          method = paste0(d.s$method, '_split'),
          message = d.s$message
        )
        d2$top_uniques <- list(NULL, NULL)
        # obtaining seq_indices is a bit more involved, but usually
        # the number of reads is not too large here and this should be fast
        # (DADA2 will separate haplotypes well given enough read depth,
        # so try_split_haplotypes will not be required)
        read_map <- bam_to_map(paste0(prefix_split, '.bam'), samtools = samtools)
        read_ids <- sapply(strsplit(names(reads), ' ', fixed = TRUE), '[', 1)
        d2$seq_indices = split(match(names(read_map), read_ids), read_map)[row.names(cons)]
        stopifnot(lengths(d2$seq_indices) == d2$n_mapped)
        # all other columns should be undefined
        d2$do_split = TRUE
        d2[, setdiff(names(d.s), names(d2))] <- NA
        d2
      } else {
        d.s
      }
    }))

    stopifnot(nrow(d) == length(check_split))
    d <- rbind(
      d[!check_split,],
      split_d
    )

    # now create new BAM file with correct reads
    if (any(d$do_split)) {
      if (all(!d$do_split)) {
        # all ASVs were split, simply move the split file
        move_bam(prefix_split, prefix)
      } else {
        # extract data for split and non-split references (using samtools view), then merge them
        prefix_merged <- paste0(prefix, '.merged')
        merge_list <- list(list(prefix, d$id[!d$do_split]),
                           list(prefix_split, d$id[d$do_split]))
        subset_combine_bam(
          out_prefix = prefix_merged,
          merge_list,
          fast = fast,
          cores = cores,
          samtools = samtools
        )
        # overwrite originals
        move_bam(prefix_merged, prefix)
        # clean up
        remove_bam(prefix_split)
      }
    }
  }

  d$do_split = NULL
  d
}

# Merges multiple BAM files, optionally selecting reference IDs from them
subset_combine_bam <- function(out_prefix,
                               sel_list,
                               cores = 1,
                               chunk_size = 400,
                               fast = FALSE,
                               write_refs = TRUE,
                               do_index = TRUE,
                               samtools = 'samtools') {
  # recursively merge large lists if necessary
  if (recursive <- length(sel_list) > chunk_size) {
    n <- max(ceiling(length(sel_list) / chunk_size), chunk_size)
    sel_chunks <- split(sel_list, ceiling(1:length(sel_list) / n))
    sel_list <- paste0(out_prefix, '_tmp_', seq_along(sel_chunks))
    for (i in seq_along(sel_chunks)) {
      subset_combine_bam(
        sel_list[i],
        sel_chunks[[i]],
        cores = cores,
        chunk_size = chunk_size,
        fast = fast,
        write_refs = write_refs,
        do_index = FALSE,
        samtools = samtools
      )
    }
  }

  # merge BAM
  bam_out <- paste0(out_prefix, '.bam')
  outfiles <- c(bam = bam_out)
  msg <- run_bash(
    c(
      samtools,
      'merge',
      '--no-PG',
      '-f', if (fast) '-1' else NULL,
      '-@', cores,
      '-o', bam_out,
      sapply(sel_list, function(x) {
        if (length(x) == 2) {
          sprintf(
            '<(samtools view --no-PG -u %s %s)',
            paste0(x[[1]], '.bam'),
            paste(x[[2]], collapse = ' ')
          )
        } else {
          stopifnot(length(x) == 1)
          paste0(x[[1]], '.bam')
        }
      })
    ),
    stderr = TRUE,
    stdout = TRUE
  )
  stopifnot(!grepl('invalid region or unknown reference', msg))
  if (length(msg) > 0 && any(grepl('coordinate sort to be lost', msg, fixed = TRUE))) {
    # resort if multiple were merged (can sometimes happen)
    # https://github.com/samtools/samtools/issues/2159
    # TODO: not entirely clear why
    merged_f <- paste0(out_prefix, '.merged_tmp.bam')
    file.rename(bam_out, merged_f)
    run_bash(
      c(
        samtools,
        'sort',
        '--no-PG',
        '-l', if (fast) 1 else 6,
        '-@', cores,
        '-o', bam_out,
        merged_f
      ),
      stderr = FALSE
    )
    file.remove(merged_f)
  }
  if (do_index) {
    run_bash(c(samtools, 'index', bam_out))
    outfiles['bam.bai'] <- paste0(out_prefix, '.bam.bai')
  }
  # merge references
  if (write_refs) {
    ref_out <- paste0(out_prefix, '.fasta')
    refs <- do.call(c, lapply(sel_list, function(x) {
      r <- read_dna(paste0(x[[1]], '.fasta'))
      if (length(x) == 2) {
        sel <- match(x[[2]], gsub(' .*', '', names(r)))
        stopifnot(!is.na(sel))
        r <- r[sel]
      }
      r
    }))
    write_dna(refs, ref_out)
    outfiles['fasta'] <- paste0(out_prefix, '.fasta')
  }

  # clean up temporary files
  if (recursive) {
    for (prefix in sel_list) {
      remove_bam(prefix, index = FALSE)
    }
  }

  outfiles
}

# Run 'samtools reheader' to change the reference names.
# Also allows removing \@PG entries, which may contain some file paths
# from the machine on which the clustering was run.
rename_bam_refs <- function(prefix,
                            out_prefix,
                            id_map,
                            ref_seq = NULL,
                            do_index = TRUE,
                            remove_pg = TRUE,
                            samtools='samtools') {
  stopifnot(!is.null(names(id_map)))
  stopifnot(!is.na(names(id_map)))
  id_trans <- setNames(paste0('@SQ\tSN:', gsub(' .*', '', id_map), '\t'),
                       paste0('@SQ\tSN:', names(id_map), '\t'))
  bam <- paste0(prefix, '.bam')
  bam_out <- paste0(out_prefix, '.bam')
  header <- run_bash(c(samtools, 'view', '-H', '--no-PG', bam), stdout = TRUE)
  if (remove_pg) {
    header <- header[!grepl('@PG\t', header, fixed=TRUE)]
  }
  header <- paste(header, collapse = '\n')
  new_header <- stringr::str_replace_all(header, stringr::fixed(id_trans))
  stopifnot(new_header != header)
  run_bash(c(samtools, 'reheader', '--no-PG', '-', bam),
           input = new_header, stdout = bam_out)
  if (is.null(ref_seq)) {
    invisible(file.copy(paste0(prefix, '.fasta'),
                        paste0(out_prefix, '.fasta'),
                        overwrite = TRUE))
  } else {
    stopifnot(!is.null(names(ref_seq)))
    stopifnot(all(names(ref_seq) %in% names(id_map)))
    ref_seq <- ref_seq[names(id_map)]
    names(ref_seq) <- id_map
    write_dna(ref_seq, paste0(out_prefix, '.fasta'))
  }
  if (do_index) {
    invisible(run_bash(c(samtools, 'index', bam_out)))
  }
}

# Map reads against references and infer a consensus sequence
ambig_consensus <- function(seqs,
                            ref_seq,
                            out_prefix,
                            consensus_threshold = 0.65,
                            consensus_by_qual = TRUE,
                            fast = FALSE,
                            homopoly_fix = FALSE,
                            min_homopoly_len = 6,
                            cores = 1,
                            minimap2 = 'minimap2',
                            samtools = 'samtools') {
  # write reference and reads to file
  ref_file <- paste0(out_prefix, '.fasta')
  reads_file <- paste0(out_prefix, '_seqs.fastq')
  stopifnot(!is.null(names(ref_seq)))
  stopifnot(!is.na(names(ref_seq)))
  write_dna(ref_seq, ref_file)
  Biostrings::writeQualityScaledXStringSet(seqs, reads_file)

  # do mapping -> consensus
  args <- c(
    reads_file,
    ref_file,
    out_prefix,
    if (isTRUE(fast)) 'true' else 'false',
    cores,
    minimap2,
    samtools,
    '-m', 'simple',
    '-c', consensus_threshold,
    # extra treatment for homopolymers: lower confidence
    if (consensus_by_qual) {
      c(
        '--use-qual',
        '--homopoly-fix',
        '--homopoly-score', '0.3',
        '--qual-calibration', ':r10.4_sup'
      )
    } else {
      NULL
    }
  )

  # read output files
  stats <- run_bash_script('ref_consensus.sh', args, stdout = TRUE)
  cons_out <- paste0(out_prefix, '_consensus.fasta')
  invisible(file.remove(reads_file))
  # invisible(file.remove(ref_file))
  out <- data.frame(
    row.names = names(ref_seq),
    consensus = as.character(read_dna(cons_out))[names(ref_seq)],
    n_mapped = parse_idxstats(stats)[names(ref_seq)],
    homopolymer_adjustments = 0,
    consensus_diffs = 0
  )
  invisible(file.remove(cons_out))

  # compare reference with consensus and fix homopolymers
  different <- ref_seq != out$consensus
  if (any(different)) {
    out$consensus_diffs[different] <- NA
    stopifnot(length(homopoly_fix) %in% c(1, length(different)))
    homopoly_fix <- different & homopoly_fix
    if (any(homopoly_fix)) {
      res <- fix_homopolymers(out$consensus[homopoly_fix],
                              ref_seq[homopoly_fix],
                              min_homopoly_len = min_homopoly_len)
      cons_ <- out$consensus
      out[homopoly_fix, c('consensus', 'homopolymer_adjustments')] <- res[, c('consensus', 'n_adjusted')]
      adjusted <- out$homopolymer_adjustments > 0
      different <- ref_seq != out$consensus
      stopifnot(out$consensus[homopoly_fix & !adjusted] == cons_[homopoly_fix & !adjusted])
      out$consensus_diffs[adjusted & !different] <- 0
    }
  }
  # calculate distance (mismatches/InDels);
  # ambiguities in the consensus can still match the sequence
  # redo reference/consensus alignment for adjusted consensus seqs.
  if (any(different)) {
    # TODO: could reuse certain alignments from homopolymer fix
    sel <- is.na(out$consensus_diffs) | out$consensus_diffs > 0
    if (any(sel)) {
      aln <- pairwise_align(out$consensus[sel], ref_seq[sel])
      out$consensus_diffs[sel] <- get_aln_stats(aln, simplify = FALSE)[, 'diffs']
      stopifnot(is.finite(out$consensus_diffs))
    }
  }

  out$consensus_ambigs <- n_ambigs(out$consensus)

  out
}


#' Attempt at "fixing" the consensus sequence in a homopolymer region
#'
#' @param cons_seq the consensus
#' @param ref_seq the reference sequence (there should be some confidence
#' that it is correct)
#' @param min_homopoly_len minimum homopolymer length to check and fix
#'
#' @returns a data frame with these columns:
#' - *consensus*: the "fixed" consensus sequence
#' - *n_adjusted* the number of "fixed" homopolymer stretches
#'
#' @details
#' Resolves N-ambiguities in the consensus if surrounded by a homopolymer stretch,
#' replacing them with the corresponding reference sequence
#' (which is assumed to be the most likely true sequence)
#'
#' Ns are usually found on the left side (or near it), as minimap2 tends to
#' left-align gaps, and if there is an ambiguous situation (base + gap), the
#' consensus becomes an N.
#'
#' 1. A consensus/reference pairwise alignment is done
#' 2. The leftmost N is identified
#' 3. The first non-N base downstream (right) of this N is assumed to be
#'    the repeated base (likely true if gaps are left-aligned).
#'    Jump to 7. if none is found.
#' 4. The longest stretch of the given base/N/gaps is searched in the aligned
#'    consensus
#' 5. If the aligned reference contains only the given base and gaps
#'    within the range identified in 4., then the consensus sequence
#'    is replaced with the aligned reference
#' 6. Go back to 2., searching for the next N
#' 7. Remove all gaps from the consensus and return this sequence
#'
fix_homopolymers <- function(cons_seq, ref_seq, min_homopoly_len = 6) {
  # function for finding homopolymer runs at a specific position
  find_run <- function(x, values, at) {
    r <- rle(values)
    cum_lengths <- cumsum(r$lengths)
    j <- which(cum_lengths >= at)[1]
    stopifnot(r$values[j] == x)
    c(
      if (j == 1) 1 else cum_lengths[j-1] + 1,
      cum_lengths[j]
    )
  }

  aln <- pairwise_align(cons_seq, ref_seq, type = 'sequence')
  aln_consensus <- as.character(aln$PatternAligned)
  aln_ref <- as.character(aln$SubjectAligned)
  # find putative homopolymer runs containing 'N' ambiguities
  # (N can be base + gap)
  n_ascii <- 78 # utf8ToInt('N')
  gap_ascii <- 45 # utf8ToInt('-')
  dna_ascii <- c(65, 67, 71, 84) # utf8ToInt('ACGT')
  do.call(rbind, lapply(seq_len(length(aln_consensus)), function(i) {
    # cat(i, ' ')
    cons <- utf8ToInt(aln_consensus[i])
    ref <- utf8ToInt(aln_ref[i])
    n_adjusted <- 0
    n_pos <- 0
    while (n_pos < length(cons)) {
      # look for first N
      n_pos <- n_pos + which(cons[(n_pos + 1):length(cons)] == n_ascii)[1]
      if (is.na(n_pos)) {
        break
      }
      # look for first non-N base downstream of first N,
      # which is assumed to be the repeated base
      rest <- tail(cons, -n_pos)
      sel_base_ascii <- rest[rest %in% dna_ascii][1]
      if (is.na(sel_base_ascii)) {
        break
      }
      # find range of the putative homopolymer run in an rle-encoded
      # consensus, where all Ns and gaps are replaced with the given base
      cons_adj <- cons
      cons_adj[cons_adj %in% c(n_ascii, gap_ascii)] <- sel_base_ascii
      rng <- find_run(sel_base_ascii, cons_adj, at = n_pos)
      # as there may be gaps, this is not the precise homopolymer length
      # but we proceed with this candidate
      check_length <- rng[2] - rng[1] + 1
      if (check_length >= min_homopoly_len) {
        stopifnot(unique(cons[rng[1]:rng[2]]) %in% c(n_ascii, sel_base_ascii, gap_ascii))
        stopifnot(all(cons_adj[rng[1]:rng[2]] == sel_base_ascii))
        # Find homopolymer range of given base in the *reference*
        # (within the range determined from the consensus).
        # The final homopolymer range can thus be smaller (due to gaps)
        # but never larger than in the consensus
        ref_sub <- ref[rng[1]:rng[2]]
        n_sel_base <- sum(ref_sub == sel_base_ascii)
        # the reference sequence must only have the given base or gaps at this location
        is_pure <- all(ref_sub == sel_base_ascii | ref_sub == gap_ascii)
        if (is_pure && n_sel_base >= min_homopoly_len) {
          stopifnot(cons[rng[1]:rng[2]] %in% c(n_ascii, gap_ascii, sel_base_ascii))
          stopifnot(length(cons[rng[1]:rng[2]]) == length(ref_sub))
          cons[rng[1]:rng[2]] <- ref_sub
          n_adjusted <- n_adjusted + 1
        }
      }
    }
    cons_out <- intToUtf8(cons[cons != gap_ascii])
    stopifnot(xor(n_adjusted == 0, gsub('-', '', aln_consensus[i], fixed=T) != cons_out))
    data.frame(consensus = cons_out, n_adjusted = n_adjusted)
  }))
}

parse_idxstats <- function(lines) {
  stats <- sapply(strsplit(lines, '\t', fixed = T), '[', c(1, 3))
  stopifnot(stats[1, ncol(stats)] == '*')
  stats <- stats[, 1:(ncol(stats) - 1), drop = F]
  setNames(as.integer(stats[2, ]), stats[1, ])
}

move_bam <- function(prefix, out_prefix) {
  for (ext in c('.fasta', '.bam', '.bam.bai')) {
    f <- paste0(prefix, ext)
    if (file.exists(f)) {
      file.rename(f, paste0(out_prefix, ext))
    }
  }
}

remove_bam <- function(prefix, out_prefix, index=TRUE) {
  invisible(file.remove(paste0(prefix, c('.fasta', '.bam', if (index) '.bam.bai'))))
}

bam_to_map <- function(bam_file, samtools = 'samtools') {
  out <- run_bash(c(samtools, 'view', bam_file), stdout = TRUE)
  out <- read.delim(textConnection(out), header = FALSE)[, c(1, 3)]
  setNames(out[[2]], out[[1]])
}

# Fixed-threshold clustering (single-linkage by default)
#
# The clustering is done with [DECIPHER::Clusterize].
cluster_fixed <- function(seqs,
                         threshold,
                         min_coverage = 0.8,
                         single_linkage = TRUE,
                         cores = 1,
                         abund_order = FALSE,
                         verbose = FALSE) {
  if (is.null(names(seqs))) {
    names(seqs) = seq_len(length(seqs))
  }
  cl <- DECIPHER::Clusterize(
    Biostrings::DNAStringSet(seqs),
    cutoff = 1 - threshold,
    minCoverage = min_coverage,
    singleLinkage = single_linkage,
    processors = cores,
    verbose = verbose
  )[names(seqs), 'cluster']
  if (abund_order) {
    # reorder cluster IDs in descending order of abundance
    match(cl, order(tabulate(cl), decreasing = TRUE))
  } else {
    # order of occurrence
    match(cl, unique(cl))
  }
}
