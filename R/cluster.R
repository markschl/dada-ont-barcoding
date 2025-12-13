

#' Auto-install packages required by the *clustering* pipeline
install_cluster_deps <- function() {
  p <- installed.packages()[,c('Package')]
  install.packages(setdiff('BiocManager', p))
  p.bioc = setdiff(c('Biostrings', 'dada2', 'DECIPHER'), p)
  if (length(p.bioc) > 0)
    BiocManager::install(p.bioc)
}


# for parallel::clusterExport
get_barcodes_export <- c(
  'get_barcodes',
  'dada2_denoise',
  'split_haplotypes',
  'cluster_fixed',
  'ambig_consensus',
  'run_bash',
  'parse_idxstats',
  'pairwise_align',
  'merge_bam',
  'align_top',
  'bam_to_map',
  'n_ambigs'
)


#' Groups raw reads with DADA2 or fixed-threshold clustering to infer the correct
#' barcode sequence.
#'
#' @param fq demultiplexed FASTQ file
#' @param out_prefix Output prefix for cluster FASTA and aligned BAM files
#' @param tmp_dir path to a temporary directory
#' @param dada_omega_a OMEGA_A parameter value(s) to use for the denoising
#'     (see ?dada2::setDadaOpt):
#'    Can be a vector of increasing values (the default), which are tried
#'    sequentially until a grouping with well-separated haplotypes and
#'    non-ambiguous consensus sequence(s) is found for the *top taxon*
#'    (*note*: others currently ignored)
#' @param omegaA_iter_threshold Only evaluate multiple `dada_omega_a` values if the
#'    read depth is not larger than the given threshold. For large samples, it is
#'    better to avoid multiple denoising rounds, and it is anyway likely that DADA2
#'    will already find all sequence variants in the first round (with lowest OMEGA_A).
#'    Haplotype variation with InDels may still be found with 'split_haplotypes',
#'    but not with DADA2 denoising at high sensitivity.
#' @param dada_min_identical Minimum number of identical sequences required to do a
#'    DADA2 denoising. Below this threshold, switch to simple fixed-threshold
#'    clustering instead (at `cluster_threshold`) and report the consensus
#'    without further attempting any haplotype splitting.
#' @param dada_min_n0 Minimum number of error-free reads
#'   (`n0` in dada-class $clustering information) needed to retain the DADA2 clustering.
#'   Switch to the fixed-threshold clustering + consensus method below this threshold.
#'   *Note*: *error-free* does not mean *identical*, as DADA2 only considers
#'   substitutions, and there can still be InDels, so `n0` is always equal or higher
#'   to the number of identical sequences.
#' @param min_haplo_depth minimum number of supporting reads needed for haplotype
#'   clusters to be included
#' @param min_haplo_freq frequency threshold to consider clusters (ASVs) as
#'    haplotypes, not noise (relative to total abundance within all haplotypes of a
#'    source tissue, as clustered with `cluster_threshold`)
#' @param max_sample_depth read a maximum of `max_sample_depth` demultiplexed reads
#'    from the input file (`fq`) for the denoising/clustering
#' @param consensus_max_depth maximum number sequences mapped against the
#'    clusters to infer the consensus sequence (if there are more, a random sample
#'    is taken)
#' @param consensus_threshold require at least the given proportion of bases
#'    to be identical at every alignment column for an unambiguous consensus call
#'    (values below the default 60% might be problematic)
#' @param consensus_by_qual Consider the read quality scores when building the
#'    consensus with [samtools consensus](https://www.htslib.org/doc/samtools-consensus.html)
#' @param cluster_threshold similarity threshold for fixed-threshold clustering
#'    of low-depth samples (see `dada_min_identical` and `dada_min_n0`)
#'    and for post-clustering of DADA2 ASVs to group potential haplotypes of the
#'    same sequenced organism together.
#' @param cluster_single_linkage do single-linkage fixed-threshold clustering;
#'    clusters will grow as far as any 2 sequences have at least `cluster_threshold`
#'    similarity. If `FALSE`, the all cluster members are compared to one centroid
#'    sequence and remain more concise.
#'    On (`TRUE`) by default, as this appears to work well for fixed-threshold
#'    clustering of low-coverage samples and for grouping of haplotypes from the
#'    same organism tissue.
#' @param split_min_identical min. number of identical sequences that needs to
#'    support *both* of top two unique sequences in a sample in order to attempt
#'    haplotype splitting (with these sequences as new references).
#' @param max_split_ratio only split an ASV into 2 parental haplotypes if the ratio
#'    of larger:smaller is not larger than `max_split_ratio`.
#'    While the haplotypes might not always have an exact 1:1 ratio, this simple
#'    constraint still enforces a certain balancing of haplotype abundances.
#'    More unbalanced ratios are only possible if DADA2 is rerun with higher sensitivity
#'    (see `dada_omega_a`).
#' @value Returns a data frame with at least the following fields:
#'    - `id`: numeric cluster (haplotype) ID, or '<num>.<split>' where '<split>' is the
#'      split haplotype number (1 or 2)
#'    - `sequence`: DADA2 ASV or most abundant unique split sequence
#'       (`NA` in case of fixed-threshold clustering)
#'    - `abundance`: number of reads used for the clustering (up to `max_sample_depth`)
#'    - `n0` from DADA2: number of error-free reads (no substitutions, there may still be InDels)
#'    - `max_identical`: max. number of identical reads
#'    - `consensus`: consensus sequence of haplotype member sequences
#'    - `consensus_ambigs`: number of ambiguous bases in the consensus (see `consensus_threshold`)
#'    - `consensus_diffs`: edit distance between haplotype sequence
#'       (ASV or dominant sequence used as reference after splitting an ASV or fixed-threshold clustering)
#'    - `method`: method from which the haplotype emerged:
#'       one of 'dada', 'dada_split', 'fixed_cluster'
get_barcodes <- function(fq,
                        out_prefix,
                        dada_omega_a = c(1e-20, 1e-10, 1e-2),
                        omegaA_iter_threshold = 1000,
                        dada_min_identical = 2,
                        dada_min_n0 = 4,
                        min_haplo_depth = 3,
                        max_sample_depth = 5000,
                        consensus_max_depth = 3000,
                        consensus_threshold = 0.65,
                        consensus_by_qual = TRUE,
                        cluster_threshold = 0.97,
                        min_haplo_freq = 0.2,
                        split_min_identical = dada_min_identical,
                        max_split_ratio = 3,
                        cores = 1,
                        minimap2 = 'minimap2',
                        samtools = 'samtools',
                        verbose = FALSE) {
  
  tmp_dir <- paste0(out_prefix, '_tmp')
  dir.create(tmp_dir, FALSE, TRUE)
  sample_name <- if (verbose) gsub('\\.fastq(\\.gz)?$', '', basename(fq))
  
  # We will need the reads with quality scores below
  # cat('read/derep '); tictoc::tic()
  reads <- Biostrings::readQualityScaledDNAStringSet(fq, nrec = max_sample_depth)
  if (length(reads) == max_sample_depth) {
    # DADA2 derepFastq does not offer a way to only read the top N reads,
    # so we need to prepare a new FASTQ file for this
    fq <- file.path(tmp_dir, '_fq_limited.fastq')
    Biostrings::writeQualityScaledXStringSet(reads, fq)
  }
  
  # DADA2: dereplicate once for multiple rounds of denoising
  dada_derep <- dada2::derepFastq(fq, qualityType = 'FastqQuality')
  # tictoc::toc()
  
  d.prev = NULL
  for (dada_attempt in seq_along(dada_omega_a)) {
    round_prefix <- file.path(tmp_dir, paste0('round_', dada_attempt))
    # Run DADA2 denoising
    # note: denoising is done in every case even if the dada_min_identical threshold
    # is not met; DADA2 is anyway very fast with low-coverage samples
    # cat('dada '); tictoc::tic()
    d <- dada2_denoise(
      dada_derep,
      max_members = consensus_max_depth,
      OMEGA_A = dada_omega_a[dada_attempt],
      # controls how many seqs. end up in ASV
      # (higher than default 1e-40, meaning that more seqs. end up as ASV members;
      # appears to be a good value for low-coverage samples)
      OMEGA_C = 1e-10
    )
    d$method = 'dada'
    d$message = ''
    # tictoc::toc()
    dada_detail <- attributes(d)[c('n_seqs', 'n_singletons')]

    # Output should never be empty with DADA2 if DETECT_SINGLETONS is on
    stopifnot(nrow(d) > 0)
    
    # low-coverage (and/or low-quality) sample?
    if (max(d$max_identical[1]) < dada_min_identical || d$n0[1] < dada_min_n0) {
      # switch to fixed-threshold clustering if there are very few identical
      # reads or or few reads without substitution errors (n0)
      # https://github.com/benjjneb/dada2/issues/943#issuecomment-585219967
      # cat('cluster '); tictoc::tic()
      cl <- cluster_fixed(reads, cluster_threshold, abund_order = TRUE)
      # consistent/understandable ID
      cl <- sprintf('group%d_seq1', cl)
      
      # for each cluster, most abundant unique read (and its count),
      # OR the first-occurring sequence (if only singletons)
      top_uniq <- lapply(split(as.character(reads), cl), function(r) {
        sort(table(r), TRUE)[1]
      })
      # Choose the reference sequence to map the reads against:
      # This is only relevant for mapping; even in case of short InDels,
      # 'samtools consensus' should be able to return the correct consensus
      ref_seq <- sapply(top_uniq, names)
      cons <- ambig_consensus(
        reads,
        ref_seq,
        out_prefix = round_prefix,
        consensus_threshold = consensus_threshold,
        consensus_by_qual = consensus_by_qual,
        fast = TRUE,
        cores = cores,
        minimap2 = minimap2,
        samtools = samtools
      )
      # cons_cmp = pairwise_align(ref_seq, cons$consensus, cores = cores)
      d <- data.frame(
        row.names = names(ref_seq),
        id = names(ref_seq),
        # here: need to parse group number, can't use numeric ID in first place
        # as IDs should be understandable in FASTA/BAM
        group = as.integer(gsub('group([0-9]+)_.*', '\\1', names(ref_seq))),
        sequence = ref_seq,
        abundance = cons$n_reads,
        consensus = cons$consensus,
        consensus_ambigs = n_ambigs(cons$consensus),
        consensus_diffs = NA,
        # cons_cmp[, 'diffs'],
        max_identical = max(unlist(top_uniq)),
        method = 'cluster_fixed',
        check.names = FALSE
      )
      d$seq_indices = split(seq_along(reads), cl)[d$id]
      
      d <- d[d$abundance >= min_haplo_depth, ]
      # tictoc::toc()
      if (verbose)
        cat(
          sprintf(
            '%s (N=%s) Fixed clusters | top n=%d | ident = %s | ambig = %s\n',
            sample_name,
            dada_detail$n_seqs,
            d$abundance[1],
            d$max_identical[1],
            d$consensus_ambigs
          ),
          file = stderr()
        )
      # Attemting another DADA2 denoising with higher OMEGA_A does not make sense,
      # as it will not yield more identical seqs. -> stop
      break
      
    } # low-coverage end
    
    # depth filter
    d <- d[d$abundance >= min_haplo_depth, ]
    
    # post-cluster to group haplotypes by (putative) source organism tissue
    # TODO: cluster by consensus? should make little difference
    d$group = cluster_fixed(d$sequence, cluster_threshold)
    
    # assign understandable IDs
    d$id = with(d, sprintf('group%d_seq%d', group, ave(group, group, FUN=seq_along)))
    stopifnot(!duplicated(d$id))
    
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
      fast = TRUE,
      cores = cores,
      minimap2 = minimap2,
      samtools = samtools
    )
    # tictoc::toc()
    # TODO: read numbers may (rarely) not be exactly the same as reads can map to the
    # "wrong" ASV (mapping not restricted) -> should we use cons$n_reads?
    d$consensus = cons$consensus
    stopifnot(names(d$consensus) == d$id)
    d$consensus_ambigs = n_ambigs(d$consensus)
    
    if (any(d$consensus_ambigs > 0)) {
      # In case of ambiguities, try further splitting
      # into 2 haplotypes
      # cat('attempt split '); tictoc::tic()
      d <- split_haplotypes(
        reads,
        d,
        prefix = round_prefix,
        max_depth = max_sample_depth,
        min_identical = split_min_identical,
        max_ratio = max_split_ratio,
        consensus_threshold = consensus_threshold,
        consensus_by_qual = consensus_by_qual,
        fast = TRUE,
        cores = cores,
        minimap2 = minimap2,
        samtools = samtools,
        verbose = verbose
      )
      # tictoc::toc()
    }
    
    # Align ASV/haplotype sequences with BAM consensus
    # (differences = mismatches/InDels, should usually be 0)
    # Ambiguities in the consensus can still match the sequence.
    cons_cmp <- pairwise_align(d$sequence, d$consensus, cores = cores)
    d$consensus_diffs = cons_cmp[, 'diffs']
    stopifnot(is.finite(d$consensus_diffs))
    
    # reorder clusters by average abundance (keep members together)
    group_abund <- ave(d$abundance, d$group, FUN=sum)
    d <- d[order(-group_abund, d$group, -d$abundance), ]

    # determine number of haplotypes
    cl1 <- d$group[1]
    cl1_abund <- d$abundance[d$group == cl1]
    nhap <- sum(cl1_abund >= min_haplo_freq * sum(cl1_abund))
    if (verbose) {
      cat(
        sprintf(
          '%s | %d hap | n = %s | n0 = %s; ident = %s | diff = %s | ambig = %s\n',
          if (dada_attempt == 1) {
            sprintf('%s (N = %d)', sample_name, dada_detail$n_seqs)
          } else {
            paste("... omegaA =", dada_omega_a[dada_attempt])
          },
          nhap,
          paste(d$abundance[1:nhap], collapse = '/'),
          paste(d$n0[1:nhap], collapse = '/'),
          paste(d$max_identical[1:nhap], collapse = '/'),
          paste(d$consensus_diffs[1:nhap], collapse = '/'),
          paste(d$consensus_ambigs[1:nhap], collapse = '/')
        ),
        file = stderr()
      )
    }
    stopifnot(length(unique(d$group[1:nhap])) == 1)
    
    # go back to previous clustering if:
    # - all haplotypes are below the abundance threshold
    # - the total reads included in all clusters are >2x lower than last time
    #   (happens for very low-depth samples)
    # TODO: more conditions?
    if (!is.null(d.prev) &&
        (nrow(d) == 0 ||
         sum(d.prev$abundance) / sum(d$abundance) > 2)) {
      if (verbose)
        cat("...back\n", file = stderr())
      d <- d.prev
      dada_attempt <- dada_attempt - 1
      break
    }
    # stop with current clustering if there are no ambiguities in consensus of *top* group
    # OR the sample depth is large
    # (in which case it is unlikely that the haplotypes are not found with the initial omegaA)
    if (nrow(d) == 0 ||
        sum(d$consensus_ambigs[1:nhap]) == 0 ||
        dada_detail$n_seqs >= omegaA_iter_threshold) {
      break
    }
    d.prev = d
    
  } # end of loop
  
  if (nrow(d) == 0) {
    return(NULL)
  }
  
  # Create final BAM files:
  # - select clusters/haplotypes above min_haplo_depth
  # - remap reads to the consensus if the consensus is not identical to the reference sequence
  #   (as it is the consensus that we will report)
  chosen_round_prefix <- file.path(tmp_dir, paste0('round_', dada_attempt))
  # We have to use the consensus if there is ambiguous bases or differences to the reference
  use_consensus <- is.na(d$consensus_diffs) | d$consensus_diffs > 0 | d$consensus_ambigs > 0
  # cat('assemble output '); tictoc::tic()
  if (all(!use_consensus)) {
    # TODO: if all refs are chosen from the source BAM file, then this is
    # not very efficient (simple copying would suffice); this may be checked in merge_bam()
    merge_bam(
      out_prefix = out_prefix,
      list(list(chosen_round_prefix, d$id)),
      cores = cores,
      samtools = samtools
    )
  } else {
    sel_reads <- reads[unlist(d$seq_indices[use_consensus])]
    ref_seq <- setNames(d$consensus[use_consensus], d$id[use_consensus])
    remap_prefix <- if (all(use_consensus)) {
      out_prefix
    } else {
      file.path(tmp_dir, 'remap_cons')
    }
    cons <- ambig_consensus(
      sel_reads,
      ref_seq,
      out_prefix = remap_prefix,
      consensus_threshold = consensus_threshold,
      consensus_by_qual = consensus_by_qual,
      fast = !all(use_consensus),
      cores = cores,
      minimap2 = minimap2,
      samtools = samtools
    )
    stopifnot(!is.na(cons$consensus))
    cons_same <- d$consensus[use_consensus] == cons$consensus
    if (any(!cons_same)) {
      # TODO: should not (and does not seem to) happen often,
      # but it is possible -> should rerun?
      d$message[use_consensus][!cons_same] <- paste0(
        d$message[use_consensus][!cons_same],
        'Re-mapping to consensus gives another consensus! '
      )
    }
    if (!all(use_consensus)) {
      merge_bam(
        out_prefix = out_prefix,
        list(
          list(remap_prefix, d$id[use_consensus]),
          list(chosen_round_prefix, d$id[!use_consensus])
        ),
        cores = cores,
        samtools = samtools
      )
    }
  }
  # tictoc::toc()
  
  # Now we delete entries of 'sequence' that are not well supported by duplicates
  d$sequence[d$max_identical <= dada_min_identical] = NA
  
  # clean up
  unlink(tmp_dir, TRUE)
  d$seq_indices = d$top_uniques = NULL
  
  # calculate homopolymer stretch length
  l <- rle(utf8ToInt(d$consensus[1]))$lengths
  d$max_homopoly_len = max(l)
  
  # flag rare haplotypes
  d$is_rare <- as.logical(ave(d$abundance, d$group, 
                              FUN=function(x) x / sum(x) < min_haplo_freq))

  attributes(d)[names(dada_detail)] = dada_detail
  attr(d, 'omega_a') = dada_omega_a[dada_attempt]
  d
}


#' Create a multiple alignment of the haplotypes from the top organism group
#' in a data frame produced by 'get_barcodes', which may further have been reordered to
#' prioritize another species in the mix if the top group was a contaminant.
align_top <- function(d,
                     out_prefix,
                     known_seq = NA,
                     cores = 1) {
  cl <- d$group[1]
  seqs <- do.call(c, lapply(which(d$group == cl & !d$is_rare), function(i) {
    d.sel <- d[i,]
    out <- if (!is.na(d.sel$sequence) && d.sel$sequence != d.sel$consensus) {
      setNames(c(d.sel$sequence, d.sel$consensus),
               paste0(d.sel$id, c('_dominant_seq', '_consensus')))
    } else {
      setNames(d.sel$consensus, d.sel$id)
    }
    names(out) <- paste(
      names(out), 
      sprintf('reads: %d | ambigs: %d | identical: %d | substitution-free: %d | %s',
              d.sel$abundance, d.sel$consensus_ambigs, d.sel$max_identical, d.sel$n0,
              d.sel$message)
    )
    out
  }))
  if (!is.null(known_seq) && !is.na(known_seq)) {
    seqs <- c(seqs, setNames(known_seq, 'known_sequence'))
  }
  # aln_seqs = seqs = Biostrings::DNAStringSet(seqs)
  if (length(seqs) > 1) {
    seqs <- Biostrings::DNAStringSet(seqs)
    aln_seqs <- DECIPHER::AlignSeqs(seqs, verbose = FALSE, processors = cores)
    Biostrings::writeXStringSet(aln_seqs, paste0(out_prefix, sprintf('_cluster%d_comparison.fasta', cl)))
  }
}

dada_learn_errors <- function(fq_paths, omega_a = 1e-20, cores = 1) {
  # To be sure we increase the NW alignment band size,
  # which is recommended for technologies with potentially some InDels (454, PacBio)
  # and also applies for Nanopore.
  dada2::learnErrors(fq_paths,
                     nbases = 1e8,
                     BAND_SIZE = 32,
                     qualityType = 'FastqQuality',
                     OMEGA_A = omega_a,
                     # this setting is hard-coded in denoise_dada2
                     OMEGA_C = 1e-10,
                     multithread = cores, 
                     verbose = 0)
  
}

#' Runs a DADA2 denoising on a FASTQ file, removes chimeras and and returns
#' an abundance-sorted data frame of ASVs with some additional information
#'
#' `singleton_threshold`: by default, singletons cannot form a new cluster, but
#' with very low-depth samples (< 5 duplicates of any sequence), DETECT_SINGLETONS
#' is turned on (see ?dada2::setDadaOpt). This increases sensitivity with InDel-rich
#' Nanopore data.
dada2_denoise <- function(x,
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
    multithread = cores,
    verbose = 0,
    ...
  )
  
  d <- dd$clustering[, c('sequence', 'abundance', 'n0', 'nunq'), drop = FALSE]
  # assign an ID, as we will reorder/subset later, but need the cluster number
  d$asv_num = 1:nrow(d)
  
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
  d$seq_indices = cl[as.character(d$asv_num)]
  
  # In addition, remember the two most abundant sequences, used later for haplotype splitting
  d$top_uniques = lapply(d$seq_indices, function(i) {
    if (length(i) > max_members) {
      # random subsampling for overly large ASVs
      i <- sample(i, max_members)
    }
    uniq_freq <- tabulate(derep$map[i])
    i <- head(order(-uniq_freq), 2)
    setNames(uniq_freq[i], names(derep$uniques)[i])
  })
  
  # max. number of identical sequences (usually identical to inferred ASV sequence)
  d$max_identical = sapply(d$top_uniques, max)
  
  # attach more information
  attr(d, 'n_seqs') = length(derep$map)
  attr(d, 'n_singletons') = sum(tabulate(derep$map) == 1)
  d$asv_num <- NULL
  d
}

#' Attempts splitting an ASV into two haplotypes.
#' Requires the output of dada2_denoise and ambig_consensus:
#' - id
#' - top_uniques
#' - consensus_ambigs
#' 
#' Splitting is done if the sum of the consensus ambiguities of the two resulting
#' sequence variants is smaller than the number of consensus ambiguities of the parent.
#'
#' Overwrites the existing BAM files and FASTA references in 'prefix'.
split_haplotypes <- function(reads,
                            d,
                            prefix,
                            max_depth = 5000,
                            min_identical = 2,
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
    # (rarely reads might switch to a different ASV)
    cons <- ambig_consensus(
      sel_reads,
      ref_seq ,
      out_prefix = prefix_split,
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
      cons <- grouped_cons[[d.s$id]]
      cons_ambigs <- n_ambigs(cons$consensus)
      # required for splitting:
      # 1) sum of split ambiguities < non-split ambiguities
      # 2) the abundance ratio is not too high
      # TODO: enforcing an 1:1 ratio with a binomial test seems too strict
      # ratio.p = binom.test(cons$n_reads[1], sum(cons$n_reads), p = 0.5)$p.value
      abund.ratio = max(cons$n_reads) / min(cons$n_reads)
      is_better <- !any(is.na(cons$consensus)) &&
        sum(cons_ambigs) < d.s$consensus_ambigs &&
        abund.ratio <= max_ratio
      if (verbose) {
        cat(
          sprintf(
            ' split %s %s: has cons: %s | ambigs %d -> %d/%d = %d | uniq %d:%d | hap %d:%d | ratio=%.2f\n',
            d.s$id,
            is_better,!any(is.na(cons$consensus)),
            d.s$consensus_ambigs,
            cons_ambigs[1],
            cons_ambigs[2],
            sum(cons_ambigs),
            d.s$top_uniques[[1]][1],
            d.s$top_uniques[[1]][2],
            cons$n_reads[1],
            cons$n_reads[2],
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
        d2 <- data.frame(
          id = row.names(cons),
          group = d.s$group,
          sequence = names(uniq),
          consensus = cons$consensus,
          consensus_ambigs = cons_ambigs,
          abundance = cons$n_reads,
          max_identical = unname(uniq),
          method = paste0(d.s$method, '_split'),
          message = d.s$message
        )
        d2$top_uniques <- list(NULL, NULL)
        # obtaining seq_indices is a bit more involved, but usually
        # the number of reads is not too large here, so this should be fast
        read_map <- bam_to_map(paste0(prefix_split, '.bam'), samtools = samtools)
        read_ids <- sapply(strsplit(names(reads), ' ', fixed = TRUE), '[', 1)
        d2$seq_indices = split(match(names(read_map), read_ids), read_map)[row.names(cons)]
        stopifnot(lengths(d2$seq_indices) == cons$n_reads)
        # all other columns should be undefined
        d2[,setdiff(names(d.s), names(d2))] <- NA
        d2$do_split = TRUE
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
        for (suffix in c('.fasta', '.bam', '.bam.bai')) {
          file.rename(paste0(prefix_split, suffix),
                      paste0(prefix, suffix))
        }
      } else {
        # extract data for split and non-split references (using samtools view), then merge them
        prefix_merged <- paste0(prefix, '.merged')
        merge_list <- list(list(prefix, d$id[!d$do_split]),
                          list(prefix_split, d$id[d$do_split]))
        merge_bam(
          out_prefix = prefix_merged,
          merge_list,
          fast = fast,
          cores = cores,
          samtools = samtools
        )
        # overwrite originals
        for (ext in c('.fasta', '.bam', '.bam.bai')) {
          file.rename(paste0(prefix_merged, ext), paste0(prefix, ext))
        }
        # clean up
        file.remove(c(paste0(
          prefix_split, c('.fasta', '.bam', '.bam.bai')
        )))
      }
    }
  }
  
  d$do_split = NULL
  d
}


merge_bam <- function(out_prefix,
                     merge_list,
                     cores = 1,
                     fast = FALSE,
                     samtools = 'samtools') {
  # merge BAM
  bam_out <- paste0(out_prefix, '.bam')
  msg <- run_bash(
    c(
      samtools,
      'merge',
      '--no-PG',
      '-f',
      if (fast) '-1' else NULL,
      '-@',
      cores,
      '-o',
      bam_out,
      sapply(merge_list, function(x) {
        sprintf(
          '<(samtools view --no-PG -uT %s %s %s)',
          paste0(x[[1]], '.fasta'),
          paste0(x[[1]], '.bam'),
          paste(x[[2]], collapse = ' ')
        )
      })
    ),
    stderr = TRUE,
    stdout = TRUE
  )
  if (length(msg) > 0 && grepl('coordinate sort to be lost', msg, fixed = TRUE)) {
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
        '-l',
        if (fast) 1 else 6,
        '-@',
        cores,
        '-o',
        bam_out,
        merged_f
      ),
      stderr = FALSE
    )
    file.remove(merged_f)
  }
  # index
  run_bash(c('samtools', 'index', bam_out))
  for (x in merge_list) {
    file.remove(paste0(x[[1]], '.fasta.fai'))
  }
  # merge references
  ref_out <- paste0(out_prefix, '.fasta')
  refs <- do.call(c, lapply(merge_list, function(x) {
    Biostrings::readDNAStringSet(paste0(x[[1]], '.fasta'))[as.character(x[[2]])]
  }))
  Biostrings::writeXStringSet(refs, ref_out)
}

ambig_consensus <- function(seqs,
                           ref_seq,
                           out_prefix,
                           consensus_threshold = 0.6,
                           consensus_by_qual = TRUE,
                           fast = FALSE,
                           cores = 1,
                           minimap2 = 'minimap2',
                           samtools = 'samtools') {
  ref_file <- paste0(out_prefix, '.fasta')
  reads_file <- paste0(out_prefix, '_seqs.fastq')
  
  stopifnot(!is.null(names(ref_seq)))
  stopifnot(!is.na(names(ref_seq)))
  
  Biostrings::writeXStringSet(Biostrings::DNAStringSet(ref_seq), ref_file)
  Biostrings::writeQualityScaledXStringSet(seqs, reads_file)
  cmd <- c(
    'scripts/ref_consensus.sh',
    reads_file,
    ref_file,
    out_prefix,
    if (isTRUE(fast)) 'true' else 'false',
    cores,
    minimap2,
    samtools,
    '-m',
    'simple',
    '-c',
    consensus_threshold,
    # extra treatment for homopolymers: lower confidence
    if (consensus_by_qual) {
      c(
        '--use-qual',
        '--homopoly-fix',
        '--homopoly-score',
        '0.3',
        '--qual-calibration',
        ':r10.4_sup'
      )
    } else {
      NULL
    }
  )
  stats <- run_bash(cmd, stdout = TRUE)
  cons_out <- paste0(out_prefix, '_consensus.fasta')
  invisible(file.remove(reads_file))
  # invisible(file.remove(ref_file))
  n_reads <- parse_idxstats(stats)
  consensus <- as.character(Biostrings::readDNAStringSet(cons_out))
  invisible(file.remove(cons_out))
  data.frame(n_reads = n_reads[names(ref_seq)], consensus = consensus[names(ref_seq)])
}

parse_idxstats <- function(lines) {
  stats <- sapply(strsplit(lines, '\t', fixed = T), '[', c(1, 3))
  stopifnot(stats[1, ncol(stats)] == '*')
  stats <- stats[, 1:(ncol(stats) - 1), drop = F]
  setNames(as.integer(stats[2, ]), stats[1, ])
}

bam_to_map <- function(bam_file, samtools = 'samtools') {
  out <- run_bash(c(samtools, 'view', bam_file), stdout = TRUE)
  out <- read.delim(textConnection(out), header = FALSE)[, c(1, 3)]
  setNames(out[[2]], out[[1]])
}

#' Fixed-threshold clustering (single-linkage by default)
cluster_fixed <- function(seqs,
                         threshold,
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
    minCoverage = 0.8,
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
