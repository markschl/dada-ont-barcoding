
#' Reads a sample table from an Excel metadata file and passes it on to [parse_sample_tab]
#'
#' @export
read_xlsx_sample_tab <- function(meta_file, sheet_name) {
  sample_tab <- openxlsx2::read_xlsx(
    meta_file,
    sheet_name,
    na.strings = c('#N/A', '')
  )
  # sample_tab <- tibble::as_tibble(sample_tab, .name_repair='unique_quiet')
  # empty cells might still have '0' (if a formula), we set these to NA
  sample_tab$sample = ifelse(!is.na(sample_tab$sample) & sample_tab$sample == '0',
                             NA,
                             sample_tab$sample)
  parse_sample_tab(sample_tab)
}

#' Parse a data frame with sample metadata
#'
#' @param sample_tab The sample table (data frame-like, see details)
#'
#' @details
#'
#' Required columns:
#'
#' - Plate
#' - *plate*: A01-H12)
#' - *amplicon*: `forward-reverse`, name as in 'primer_index' column of primer table;
#'   see [parse_primer_tab]
#' - *indexes*: `forward-reverse`, name as in 'primer_index' column of primer table;
#'   see [parse_primer_tab]
#' - *sample*: sample name, should not be duplicated
#' - *sample_type*: Sample type such as 'negative control', will show up in reports;
#'   can be left empty
#'
#' Optional columns:
#'
#' - *taxon*: Known or suspected taxon (any rank), will be compared to
#'   the sequence-based identification and used for detecting contamination
#' - *known sequence*: Already known sequence (if any), e.g. from a previous Sanger
#'   sequencing; will be compared with the Nanopore-derived sequence
#'
#' @returns a valid sample metadata table (data frame)
#'   with sample names de-duplicated (if necessary)
#'
#' @details
#' ## Note on amplicon multiplexing
#'
#' Primers are searched in the order that amplicons appear in the sample table
#' With *nested amplicons*, the shorter one should be placed at the end.
#' Also place amplicons with *little data* before other amplicons with similar
#' primers to make sure that they don't get "swallowed" in case of unspecific
#' primer or sample index matching.
#'
#' @export
parse_sample_tab <- function(sample_tab) {
  # normalize & validate headers
  names(sample_tab) <- gsub('[- .]+', '_', tolower(names(sample_tab)))
  required_cols <- c('plate', 'well', 'amplicon', 'indexes', 'sample', 'sample_type')
  missing_cols <- setdiff(required_cols, names(sample_tab))
  if (length(missing_cols) > 0) {
    stop('The sample metadata table is missing some columns: ',
         paste(missing_cols, collapse=', '))
  }

  # remove barcodes without sample name
  sample_tab <- sample_tab[!is.na(sample_tab$sample),]

  # Check for duplicates

  sample_tab$unique_sample = sample_tab$sample
  # sample_tab <- sample_tab %>%
  #   relocate(unique_sample, .after='sample')
  dup_s <- unique(sample_tab$unique_sample[duplicated(sample_tab$unique_sample)])
  if (length(dup_s[dup_s != '0']) > 0) {
    warning(paste0('The following samples are duplicated: ',
                   paste(dup_s, collapse=', '),
                   '. The plate/coordinate was added to the name.'))
    sel <- sample_tab$unique_sample %in% dup_s
    sample_tab$unique_sample[sel] = with(sample_tab[sel, ], paste(unique_sample, plate, well, sep='_'))
  }

  # check known sequence (enforce IUPAC alphabet)
  if (!is.null(sample_tab$known_sequence) && any(!is.na(sample_tab$known_sequence))) {
    # do some simple cleaning
    sample_tab$known_sequence <- gsub('[\r\n\t -.]', '', sample_tab$known_sequence)
    i <- charToRaw(paste0(Biostrings::DNA_ALPHABET, collapse=''))
    i <- i[i >= 65]
    sel <- !is.na(sample_tab$known_sequence)
    is_valid <- sapply(sample_tab$known_sequence[sel] ,
                       function(s) all(charToRaw(s) %in% i))
    if (!all(is_valid)) {
      stop("Some samples have an invalid 'known seqence': ",
           paste(seq_tab$sample[sel][!is_valid], collapse=', '))
    }
  }

  # convert amplicon to factor (order matters)
  amplicons <- unique(sample_tab$amplicon)
  names(amplicons) = amplicons  # useful for lapply -> bind_rows
  sample_tab$amplicon = factor(sample_tab$amplicon, amplicons)
  sample_tab
}

#' @export
read_xlsx_primer_tab <- function(meta_file, sheet_name, ...) {
  primer_tab <- openxlsx2::read_xlsx(meta_file, sheet_name, na.strings = c(''))
  parse_primer_tab(primer_tab, ...)
}

#' Parse a data frame with primer information
#'
#' @param primer_tab Primer table (data frame-like; see details)
#' @param amplicons Character vector of amplicons pooled in the library;
#' their names should be in the form `forward-reverse` (primer names should not contain dashes)
#'
#' @details
#'
#' Columns of `primer_tab`:
#'
#' - *primer-index*: Primer name and sample index name, delimited by a dash;
#'   the primer names **must not contain dashes** and both *primer* and *index*
#'   should correspond to the names in the *amplicon* column of the sample table
#'   (see [parse_sample_tab]).
#' - *index_seq* sample index sequence
#' - *primer_seq* primer sequence (repeated)
#'
#' If `amplicon` is not provided, the forward/reverse primers are assumed to be
#' in order. `amplicon` may be obtained from the *amplicon* column of the sample
#' table (see [parse_sample_tab]).
#'
#' @returns a named list of amplicon metadata (named by amplicons = `forward-reverse`).
#' Each amplicon-specific item is another list with *two entries* (for forward and reverse),
#' which are again lists with following entries:
#'
#'  - *primer*: named vector of primer sequences
#'  - *index*: named vector of index sequences
#'  - *index_len*: the length of the index sequences (no variable length currently allowed,
#'    but lengths may vary between forward and reverse sample indexes)
#'
#' @export
parse_primer_tab <- function(primer_tab, amplicons = NULL) {
  # normalize & validate headers
  names(primer_tab) <- gsub('[- .]+', '_', tolower(names(primer_tab)))
  required_cols <- c('primer_index', 'index_seq', 'primer_seq')
  missing_cols <- setdiff(required_cols, names(primer_tab))
  if (length(missing_cols) > 0) {
    stop('The primers table is missing some columns: ',
         paste(missing_cols, collapse=', '))
  }
  primer_tab$primer_name <- gsub('-.*', '', primer_tab$primer_index)
  stopifnot(!duplicated(primer_tab$primer_index))

  # validate amplicons and primer names
  primers <- split(primer_tab, primer_tab$primer_name)
  if (is.null(amplicons)) {
    stopifnot(length(primers) == 2)
    amplicons <- unique(primer_tab$primer_name)
  }
  amp_primer_names <- strsplit(amplicons, '-', fixed=TRUE)
  names(amp_primer_names) <- amplicons

  primers_only <- setdiff(names(primers), unlist(amp_primer_names))
  amp_only <- setdiff(unlist(amp_primer_names), names(primers))
  if (length(primers_only) > 0) {
    stop("Some primer name(s) found only in the primers table sheet, but not ",
         "in the amplicons (samples table): ",
         paste(primers_only, collapse=', '))
  }
  if (length(amp_only) > 0) {
    stop("Some primer name(s) found only in the amplicon names (samples table), ",
         "but not in the primers table: ",
         paste(amp_only, collapse=', '))
  }

  # collect per-amplicon primer/index information
  lapply(amp_primer_names, function(p) {
    if (length(p) != 2) {
      stop("The 'amplicon' column in the samples list needs to be in the form 'forward-reverse' ",
           "where the forward/reverse primer name has to match the primer names in ",
           "the 'primer-index' column of the primers table (not dashes allowed in primer names).")
    }
    comb <- list(forward=primers[[p[1]]], reverse=primers[[p[2]]])
    lapply(comb, function(d) {
      primer <- unique(d$primer_seq)
      if (length(primer) != 1) {
        stop(sprintf("The primer '%s' does not have the same sequence in all rows of the 'primer' column",
                     d$primer_name[1]))
      }
      stopifnot(!duplicated(d$index_seq))
      index_len <- unique(nchar(d$index_seq))
      if (length(index_len) != 1) {
        # TODO: in theory we could support this (just affects how we trim the sequences before the primer)
        stop(sprintf("The barcodes for primer '%s' should all have the same length", primer))
      }
      list(
        primer = setNames(primer, d$primer_name[1]),
        index = setNames(d$index_seq, d$primer_index),
        index_len = index_len
      )
    })
  })
}

#' Search and remove primers and sample indexes
#'
#' @param fq_paths Character vector of FASTQ file paths (compressed or uncompressed)
#' @param amplicon_primers Amplicon specification, as returned by [parse_primer_tab]
#' @param sample_tab Sample table as returned by [parse_sample_tab]
#' @param out_dir Output directory containing the demultiplexed files
#' (`fprimer_fidx-rprimer_ridx.fastq.gz`)
#' @param primer_max_err Maximum allowed error rate in primers
#' (edit distance = substitutions/InDels). Default: 0.2 (20%)
#' @param idx_max_diffs Maximum allowed differences in forward/reverse sample
#' index sequences; don't set too high to avoid tag switching ("index hopping")
#' resulting in a "background noise" of unspecific sequences. Default: 0 differences
#' @param min_barcode_length Remove trimmed sequences shorter than `min_barcode_length`;
#' make sure that the value is shorter than your expected minimal barcode length
#' for *all* amplicons (default: 50 bp)
#' @param error_threshold Maximum sequencing errors allowed per sequence
#' (as estimated from quality scores). This threshold impacts how many sequences
#' pass the filter in (see [do_demux]). Setting it too low may remove too many
#' sequences, but setting it too high may lead to noisy/errorneous results.
#' The default value of 2.5 seems to work well with R10.4 data and ~0.5-1.5 kb amplicons,
#' removing about half of the trimmed sequences.
#  The value may still be adjusted based on how Figures 4 and 5 in the HTML report.
#'
#' @details
#' This function calls [do_primer_search] and [do_demux] sequentially.
#'
#' Requires [seqtool](https://github.com/markschl/seqtool) v0.4 or higher
#' installed either system-wide (in `PATH`) or locally
#' (provide path with `set_program_path('seqtool', 'path/to/st')`).
#'
#' @export
do_trim_demux <- function(fq_paths,
                          amplicon_primers,
                          sample_tab,
                          out_dir,
                          primer_max_err = 0.2,
                          idx_max_diffs = 0,
                          min_barcode_length = 50,
                          error_threshold = 2.5,
                          cores = 1,
                          keep_trimmed = FALSE) {
  trim_dir <- file.path(out_dir, '_trim')
  trim <- do_primer_search(
    fq_paths,
    amplicon_primers,
    trim_dir,
    primer_max_err = primer_max_err,
    idx_max_diffs = idx_max_diffs,
    min_barcode_length = min_barcode_length,
    cores = cores
  )
  n <- trim$stats$counts
  if (sum(n$count[n$valid]) == 0) {
    stop('No reads with primer/sample index found!')
  }
  seq_tab <- do_demux(trim$trimmed_fq, out_dir, sample_tab, error_threshold = error_threshold)
  if (!keep_trimmed) {
    unlink(trim_dir, TRUE)
  }
  list(seq_tab = seq_tab, trim_stats = trim$stats)
}

#' Search and remove primer sequences
#'
#' Requires [seqtool](https://github.com/markschl/seqtool) v0.4 or higher
#'
#' @export
do_primer_search <- function(fq_paths,
                              amplicon_primers,
                              out_dir,
                              primer_max_err = 0.2,
                              idx_max_diffs = 0,
                              min_barcode_length = 50,
                              cores = 1) {
  amplicons <- names(amplicon_primers)
  seqtool <- get_program('st', long_name = 'seqtool')
  for (amplicon in amplicons) {
    amp_search_dir <- file.path(out_dir, amplicon)
    demux_fq <- file.path(amp_search_dir, 'trimmed.fastq.zst')

    if (length(fq_paths) > 0) {
      cat(amplicon, sep='\n', file=stderr())
      unlink(amp_search_dir)  # clean up old files
      dir.create(amp_search_dir, FALSE, TRUE)
      amp_pr <- amplicon_primers[[amplicon]]
      seq_prefix <- file.path(amp_search_dir, '')
      primer_paths <- paste0(seq_prefix, c('fwd', 'rev'), '_primers.fasta')
      idx_paths <- paste0(seq_prefix, c('fwd', 'rev'), '_idx.fasta')
      for (i in 1:2) {
        write_dna(amp_pr[[i]]$primer, primer_paths[i])
        write_dna(amp_pr[[i]]$index, idx_paths[i])
      }
      search_opts <- c(
        '-p', primer_max_err,
        '-b', idx_max_diffs,
        '-l', min_barcode_length,
        '-o', amp_search_dir,
        '-c', 'zst',
        '-t', cores,
        '-s', seqtool
      )
      run_bash(c('scripts/find-primers.sh', search_opts, seq_prefix, fq_paths))
    }
    # input for next amplicon search
    fq_names <- paste0('no_',
                       c('primer.fwd', 'primer.rev', 'index.fwd', 'index.rev'),
                       '.fastq.zst')
    fq_paths <- file.path(amp_search_dir, fq_names)
  }

  # read stats files
  # primer start position
  pos <- do.call(rbind, lapply(amplicons, function(amplicon) {
    d <- read.delim(file.path(out_dir, amplicon, 'primer_pos.tsv'),
                  colClasses=c('character', 'integer', 'integer'))
    d$amplicon <- amplicon
    d$dir <- factor(d$dir, c('fwd', 'rev'), c('forward', 'reverse'))
    aggregate(count ~ amplicon + dir + pos, data=d, sum)
  }))

  # amplicon length
  category_trans <- c(trimmed='regular', concatenated='concatenated products')
  l <- do.call(rbind, lapply(names(amplicon_primers), function(amplicon) {
    d <- read.delim(file.path(out_dir, amplicon, 'length_stats.tsv'),
                    colClasses=c('character', 'integer', 'integer'))
    d$amplicon <- amplicon
    d
  }))
  l$category <- trimws(gsub('\\.fastq\\.zst$', '', l$file))
  l$category = factor(category_trans[l$category], category_trans)
  l$file <- NULL

  # primer/index counts
  n <- do.call(rbind, lapply(seq_along(amplicons), function(i) {
    d <- read.delim(file.path(out_dir, amplicons[i], 'trim_counts.tsv'),
               colClasses=c('character', 'integer'))
    d$file <- trimws(gsub('\\.fastq\\.zst$', '', d$file))
    d$valid <- d$file == 'trimmed'
    if (i != length(amplicons)) {
      # all data used in subsequent primer searches should be excluded
      d <- d[!startsWith(d$file, 'no_'),]
    }
    d$amplicon <- amplicons[i]
    d
  }))
  s <- stringr::str_split_fixed(n$file, stringr::fixed('.'), 2)
  n$direction <- s[,2]
  n$category = ifelse(
    n$valid,
    'valid amplicon',
    gsub('_', ' ', s[,1])
  )
  n$file <- NULL

  # quality
  q <- do.call(rbind, lapply(seq_along(amplicons), function(i) {
    d <- read.delim(file.path(out_dir, amplicons[i], 'qual_stats.tsv'),
               colClasses=c('character', rep('integer', 7)),
               na.strings = 'undefined')
    d$amplicon <- amplicons[i]
    d$exp_err = as.integer(gsub(r'{\(\d+, (\d+)\]}', '\\1', d$exp_err))
    d$primer_mis = d$f_primer_mis + d$r_primer_mis
    d$idx_mis = d$f_idx_mis + d$r_idx_mis
    if (i != length(amplicons)) {
      # all data used in the subsequence primer searches should be excluded
      d <- d[!is.na(d$primer_mis) & !is.na(d$idx_mis),]
    }
    d
  }))

  list(
    trimmed_fq = file.path(out_dir, amplicons, 'trimmed.fastq.zst'),
    stats = list(
      position = pos, amplicon_len = l, counts = n, quality = q
    )
  )
}

#' Group sequences by sample index combination
#'
#' The provided FASTQ file should contain the sample indexes in the sequence headers
#' as follows: `<id> fi=<fwd-idx-name> ri=<rev-idx-name>` (the position in the header
#' does not matter). Such FASTQ files are produced by [do_primer_search].
#'
#' Requires [seqtool](https://github.com/markschl/seqtool) v0.4 or higher
#'
#' @returns `sample_tab` with a `reads_path` column and additional rows
#' for samples found in the reads but not in the sample table
#'
#' @details
#' The read files are named as follows: `fprimer-fidx--rprimer-ridx.fastq.gz`
#'
#' @export
do_demux <- function(primer_search_fq,
                      out_dir,
                      sample_tab,
                      error_threshold = 2.5,
                      min_barcode_length = 50) {

  stopifnot(c('amplicon', 'indexes', 'sample', 'sample_type', 'taxon', 'known_sequence') %in% names(sample_tab))
  seqtool <- get_program('st', long_name = 'seqtool')
  run_bash(c('scripts/filter-split.sh',
             error_threshold,
             # **note**: currently, we have only one threshold both in 'find-primers' and 'demux'
             min_barcode_length,
             out_dir,
             seqtool,
             primer_search_fq))
  trimmed_fq <- file.path(out_dir, list.files(out_dir, pattern='.fastq.gz$'))

  # initialize the sequence table
  sample_tab$i_ <- seq_len(nrow(sample_tab))
  trimmed_fq <- sort(trimmed_fq)
  t <- data.frame(
    reads_path = trimmed_fq,
    primer_indexes = gsub('\\.fastq.gz$', '', basename(trimmed_fq)),
    i_ = nrow(sample_tab) + seq_len(length(trimmed_fq))
  )
  t <- cbind(t, stringr::str_match(
    t$primer_indexes, '(?<fpr>[^-]+)-(?<fidx>[^-]+)--(?<rpr>[^-]+)-(?<ridx>[^-]+)'
  ))
  t$amplicon <- paste(t$fpr, t$rpr, sep='-')
  t$indexes = paste(t$fidx, t$ridx, sep='-')
  t <- merge(
    sample_tab,
    t[c('reads_path', 'primer_indexes', 'amplicon', 'indexes', 'i_')],
    by = c('amplicon', 'indexes'),
    all = TRUE
  )
  t <- t[order(t$i_.x, t$i_.y),]
  t$i_.x <- t$i_.y <- NULL
  p <- stringr::str_split_fixed(t$amplicon, '-', 2)
  i <- stringr::str_split_fixed(t$indexes, '-', 2)
  t$primer_indexes <- sprintf('%s-%s--%s-%s', p[,1], i[,1], p[,2], i[,2])
  if (length(unique(na.omit(t$amplicon))) > 0) {
    t$indexes <- t$primer_indexes
  }
  t$category <- ifelse(is.na(t$sample_type) | t$sample_type == '',
                       'sample with reads',
                       as.character(t$`sample type`))
  t$category[is.na(t$reads_path)] <- 'no reads'
  t$category[is.na(t$sample)] <- 'unused index combination'
  t
}

#' Infer barcode sequences for all samples
#'
#' High-level function calling [infer_barcodes] on every sample,
#' Parallel processing possible with `cores` > 1. Custom parallel processing
#' frameworks can be plugged in by providing `parallel_lapply_fn`.
#'
#' @param seq_tab Sequence metadata table returned by [do_trim_demux] or [do_demux]
#' @param dada_err Result of [dada_learn_errors]
#' @param aln_out Optional output directory for BAM alignments (none saved if `NULL`)
#' @param tmp_dir Optional temporary directory (e.g. [RAM drive](https://en.wikipedia.org/wiki/RAM_drive)
#'    such as `/run/user/1000/DadaNanoBC` on Linux, for faster processing)
#'
#' @returns Input table (`seq_tab`) with a new column `clustering`, which is a list
#' of data frames as returned by [infer_barcodes]
#'
#' @export
do_infer_all_barcodes <- function(seq_tab,
                              dada_err,
                              aln_out = NULL,
                              tmp_dir = NULL,
                              parallel_lapply_fn = NULL,
                              ...,
                              cores = 1) {
  # for parallel::clusterExport
  cluster_export <- c('get_program',
                      'compare_seqs',
                      'write_dna',
                      .infer_barcodes_export)
  idx <- seq_len(nrow(seq_tab))
  stopifnot(!is.null(seq_tab$indexes))
  stopifnot(!is.na(seq_tab$indexes))
  seq_tab$clustering = vector(mode = 'list', nrow(seq_tab))
  names(seq_tab$clustering) <- seq_tab$indexes
  names(idx) <- seq_tab$indexes
  sel_comb <- seq_tab$indexes[!is.na(seq_tab$reads_path)]
  if (length(sel_comb) > 0) {
    idx <- idx[sel_comb]
    batch_size <- max(1, min(12, ceiling(length(idx) / cores / 10)))
    idx_batches <- split(idx, ceiling(1:length(idx) / batch_size))
    process_fn <- function(indexes) {
      lapply(indexes, function(i) {
        # i=which(seq_tab$indexes=='ITS5_bc485-ITS4_bc082')
        fq <- seq_tab$reads_path[i]
        known_seq <- seq_tab$known_sequence[i]
        # cat(fq, '\n')
        idx <- seq_tab$indexes[i]
        if (!is.null(aln_out)) {
          alignment_prefix <- file.path(aln_out, idx, idx)
          dir.create(dirname(alignment_prefix), FALSE, TRUE)
        } else {
          alignment_prefix <- NULL
        }
        d <- infer_barcodes(
          fq,
          dada_err,
          alignment_prefix = alignment_prefix,
          id_prefix = idx,
          tmp_dir = tmp_dir,
          ...
        )
        if (!is.null(d)) {
          bam_out <- if (!is.null(alignment_prefix)) {
            paste0(alignment_prefix, '_seq_comparison.bam')
          }
          d <- compare_seqs(d,
                            bam_out = bam_out,
                            tmp_dir = tmp_dir,
                            known_seq = known_seq)
        }
        d
      })
    }
    res <- if (is.null(parallel_lapply_fn)) {
      process_parallel(idx_batches,
                       process_fn,
                       cores = cores,
                       export = cluster_export)
    } else {
      parallel_lapply_fn(idx_batches, process_fn)
    }
    res <- unlist(unname(res), recursive = FALSE)
    stopifnot(sel_comb == names(res))
    seq_tab$clustering[sel_comb] <- res
    seq_tab <- propagate_data(seq_tab)
  }
  seq_tab
}


#' Combine all BAM alignment files into a single one
#'
#' @param seq_tab Sequence table as returned by [do_infer_all_barcodes] or downstream functions
#' @param aln_dir Alignment directory (as provided to [do_infer_all_barcodes])
#' @param out_prefix Output prefix for '.bam', '.bam.bai' and '.fasta' files
#'
#' @returns Named list of output files
#'
#' @export
do_combine_alignments <- function(seq_tab, aln_dir, out_prefix, top_only = FALSE) {
  # combine mapped reads
  idx_path <- file.path(aln_dir, seq_tab$indexes, seq_tab$indexes)
  sel_i <- which(!is.na(seq_tab$omega_a))
  sel_list <- if (top_only) {
    lapply(sel_i, function(j) {
      d <- seq_tab$clustering[[j]]
      sel <- d$taxon_num == d$taxon_num[1] & !d$is_rare
      if (d$taxon_num[1] != 1) {
        # list was reordered, which suggests contamination:
        # also include top taxon
        max_abund <- ave(d$n_mapped, d$taxon_num, FUN = max)
        sel <- sel | max_abund == max(max_abund) & !d$is_rare
      }
      list(idx_path[j], d$full_id[sel])
    })
  } else {
    idx_path[sel_i]
  }
  samtools <- get_program('samtools')
  outfiles <- subset_combine_bam(out_prefix, sel_list, samtools = samtools)
  # combine sequence comparisons
  sel <- !is.na(seq_tab$has_seq_comparison) &
    seq_tab$has_seq_comparison
  cmp_out <- subset_combine_bam(
    paste0(out_prefix, '_seq_comparison'),
    paste0(idx_path, '_seq_comparison')[sel],
    write_refs = FALSE,
    samtools = samtools
  )
  c(outfiles, cmp_out)
}


#' Propagate information from nested cluster tables to the main seq_tab
propagate_data <- function(seq_tab, extra_seq_cols = NULL) {
  # Propagate sample-level attributes to seq_tab
  attr_cols <- c('n_reads', 'n_singletons', 'omega_a', 'has_seq_comparison')
  for (a in attr_cols) {
    seq_tab[[a]] = sapply(seq_tab$clustering, function(cl)
      attr(cl, a) %||% NA)
  }
  seq_tab$n_reads[is.na(seq_tab$n_reads)] = 0
  seq_tab$singleton_frac = with(seq_tab, n_singletons / n_reads)

  # Propagate extra information from the top taxon to seq_tab
  summary_cols <- list(list(c('n_mapped', 'max_identical', 'n0'), sum), list(c('max_homopoly_len'), max))
  for (sc in summary_cols) {
    for (col in sc[[1]]) {
      seq_tab[[paste0('top_', col)]] = sapply(seq_tab$clustering, function(d)
        if (!is.null(d)) {
          d <- d[d$taxon_num == d$taxon_num[1], ]
          sc[[2]](d[[col]])
        } else
          NA)
    }
  }

  seq_cols <- c(
    'sequence',
    'consensus',
    'consensus_diffs',
    'consensus_ambigs',
    'homopolymer_adjustments',
    'known_seq_diffs',
    'method',
    'message',
    extra_seq_cols
  )
  for (col in seq_cols) {
    seq_tab[[paste0('top_seq_', col)]] = sapply(seq_tab$clustering, function(d) {
      d[[col]][1] %||% NA
    })
  }

  # more information
  seq_tab$top_n_reads = sapply(seq_tab$clustering, function(d) {
    if (is.null(d))
      0
    else
      sum(d$taxon_num == d$taxon_num[1] & !d$is_rare)
  })
  seq_tab$target_cluster_frac = with(seq_tab, top_n_mapped / n_reads)

  seq_tab
}


#' Run the taxonomy assignment and compare with provided taxonomic names of specimens
#'
#' @param seq_tab Sequence table as returned by [do_infer_all_barcodes] or downstream functions
#' @param db_file UTAX-formatted sequence database (optionally compressed);
#' see [load_taxdb]
#' @param gbif_cache_file Path to cache file for saving GBIF taxonomy lookup results
#' (may be shared across datasets) (see [get_gbif_taxa])
#' @param tmp_dir optional temporary directory for saving intermediate sequence files
#' @param confidence_threshold SINTAX bootstrap threshold (default: 0.8 / 80%)
#' @param summary_ranks optional named character vector with ranks that should
#' show up in the lineages (default: determine automatically)
#' @param known_contaminants nested list specifying known contaminant names,
#' such as `list(family=c('Aspergillaceae', ...), genus=c(...))`
#' @param likely_kingdom kingdom name providing some guidance for the GBIF
#' name search (see [get_gbif_taxa])
#' @param contam_rank_delta Require at least N additional ranks to be matching
#' between the provided (e.g. morphology-based) and the auto-assigned
#' sequence-based taxonomic lineages for a taxon to be determined as the
#' "correct" taxon in the mix, and the top taxon being down-ranked
#'
#' @returns
#' `seq_tab`, with the nested data frames in `clustering` **reordered** as such
#' that dominant taxa flagged as contamination are moved to the bottom.
#' The following new columns are added:
#' - *unique_id*: unique sequence id in the form: `u1`, `u2`, etc.
#' - *taxon*: auto-inferred taxon name
#'   (*may not be correct, confirm with BLAST or other identification methods*)
#' - *short_lineage*: taxonomic lineage of the auto-inferred taxon
#'   (provided or auto-inferred `summary_ranks` shown)
#' - *is_contaminant*: (logical) `TRUE` for taxa recognized as contamination
#' - *matching_ranks*: Number of matching taxonomic ranks in the comparison
#'   between the provided and auto-assigned sequence-based taxonomy
#'   (usually: kingdom/domain, phylum, class, order, family, genus, species).
#'   The comparison can only be done up to the highest defined rank in either
#'   lineage.
#' - *mismatching_ranks*: Number of non-matching ranks in the comparison
#'   between the provided and auto-assigned sequence-based taxonomy
#' - *unspecific*: `TRUE` for all taxa that are not at the top in the clusters table
#'   after the contamination ranking
#'
#' The main `seq_tab` table further contains some columns with data from the
#' top sequence: *top_seq_matching_ranks*, *top_seq_mismatching_ranks*,
#' *top_seq_taxon*, *top_seq_short_lineage*
#'
#' @details
#'
#' Sequences are flagged as contamination (*is_contaminant*) if:
#' - there exists a less abundant taxon with the correct kingdom
#'   (compared to *taxon*, GBIF name comparison guided by `likely_kingdom`)
#' - a taxon has has at least `contam_rank_delta` more consistent (matching)
#'   taxonomic ranks in the GBIF lineage than the top taxon
#'   (explanation for *matching_ranks* above); undefined ranks are excluded
#'   from the comparison
#'
#' @export
do_assign_compare_taxonomy <- function(seq_tab,
                                       db_file,
                                       gbif_cache_file,
                                       tmp_dir = NULL,
                                       confidence_threshold = 0.8,
                                       summary_ranks = NULL,
                                       known_contaminants = NULL,
                                       likely_kingdom = NULL,
                                       contam_rank_delta = 3,
                                       cores = 1) {
  tax <- do_assign_taxonomy(
    seq_tab,
    db_file,
    tmp_dir = tmp_dir,
    confidence_threshold = confidence_threshold,
    summary_ranks = summary_ranks,
    cores = cores
  )
  seq_tab <- do_compare_taxonomy(
    tax$seq_tab,
    tax$seq_lineages,
    gbif_cache_file,
    known_contaminants = known_contaminants,
    likely_kingdom = likely_kingdom,
    contam_rank_delta = contam_rank_delta
  )
  attr(seq_tab, 'summary_ranks') <- tax$summary_ranks
  seq_tab
}

#' Auto-assign taxonomic names/lineages
#'
#' @export
do_assign_taxonomy <- function(seq_tab,
                               db_file,
                               tmp_dir = NULL,
                               confidence_threshold = 0.8,
                               summary_ranks = NULL,
                               cores = 1) {
  # De-duplicate sequences
  stopifnot(!is.null(seq_tab$clustering))
  unique_seqs <- sort(unique(unlist(
    lapply(seq_tab$clustering, function(d)
      d$consensus)
  )))
  stopifnot(!is.null(unique_seqs))
  names(unique_seqs) = paste0('u', seq_along(unique_seqs))
  unique_map <- setNames(names(unique_seqs), unique_seqs)

  # assign using SINTAX
  seqs_tmp <- tempfile('seqs', tmpdir = tmp_dir %||% tempdir(), fileext =
                         '.fasta')
  write_dna(unique_seqs, seqs_tmp)
  seq_lineages <- assign_taxonomy_sintax(
    seqs_tmp,
    db_file,
    confidence_threshold = confidence_threshold,
    threads = cores
  )
  seq_lineages <- as.matrix(seq_lineages)
  stopifnot(length(setdiff(
    names(unique_seqs), rownames(seq_lineages)
  )) == 0)
  file.remove(seqs_tmp)

  # summarize higher ranks at appropriate level
  # (for reports)
  if (is.null(summary_ranks)) {
    top_seq <- na.omit(sapply(seq_tab$clustering, function(t)
      t$consensus[1] %||% NA))
    l <- seq_lineages[match(unique_map[top_seq], rownames(seq_lineages)), , drop = FALSE]
    n.names = apply(l, 2, function(x)
      sum(cumsum(prop.table(
        sort(table(x), TRUE)
      )) >= 0.7))
    def.frac = apply(l, 2, function(x)
      mean(!is.na(x)))
    # Determine the last rank with <= 20 taxa and also >= 40% lineages need to be non-NA.
    # All lower ranks with more taxa/more NAs will be ignored
    i <- tail(which(n.names <= 20 & def.frac >= 0.4), 1)
    if (length(i) == 0) {
      i <- 1
    }
    l <- l[, 1:i, drop = FALSE]
    # The last column in 'l' will have the reported name
    # Thinking about plots where rare names will be lumped together in an 'Other' category,
    # we don't want to include ranks that would all have the same name for those
    # top ~ 10 taxa
    abundant.taxa = names(head(sort(table(l[, ncol(l)]), TRUE), 10))
    n.names.top = apply(l[l[, ncol(l)] %in% abundant.taxa, , drop = FALSE], 2, function(x)
      length(unique(na.omit(x))))
    # Select 2nd rank (usually phylum) + max. 2 ranks including the one selected above
    summary_rank_i <- sort(unique(c(
      min(2, length(n.names.top)), tail(which(n.names.top > 1), 2), length(n.names.top)
    )))
    summary_ranks <- colnames(seq_lineages)[summary_rank_i]
  } else {
    summary_ranks <- intersect(summary_ranks, colnames(seq_lineages))
    stopifnot(length(summary_ranks) > 0)
  }

  # taxa names
  seq_lineages <- cbind(name = make_taxon_name(seq_lineages), seq_lineages)

  # add unique seq codes, taxa names and short lineages to seq_tab
  has_seqs <- !sapply(seq_tab$clustering, is.null)
  seq_tab$clustering[has_seqs] = lapply(seq_tab$clustering[has_seqs], function(d) {
    d$unique_id <- unique_map[d$consensus]
    d$taxon = seq_lineages[d$unique_id, 'name']
    d$short_lineage = apply(seq_lineages[d$unique_id, summary_ranks, drop = FALSE], 1, function(x) {
      if (all(is.na(x)))
        NA
      else
        paste(na.omit(x), collapse = ' / ')
    })
    d
  })

  list(
    seq_tab = seq_tab,
    seq_lineages = seq_lineages,
    summary_ranks = summary_ranks
  )
}


#' @export
do_compare_taxonomy <- function(seq_tab,
                                seq_lineages,
                                gbif_cache_file,
                                known_contaminants = NULL,
                                likely_kingdom = NULL,
                                contam_rank_delta = 3) {
  # blacklist contaminants
  has_seqs <- !sapply(seq_tab$clustering, is.null)
  seq_tab$clustering[has_seqs] = lapply(seq_tab$clustering[has_seqs], function(d) {
    stopifnot(!is.null(d$unique_id))
    lineages <- seq_lineages[d$unique_id, , drop = F]
    # assign to contaminant if any member of a taxon has a known contaminant name
    d$is_contaminant = FALSE
    for (rank in setdiff(names(known_contaminants), colnames(lineages))) {
      taxa <- known_contaminants[[rank]]
      as.logical(ave(lineages[, rank] %in% taxa, d$taxon_num, FUN = max))
    }
    d
  })

  # try to find other contamination based on name matching
  if (!('taxon' %in% names(seq_tab))) {
    seq_tab$taxon = NA
  }

  taxon <- unique(na.omit(seq_tab$taxon))
  if (length(taxon) > 0) {
    # get lineages from GBIF
    likely_kingdom <- as.character(likely_kingdom %||% NA)
    seq_lineages_gbif <- get_gbif_taxa(
      seq_lineages,
      cache_file = gbif_cache_file,
      likely_kingdom = likely_kingdom,
      verbose = TRUE
    )
    row.names(seq_lineages_gbif) <- row.names(seq_lineages)
    provided_lineages_gbif <- get_gbif_taxa(
      taxon,
      cache_file = gbif_cache_file,
      likely_kingdom = likely_kingdom,
      verbose = TRUE
    )
    row.names(provided_lineages_gbif) <- taxon

    has_taxon <- has_seqs & !is.na(seq_tab$taxon)
    seq_tab$clustering[has_taxon] <- lapply(which(has_taxon), function(i) {
      # cat(i, '\n')
      d <- seq_tab$clustering[[i]]
      taxon <- seq_tab$taxon[i]
      # sort by taxon and abundance (in case it has been reordered above)
      max_abund <- ave(d$n_mapped, d$taxon_num, FUN = max)
      d <- d[order(-max_abund, d$taxon_num, -d$n_mapped), ]
      # calculate taxa overlap (proportion of matching names in the lineage)
      provided_lineage <- provided_lineages_gbif[taxon, ]
      seq_lineages <- seq_lineages_gbif[d$unique_id, , drop = FALSE]
      lineage_cmp <- t(apply(seq_lineages, 1, function(l)
        l == provided_lineage))
      d$matching_ranks = rowSums(lineage_cmp, na.rm = TRUE)
      d$mismatching_ranks = rowSums(!lineage_cmp, na.rm = TRUE)
      if (all(d$taxon_num == d$taxon_num[1])) {
        return(d)
      }
      # when doing pairwise comparisons of all seqs to the top one, make sure
      # that both lineages have the same length
      # (for the longer lineage, set, values beyond the shorter lineage to NA)
      lineage_cmp <- ifelse(is.na(lineage_cmp[rep(1, nrow(lineage_cmp)), ]) |
                              is.na(lineage_cmp),
                            NA,
                            lineage_cmp)
      rownames(lineage_cmp) <- d$taxon
      matching_ranks_adj <- rowSums(lineage_cmp, na.rm = TRUE)
      # average the mismatch score by taxon (we only move whole clusters to the front)
      matching_ranks_adj <- ave(matching_ranks_adj, d$taxon_num)
      # by how many ranks does the overlap improve if choosing another variant?
      match_delta <- matching_ranks_adj - matching_ranks_adj[1]
      # flag already assigned (known) contaminants
      match_delta[d$is_contaminant] = 0
      # ignore if < 3 ranks improvement over top taxon
      match_delta[match_delta < contam_rank_delta] = 0
      # also look at kingdom mismatches
      if (any(lineage_cmp[, 1], na.rm = TRUE) &
          !all(lineage_cmp[, 1], na.rm = TRUE)) {
        match_delta[!is.na(lineage_cmp[, 1]) & lineage_cmp[, 1]] <- Inf
      }
      if (any(match_delta > 0)) {
        best_i <- which.max(match_delta)
        best_i <- which(d$taxon_num == d$taxon_num[best_i])[1]
        if (any(d$taxon_num != d$taxon_num[1]) && best_i > 1) {
          # flag all groups up to the best one
          d$is_contaminant[1:(best_i - 1)] = TRUE
        }
      }
      d
    })
  }

  # move contaminant taxa down to be located *after* the most abundant
  # non-contaminant taxon
  seq_tab$has_contamination = sapply(seq_tab$clustering, function(d)
    isTRUE(d$is_contaminant[1]))
  seq_tab$clustering[has_seqs] = lapply(which(has_seqs), function(i) {
    # cat(i, ' ')
    d <- seq_tab$clustering[[i]]
    # also sort by taxon and abundance (in case it has been reordered above)
    taxon_abund <- ave(ifelse(d$is_rare, 0, d$n_mapped), d$taxon_num, FUN =
                         sum)
    d <- d[order(d$is_contaminant, -taxon_abund, d$taxon_num, -d$n_mapped), ]

    # add 'unspecific'
    d$unspecific = d$taxon_num != d$taxon_num[1] | d$is_rare
    d
  })

  # check group over-abundance and suspiciously similar taxonomy
  seq_tab$clustering[has_seqs] = lapply(seq_tab$clustering[has_seqs], function(d) {
    # check group1 over-abundance
    grp <- unique(d$taxon_num)
    if (length(grp) > 1) {
      grp_abund <- ave(d$n_mapped, d$taxon_num, FUN = sum)
      is_grp1 <- d$taxon_num == grp[1]
      abund_ratio <- grp_abund[is_grp1][1] / grp_abund
      grp2_ratio <- abund_ratio[d$taxon_num == grp[2]][1]
      if (grp2_ratio < 3 && grp2_ratio >= 1) {
        d$message[is_grp1] <- paste(d$message[is_grp1], '< 3 times overabundant;')
      }
      # suspiciously similar within 8x abundance range?
      lineages <- seq_lineages[d$unique_id, , drop = F]
      lineage_cmp <- t(apply(lineages, 1, function(l)
        l == lineages[1, ]))
      rank_match <- ave(rowSums(lineage_cmp, na.rm = TRUE), d$taxon_num)
      rank_match_delta <- rank_match[1] - rank_match
      rank_mismatch <- ave(rowSums(!lineage_cmp, na.rm = TRUE), d$taxon_num)
      if (any(!is_grp1 &
              abund_ratio < 8 &
              rank_mismatch < 2 & rank_match_delta < 2)) {
        d$message[is_grp1] <- paste(d$message[is_grp1], 'Related taxa in mix;')
      }
    }
    d
  })

  propagate_data(
    seq_tab,
    extra_seq_cols = c(
      'matching_ranks',
      'mismatching_ranks',
      'taxon',
      'short_lineage'
    )
  )
}

process_parallel <- function(data,
                             func,
                             cores = 1,
                             export = NULL) {
  if (cores > 1) {
    cl <- parallel::makeCluster(cores)
    if (!is.null(export))
      parallel::clusterExport(cl = cl, export)
    tryCatch(
      parallel::parLapply(cl = cl, data, func),
      finally = parallel::stopCluster(cl)
    )
  } else {
    lapply(data, func)
  }
}
