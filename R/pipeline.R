
setup_pipeline <- function() {
  p <- installed.packages()[,c('Package')]
  p <- setdiff(c('openxlsx2', 'stringr', 'pbapply'), p)
  if (length(p) > 0) {
    install.packages(p)
  }
  install_cluster_deps()
  install_sequtil_deps()
  install_taxonomy_deps()
  require_bash()
}

read_xlsx_sample_tab <- function(meta_file, sheet_name) {
  sample_tab <- openxlsx2::read_xlsx(
    meta_file,
    sheet_name,
    na.strings = c('#N/A', '')
  )
  sample_tab <- as_tibble(sample_tab, .name_repair='unique_quiet') %>% 
    # empty cells might still have '0' (if a formula), we set these to NA
    mutate(sample = ifelse(!is.na(sample) & sample == '0', NA, sample))
  parse_sample_tab(sample_tab)
}

parse_sample_tab <- function(sample_tab) {
  # normalize headers
  names(sample_tab) <- gsub('[- ]+', '_', names(sample_tab))
  
  # check for invalid characters
  new_sample <- gsub('\\s', '_', iconv(sample_tab$sample, 'utf-8', 'ascii'), perl=T)
  sel <- !is.na(new_sample) & new_sample != sample_tab$sample
  if (any(sel)) {
    warning(paste0('The following samples had spaces or non-ASCII characters and were adjusted: ',
                   paste(new_sample[sel], collapse=', ')))
    sample_tab$sample = new_sample
  }
  
  # remove barcodes without sample name
  sample_tab <- sample_tab[!is.na(sample_tab$sample),]
  
  # Check for duplicates
  
  sample_tab$unique_sample = sample_tab$sample
  sample_tab <- sample_tab %>% 
    relocate(unique_sample, .after='sample')
  dup_s <- unique(sample_tab$unique_sample[duplicated(sample_tab$unique_sample)])
  if (length(dup_s[dup_s != '0']) > 0) {
    warning(paste0('The following samples are duplicated: ',
                   paste(dup_s, collapse=', '),
                   '. We simply append the plate/coordinate to the name.'))
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

read_xlsx_primer_tab <- function(meta_file, sheet_name, ...) {
  primer_tab <- openxlsx2::read_xlsx(meta_file, sheet_name, na.strings = c(''))
  primer_tab$primer_name <- gsub('_.*', '', primer_tab$primer_barcode)
  stopifnot(!duplicated(primer_tab$primer_barcode))
  parse_primer_tab(primer_tab, ...)
}

parse_primer_tab <- function(primer_tab, amplicons = NULL) {
  
  if (is.null(amplicons)) {
    stopifnot(length(primers) == 2)
    amplicons <- unique(primers$primer_name)
  }
  
  primers <- split(primer_tab, primer_tab$primer_name)
  
  amp_primer_names <- strsplit(amplicons, '_', fixed=TRUE)
  names(amp_primer_names) <- amplicons
  
  primers_only <- setdiff(names(primers), unlist(amp_primer_names))
  amp_only <- setdiff(unlist(amp_primer_names), names(primers))
  if (length(primers_only) > 0) {
    stop("Some primer name(s) found only in the 'primer' sheet, but not ",
         "in the 'sample_list' sheet: ",
         paste(primers_only, collapse=', '))
  }
  if (length(amp_only) > 0) {
    stop("Some primer name(s) found only in the 'sample_list' sheet, but not ",
         "in the 'primer' sheet: ",
         paste(amp_only, collapse=', '))
  }
  
  # collect per-amplicon primer/index information
  amplicon_primers <- lapply(amp_primer_names, function(p) {
    if (length(p) != 2) {
      stop("The 'amplicon' column in 'sample_list' needs to be in the form 'forward_reverse' ",
           "where the forward/reverse primer name has to match the primer names in ",
           "the 'primer_barcode' column of the 'primers' sheet.")
    }
    # TODO: direction not actually needed
    comb <- list(forward=primers[[p[1]]], reverse=primers[[p[2]]])
    lapply(comb, function(d) {
      primer <- unique(d$primer)
      if (length(unique(d$primer)) != 1) {
        stop(sprintf("The primer '%s' does not have the same sequence in all rows of the 'primer' column",
                     d$primer_name[1]))
      }
      stopifnot(!duplicated(d$barcode))
      idx_len <- unique(nchar(d$barcode))
      if (length(idx_len) != 1) {
        # TODO: in theory we could support this (just affects how we trim the sequences before the primer)
        stop(sprintf("The barcodes for primer '%s' should all have the same length", primer))
      }
      list(
        primer = setNames(primer, d$primer_name[1]),
        barcode = setNames(d$barcode, d$primer_barcode),
        barcode_len = idx_len
      )
    })
  })
}

run_primer_search <- function(fq_paths,
                              amplicon_primers,
                              search_dir,
                              primer_max_err = 2.5,
                              idx_max_diffs = 0,
                              min_barcode_length = 50,
                              overwrite = FALSE,
                              ...) {
  amplicons <- names(amplicon_primers)
  for (amplicon in amplicons) {
    amp_search_dir <- file.path(search_dir, amplicon)
    demux_fq <- file.path(amp_search_dir, 'trimmed.fastq.zst')
    
    if (overwrite || !file.exists(demux_fq) && length(fq_paths) > 0) {
      cat(amplicon, sep='\n', file=stderr())
      unlink(amp_search_dir)  # clean up old files
      dir.create(amp_search_dir, FALSE, TRUE)
      amp_pr <- amplicon_primers[[amplicon]]
      seq_prefix <- file.path(amp_search_dir, '')
      primer_paths <- paste0(seq_prefix, c('fwd', 'rev'), '_primers.fasta')
      idx_paths <- paste0(seq_prefix, c('fwd', 'rev'), '_idx.fasta')
      for (i in 1:2) {
        write_dna(amp_pr[[i]]$primer, primer_paths[i])
        write_dna(amp_pr[[i]]$barcode, idx_paths[i])
      }
      search_opts <- c(
        '-p', primer_max_err,
        '-b', idx_max_diffs,
        '-l', min_barcode_length,
        '-o', amp_search_dir,
        '-c', 'zst',
        '-t', cores,
        '-s', .programs$seqtool
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
  # TODO: move these files to better location?
  
  # primer start position
  pos <- do.call(rbind, lapply(amplicons, function(amplicon) {
    d <- read.delim(file.path(search_dir, amplicon, 'primer_pos.tsv'),
                  colClasses=c('character', 'integer', 'integer'))
    d$amplicon <- amplicon
    d$dir <- factor(d$dir, c('fwd', 'rev'), c('forward', 'reverse'))
    aggregate(count ~ amplicon + dir + pos, data=d, sum)
  }))
  
  # amplicon length
  category_trans <- c(trimmed='regular', concatenated='concatenated products')
  l <- do.call(rbind, lapply(names(amplicon_primers), function(amplicon) {
    d <- read.delim(file.path(search_dir, amplicon, 'length_stats.tsv'),
                    colClasses=c('character', 'integer', 'integer'))
    d$amplicon <- amplicon
    d
  }))
  l$category <- trimws(gsub('\\.fastq\\.zst$', '', l$file))
  l$category = factor(category_trans[l$category], category_trans)
  l$file <- NULL

  # primer/index counts
  n <- do.call(rbind, lapply(seq_along(amplicons), function(i) {
    d <- read.delim(file.path(search_dir, amplicons[i], 'trim_counts.tsv'),
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
    d <- read.delim(file.path(search_dir, amplicons[i], 'qual_stats.tsv'),
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

  list(position = pos, amplicon_len = l, counts = n, quality = q)
}


run_demux <- function(primer_search_fq, 
                      outdir, 
                      sample_tab,
                      overwrite_threshold = 0.1,
                      force_rerun = FALSE,
                      error_threshold = 2.5, 
                      min_barcode_length = 50,
                      ...) {
  
  names(sample_tab) <- tolower(names(sample_tab))
  stopifnot(c('amplicon', 'indexes', 'sample', 'sample_type', 'morpho_taxon', 'known_sequence') %in% names(sample_tab))
  
  list_demux_files <- function(demux_dir) {
    file.path(demux_dir, list.files(demux_dir, pattern='.fastq.gz$'))  
  }
  
  trimmed_fq <- list_demux_files(outdir)
  
  # run if <10% of sample files found
  # (cannot assume 100% as some index combinations may have no reads)
  if (force_rerun || length(trimmed_fq) < overwrite_threshold * nrow(sample_tab)) {
    run_bash(c('scripts/filter-split.sh',
               error_threshold,
               # **note**: currently, we have only one threshold both in 'find-primers' and 'demux'
               min_barcode_length,
               outdir,
               .programs$seqtool,
               primer_search_fq))
    trimmed_fq <- list_demux_files(outdir)
  }
  
  # initialize the sequence table
  sample_tab$i_ <- seq_len(nrow(sample_tab))
  t <- data.frame(
    reads_path = sort(trimmed_fq),
    primer_indexes = gsub('\\.fastq.gz$', '', basename(trimmed_fq)),
    i_ = nrow(sample_tab) + seq_len(length(trimmed_fq))
  )
  t <- cbind(t, stringr::str_match(
    t$primer_indexes, '(?<fpr>[^_]+)_(?<fidx>[^-]+)-(?<rpr>[^_]+)_(?<ridx>[^-]+)'
  ))
  t$amplicon <- paste(t$fpr, t$rpr, sep='_')
  t$indexes = paste(t$fidx, t$ridx, sep='-')
  t <- merge(
    sample_tab,
    t[c('reads_path', 'primer_indexes', 'amplicon', 'indexes', 'i_')],
    by = c('amplicon', 'indexes'),
    all = TRUE
  )
  t <- t[order(t$i_.x, t$i_.y),]
  t$i_.x <- t$i_.y <- NULL
  p <- stringr::str_split_fixed(t$amplicon, '_', 2)
  i <- stringr::str_split_fixed(t$indexes, '-', 2)
  t$primer_indexes <- sprintf('%s_%s-%s_%s', p[,1], i[,1], p[,2], i[,2])
  t$category <- ifelse(is.na(t$sample_type), 'sample with reads', as.character(t$`sample type`))
  t$category[is.na(t$reads_path)] <- 'no reads'
  t$category[is.na(t$sample)] <- 'unused index combination'
  t
}

# run clustering
cluster_fq <- function(fq, known_seq, aln_out, opts = NULL, tmp_dir = NULL) {
  # cat(fq, '\n')
  out_prefix <- get_aln_prefix(fq, aln_out)
  args <- modifyList(opts %||% list(), 
                    list(fq = fq, out_prefix = out_prefix,
                         tmp_dir = tmp_dir,
                         minimap2 = .programs$minimap2,
                         samtools = .programs$samtools))
  d <- do.call(get_barcodes, args)
  if (!is.null(d)) {
    align_top(d, out_prefix, known_seq = known_seq)
  }
  d
}

get_aln_prefix <- function(fq, aln_out) {
  out_name <- gsub('\\.[^\\.]+', '', basename(fq))
  out_dir <- file.path(aln_out, out_name)
  dir.create(out_dir, FALSE, TRUE)
  file.path(out_dir, out_name)
}

run_learn_errors <- function(seq_tab, 
                             cache_file,
                             ...,
                             overwrite = FALSE) {
  run_or_load(cache_file, function() {
    dada_learn_errors(sample(na.omit(seq_tab$reads_path)), ...)
  }, rerun = overwrite)
}


get_amplicon_opts <- function(opts, amplicons) {
  opts <- opts %||% list()
  out <- opts[amplicons]
  names(out) <- amplicons
  opts[amplicons] <- NULL
  for (amplicon in amplicons) {
    out[[amplicon]] <- modifyList(opts, out[[amplicon]] %||% list())
  }
  out
}


run_clustering <- function(seq_tab,
                           dada_err,
                           cache_file,
                           aln_out, 
                           opts = NULL,
                           tmp_dir = NULL,
                           cores = 1,
                           overwrite = FALSE) {
  # for parallel::clusterExport
  cluster_export <- c('cluster_fq', 'opts', 'tmp_dir',
                     '.programs', 'dada_err', 'aln_out', 'get_aln_prefix',
                     'seq_tab', get_barcodes_export)
  idx <- seq_len(nrow(seq_tab))
  stopifnot(!is.null(seq_tab$primer_indexes))
  stopifnot(!is.na(seq_tab$primer_indexes))
  seq_tab$clustering = vector(mode='list', nrow(seq_tab))
  names(seq_tab$clustering) <- seq_tab$primer_indexes
  names(idx) <- seq_tab$primer_indexes
  sel_comb <- seq_tab$primer_indexes[!is.na(seq_tab$reads_path)]
  if (length(sel_comb) > 0) {
    res <- run_or_load(cache_file, function() {
      parallel_lapply(idx[sel_comb], function(i) {
        # i=which(seq_tab$primer_indexes=='...')
        fq <- seq_tab$reads_path[i]
        known_seq <- seq_tab$known_sequence[i]
        cluster_fq(fq, known_seq, aln_out, opts = opts, tmp_dir = tmp_dir)
      }, cores=cores, export=cluster_export)
    }, rerun = overwrite)
    if (length(intersect(sel_comb, names(res))) != length(sel_comb)) {
      stop(sprintf('It seems that %s is outdated, please delete', cache_file))
    }
    seq_tab$clustering[sel_comb] <- res[sel_comb]
    seq_tab <- propagate_data(seq_tab)
  }
  seq_tab
}


#' Propagate information from nested cluster tables to the main seq_tab
propagate_data <- function(seq_tab, haplo_extra_cols=NULL) {
  # Propagate sample-level attributes to seq_tab
  attr_cols <- c('n_reads', 'n_singletons', 'omega_a')
  for (a in attr_cols) {
    seq_tab[[a]] = sapply(seq_tab$clustering, function(cl) attr(cl, a) %||% NA)
  }
  seq_tab$n_reads[is.na(seq_tab$n_reads)] = 0
  seq_tab$singleton_frac = with(seq_tab, n_singletons / n_reads)
  
  # Propagate extra information from the top taxon to seq_tab
  summary_cols <- list(
    list(c('abundance', 'max_identical', 'n0'), sum),
    list(c('max_homopoly_len'), max)
  )
  for (sc in summary_cols) {
    for (col in sc[[1]]) {
      seq_tab[[paste0('top_', col)]] = sapply(seq_tab$clustering, function(d) 
        if (!is.null(d)) {
          d <- filter(d, taxon_num == taxon_num[1])
          sc[[2]](d[[col]])
        } else NA
      )
    }
  }
  
  haplo_cols <- c('consensus_diffs', 'consensus_ambigs',
                 'sequence', 'consensus',
                 'method',
                 haplo_extra_cols)
  for (col in haplo_cols) {
    seq_tab[[paste0('top_seq_', col)]] = sapply(seq_tab$clustering, function(d) {
      d[[col]][1] %||% NA
    })
  }
  
  # more information
  seq_tab$top_n_reads = sapply(seq_tab$clustering, function(d) {
    if (is.null(d)) 0 else sum(d$taxon_num == d$taxon_num[1] & !d$is_rare)
  })
  seq_tab$target_cluster_frac = with(seq_tab, top_abundance / n_reads)
  
  seq_tab
}

load_taxdb <- function(db_dir, db_type, db_urls, db_tax_urls=NULL, ...) {
  db_hash <- tools::md5sum(bytes=charToRaw(paste(sort(db_urls), collapse='')))
  taxdb_file <- file.path(db_dir, paste0(db_hash, '.fasta.gz'))
  if (!file.exists(taxdb_file)) {
    message('Downloading taxonomy database: ', paste(db_urls, collapse=', '))
    stopifnot(db_type %in% c('unite', 'qiime_fasta', 'qiime_qza'))
    load_fn <- get(paste0('load_', db_type))
    load_fn(db_urls, taxdb_file, tax_urls=db_tax_urls)
  }
  taxdb_file
}


run_taxonomy_assignment <- function(seq_tab,
                                    prefix,
                                    taxdb_file,
                                    gbif_taxdb_dir = 'taxdb',
                                    known_contaminants = NULL,
                                    confidence_threshold = 0.8,
                                    kingdom = NULL,
                                    contam_rank_delta = 3,
                                    summary_ranks = NULL,
                                    ...) {
  # Since there can be duplicates, we dereplicate before doing further analyses.
  unique_seqs <- sort(unique(unlist(
    lapply(seq_tab$clustering, function(d)
      d$consensus)
  )))
  names(unique_seqs) = paste0('u', seq_along(unique_seqs))
  unique_map <- setNames(names(unique_seqs), unique_seqs)
  unique_seqs_file <- paste0(prefix, '_unique_seqs.fasta')
  seqs_updated <- !file.exists(unique_seqs_file) ||
    !isTRUE(all.equal(as.character(read_dna(unique_seqs_file))[names(unique_seqs)], unique_seqs))
  if (seqs_updated) {
    write_dna(unique_seqs, unique_seqs_file)
  }
  
  # assign using SINTAX
  auto_lineages <- run_or_read(paste0(prefix, '_tab.tsv'), function() {
    assign_taxonomy(unique_seqs,
                    taxdb_file,
                    threshold = confidence_threshold,
                    vsearch = .programs$vsearch,
                    threads = cores)
  }, rerun = seqs_updated)
  auto_lineages <- as.matrix(auto_lineages)
  stopifnot(length(setdiff(names(unique_seqs), rownames(auto_lineages))) == 0)
  
  # summarize higher ranks at appropriate level
  if (is.null(summary_ranks)) {
    top_seq <- na.omit(sapply(seq_tab$clustering, function(t) t$consensus[1] %||% NA))
    l <- auto_lineages[match(unique_map[top_seq], rownames(auto_lineages)), , drop = FALSE]
    n.names = apply(l, 2, function(x) sum(cumsum(prop.table(sort(table(x), TRUE))) >= 0.7))
    def.frac = apply(l, 2, function(x) mean(!is.na(x)))
    i <- tail(which(n.names <= 20 & def.frac >= 0.4), 1)
    if (length(i) == 0) {
      i <- 1
    }
    l <- l[, 1:i, drop = FALSE]
    abundant.taxa = names(head(sort(table(l[, ncol(l)]), T), 10))
    n.unique = apply(l[l[, ncol(l)] %in% abundant.taxa, , drop=FALSE],
                     2,
                     function(x) length(unique(na.omit(x))))
    summary_ranks <- colnames(auto_lineages)[sort(unique(c(2, tail(which(n.unique > 1), 2))))]
    saveRDS(summary_ranks, paste0(prefix, '_summary_ranks.rds'))
  }

  # create taxa names
  for (rank in c('genus', 'species')) {
    if (!(rank %in% colnames(auto_lineages))) {
      auto_lineages <- cbind(auto_lineages, NA)
      colnames(auto_lineages)[ncol(auto_lineages)] <- rank
    }
  }
  highest_taxon <- apply(auto_lineages, 1, function(l) tail(na.omit(l), 1)[1])
  taxname <- ifelse(is.na(highest_taxon), 'Unknown', paste('Unknown', highest_taxon))
  has_spec <- !is.na(auto_lineages[, 'species'])
  gen_only <- !is.na(auto_lineages[, 'genus']) & !has_spec
  taxname[has_spec] <- auto_lineages[has_spec, 'species']
  taxname[gen_only] <- paste(auto_lineages[gen_only & !has_spec, 'genus'], 'sp.')
  auto_lineages <- cbind(name = taxname, auto_lineages)

  # taxa names and short lineages add to seq_tab
  has_seqs <- seq_tab$n_reads > 0
  seq_tab$clustering[has_seqs] = lapply(seq_tab$clustering[has_seqs], function(d) {
    unique_id <- unique_map[d$consensus]
    d$taxon = auto_lineages[unique_id, 'name']
    d$short_lineage = apply(
      auto_lineages[unique_id, summary_ranks, drop = FALSE],
      1,
      function(x) {
        if (all(is.na(x))) NA else paste(na.omit(x), collapse = ' / ')
      }
    )
    d
  })
  
  # blacklist contaminants
  seq_tab$clustering[has_seqs] = lapply(seq_tab$clustering[has_seqs], function(d) {
    unique_id <- unique_map[d$consensus]
    stopifnot(!is.na(unique_id))
    lineages <- auto_lineages[unique_id, , drop = F]
    # assign to contaminant if any member of a taxon has a known contaminant genus
    d$is_contaminant = FALSE
    for (rank in setdiff(names(known_contaminants), colnames(lineages))) {
      taxa <- known_contaminants[[rank]]
      as.logical(ave(lineages[, rank] %in% taxa, d$taxon_num, FUN = max))
    }
    d
  })
  
  # get lineages from GBIF (used below)
  gbif_db <- file.path(gbif_taxdb_dir, sprintf('gbif_tax_%s.tsv', kingdom %||% ''))
  auto_lineages_gbif <- gbif_taxa(
    auto_lineages,
    cache_file = gbif_db,
    likely_kingdom = kingdom,
    verbose = TRUE
  )
  row.names(auto_lineages_gbif) <- row.names(auto_lineages)
  
  # try to find other contamination based on name matching
  if (!('morpho_taxon' %in% names(seq_tab))) {
    seq_tab$morpho_taxon = NA
  }
  
  morpho_taxon <- unique(na.omit(seq_tab$morpho_taxon))
  if (length(morpho_taxon) > 0) {
    # match against GBIF
    kingdom <- kingdom %||% NA
    morpho_lineages_gbif <- gbif_taxa(
      morpho_taxon,
      cache_file = gbif_db,
      likely_kingdom = kingdom,
      verbose = TRUE
    )
    row.names(morpho_lineages_gbif) <- morpho_taxon
    
    has_morphospec <- has_seqs & !is.na(seq_tab$morpho_taxon)
    seq_tab$clustering[has_morphospec] <- lapply(which(has_morphospec), function(i) {
      # cat(i, '\n')
      d <- seq_tab$clustering[[i]]
      taxon <- seq_tab$morpho_taxon[i]
      # sort by taxon and abundance (in case it has been reordered above)
      max_abund <- ave(d$abundance, d$taxon_num, FUN = max)
      d <- d[order(-max_abund, d$taxon_num, -d$abundance), ]
      # calculate taxa overlap (proportion of matching names in the lineage)
      morpho_lineage <- morpho_lineages_gbif[taxon, ]
      unique_id <- unique_map[d$consensus]
      stopifnot(!is.na(unique_id))
      seq_lineages <- auto_lineages_gbif[unique_id,, drop = FALSE]
      lineage_cmp <- t(apply(seq_lineages, 1, function(l) l == morpho_lineage))
      d$matching_ranks = rowSums(lineage_cmp, na.rm = TRUE)
      d$mismatching_ranks = rowSums(!lineage_cmp, na.rm = TRUE)
      if (all(d$taxon_num == d$taxon_num[1])) {
        return(d)
      }
      # when doing pairwise comparisons of all seqs to the top one, make sure
      # that both lineages have the same length
      # (for the longer lineage, set, values beyond the shorter lineage to NA)
      lineage_cmp <- ifelse(
        is.na(lineage_cmp[rep(1, nrow(lineage_cmp)),]) | is.na(lineage_cmp), 
        NA, 
        lineage_cmp
      )
      rownames(lineage_cmp) <- d$taxon
      matching_ranks_adj <- rowSums(lineage_cmp, na.rm=TRUE)
      # average the mismatch score by taxon (we only move whole clusters to the front)
      matching_ranks_adj <- ave(matching_ranks_adj, d$taxon_num)
      # by how many ranks does the overlap improve if choosing another variant?
      match_delta <- matching_ranks_adj - matching_ranks_adj[1]
      # flag already assigned (known) contaminants
      match_delta[d$is_contaminant] = 0
      # ignore if < 3 ranks improvement over top taxon
      match_delta[match_delta < contam_rank_delta] = 0
      # also look at kingdom mismatches
      if (any(lineage_cmp[,1], na.rm=TRUE) & !all(lineage_cmp[,1], na.rm=TRUE)) {
        match_delta[!is.na(lineage_cmp[,1]) & lineage_cmp[,1]] <- Inf
      }
      if (any(match_delta > 0)) {
        best_i <- which.max(match_delta)
        best_i <- which(d$taxon_num == d$taxon_num[best_i])[1]
        if (any(d$taxon_num != d$taxon_num[1]) && best_i > 1) {
          # flag all groups up to the best one
          d$is_contaminant[1:(best_i-1)] = TRUE
        }
      }
      d
    })
  }
  
  # move contaminant compound clusters down to be located *after*
  # the most abundant non-contaminant taxon,
  # and after this, do some more checks
  seq_tab$has_contamination = sapply(seq_tab$clustering, function(d) isTRUE(d$is_contaminant[1]))
  seq_tab$clustering[has_seqs] = lapply(which(has_seqs), function(i) {
    # cat(i, ' ')
    d <- seq_tab$clustering[[i]]
    # also sort by taxon and abundance (in case it has been reordered above)
    taxon_abund <- ave(ifelse(d$is_rare, 0, d$abundance), d$taxon_num, FUN=sum)
    d <- d[order(d$is_contaminant, -taxon_abund, d$taxon_num, -d$abundance), ]

    # add 'unspecific'
    d$unspecific = d$taxon_num != d$taxon_num[1] | d$is_rare
    
    # check group1 over-abundance
    grp <- unique(d$taxon_num)
    if (length(grp) > 1) {
      grp_abund <- ave(d$abundance, d$taxon_num, FUN=sum)
      is_grp1 <- d$taxon_num == grp[1]
      abund_ratio <- grp_abund[is_grp1][1] / grp_abund
      if (abund_ratio[d$taxon_num == grp[2]][1] < 3) {
        d$message[is_grp1] <- paste(
          d$message[is_grp1],
          '< 3 times overabundant;'
        )
      }
      # also check if taxonomy is suspiciously similar
      # within 8x abundance range
      unique_id <- unique_map[d$consensus]
      seq_lineages <- auto_lineages_gbif[unique_id,, drop = FALSE]
      lineage_cmp <- t(apply(seq_lineages, 1, function(l) l == seq_lineages[1,]))
      rank_match <- ave(rowSums(lineage_cmp, na.rm = TRUE), d$taxon_num)
      rank_match_delta <- rank_match[1] - rank_match
      rank_mismatch <- ave(rowSums(!lineage_cmp, na.rm = TRUE), d$taxon_num)
      if (any(!is_grp1 & abund_ratio < 8 & rank_mismatch < 2 & rank_match_delta < 2)) {
        d$message[is_grp1] <- paste(
          d$message[is_grp1],
          'Related taxa in mix;'
        )
      }
    }
    d
  })
  seq_tab
}

propagate_tax_data <- function(seq_tab) {
  propagate_data(
    seq_tab,
    haplo_extra_cols <- c(
      'matching_ranks',
      'mismatching_ranks',
      'taxon',
      'short_lineage',
      'message'
    )
  )  
}

compare_known_seqs <- function(seq_tab, aln_out, cores = 1) {
  if (!('known sequence' %in% names(seq_tab))) {
    seq_tab$known_sequence = NA
  }
  
  has_known_seqs <- seq_tab$n_reads > 0 & !is.na(seq_tab$known_sequence)
  if (any(has_known_seqs)) {
    seq_tab$clustering[has_known_seqs] = pbapply::pblapply(which(has_known_seqs), function(i) {
      d <- seq_tab$clustering[[i]]
      d$known_seq_diffs <- diffs <- match_both_orient(
        d$consensus,
        rep(seq_tab$known_sequence[i], nrow(d)),
        cores=cores
      )
      d$known_seq_dir <- as.integer(names(diffs))
      d
    })
  }
  
  seq_tab$known_seq_diffs = sapply(seq_len(nrow(seq_tab)), function(i) {
    d <- seq_tab$clustering[[i]]
    if (!is.null(d)) {
      sel <- !d$is_contaminant & !d$is_rare
      sd <- d$known_seq_diffs[sel]
      if (length(sd) > 0) {
        best_i <- which.min(sd)
        # set correct orientation if necessary
        if (d$known_seq_dir[best_i] == 2) {
          seq_tab$known_sequence[i] <- rev_complement(seq_tab$known_sequence[i])
        }
        if (d$taxon_num[sel][best_i] == d$taxon_num[1]) {
          if (sd[best_i] > 0) {
            out_prefix <- get_aln_prefix(seq_tab$reads_path[i], aln_out)
            align_top(d, out_prefix, known_seq = seq_tab$known_sequence[i])
          }
        } else {
          is_g1 <- d$taxon_num == d$taxon_num[1]
          d$message[is_g1] <- paste(d$message[is_g1], 'Known seq. matches other taxon in mix;')
        }
        return(sd[best_i])
      }
    }
    NA
  })
  
  seq_tab
}

#' Check installed programs
get_programs <- function(bin_dir='bin') {
  required_programs <- c(seqtool='st', samtools='samtools', minimap2='minimap2', vsearch='vsearch')
  out <- lapply(names(required_programs), function(program) {
    bin_name <- required_programs[program]
    if (bash_cmd_succeeds(c(bin_name, '--version'))) {
      bin_name
    } else if (bash_cmd_succeeds(c(file.path(bin_dir, bin_name), '--version'))) {
      file.path(bin_dir, bin_name)
    } else {
      stop(
        sprintf("The program '%s' was not found. ", bin_name),
        "It needs to be installed either system-wide and visible to R (in $PATH), ",
        "or installed in the 'bin' subdirectory. ",
        "Consider installing by runing 'scripts/install-deps.sh' in the (Bash) console. ",
        "For more information, see..."
        # TODO: reference tutorial
      )
    }
  })
  names(out) = names(required_programs)
  out
}

require_bash <- function() {
  if (Sys.which('bash') == '') {
    stop('The Bash shell is missing. On Windows, WSL2 is needed (https://learn.microsoft.com/en-us/windows/wsl/install)')
  }
}

#' Run a command in Bash, which should work on Linux, OS-X and Windows WSL
#' **Note**: assuming **no** spaces in arguments
run_bash <- function(cmd, ...) {
  cmd <- paste(cmd, collapse=' ')
  system2('bash', c('-c', shQuote(cmd)), ...)
}

bash_cmd_succeeds <- function(cmd) {
  tryCatch(is.null(suppressWarnings(attr(run_bash(cmd, stderr=TRUE), 'status'))),
           error=function(e) FALSE)
}


parallel_lapply <- function(data, func, cores=1, export=NULL) {
  if (cores > 1) {
    cl <- parallel::makeCluster(cores)
    if (!is.null(export))
      parallel::clusterExport(cl=cl, export)
  } else {
    cl <- NULL
  }
  tryCatch(pbapply::pblapply(data, func, cl=cl),
           finally=if (!is.null(cl)) parallel::stopCluster(cl))
}

run_or_load <- function(rds_file, run_fn, rerun=FALSE) {
  if (!file.exists(rds_file) || rerun) {
    res <- run_fn()
    saveRDS(res, rds_file)
    res
  } else {
    readRDS(rds_file)
  }
}

run_or_read <- function(tsv_file, run_fn, rerun=FALSE) {
  if (!file.exists(tsv_file) || rerun) {
    res <- run_fn()
    write.table(res, tsv_file, sep='\t', quote=F, na='')
    res
  } else {
    read.delim(tsv_file, na.strings='')
  }
}





#### Setup globals ############################################

.programs <- get_programs()

###############################################################

