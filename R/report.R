

create_report <- function(seq_tab,
                          outfile,
                          low_abund_threshold = 20,
                          n_curate = 4,
                          max_seqs = 6,
                          ...) {
  
  seq_tab_def <- seq_tab[!is.na(seq_tab$sample),]

  # TODO: should issues only be reported for most abundant seq or for all?
  #   (right now: for all)
  seq_tab_def$issues = sapply(1:nrow(seq_tab_def), function(i) {
    # cat(i, ' ')
    d <- seq_tab_def$clustering[[i]]
    d.flt = if (!is.null(d) && nrow(d) > 0) d[!d$unspecific,] else NULL
    out <- c()
    if (!is.null(d.flt) && nrow(d.flt) > 0) {
      if (nrow(d.flt) > 4) {
        out <- c(out, 'many-variants')
      }
      contam <- d$abundance[1] != max(d$abundance)
      if (contam) {
        out <- c(out, 'contamination')
      }
      if (!is.null(d.flt$known_seq_diffs) && all(!is.na(d.flt$known_seq_diffs)) &&
          any(d.flt$known_seq_diffs > 0)) {
        # TODO: threshold?
        if (contam && any(d$known_seq_diffs[d$is_contaminant] < 2)) {
          out <- c(out, 'known-seq-contamination')
        } else {
          out <- c(out, 'known-seq-diffs')
        }
      }
      if (!is.na(d.flt$consensus_diffs[1]) && d.flt$consensus_diffs[1] > 0) {
        out <- c(out, 'consensus-diffs')
      } else if (all(!is.na(d.flt$consensus_diffs)) && any(d.flt$consensus_diffs > 0)) {
        out <- c(out, '[consensus-diffs]')
      }
      if (d.flt$consensus_ambigs[1] > 0) {
        out <- c(out, 'ambig-consensus')
      } else if (any(d.flt$consensus_ambigs > 0)) {
        out <- c(out, '[ambig-consensus]')
      }
      # TODO: -> ?
      if (d.flt$abundance[1] < low_abund_threshold) {
        out <- c(out, 'low-coverage')
      } 
      if (!is.null(d.flt$mismatching_ranks)) {
        max_mm <- max(d.flt$mismatching_ranks)
        if (max_mm >= 3) {
          out <- c(out, 'tax-mismatch')
        } else if (max_mm == 2) {
          out <- c(out, '[tax-mismatch]')
        }
      }
      if (!is.null(d.flt$message) && !is.na(d.flt$message[1]) && any(d.flt$message != '')) {
        out <- c(out, paste(unique(d.flt$message[d.flt$message != '']), collapse=', '))
      }
    } else {
      out <- c(out, 'low-coverage')
    }
    
    paste(out, collapse=', ')
  })
  
  
  library(openxlsx2)
  
  # general information
  
  max_n_seq <- max(seq_tab_def$top_n_seqs, na.rm=T)
  n_clust <- min(max_n_seq, max_seqs %||% 6)
  n_curate <- min(max_n_seq, n_curate %||% 2)
  
  # filtered data
  
  # TODO: where does amplicon order get lost in sample_tab -> seq_tab?
  seq_tab_def <- seq_tab_def[order(seq_tab_def$amplicon, seq_tab_def$plate, seq_tab_def$well),]
  seq_tab_def$top_seq_taxon <- gsub('_', ' ', seq_tab_def$top_seq_taxon)
  seq_tab_def$`matching ranks` <- with(seq_tab_def,
                                      ifelse(is.na(top_seq_matching_ranks), NA,
                                      paste(top_seq_matching_ranks, 
                                            top_seq_matching_ranks + top_seq_mismatching_ranks, sep='/')))
  seq_tab_def$final_comment <- ifelse(
    !is.na(seq_tab_def$top_seq_mismatching_ranks) & seq_tab_def$top_seq_mismatching_ranks >= 3,
    'unexpected taxon',
    NA
  )
  # top_abundance < low_abund_threshold ~ 'low coverage',
  # top_consensus_diffs > 0 | !is.na(known_seq_diffs) & known_seq_diffs > 0 & !contaminated ~ 'uncertain sequence',
  # top_consensus_ambigs > 0 ~ 'seq. variation',
  seq_tab_def$link = '→ data'

  fixed_cols <- c(
    'amplicon', 'plate', 'well',
    'indexes', 'data link'='link', 'sample', 'sample type',
    'issues', '# seqs'='top_n_seqs', '# ambig'='top_seq_consensus_ambigs',
    '# seqs'='n_seqs', '% assigned'='target_cluster_frac',
    'morphospecies', 'auto-assigned'='top_seq_taxon', 'matching ranks',
    'seq'='known sequence', 'diffs'='known_seq_diffs',
    'comment'='final_comment'
  )
  
  names(fixed_cols) = ifelse(names(fixed_cols) == '', fixed_cols, names(fixed_cols))
  
  col_widths <- c(
    12, 5,
    5, 12, 8, 12, 10,
    25, 8, 8, 9, 8,
    22, 22, 8,
    5, 5,
    18
  )
  
  
  # more colums...
  
  # read coverage
  
  abund.m = t(simplify2array(lapply(seq_tab_def$clustering, function(d) {
    abund <- if (is.null(d) || nrow(d) == 0) NA else d$abundance[!d$unspecific]
    abund[1:n_clust]
  }), except=NA))
  colnames(abund.m) = paste0('cov_seq', 1:n_clust)
  
  # mismatches (flags)
  
  ambig.m = t(simplify2array(lapply(seq_tab_def$clustering, function(d) {
    a <- if (is.null(d$consensus_ambigs)) NA else d$consensus_ambigs[!d$unspecific]
    a[1:n_clust]
  }), except=NA))
  flags.m = ifelse(ambig.m > 0, 'a', NA)
  colnames(flags.m) = paste0('flags_seq', 1:n_clust)
  
  # (curated) sequences
  # with associated taxonomy
  d <- unique(do.call(rbind,
    lapply(seq_tab_def$clustering, function(d) cbind(d$consensus, d$taxon))
  ))
  rownames(d) <- paste0('u', seq_len(nrow(d)))
  unique_seqs <- setNames(names(d[,1]), d[,1])
  unique_seq_taxa <- d[,2]
  stopifnot(!duplicated(unique_seqs))
  uniq.m = as.data.frame(t(simplify2array(lapply(seq_tab_def$clustering, function(d) {
    seq <- if (is.null(d$consensus)) NA else d$consensus[!d$unspecific]
    unique_seqs[seq][1:n_curate]
  }), except=NA)))
  colnames(uniq.m) = paste0('unq', 1:n_curate)
  curseq.m = curfasta.m = matrix(NA, nrow=nrow(seq_tab_def), ncol=n_curate)
  colnames(curseq.m) = paste0('seq', 1:n_curate)
  colnames(curfasta.m) = paste0('FA', 1:n_curate)
  
  # combine all
  
  out <- cbind(seq_tab_def[,fixed_cols], 
              curseq.m, curfasta.m, uniq.m,
              flags.m, abund.m)
  
  rest_i <- (length(fixed_cols) + 1):ncol(out)
  cols <- c(
    fixed_cols,
    setNames(colnames(out)[rest_i], gsub('.+?_(.*)', '\\1', colnames(out)[rest_i]))
  )
  
  names(out) = names(cols)
  
  
  # top headers
  
  headers <- list(
    Overview = list(rng=c('issues', 'target_cluster_frac'), bg='#e5f5e0'),
    Taxonomy = list(rng=c('morphospecies', 'matching ranks'), bg='#deebf7'),
    `Known sequences` = list(rng=c('known sequence', 'known_seq_diffs'), bg='#fee0d2'),
    `Curation` = list(rng = c('final_comment', paste0('FA', n_curate)), bg='#ffedd5'),
    `Uniques` = list(rng = c('unq1', paste0('unq', n_curate)), bg='#efedf5'),
    `[Flags]` = list(rng=length(fixed_cols) + 3*n_curate + c(1, n_clust)),
    `Read coverage` = list(rng = length(fixed_cols) + 3*n_curate + n_clust + c(1, n_clust))
  )
  
  # header tooltips
  
  info <- list(
    'indexes' = 'Forward-reverse sample index combination',
    'top_n_seqs' = paste(
      'Number of haplotypes/polymorphic sequence variants for the top taxon'
    ),
    'top_seq_consensus_ambigs' = paste(
      '# of ambiguous bases in the alignment of the top sequence variant',
      'The presence of ambiguities indicates that not all haplotypes/polymorphisms',
      'were resolved by the clustering.',
      'Click on "→ data" to open a folder with BAM files,',
      'which can be inspected in a sequence viewer'
    ),
    'n_seqs' = paste(
      'Number of analyzed raw Nanopore reads'
    ),
    'target_cluster_frac' = paste(
      'Proportion of sequences that belong to the top taxon.',
      'Proportions < 50% may due to contamination or low sequence quality,',
      'or can also appear with low read depths (<500)'
    ),
    'issues' = 
      paste('Comma-separated list of issues. Issues [in brackets] are not necessarily ',
            'critical problems (may be false warnings), but may be inspected nonetheless\r\n') +
      fmt_txt('ambig-consensus:', bold=T) + paste(
        'At least one of the top taxon consensus sequences have ambiguities,',
        'suggesting either unresolved sequence variation, or sequencing errors',
        'such as long homopolmyers.',
        'To investigate, click on "→ data" to open a folder with alignment files.\r\n') +
      fmt_txt('[consensus-diffs]:', bold=T) + paste(
        'At least one of the top taxon consensus sequences does not match the',
        'dominant unique sequence. This is usually not a real issue.\r\n') +
      fmt_txt('tax-mismatch:', bold=T) + paste(
        'Inconsistency found between the name in the morphospecies column',
        'and the auto-assigned taxon (see also "matching ranks").',
        'Strong inconsistencies might also indicate a preparation error or',
        'contamination, although contamination is usually automatically recognized\r\n') +
      fmt_txt('contamination:', bold=T) + paste(
        'The most abundant taxon is suspected to be a contamination',
        'and was therefore removed.\r\n') +
      fmt_txt('known-seq-diffs:', bold=T) + paste(
        'There is at least one mismatch between the top taxon and the provided',
        '"known" sequence.',
        'To investigate, click on "→ data" to open a folder containing ...seq_comparison.fasta\r\n') +
      fmt_txt('known-seq-contamination:', bold=T) + paste(
        'The provided "known" sequence appears to be a contamination itself\r\n') +
      fmt_txt('many-variants:', bold=T) + paste(
        'Strong sequence variability for the top taxon',
        '(> 4 abundant variants above the frequency threshold).',
        'Maybe multiple individuals are present in the DNA mix?\r\n') +
      fmt_txt('low-coverage:', bold=T) + paste(
        sprintf('< %s Nanopore sequences support the top taxon;', low_abund_threshold),
        'Low-coverage samples may have more ambiguities, and there is a higher',
        'risk for errors (possibly investigate the BAM alignments)\r\n'),
    'top_seq_taxon' = paste(
      'Auto-assigned name for the top selected taxon'
    ),
    'matching ranks' = paste(
      'This column contains the number of matching/total taxonomic ranks',
      'from the comparison of the morphological and automatic sequence-based identifications'
    ),
    'known_seq_diffs' = paste(
      'Number of differences (substitutions/InDels) between the provided "known" sequence',
      'and the top selected sequence'
    ),
    'final_comment' = paste(
      'Comment that indicates why a certain sequence was chosen (or not).',
      'If there are >2 mismatching ranks in the comparison between morphological',
      'and sequence-based identification (see "matching ranks"), "unexpected taxon"',
      'is automatically pre-filled, and all sequences in the curation section are hidden.'
    ),
    'seq1' = paste(
      'Finally chosen sequence(s) (autofilled unless there are problems).',
      'This can be edited by manually filling another sequence from the "details"',
      'workbook or from an alignment, etc.'
    ),
    'FA1' = paste(
      'Sequences in FASTA format (can be copied)'
    )
  )
  
  row_i <- 2+1:nrow(out)
  row_i_all <- 1:(nrow(out)+2)
  
  ref <- function(col, row, lock=NA) {
    stopifnot(!is.na(col))
    stopifnot(!is.na(row))
    paste0(if (!is.na(lock) && lock == 'col') '$' else '',
           int2col(col),
           if (!is.na(lock) && lock == 'row') '$' else '',
           row)
  }
  
  wb <- wb_workbook()
  wb$add_worksheet('samples')
  
  # data and headers
  
  wb$add_data(x=out, na.strings='', start_row = 2, with_filter = T)
  wb$freeze_pane(first_active_row = 3)
  wb$add_font(dims = wb_dims(rows = 2, cols = 1:ncol(out)),
              bold = TRUE)
  wb$add_cell_style(dims = wb_dims(rows = 2, cols = 1:ncol(out)), wrap_text = T)
  for (h in names(headers)) {
    d <- headers[[h]]
    if ('character' %in% class(d$rng)) {
      d$rng = match(d$rng, cols)
    }
    cols_ <- d$rng[1]:d$rng[2]
    wb$add_data(x=h, col_names=F, start_col=cols_[1])
    wb$add_font(dims = wb_dims(cols=cols_[1], rows=1), bold = T, size = 13)
    wb$merge_cells(dims = wb_dims(rows=1, cols=cols_))
    if (!is.null(d$bg)) {
      wb$add_fill(dims = wb_dims(rows=row_i_all, cols=cols_),
                  color = wb_color(d$bg))
    }
  }
  
  for (h in names(info)) {
    wb$add_comment(dims = wb_dims(cols = which(cols == h), rows = 2), 
                   comment = wb_comment(text = info[[h]], author = 'info'))
  }
  
  wb$set_col_widths(
    cols = 1:ncol(out),
    widths = c(col_widths, rep(5, 3*n_curate), rep(8, 2*n_clust)),
    hidden = c(rep(F, length(fixed_cols)),
               rep(F, 2*n_curate),
               # show first two "uniques" columns
               c(F, F, rep(T, max(0, n_curate-2)))[1:n_curate],
               rep(T, n_clust), 
               rep(F, n_clust))
  )
  
  
  # curation output
  
  
  fa_desc_formula <- sprintf(
    '%s&"-"&%s&" "&%s&" / "&%s', 
    ref(match('plate', cols), row_i, lock='col'),
    ref(match('well', cols), row_i, lock='col'),
    ref(match('morphospecies', cols), row_i, lock='col'),
    ref(match('top_seq_taxon', cols), row_i, lock='col')
  )

  
  # highlight ambigs
  
  wb$add_dxfs_style(
    name = 'ambig',
    font_color = wb_color(hex = "bb34bd")
  )
  ranges <- list(
    length(fixed_cols) + c(1, n_curate), 
    length(fixed_cols) + c(1, n_clust)
  )
  for (i in 0:2) {
    wb$add_conditional_formatting(
      dims = wb_dims(rows=row_i, cols=(length(fixed_cols) + i*n_curate) + c(1, n_curate)),
      style = 'ambig',
      rule = sprintf('%s3="a"', int2col(length(fixed_cols) + 3 * n_curate + 1)),
    )
  }
  
  
  wb$add_dxfs_style(
    name = 'edited',
    font_color = wb_color(hex = "0108ff")
  )
  
  for (i in 1:n_curate) {
    u <- ref(match(paste0('unq', i), cols), row_i)
    wb$add_formula(
      x=sprintf(
        'IF(AND(ISERROR(MATCH(%s, _internal!$A2:$A20)), %s <> ""), VLOOKUP(%s, sequences!$A:$B, 2, 0), "")',
        ref(match('final_comment', cols), row_i),
        u, u
      ),
      dims=wb_dims(cols=match(paste0('seq', i), cols), rows=row_i)
    )
    
    wb$add_formula(
      x=sprintf(
        'IF(%s = "", "", ">"&%s&" "&%s&CHAR(13)&%s)', 
        ref(match(paste0('seq', i), cols), row_i),
        ref(match('sample', cols), row_i, lock='col'),
        fa_desc_formula,
        ref(match(paste0('seq', i), cols), row_i)
      ),
      dims=wb_dims(cols=match(paste0('FA', i), cols), rows=row_i)
    )
    
    for (prefix in c('seq', 'FA')) {
      wb$add_conditional_formatting(
        dims = wb_dims(rows=row_i, cols=match(paste0(prefix, i), cols)),
        style = 'edited',
        rule = sprintf('%s <> VLOOKUP(%s, sequences!$A:$B, 2, 0)', 
                       ref(match(paste0('seq', i), cols), min(row_i), lock='col'),
                       ref(match(paste0('unq', i), cols), min(row_i), lock='col'))
      )
    }
  }
  
  # hyperlinks
  
  wb$add_hyperlink(dims = wb_dims(cols = match('link', cols), rows = row_i),
                   target = paste0('./alignments/', seq_tab_def$primer_indexes))
  
  # formatting
  
  # wb$add_numfmt(dims=paste0('M2:M', nrow(d)), numfmt='0%')
  haplo_format <- function(rows, cols) {
    wb$add_conditional_formatting(
      dims = wb_dims(rows=rows, cols=cols),
      style = c('#fafff8', '#FFEB84', 'orange'),
      rule = c(1, 2, 5),
      type = 'colorScale'
    )
  }
  ambig_format <- function(rows, cols) {
    wb$add_conditional_formatting(
      dims = wb_dims(rows=rows, cols=cols),
      style = c('#fafff8', '#FFEB84', 'orange'),
      rule = c(1, 2, 5),
      type = 'colorScale'
    )
  }
  abund_format <- function(rows, cols) {
    wb$add_conditional_formatting(
      dims = wb_dims(rows=rows, cols=cols),
      style = c('#F8696B', '#FFEB84', '#7adb93'),
      rule = c(0, 10, 100),
      type = 'colorScale'
    )
  }
  
  haplo_format(row_i, match('top_n_seqs', cols))
  ambig_format(row_i, match('top_seq_consensus_ambigs', cols))
  abund_format(row_i, length(fixed_cols) + 3*n_curate + n_clust + 1:n_clust)
  
  wb$add_conditional_formatting(
    dims = wb_dims(rows=row_i, cols=match('n_seqs', cols)),
    style = c('#fc8d59', '#ffffbf', '#99d594'),
    rule = c(low_abund_threshold, 200, 1000),
    type = 'colorScale'
  )
  wb$add_conditional_formatting(
    dims = wb_dims(rows=row_i, cols=match('target_cluster_frac', cols)),
    style = c('#d95f0e', '#FFEB84', 'white'),
    rule = c(0, 0.3, 1),
    type = 'colorScale'
  )
  
  wb$add_numfmt(
    dims = wb_dims(row_i, match('target_cluster_frac', cols)),
    numfmt = '0 %'
  )
  
  # taxonomy issues
  
  lv <- sprintf('%s-suspicious-taxon', c('slightly', 'quite', 'very'))
  wb$add_dxfs_style(
    name = lv[1],
    font_color = wb_color(hex = "987341"),
    bg_fill = wb_color(hex = "e9f3fb")
  )
  wb$add_dxfs_style(
    name = lv[2],
    font_color = wb_color(hex = "ce9a11"),
    bg_fill = wb_color(hex = "e9f3fb")
  )
  wb$add_dxfs_style(
    name = lv[3],
    font_color = wb_color(hex = "f8696b"),
    bg_fill = wb_color(hex = "e9f3fb"),
    text_bold = T
  )
  tax_format <- function(rows, match_col, colrng) {
    ref1 <- ref(match_col, min(rows), lock='col')
    wb$add_conditional_formatting(
      dims = wb_dims(rows=rows, cols=colrng),
      style = lv[1],
      rule = sprintf('RIGHT(%s,1)-LEFT(%s,1)=1', ref1, ref1),
    )
    wb$add_conditional_formatting(
      dims = wb_dims(rows=rows, cols=colrng),
      style = lv[2],
      rule = sprintf('RIGHT(%s,1)-LEFT(%s,1)=2', ref1, ref1),
    )
    wb$add_conditional_formatting(
      dims = wb_dims(rows=rows, cols=colrng),
      style = lv[3],
      rule = sprintf('RIGHT(%s,1)-LEFT(%s,1)>=3', ref1, ref1),
    )
  }
  
  colrng <- match(c('top_seq_taxon', 'matching ranks'), cols)
  tax_format(row_i, colrng[2], colrng)
  
  
  # Details sheet
  
  dcols <- c(
    'indexes', 'data link'='link', 'name'='id', 'unique id'='unique_id',
    'taxon group'='group', 'sequence', 'FASTA'='fasta',
    '# reads'='abundance', '# cons. diffs'='consensus_diffs',  '# ambig'='consensus_ambigs',
    'auto-taxon'='taxon', 'matching ranks',
    'contaminant'='is_contaminant', 'unspecific'
  )
  names(dcols) = ifelse(names(dcols) == '', dcols, names(dcols))
  
  dcol_widths <- c(
    12, 8, 12, 10,
    10, 15, 15,
    8, 8, 8,
    22, 8,
    10, 10
  )
  
  l <- lapply(1:nrow(seq_tab_def), function(i) {
    d <- seq_tab_def$clustering[[i]]
    if (!is.null(d)) {
      d$indexes <- seq_tab_def$indexes[i]
      spec <- !d$unspecific
      d$link <- '→ data'
      d$fasta <- NA
      if (!is.null(d$matching_ranks) && any(!is.na(d$matching_ranks))) {
        d$`matching ranks` = with(d, paste(matching_ranks, matching_ranks + mismatching_ranks, sep='/'))
      } else {
        d$`matching ranks` <- NA
      }
      d$unique_id <- unique_seqs[d$consensus]
      d$consensus<- ''
      d <- d[dcols]
      names(d) <- names(dcols)
      # adding extra (needed for data links) info to thend
      cbind(d, primer_indexes = seq_tab_def$primer_indexes[i])
    }
  })
  d <- do.call(rbind, l)

  lgl <- sapply(d, is.logical)
  d[lgl] = lapply(d[lgl], function(l) ifelse(l, 'x', ''))
  
  wb$add_worksheet('details')
  wb$add_data(x=d[names(dcols)], na.strings='', with_filter = T)
  wb$freeze_pane(first_active_row = 2)
  wb$add_font(dims = wb_dims(rows = 1, cols = 1:ncol(d)), bold = T)
  wb$add_cell_style(dims = wb_dims(rows = 1, cols = 1:ncol(d)), wrap_text = T)
  wb$set_col_widths(
    cols = 1:length(dcol_widths),
    widths = dcol_widths
  )
  
  drow_i <- 1 + 1:nrow(d)
  uref <- sprintf('$%s%d', int2col(match('unique_id', dcols)), drow_i)
  wb$add_formula(
    x=sprintf('IFERROR(VLOOKUP(%s, sequences!$A:$B, 2, 0), "")', uref),
    dims=wb_dims(drow_i, match('sequence', dcols))
  )
  wb$add_formula(
    x=sprintf('IFERROR(VLOOKUP(%s, sequences!$A:$C, 3, 0), "")', uref),
    dims=wb_dims(drow_i, match('taxon', dcols))
  )

  seq_i <- match('sequence', dcols)
  wb$add_formula(
    x=sprintf('IF(%s = "", "", ">"&%s&"-"&%s&" "&%s&" "&%s&CHAR(13)&SUBSTITUTE(%s, "-", ""))', 
              ref(seq_i, drow_i),
              ref(match('indexes', dcols), drow_i),
              ref(match('unique_id', dcols), drow_i),
              ref(match('id', dcols), drow_i),
              ref(match('taxon', dcols), drow_i),
              ref(seq_i, drow_i)),
    dims=wb_dims(drow_i, match('fasta', dcols))
  )
  
  wb$add_hyperlink(dims = wb_dims(cols = match('link', dcols), rows = drow_i),
                   target = paste0('./alignments/', d$primer_indexes))

  # formatting
  
  ambig_format(drow_i, match('consensus_ambigs', dcols))
  ambig_format(drow_i, match('consensus_diffs', dcols))
  abund_format(drow_i, match('abundance', dcols))
  
  colrng <- match(c('taxon', 'matching ranks'), dcols)
  tax_format(drow_i, colrng[2], colrng)
  
  # Sequences sheet
  
  wb$add_worksheet('sequences')
  d <- data.frame(
    unique_name = unique_seqs, 
    sequence = names(unique_seqs),
    taxon = unique_seq_taxa[unique_seqs]
  )
  wb$add_data(x=d)
  
  # Internals
  
  wb$add_worksheet('_internal')
  wb$add_data(x=data.frame(auto_comments = unique(na.omit(out$comment))))
  
  # Save
  
  wb$save(outfile)
}