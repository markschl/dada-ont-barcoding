
install_taxonomy_deps <- function() {
  p <- installed.packages()[,c('Package')]
  p <- setdiff(c('rbgif', 'stringr'), p)
  if (length(p) > 0) {
    install.packages(p)
  }
}

.tax_ranks = c('kingdom',
               'phylum',
               'class',
               'order',
               'family',
               'genus',
               'species')

load_utax_unite <- function(url, outfile, required_ranks = 5) {
  download.file(url, 'unite.tar.gz', timeout = 600, quiet = TRUE)
  untar('unite.tar.gz', exdir = 'unite')
  fasta <- list.files('unite', pattern = '.fasta')
  fasta <- fasta[grepl('dynamic', fasta)][1]
  tax <- list.files('unite', pattern = '.txt')
  tax <- tax[grepl('dynamic', tax)]
  unite <- Biostrings::readDNAStringSet(file.path('unite', fasta))
  # clean the taxonomy
  tax <- read.delim(file.path('unite', tax))
  lineages <- stringr::str_split_fixed(tax[, 2], ';', 8)[, 1:7]
  lineages <- apply(lineages, 2, function(t)
    gsub('([a-z])__+.+?_Incertae_sedis', '\\1__', t))
  lineages[, 7] = gsub('s__[^;]*_sp\\.?$', 's__', lineages[, 7])
  # convert to UTAX format
  # removing lineages not defined at the family level
  # (as those appear to be so poorly annotated that they confuse the classifier)
  names(unite) = qiime_to_utax(tax$Feature.ID, lineages, required_ranks =
                                 required_ranks)
  # remove seqs. with these empty lineages from the database
  unite <- unite[!endsWith(names(unite), 'tax=')]
  dir.create(dirname(outfile), FALSE, TRUE)
  writeXStringSet(unite, outfile, compress = endsWith(outfile, '.gz'))
  unlink('unite', TRUE)
  invisible(file.remove('unite.tar.gz'))
}

qiime_to_utax <- function(ids, lineages, required_ranks = 0) {
  lineages <- apply(lineages, 2, function(t)
    gsub('[:,]', '_', t)) # escape reserved chars
  lineages[grepl('^[a-z]__$', lineages, perl = T)] = NA
  lineages <- apply(lineages, 2, function(t)
    gsub('([a-z])__(.*)', '\\1:\\2', t))
  def_rank <- apply(lineages, 1, function(l)
    tail(which(!is.na(l)), 1)[1])
  lineages[is.na(def_rank) | def_rank < required_ranks, ] = NA
  sprintf('%s;tax=%s', ids, apply(lineages, 1, function(l)
    paste(na.omit(l), collapse = ',')))
}

assign_taxonomy <- function(seqs,
                           utax_db,
                           threshold = 0.8,
                           threads = 1,
                           verbose = FALSE,
                           ...) {
  seq_file <- '_sintax.seqs.fasta' # TODO: where to place temp. file?
  stopifnot(!is.null(names(seqs)))
  Biostrings::writeXStringSet(Biostrings::DNAStringSet(seqs), seq_file)
  cat(make_fasta(seqs), sep = '', file = seq_file)
  out <- run_bash(
    c(
      'vsearch',
      '-sintax', seq_file,
      '-db', utax_db,
      '-tabbedout', '-',
      '-strand', 'both',
      '-sintax_cutoff', threshold,
      if (verbose) NULL else '-quiet',
      '-threads', threads
    ),
    stdout = TRUE
  )
  tax <- read.delim(
    textConnection(out),
    header = FALSE,
    col.names = c('id', 'details', 'strand', 'lineage')
  )
  file.remove(seq_file)
  lineages <- parse_utax_lineages(tax$lineage, ...)
  rownames(lineages) = tax$id
  lineages
}

parse_utax_lineages <- function(lineages, ranks = .tax_ranks) {
  lineages <- strsplit(lineages, ',', fixed = TRUE)
  rank_chars <- substr(ranks, 1, 1)
  lineages <- t(simplify2array(lapply(lineages, function(l) {
    l <- strsplit(l, ':', fixed = T)
    if (length(l) == 0)
      return(rep(NA, length(rank_chars)))
    r <- sapply(l, '[', 1)
    n <- sapply(l, '[', 2)
    stopifnot(sapply(l, length) == 2)
    stopifnot(r %in% rank_chars)
    n[match(rank_chars, r)]
  }), except = NA))
  colnames(lineages) = ranks
  lineages[, 'species'] = gsub('_', ' ', lineages[, 'species'])
  lineages[lineages == ''] = NA
  lineages
}

def_rank <- function(tax) {
  apply(tax, 1, function(x) max(0, which(!is.na(x))))
}

normalize_taxa <- function(taxa) {
  taxa_norm <- gsub('[\r\n"\']', '', trimws(taxa))
  taxa_norm <- gsub('[\t]', ' ', taxa_norm)
  taxa_norm <- gsub(' spec\\.?$', ' sp.', taxa_norm)
  taxa_norm <- gsub('^(cf\\.?|[Uu]nknown) ', '', taxa_norm)
}

gbif_taxa <- function(taxa,
                     cache_file,
                     known_kingdom = NA,
                     likely_kingdom = NA,
                     verbose = FALSE) {
  stopifnot(!is.null(known_kingdom))
  stopifnot(!is.null(likely_kingdom))
  ranks <- .tax_ranks
  lineages <- NULL
  taxa_norm <- normalize_taxa(taxa)
  if (file.exists(cache_file)) {
    lineages <- read.delim(cache_file, na = '')
    row.names(lineages) <- lineages$taxon
    lineages$taxon <- NULL
    query_taxa <- na.omit(setdiff(taxa_norm, row.names(lineages)))
  } else {
    query_taxa <- unique(na.omit(taxa_norm))
  }
  if (length(query_taxa) > 0) {
    if (verbose)
      message('Matching ', length(query_taxa), ' taxa against GBIF backbone')
    
    # species match
    # also with alternative kingdom (narrowing down the kingdom is sometimes successful)
    gbif_match <- function(q, ...) {
      d <- data.frame(name = q, ...)
      u <- unique(d)
      t <- as.data.frame(suppressMessages(rgbif::name_backbone_checklist(u)))
      t <- t[match(1:nrow(u), t$verbatim_index),]
      stopifnot(t$verbatim_name == u$name)
      for (r in setdiff(ranks, names(t))) {
        t[[r]] = NA_character_
      }
      m <- match(d$name, u$name)
      stopifnot(!is.na(m))
      t <- t[m, ranks]
      row.names(t) <- NULL
      t
    }
    
    try_kingdoms <- as.character(unique(c(known_kingdom, likely_kingdom)))
    tax.out <- data.frame(lapply(setNames(ranks, ranks), function(c) rep(NA_character_, length(query_taxa))))
    for (kingdom in try_kingdoms) {
      sel <- is.na(tax.out$species)
      if (any(sel)) {
        tax <- gbif_match(query_taxa[sel], kingdom = kingdom)
        r <- def_rank(tax)
        is_better <- r > def_rank(tax.out[sel,]) & r > 1
        tax.out[sel,][is_better,] <- tax[is_better,]
      }
    }
    # for all bad hits, try searching genus names,
    query_genera <- sapply(strsplit(query_taxa, ' ', fixed=TRUE), '[', 1)
    for (kingdom in try_kingdoms) {
      sel <- is.na(tax.out$genus) &
        query_taxa != query_genera &
        grepl('^[A-Z]', query_genera) & 
        nchar(query_genera) >= 5
      if (any(sel)) {
        tax <- gbif_match(query_genera[sel], genus = query_genera[sel], kingdom = kingdom)
        r <- def_rank(tax.out[sel,])
        is_better <- def_rank(tax) > r & r > 1
        tax.out[sel,][is_better,] <- tax[is_better,]
      }
    }
    # TODO: names may still not be found
    # (https://docs.ropensci.org/rgbif/articles/taxonomic_names.html#too-many-choices-problem)
    
    # remove "false positives":
    # we know that 'Fungus indet.' does not match the correct rank in GBIF,
    # so correct for that
    sel <- grepl('Fungus\\b', query_taxa)
    if (any(sel)) {
      tax.out[sel, 'kingdom'] = 'Fungi'
      tax.out[sel, setdiff(ranks, 'kingdom')] = NA
    }
    
    row.names(tax.out) <- query_taxa
    lineages <- rbind(lineages, tax.out)
    write.table(cbind(taxon = row.names(lineages), lineages), 
                cache_file, sep = '\t', na = '', quote = FALSE, row.names = FALSE)
  }
  stopifnot(!is.na(match(taxa_norm, row.names(lineages))))
  lineages[taxa_norm,]
}
