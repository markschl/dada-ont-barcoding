
.valid_ranks <- c('domain',
                  'kingdom',
                  'phylum',
                  'class',
                  'order',
                  'family',
                  'genus',
                  'species')

#' Download a taxonomic database and save it in UTAX format
#'
#' @param db_url Character vector with one or more file URLs; sequences from all
#'   URLs are combined into one final database.
#' @param db_type Database type (character) (see details)
#' @param db_file Output file (will be GZIPped FASTA, `.fasta.gz`)
#' @param db_tax_url (optional) URL(s) pointing to file(s) containing the taxonomy for the
#'   sequences in the files at `db_url` (currently implemented for *qiime_qza* format)
#'
#' @details
#'
#' Currently, three database types are implemented (more may follow):
#'
#' - *"unite"*: UNITE ITS database (https://unite.ut.ee/repository.php):
#'   choose a dataset (e.g. [all Eukaryotes](https://doi.plutof.ut.ee/doi/10.15156/BIO/3301243)),
#'   click on the file at *Downloads* and copy the actual file URL once the file
#'   starts downloading in the browser.
#'   The SH sequence set with a *dynamic clustering threshold* is chosen from the archive
#'   and imported.
#' - *"qiime_fasta"*: FASTA file with lineages in headers in the form
#'   `<seqid> rank1__name; rank2__name; ...`
#'   whereby `rank1` and `rank2` are one-letter codes such as `k` (for kingdom), `g` (for genus), etc.
#'   If only sequences are in the file, specify the URL of a tab-delimited file
#'   mapping sequence IDs to lineages.
#'   One example of database providing such files [GTDB](https://gtdb.ecogenomic.org)
#'   (e.g. [here](https://data.gtdb.aau.ecogenomic.org/releases/release226/226.0/genomic_files_reps))
#' - *"qiime_qza"*: [QIIME2 QZA](https://amplicon-docs.qiime2.org/en/stable/explanations/archives.html) format,
#'   which is essentially a ZIP file containing a *qiime_fasta*-formatted FASTA
#'   file, while another QZA file should contain the tab-delimited taxonomy (provide with `db_tax_url`).
#'
#' This function may be expanded and made more flexible in the future
#'
#' @export
load_taxdb <- function(db_url,
                       db_type = c('unite', 'qiime_fasta', 'qiime_qza'),
                       db_file, db_tax_url=NULL, ...) {
  db_type <- match.arg(db_type)
  load_fn <- get(paste0('load_', db_type))
  load_fn(db_url, db_file, tax_urls=db_tax_url)
}

load_unite <- function(url, outfile, ...) {
  stopifnot(length(url) == 1)  # multiple files not implemented
  download.file(url, 'unite.tar.gz', timeout = 1200, quiet = TRUE)
  untar('unite.tar.gz', exdir = 'unite')
  fasta <- list.files('unite', pattern = '.fasta')
  fasta <- fasta[grepl('dynamic', fasta)][1]
  tax <- list.files('unite', pattern = '.txt')
  tax <- tax[grepl('dynamic', tax)]
  seqs <- read_dna(file.path('unite', fasta))
  tax <- read.delim(file.path('unite', tax))
  tax <- setNames(tax[,2], tax[,1])
  tax <- tax[names(seqs)]
  tax <- gsub(';sh__.+', '', tax)
  stopifnot(!is.na(tax))
  convert_write_utax(seqs, outfile, taxonomy=tax, ...)
  unlink('unite', TRUE)
  invisible(file.remove('unite.tar.gz'))
}

load_qiime_fasta <- function(urls, outfile, tax_urls=NULL, ...) {
  stopifnot(is.null(tax_urls))  # taxonomy expected in FASTA file
  tmp_db <- paste0(outfile, '.tmp.gz')
  seqs <- do.call(c, lapply(urls, function(url) {
    download.file(url, tmp_db, timeout = 600, quiet = TRUE)
    read_dna(tmp_db)
  }))
  convert_write_utax(seqs, outfile, ...)
  invisible(file.remove(tmp_db))
}

load_qiime_qza <- function(urls, outfile, tax_urls=NULL, ...) {
  tmp_db_dir <- paste0(outfile, '.taxdb_tmp')
  tmp_taxdb_dir <- paste0(outfile, '.taxdb_tmp_tax')
  seqs <- do.call(c, lapply(urls, function(url) {
    read_dna(load_qza(url, tmp_db_dir, 'dna-sequences.fasta'))
  }))
  taxonomy <- if (!is.null(tax_urls)) {
    tax <- do.call(rbind, lapply(tax_urls, function(url) {
      read.delim(load_qza(url, tmp_taxdb_dir, 'taxonomy.tsv'))
    }))
    tax <- setNames(tax[,2], tax[,1])
    tax <- tax[names(seqs)]
    stopifnot(!is.na(tax))
    tax
  } else {
    NULL
  }
  convert_write_utax(seqs, outfile, taxonomy=taxonomy, ...)
  unlink(tmp_db_dir, TRUE)
  if (!is.null(tax_url)) {
    unlink(tmp_taxdb_dir, TRUE)
  }
}

load_qza <- function(url, outdir, fname) {
  qza <- paste0(outdir, '.qza')
  download.file(url, qza, timeout = 600, quiet = TRUE)
  unlink(outdir, TRUE)
  dir.create(outdir, recursive = TRUE)
  unzip(qza, exdir = outdir)
  file.remove(qza)
  file.path(outdir, list.files(outdir)[1], 'data', fname)
}


# Convert QIIME-style lineages in sequence names (separated from sequence IDs by space)
# to UTAX-style taxonomic lineages
convert_utax <- function(seqs,
                         taxonomy = NULL,
                         unknown_pat = '(undefined|unknown|incertae[_ ]sedis)',
                         sp_pat = r'{[ _]+(sp|spec)\.?$}',
                         # by default, we require at least an ~order name
                         required_ranks = 4,
                         ...) {

  if (is.null(taxonomy)) {
    s <- stringr::str_split_fixed(names(seqs), '\\s+', 2)
    names(seqs) <- s[,1]
    taxonomy <- s[,2]
  }
  lineages <- strsplit(taxonomy, ';', fixed=TRUE)
  lineages <- lapply(lineages, function(l) {
    s <- stringr::str_split(l, stringr::fixed('__'), 2)
    if (any(lengths(s) == 1)) {
      stop("Lineage not in form 'rank__Taxon' where 'rank' is a 1-letter rank code: ",
           paste(l, collapse=';'))
    }
    s <- do.call(rbind, s)
    setNames(trimws(s[,2]), trimws(s[,1]))
  })
  rank_codes <- unique(unlist(lapply(lineages, names)))
  o <- do.call(rbind, lapply(lineages, function(l) match(rank_codes, names(l))))
  rank_codes <- rank_codes[order(colMeans(o))]
  valid_codes <- substr(.valid_ranks, 1, 1)
  invalid_codes <- setdiff(rank_codes, valid_codes)
  if (length(invalid_codes) > 0) {
    stop('Unknown/invalid rank codes in taxonomy: ', paste(invalid_codes, collapse=', '))
  }
  lineages <- do.call(rbind, lapply(lineages, '[', rank_codes))
  colnames(lineages) <- valid_codes[match(rank_codes, valid_codes)]
  # remove GTDB-style annotations from end
  lineages[,ncol(lineages)] <- gsub(' +\\[\\w+=.*?\\]', '', lineages[,ncol(lineages)])
  # set unknown names to NA
  lineages[grepl(unknown_pat, lineages, ignore.case=TRUE, perl=TRUE)] <- NA
  if ('s' %in% colnames(lineages)) {
    lineages[grepl(sp_pat, lineages[,'s'], ignore.case=TRUE, perl=TRUE), 's'] <- NA
  }
  lineages[lineages == ''] <- NA

  # rank level filter
  sel <- rank_level(lineages) >= required_ranks
  lineages <- lineages[sel, , drop=FALSE]
  seqs <- seqs[sel]

  # convert to UTAX
  # escape reserved chars
  lineages <- apply(lineages, 2, function(t) gsub('[:,\\s]', '_', t, perl=TRUE))
  lineages <- apply(lineages, 1, function(l) {
    l <- na.omit(l)
    paste(paste(names(l), l, sep=':'), collapse=',')
  })
  names(seqs) <- sprintf('%s;tax=%s', names(seqs), lineages)
  seqs
}

convert_write_utax <- function(seqs, outfile, ...) {
  seqs <- convert_utax(seqs, ...)
  if (!dir.exists(dirname(outfile))) {
    dir.create(dirname(outfile), FALSE, TRUE)
  }
  write_dna(seqs, outfile, compress = endsWith(outfile, '.gz'))
}

write_taxdb <- function(seqs, outfile) {
  write_dna(utax_seqs, compress = endsWith(outfile, '.gz'))
}

#' Run taxonomy assignment (SINTAX algorithm)
#'
#' This requires [VSEARCH](https://github.com/torognes/vsearch) to be installed.
#'
#' @export
assign_taxonomy_sintax <- function(seq_file,
                                   utax_db,
                                   confidence_threshold = 0.8,
                                   tmp_prefix = NULL,
                                   threads = 1) {
  vsearch <- get_program('vsearch')
  out <- run_bash(
    c(
      vsearch,
      '-sintax', seq_file,
      '-db', utax_db,
      '-tabbedout', '-',
      '-strand', 'both',
      '-sintax_cutoff', confidence_threshold,
      '-quiet',
      '-threads', threads
    ),
    stdout = TRUE
  )
  tax <- read.delim(
    textConnection(out),
    header = FALSE,
    col.names = c('id', 'details', 'strand', 'lineage')
  )
  lineages <- parse_utax_lineages(tax$lineage)
  rownames(lineages) = tax$id
  lineages
}


parse_utax_lineages <- function(lineages) {
  lineages <- strsplit(as.character(lineages), ',', fixed = TRUE)
  lineages <- lapply(lineages, function(l) {
    s <- stringr::str_split_fixed(l, stringr::fixed(':'), 2)
    setNames(s[,2], s[,1])
  })
  rank_codes <- unique(unlist(lapply(lineages, names)))
  if (length(rank_codes) == 1 && is.na(rank_codes)) {
    # all-NA -> return just kingdom = NA
    rank_codes <- 'k'
  }
  o <- do.call(rbind, lapply(lineages, function(l) match(rank_codes, names(l))))
  rank_codes <- rank_codes[order(colMeans(o))]
  valid_codes <- substr(.valid_ranks, 1, 1)
  invalid_codes <- setdiff(rank_codes, valid_codes)
  if (length(invalid_codes) > 0) {
    stop('Unknown/invalid rank codes in assigned taxonomy: ', paste(invalid_codes, collapse=', '))
  }
  lineages <- do.call(rbind, lapply(lineages, '[', rank_codes))
  colnames(lineages) <- .valid_ranks[match(rank_codes, valid_codes)]
  if ('species' %in% colnames(lineages)) {
    lineages[, 'species'] = gsub('_', ' ', lineages[, 'species'])
  }
  lineages[lineages == ''] = NA
  lineages
}

make_taxon_name <- function(lineages) {
  for (rank in c('genus', 'species')) {
    if (!(rank %in% colnames(lineages))) {
      lineages <- cbind(lineages, NA)
      colnames(lineages)[ncol(lineages)] <- rank
    }
  }
  highest_taxon <- apply(lineages, 1, function(l) tail(na.omit(l), 1)[1])
  taxname <- ifelse(is.na(highest_taxon), 'Unknown', paste('Unknown', highest_taxon))
  has_spec <- !is.na(lineages[, 'species'])
  gen_only <- !is.na(lineages[, 'genus']) & !has_spec
  taxname[has_spec] <- lineages[has_spec, 'species']
  taxname[gen_only] <- paste(lineages[gen_only & !has_spec, 'genus'], 'sp.')
  taxname
}

rank_level <- function(tax) {
  apply(tax, 1, function(x) max(0, which(!is.na(x))))
}

normalize_taxa <- function(taxa) {
  t <- gsub('[\r\n"\']', '', trimws(taxa))
  t <- gsub('[\t]', ' ', t)
  t <- gsub('[|/]', '_', t)
  t <- gsub(' +(spec|sp|indet)\\.?$', '', t)
  t <- gsub('^(cf\\.?|[Uu]nknown) ', '', t)
  gsub('Fungus\\b', 'Fungi', t)
}

#' Match taxonomic names/lineages against the GBIF taxonomy
#'
#' This function uses the GBIF web API
#' (through [rgbif](https://docs.ropensci.org/rgbif/reference/name_backbone_checklist.html))
#' to retrieve the currently accepted taxon name and lineage
#' according to the GBIF
#' [backbone taxonomy](https://www.gbif.org/dataset/d7dddbf4-2cf0-4f39-9b2a-bb099caae36c).
#'
#' @param taxa A character vector of taxonomic names or a matrix/data frame with lineages;
#' columns must have rank names and the first column should be *name*, containing
#' taxonomic names such as "Russula sp." or "Bacteria". The lineages don't need
#' to be complete (NA or only few ranks allowed), but their presence may help
#' with the matching.'
#' @param cache_file Path where past searches are stored in a tab-separated text format.
#' This avoids havoing to query the GBIF API repeatedly for thousands of taxa.
#' Delete the file (or individual lines) to re-query the GBIF API for the taxa.
#' @param known_kingdom provide a known kingdom as alternative to providing the
#' kingdom in a matrix of lineages. Supplying the kingdom may help to resolve
#' ambiguous matches.
#' @param likely_kingdom If the kingdom is not known for sure (but very likely),
#' it can be provided here. The GBIF API will be searched without and with kingdom,
#' and the result yielding a more complete lineage is returned
#'
#' @returns A matrix with lineages; rows are in the same order as in the `taxa`
#' vector or matrix
#'
#' @export
get_gbif_taxa <- function(taxa,
                          cache_file,
                          known_kingdom = NA,
                          likely_kingdom = NA,
                          verbose = FALSE) {
  gbif_ranks <- c('kingdom',
                  'phylum',
                  'class',
                  'order',
                  'family',
                  'genus',
                  'species')

  taxa_key <- function(d) {
    apply(d, 1, function(l) paste(gsub('[|\\s]', '_', l, perl=TRUE), collapse='|'))
  }

  gbif_match <- function(d) {
    stopifnot(inherits(d, 'data.frame'))
    stopifnot('name' %in% names(d))
    # TODO: what would be the best way to map uniques back to original position?
    u <- unique(d)
    if (nrow(u) != nrow(d)) {
      m <- match(taxa_key(d), taxa_key(u))
    } else {
      m <- seq_len(nrow(d))
    }
    t <- as.data.frame(suppressMessages(rgbif::name_backbone_checklist(d)))
    t <- t[match(1:nrow(d), t$verbatim_index), ]
    stopifnot(t$verbatim_name == d$name)
    for (r in setdiff(gbif_ranks, names(t))) {
      t[[r]] = NA_character_
    }
    t <- t[gbif_ranks]
    row.names(t) <- NULL
    t[m,]
  }

  stopifnot(!is.null(known_kingdom))
  stopifnot(!is.null(likely_kingdom))

  if (inherits(taxa, 'data.frame') || inherits(taxa, 'matrix')) {
    taxa <- as.data.frame(taxa)
    stopifnot('name' %in% names(taxa))
  } else {
    stopifnot(is.null(dim(taxa)))
    taxa <- data.frame(name = taxa)
  }
  names(taxa)[names(taxa) == 'domain'] <- 'kingdom'
  for (r in setdiff(gbif_ranks, names(taxa))) {
    taxa[[r]] <- NA_character_
  }
  taxa <- taxa[c('name', gbif_ranks)]
  taxa$name <- normalize_taxa(taxa$name)
  lineages <- NULL
  taxa <- data.frame(key = taxa_key(taxa), taxa)
  unique_taxa <- unique(taxa)
  stopifnot(!duplicated(unique_taxa$key))
  row.names(unique_taxa) <- unique_taxa$key
  unique_taxa$key <- NULL
  if (file.exists(cache_file)) {
    lineages <- read.delim(cache_file, na = '')
    row.names(lineages) <- lineages$taxon
    lineages$taxon <- NULL
    query_taxa <- unique_taxa[setdiff(row.names(unique_taxa), row.names(lineages)),]
  } else {
    query_taxa <- unique_taxa
  }
  if (nrow(query_taxa) > 0) {
    if (verbose)
      message('Matching ', nrow(query_taxa), ' taxa against GBIF backbone')

    # try without/with alternative kingdom
    # (narrowing down the putative kingdom is sometimes successful)
    try_kingdoms <- as.character(unique(c(known_kingdom, likely_kingdom)))
    tax.out <- data.frame(lapply(setNames(gbif_ranks, gbif_ranks),
                                 function(c) rep(NA_character_, nrow(query_taxa))),
                          row.names = row.names(query_taxa))
    for (kingdom in try_kingdoms) {
      sel <- is.na(tax.out$species)
      if (any(sel)) {
        d <- query_taxa[sel,]
        d$kingdom <- d$kingdom[is.na(d$kingdom)] <- kingdom
        tax <- gbif_match(d)
        r <- rank_level(tax)
        is_better <- r > rank_level(tax.out[sel,]) & (is.na(kingdom) | d$kingdom == kingdom | r > 1)
        tax.out[sel,][is_better,] <- tax[is_better,]
      }
    }
    # for all bad hits WITHOUT phylum-family, try searching genus names
    # inferred from the species
    if (all(sapply(query_taxa[gbif_ranks[2:5]], function(x) all(is.na(x))))) {
      sel <- is.na(query_taxa$genus) & is.na(tax.out$species) & is.na(is.na(tax.out$genus))
      s <- strsplit(query_taxa$name, ' ', fixed=TRUE)
      word1 <- sapply(s, '[[', 1)
      sel[lengths(s) == 1 | nchar(word1) < 5] <- FALSE
      query_taxa$genus[sel] <- word1[sel]
      for (kingdom in try_kingdoms) {
        if (any(sel)) {
          d <- query_taxa[sel,]
          d$kingdom <- d$kingdom[is.na(d$kingdom)] <- kingdom
          tax <- gbif_match(d)
          r <- rank_level(tax)
          is_better <- r > rank_level(tax.out[sel,]) & (is.na(kingdom) | d$kingdom == kingdom | r > 1)
          tax.out[sel,][is_better,] <- tax[is_better,]
          sel[is_better] <- FALSE
        }
      }
    }
    # *note*: names may still not be found
    # (https://docs.ropensci.org/rgbif/articles/taxonomic_names.html#too-many-choices-problem)
    lineages <- rbind(lineages, tax.out)
    write.table(cbind(taxon = row.names(lineages), lineages),
                cache_file, sep = '\t', na = '', quote = FALSE, row.names = FALSE)
  }
  stopifnot(!is.na(match(taxa$key, row.names(lineages))))
  lineages[taxa$key,]
}
