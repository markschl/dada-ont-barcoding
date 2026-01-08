

#### Helper functions ##########################################################


init_config <- function(analysis_dir) {
  config <- yaml::read_yaml(file.path(analysis_dir, 'config.yaml'))
  # some minimal configuration defaults (needed in report, etc.)
  # TODO: these defaults are the same in the function args, some redundancy here
  minimal_config_defaults <- list(
    reads_path = file.path(analysis_dir, 'reads.fastq.gz'),
    analysis_dir = analysis_dir,
    tmp_dir = file.path(analysis_dir, 'tmp'),
    taxdb_dir = 'taxdb',
    demultiplex = list(
      error_threshold = 2.5,
      primer_max_err = 0.2,
      idx_max_diffs = 0,
      min_barcode_length = 50,
      keep_trimmed = FALSE
    ),
    cluster = list(
      max_sample_depth = 5000,
      consensus_threshold = 0.65,
      fixed_cluster_threshold = 0.97,
      min_variant_freq = 0.2,
      dada_min_identical = 2,
      dada_min_n0 = 4
    ),
    taxonomy = list(
      confidence_threshold = 0.8,
      contam_rank_delta = 3
    ),
    output = list(
      report = list(
        low_abund_threshold = 20
      ),
      alignments = 'top_combined'
    )
  )
  config <- modifyList(minimal_config_defaults, config)
  stopifnot(file.exists(config$reads_path))
  # init directories
  dir.create(analysis_dir, FALSE)
  dir.create(config$taxdb_dir, FALSE)
  # output
  stopifnot(config$output$alignments %in% c('separate', 'top_combined', 'combined'))
  config$alignment_dir <- if (length(config$output$alignments) > 0) {
    file.path(ifelse('separate' %in% config$output$alignments, analysis_dir, config$tmp_dir),
              'separate_alignments')
  }
  # read metadata
  meta_file <- file.path(analysis_dir, 'meta.xlsx')
  config$sample_tab <- read_xlsx_sample_tab(meta_file, 'sample_list')
  config$amplicons <- read_xlsx_primer_tab(meta_file, 'primers',
                                           amplicons=levels(config$sample_tab$amplicon))
  # taxonomy assignment options
  stopifnot(!is.null(config$taxonomy))
  config$amplicon_taxdb <- init_nested_opts(config$taxonomy, names(config$amplicons),
                                            description = 'taxonomy')
  config
}

#' Propagates global options to nested sub-sections overriding them
#' (used for taxonomy databases)
init_nested_opts <- function(opts, nested_keys, description = NULL) {
  opts <- opts %||% list()
  global_keys <- setdiff(names(opts), nested_keys)
  out <- list()
  for (key in nested_keys) {
    out[[key]] <- modifyList(opts[global_keys], opts[[key]] %||% list())
  }
  out
}

run_with_cfg <- function(fn, ..., cfg = NULL, inner_fn = NULL, exclude = NULL) {
  argspec <- names(formals(fn))
  if (!is.null(inner_fn)) {
    inner_argspec <- names(formals(inner_fn))
    stopifnot(!('...' %in% inner_argspec))
    argspec <- c(argspec, inner_argspec)
  }
  args <- list(...)
  if (is.null(names(args)) && length(args) > 0)
    names(args) <- NA
  is_positional <- names(args) == '' | is.na(names(args))
  names(args)[is_positional] <- argspec[seq_len(sum(is_positional))]
  args <- modifyList(args %||% list(), cfg %||% list())
  used <- intersect(names(args), argspec)
  unused <- setdiff(names(args), c(argspec, exclude))
  if (length(unused) > 0) {
    warning(
      'The following config entries are not used by ', deparse(substitute(fn)),
      ": ", paste(unused, collapse = ", ")
    )
  }
  do.call(fn, args[used])
}

ensure_taxdb <- function(taxdb_dir, cfg) {
  stopifnot(!is.null(cfg$db_url))
  stopifnot(length(cfg$db_url) >= 1)
  stopifnot(!is.null(cfg$db_type))
  all_urls <- sort(c(cfg$db_url, cfg$db_tax_url))
  db_hash <- tools::md5sum(bytes=charToRaw(paste(all_urls, collapse='')))
  db_file <- file.path(taxdb_dir, paste0(db_hash, '.fasta.gz'))
  if (!file.exists(db_file)) {
    load_taxdb(db_url=cfg$db_url, db_type=cfg$db_type, db_file=db_file, db_tax_url=cfg$db_tax_url)
  }
  db_file
}

assign_taxonomy_amplicons <- function(seq_tab,
                                      taxdb_cfg,
                                      taxdb_dir,
                                      tmp_dir = NULL,
                                      cores = 1) {
  amp_summary_ranks <- setNames(vector('list', 3), names(taxdb_cfg))
  amplicons <- unique(seq_tab$amplicon)
  for (amplicon in amplicons) {
    cat(amplicon, '\n')
    sel <- seq_tab$amplicon == amplicon
    cfg <- taxdb_cfg[[amplicon]]
    db_file <- ensure_taxdb(taxdb_dir, cfg)
    gbif_cache_file <- file.path(taxdb_dir,
                                 paste0('gbif_taxa_', cfg$likely_kingdom %||% '', '.tsv'))
    ret <- run_with_cfg(
      do_assign_compare_taxonomy,
      seq_tab = seq_tab[sel, ],
      gbif_cache_file = gbif_cache_file,
      db_file = db_file,
      cores = cores,
      cfg = cfg,
      exclude = amplicons
    )
    amp_summary_ranks[[amplicon]] <- attr(ret, 'summary_ranks')
    seq_tab[, setdiff(names(ret), names(seq_tab))] <- NA
    seq_tab[sel, ] <- ret
  }
  list(seq_tab = seq_tab,
       amp_summary_ranks = setNames(amp_summary_ranks, names(taxdb_cfg)))
}

recluster_contaminated <- function(seq_tab,
                                   dada_err,
                                   cluster_cfg,
                                   taxdb_cfg,
                                   taxdb_dir,
                                   max_contam_sample_depth = 1e7,
                                   amp_summary_ranks = NULL,
                                   tmp_dir = NULL,
                                   cores = 1) {
  max_sample_depth <- cluster_cfg$max_sample_depth %||% max_contam_sample_depth
  do_recluster <- seq_tab$has_contamination & !is.na(seq_tab$n_reads) & seq_tab$n_reads >= max_sample_depth
  if (any(do_recluster)) {
    cluster_cfg <- modifyList(cluster_cfg %||% list(),
                              list(max_sample_depth = max_contam_sample_depth))
    for (amplicon in unique(seq_tab$amplicon[do_recluster])) {
      sel <- do_recluster & seq_tab$amplicon == amplicon
      amp_seqtab <- run_with_cfg(
        do_infer_all_barcodes,
        seq_tab[sel, ],
        dada_err,
        aln_out = config$alignment_dir,
        tmp_dir = tmp_dir,
        cores = cores,
        cfg = cluster_cfg
      )
      taxdb_cfg$summary_ranks <- amp_summary_ranks[[amplicon]]
      tax <- assign_taxonomy_amplicons(amp_seqtab,
                                       taxdb_cfg,
                                       taxdb_dir,
                                       tmp_dir = tmp_dir,
                                       cores = cores)
      stopifnot(seq_tab$indexes[sel] == tax$seq_tab$indexes)
      seq_tab[sel, ] <- tax$seq_tab
    }
  }
  seq_tab
}


#### Pipeline definition #######################################################

# obtain working directory, number of workers and program paths from environment vars
Sys.setenv('DadaNanoBC_ANALYSIS_DIR' = 'analysis/20251209_ONTbarcoding_16S_ITS-test')
Sys.setenv('DadaNanoBC_WORKERS' = 8)

analysis_dir <- Sys.getenv('DadaNanoBC_ANALYSIS_DIR', 'analysis')
n_workers <- as.integer(Sys.getenv('DadaNanoBC_WORKERS', '4'))
for (dep in c('seqtool', 'samtools', 'minimap2', 'vsearch')) {
  path <- Sys.getenv(paste0('DadaNanoBC_', dep))
  if (nzchar(path) > 0) {
    options(setNames(path, paste0('DadaNanoBC.', dep)))
  }
}

# install pipeline packages if needed
missing.pkgs <- setdiff(c('targets', 'tarchetypes', 'crew'),
                        installed.packages()[, 'Package'])
if (length(missing.pkgs) > 0) {
  install.packages(missing.pkgs)
}

library(targets)
library(tarchetypes)
library(crew)

single_ctrl <- crew_controller_local('full', workers = 1, seconds_idle = 900)
multi_ctrl <- crew_controller_local('parallel', workers = n_workers, seconds_idle = 900)

pkgs <- c('DadaNanoBC',
          'yaml',
          # for HTML report:
          'rmarkdown',
          'bookdown',
          'pander',
          'dplyr',
          'ggplot2',
          'patchwork',
          'xfun')

tar_option_set(
  packages = pkgs,
  controller = crew_controller_group(single_ctrl, multi_ctrl),
  resources = tar_resources(crew = tar_resources_crew(controller = 'full')),
  seed = 42
)

list(
  # Read sample and primer information and initialize configuration
  tar_target(
    config,
    init_config(analysis_dir)
  ),
  # Do primer search and demultiplexing
  tar_target(
    trim_demux,
    run_with_cfg(
      do_trim_demux,
      config$reads_path,
      config$amplicons,
      config$sample_tab,
      out_dir = file.path(config$tmp_dir, 'demux'),
      cores = n_workers,
      cfg = config$demultiplex
    )
  ),
  # DADA2 learnErrors() on a (pseudo-)random selection
  tar_target(
    dada_err,
    dada_learn_errors(
      sample(na.omit(trim_demux$seq_tab$reads_path)),
      omega_a = config$cluster$omega_a[1] %||% 1e-20,
      cores = n_workers,
      nbases = 1e6
    )
  ),
  # Do clustering with parallelization handled by 'targets'
  # (progress shown)
  tar_group_size(
    group_seqtab,
    trim_demux$seq_tab,
    8  # TODO: optimal batch size?
  ),
  tar_target(
    cluster_seqtab,
    run_with_cfg(
      do_infer_all_barcodes,
      group_seqtab,
      dada_err,
      aln_out = config$alignment_dir,
      tmp_dir = config$cluster$tmp_dir,
      cores = 1,
      inner_fn = infer_barcodes,
      cfg = config$cluster
    ),
    pattern = map(group_seqtab),
    resources = tar_resources(crew = tar_resources_crew(controller = 'parallel'))
  ),
  tar_target(
    cluster_seqtab_combined,
    bind_rows(cluster_seqtab) %>% select(-tar_group)
  ),
  # Taxonomy assignment
  tar_target(
    taxonomy,
    assign_taxonomy_amplicons(
      cluster_seqtab_combined,
      config$amplicon_taxdb,
      config$taxdb_dir,
      tmp_dir = config$tmp_dir,
      cores = n_workers
    )
  ),
  # Re-cluster contaminated samples with more reads if necessary
  tar_target(
    recluster_contam_seqtab,
    recluster_contaminated(
      taxonomy$seq_tab,
      dada_err,
      cluster_cfg = config$cluster,
      taxdb_cfg = config$amplicon_taxdb,
      taxdb_dir = config$taxdb_dir,
      max_contam_sample_depth = config$cluster$max_contam_sample_depth,
      amp_summary_ranks = taxonomy$amp_summary_ranks,
      tmp_dir = config$tmp_dir,
      cores = n_workers
    )
  ),
  # output files
  tar_render(
    rmd_report,
    'analysis.Rmd',
    file.path(analysis_dir, 'report.html')
  ),
  tar_target(xl_report, {
    out <- file.path(analysis_dir, 'report.xlsx')
    run_with_cfg(create_report, recluster_contam_seqtab, out, cfg = config$report)
    out
  }, format = 'file'),
  tar_target(
    combined_bam_top,
    if ('top_combined' %in% config$output$alignments)
      do_combine_alignments(
        recluster_contam_seqtab,
        config$alignment_dir,
        outdir = file.path(analysis_dir, 'top_alignments'),
        top_only = TRUE
      ),
    format = 'file'
  ),
  tar_target(
    combined_bam,
    if ('combined' %in% config$output$alignments)
      do_combine_alignments(
        recluster_contam_seqtab,
        config$alignment_dir,
        outdir = file.path(analysis_dir, 'alignments'),
        top_only = FALSE
      ),
    format = 'file'
  )
)
