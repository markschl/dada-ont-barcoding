

#### Helper functions ##########################################################


init_config <- function(config_file, analysis_dir) {
  config <- yaml::read_yaml(config_file)
  config$analysis_dir <- analysis_dir
  # some minimal configuration defaults (needed in report, etc.)
  # TODO: these defaults are the same in the function args, some redundancy here
  minimal_config_defaults <- list(
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
      report = list(low_abund_threshold = 20),
      alignments = 'separate'
    )
  )
  config <- modifyList(minimal_config_defaults, config)
  # output
  stopifnot(config$output$alignments %in% c('separate', 'top_combined', 'combined'))
  if (length(config$output$alignments) > 0) {
    if ('separate' %in% config$output$alignments) {
      config$alignment_dir <- config$output$report$bam_dir <- file.path(analysis_dir, 'separate_alignments')
    } else {
      config$alignment_dir <- file.path(tmp_dir, 'separate_alignments')
    }
  }
  # taxonomy assignment options
  stopifnot(!is.null(config$taxonomy))
  config$taxdb <- init_nested_opts(config$taxonomy, names(config$amplicons), description = 'taxonomy')
  config
}

.config_sections <- c(
  'demultiplex',
  'cluster',
  'taxdb',
  'output',
  'alignment_dir'
)

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

run_with_cfg <- function(fn,
                         ...,
                         cfg = NULL,
                         inner_fn = NULL,
                         exclude = NULL) {
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
      'The following config entries are not used by ',
      deparse(substitute(fn)),
      ": ",
      paste(unused, collapse = ", ")
    )
  }
  do.call(fn, args[used])
}

get_taxdb <- function(taxdb_dir, cfg) {
  stopifnot(!is.null(cfg$db_url), length(cfg$db_url) >= 1,
            !is.null(cfg$db_type))
  all_urls <- sort(c(cfg$db_url, cfg$db_tax_url))
  db_hash <- tools::md5sum(bytes = charToRaw(paste(all_urls, collapse = '')))
  db_file <- file.path(taxdb_dir, paste0(db_hash, '.fasta.gz'))
  if (!file.exists(db_file)) {
    load_taxdb(
      db_url = cfg$db_url,
      db_type = cfg$db_type,
      db_file = db_file,
      db_tax_url = cfg$db_tax_url
    )
  }
  db_file
}

assign_taxonomy_amplicons <- function(seq_tab, taxdb_cfg, taxdb_dir, cores = 1) {
  amplicons <- unique(seq_tab$amplicon)
  amp_summary_ranks <- setNames(vector('list', length(amplicons)), amplicons)
  for (amplicon in amplicons) {
    cat(amplicon, '\n')
    sel <- seq_tab$amplicon == amplicon
    cfg <- taxdb_cfg[[amplicon]]
    db_file <- get_taxdb(taxdb_dir, cfg)
    gbif_cache_file <- file.path(taxdb_dir,
                                 paste0('gbif_taxa_', cfg$likely_kingdom %||% '', '.tsv'))
    ret <- run_with_cfg(
      do_assign_compare_taxonomy,
      seq_tab = seq_tab[sel, ],
      gbif_cache_file = gbif_cache_file,
      db_file = db_file,
      cores = cores,
      cfg = cfg,
      exclude = c(amplicons, 'db_url', 'db_type', 'db_tax_url')
    )
    amp_summary_ranks[[amplicon]] <- base::attr(ret, 'summary_ranks')
    seq_tab[, setdiff(names(ret), names(seq_tab))] <- NA
    seq_tab[sel, ] <- ret
  }
  list(seq_tab = seq_tab,
       amp_summary_ranks = setNames(amp_summary_ranks, amplicons))
}

recluster_contaminated <- function(seq_tab,
                                   dada_err,
                                   cluster_cfg,
                                   taxdb_cfg,
                                   taxdb_dir,
                                   aln_dir,
                                   max_sample_depth,
                                   max_contam_sample_depth = NULL,
                                   amp_summary_ranks = NULL,
                                   cores = 1) {
  do_recluster <- seq_tab$has_contamination &
    !is.na(seq_tab$n_reads) & seq_tab$n_reads >= max_sample_depth
  if (any(do_recluster)) {
    cluster_cfg <- modifyList(cluster_cfg %||% list(),
                              list(max_sample_depth = max_contam_sample_depth))
    for (amplicon in unique(seq_tab$amplicon[do_recluster])) {
      sel <- do_recluster & seq_tab$amplicon == amplicon
      amp_seqtab <- run_with_cfg(
        do_infer_all_barcodes,
        seq_tab[sel, ],
        dada_err,
        aln_out = aln_dir,
        cores = cores,
        cfg = cluster_cfg,
        inner_fn = infer_barcode,
        exclude = 'max_contam_sample_depth'
      )
      taxdb_cfg$summary_ranks <- amp_summary_ranks[[amplicon]]
      tax <- assign_taxonomy_amplicons(amp_seqtab, taxdb_cfg, taxdb_dir, cores = cores)
      stopifnot(seq_tab$indexes[sel] == tax$seq_tab$indexes)
      seq_tab[sel, ] <- tax$seq_tab
    }
  }
  seq_tab
}


#### Pipeline definition #######################################################

# settings
# modify with Sys.getenv('DadaNanoBC_ANALYSIS_DIR' = ...)

analysis_dir <- Sys.getenv('DadaNanoBC_ANALYSIS_DIR', 'analysis')
n_workers <- as.integer(Sys.getenv('DadaNanoBC_WORKERS', max(parallel::detectCores(), 8)))
tmp_dir = file.path(analysis_dir, 'tmp')
taxdb_dir = 'taxdb'
dir.create(taxdb_dir, FALSE, TRUE)

message(sprintf(
  "Running pipeline in '%s' with %d workers",
  analysis_dir,
  n_workers
))


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

pkgs <- c(
  'DadaNanoBC',
  'yaml',
  # for HTML report:
  'rmarkdown',
  'bookdown',
  'pander',
  'dplyr',
  'ggplot2',
  'patchwork',
  'xfun'
)

tar_option_set(
  packages = pkgs,
  controller = crew_controller_group(single_ctrl, multi_ctrl),
  resources = tar_resources(crew = tar_resources_crew(controller = 'full')),
  seed = 42
)

unlist(list(
  # Read sample and primer information and initialize configuration
  tar_target(
    config_file,
    file.path(analysis_dir, 'config.yaml'),
    format = 'file'
  ),
  tar_target(
    config,
    init_config(config_file, analysis_dir)
  ),
  lapply(
    .config_sections,
    function(section) {
      tar_target_raw(
        paste0('config_', section),
        substitute(config[[sec]], list(sec = section))
      )
    }
  ),
  tar_target(
    meta_file,
    file.path(analysis_dir, 'meta.xlsx'),
    format = 'file'
  ),
  tar_target(
    sample_tab,
    read_xlsx_sample_tab(meta_file, 'sample_list')
  ),
  tar_target(
    amplicons,
    read_xlsx_primer_tab(meta_file, 'primers',
                         amplicons = levels(sample_tab$amplicon))
  ),
  # Do primer search and demultiplexing
  tar_target(
    reads_file,
    file.path(analysis_dir, 'reads.fastq.gz'),
    format = 'file'
  ),
  tar_target(
    trim_demux,
    run_with_cfg(
      do_trim_demux,
      reads_file,
      amplicons,
      sample_tab,
      out_dir = file.path(tmp_dir, 'demux'),
      cores = n_workers,
      cfg = config_demultiplex
    )
  ),
  # DADA2 learnErrors() on a (pseudo-)random selection
  tar_target(
    dada_err,
    dada_learn_errors(
      sample(na.omit(trim_demux$seq_tab$reads_path)),
      omega_a = config_cluster$omega_a[1] %||% 1e-20,
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
      aln_out = config_alignment_dir,
      cores = 1,
      inner_fn = infer_barcode,
      cfg = config_cluster,
      exclude = 'max_contam_sample_depth'
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
      config_taxdb,
      taxdb_dir,
      cores = n_workers
    )
  ),
  # Re-cluster contaminated samples with more reads if necessary
  tar_target(
    recluster_contam,
    recluster_contaminated(
      taxonomy$seq_tab,
      dada_err,
      cluster_cfg = config_cluster,
      taxdb_cfg = config_taxdb,
      taxdb_dir = taxdb_dir,
      aln_dir = config_alignment_dir,
      max_sample_depth = config_cluster$max_sample_depth,
      max_contam_sample_depth = config_cluster$max_contam_sample_depth,
      amp_summary_ranks = taxonomy$amp_summary_ranks,
      cores = n_workers
    )
  ),
  # output files
  tar_render(
    rmd_report,
    system.file('report.Rmd', package = 'DadaNanoBC'),
    file.path(normalizePath(analysis_dir), 'report.html')
  ),
  tar_target(xl_report, {
    out <- file.path(analysis_dir, 'report.xlsx')
    run_with_cfg(create_excel_report,
                 recluster_contam,
                 out,
                 cfg = config_output$report)
    out
  }, format = 'file'),
  tar_target(
    combined_bam_top,
    if ('top_combined' %in% config_output$alignments)
      do_combine_alignments(
        recluster_contam,
        config_alignment_dir,
        outdir = file.path(analysis_dir, 'top_alignments'),
        top_only = TRUE
      ),
    format = 'file'
  ),
  tar_target(
    combined_bam,
    if ('combined' %in% config_output$alignments)
      do_combine_alignments(
        recluster_contam,
        config_alignment_dir,
        outdir = file.path(analysis_dir, 'alignments'),
        top_only = FALSE
      ),
    format = 'file'
  )
), recursive = FALSE)
