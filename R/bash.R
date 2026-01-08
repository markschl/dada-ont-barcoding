
.env <- new.env(parent = emptyenv())
.env$bash_present <- NULL
.env$programs <- list()

require_bash <- function() {
  if (Sys.which('bash') == '') {
    stop('The Bash shell was not found. On Windows, WSL2 is needed (https://learn.microsoft.com/en-us/windows/wsl/install)')
  }
  .env$bash_present <- TRUE
}

# Run a command in Bash, which should work on Linux, OS-X and Windows WSL
# **Note**: assuming **no** spaces in arguments (no spaces in paths!)
run_bash <- function(cmd, stdout = '', stderr = '', ...) {
  if (is.null(.env$bash_present)) {
    require_bash()
  }
  cmd <- paste(cmd, collapse=' ')
  full_cmd <- paste0('set -euo pipefail; ', cmd)
  rv <- system2('bash', c('-c', shQuote(full_cmd)), stdout = stdout, stderr = stderr, ...)
  code <- if (isTRUE(stdout) || isTRUE(stderr)) {
    attr(rv, 'status') %||% 0
  } else {
    rv
  }
  if (code != 0) {
    stop('Command failed: ', cmd)
  }
  rv
}

run_bash_script <- function(script, args, ...) {
  script <- system.file('bash_scripts', script, package = 'DadaNanoBC')
  stopifnot(nzchar(script))
  run_bash(c(script, args), ...)
}

#' Set the path to one or more programs
#'
#' @examples
#' set_program_path(samtools = '/path/to/samtools', minimap2 = '/path/to/minimap2')
#'
#' @export
set_program_path <- function(...) {
  dots <- list(...)
  if (length(dots) == 0)
    return()
  names(dots) <- paste0('DadaNanoBC', names(dots))
  do.call(options, dots)
  invisible(NULL)
}

find_program <- function(program, full_name = NULL) {
  path <- .env$programs[[program]]
  if (!is.null(path)) {
    return(path)
  }
  full_name <- full_name %||% program
  path <- getOption(
    paste0('DadaNanoBC.', full_name),
    Sys.which(program)
  )
  if (!nzchar(path)) {
    path <- Sys.getenv(paste0('DadaNanoBC_', full_name))
    if (!nzchar(path)) {
      return(NULL)
    }
  }
  .env$programs[[program]] <- path
  path
}

get_program <- function(program, full_name=NULL) {
  path <- find_program(program, full_name = full_name)
  if (is.null(path))
    stop(sprintf(paste(
      "The program '%s' was not found.",
      "It needs to be installed either system-wide and visible to R (in $PATH), ",
      "or its path can be registered",
      "(run DadaNanoBC::check_system_requirements() for more information)."
    )), full_name)
  path
}

#' Check external software requirements and report missing programs
#'
#' @export
check_system_requirements <- function() {
  programs <- list(
    samtools = list(name = 'samtools', url = 'https://www.htslib.org'),
    minimap2 = list(name = 'minimap2', url = 'https://lh3.github.io/minimap2'),
    VSEARCH = list(
      name = 'vsearch',
      description = 'VSEARCH',
      url = 'https://github.com/torognes/vsearch'
    ),
    seqtool = list(
      name = 'st',
      url = 'https://github.com/markschl/seqtool',
      full_name = 'seqtool',
      description = 'seqtool (v0.4 or higher)',
      download = 'https://github.com/markschl/seqtool/releases'
    )
  )

  paths <- lapply(programs, function(p) find_program(p$name, full_name = p$full_name))
  missing <- sapply(paths, is.null)

  if (any(!missing)) {
    cat('Found:', sapply(names(paths)[!missing], function(p) {
      info <- programs[[p]]
      sprintf('%s: %s', info$description %||% info$full_name %||% p, paths[p])
    }), sep = '\n')
  }

  if (any(missing)) {
    cat(
      '',
      'Not found:',
      sapply(names(paths)[missing], function(p) {
        info <- programs[[p]]
        sprintf(
          '%s: %s%s',
          info$full_name %||% p,
          info$url,
          if (!is.null(info$download))
            sprintf(' (download from %s)', info$download)
          else
            ''
        )
      }),
      "",
      "All software needs to be visible to R (in $PATH), or run Sys.setenv(DadaNanoBC_<tool> = 'path/to/tool'). ",
      "`set_program_path('program', 'path/to/program') works as well within the same R instance ",
      "(does not work in `targets` pipelines)",
      sep = '\n'
    )
    stop('Some programs not found!')
  }
}
