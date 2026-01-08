
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

find_program <- function(program) {
  path <- .env$programs[[program]]
  if (!is.null(path)) {
    return(path)
  }
  path <- .env$programs[[program]] <- getOption(
    paste0('DadaNanoBC.', program),
    Sys.which(program)
  )
  if (!nzchar(path))
    return(NULL)
  path
}

get_program <- function(program, long_name=NULL) {
  path <- find_program(program)
  if (is.null(path))
    stop(sprintf(paste(
      "The program '%s' was not found.",
      "It needs to be installed either system-wide and visible to R (in $PATH),",
      "or its path can be registered with `set_program_path('%s', 'path/to/%s')`."
      # "Consider installing by runing 'scripts/install-deps.sh' in the (Bash) console. ",
      # "For more information, see..."
      # # TODO: reference tutorial
    )), long_name, program, program)
  path
}

check_system_requirements <- function() {
  programs <- list(
    samtools = list(name = 'samtools', url = 'https://www.htslib.org'),
    minimap2 = list(name = 'minimap2', url = 'https://lh3.github.io/minimap2'),
    VSEARCH = list(name = 'VSEARCH', url = 'https://github.com/torognes/vsearch'),
    seqtool = list(
      name = 'st',
      url = 'https://github.com/markschl/seqtool',
      long_name = 'seqtool',
      download = 'https://github.com/markschl/seqtool/releases'
    )
  )

  paths <- lapply(programs, function(p)
    find_program(p$name))
  missing <- sapply(paths, is.null)

  if (any(!missing)) {
    cat('Found:', sapply(names(paths)[!missing], function(p) {
      info <- programs[[p]]
      sprintf('%s: %s', info$long_name %||% p, paths[p])
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
          info$long_name %||% p,
          info$url,
          if (!is.null(info$download))
            sprintf(' (download from %s)', info$download)
          else
            ''
        )
      }),
      "",
      "Missing software needs to be installed either system-wide and visible to R (in $PATH), ",
      "or its path can be registered with `set_program_path('program', 'path/to/program')",
      sep = '\n'
    )
    stop('Some programs not found!')
  }
}
