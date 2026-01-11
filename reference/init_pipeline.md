# Scaffold the targets Pipeline

Copies the pipeline code (`_targets.R`) to current directory (or `path`)
and initializes an analysis directory (without overwriting files).

## Usage

``` r
init_pipeline(path = ".", bash = FALSE, analysis_dir = "analysis")
```

## Arguments

- path:

  Where to create the `_targets.R` and other files

- bash:

  (logical) if `TRUE`, initialize in "Bash mode", which means that an
  `infer_barcodes` Bash script is copied along with `_targets.R`

- analysis_dir:

  Analysis directory
