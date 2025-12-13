#!/usr/bin/env Rscript

run_name = commandArgs(trailingOnly=T)
stopifnot(length(run_name) == 1)

rmarkdown::render(
  'analysis.Rmd', 
  output_format = 'bookdown::html_document2',
  output_dir = file.path('analysis', run_name),
  params = list(run_name = run_name)
)
