# Learn error rates

Calls
[dada2::learnErrors](https://rdrr.io/pkg/dada2/man/learnErrors.html)
with settings adjusted for Nanopore data and the clustering procedure
([infer_barcode](https://markschl.github.io/DadaNanoBC/reference/infer_barcode.md))

## Usage

``` r
dada_learn_errors(fq_paths, omega_a = 1e-20, cores = 1, ...)
```
