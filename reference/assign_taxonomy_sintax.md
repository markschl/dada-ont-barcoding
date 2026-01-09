# Run taxonomy assignment (SINTAX algorithm)

This requires [VSEARCH](https://github.com/torognes/vsearch) to be
installed.

## Usage

``` r
assign_taxonomy_sintax(
  seq_file,
  utax_db,
  confidence_threshold = 0.8,
  tmp_prefix = NULL,
  threads = 1
)
```
