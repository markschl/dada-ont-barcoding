# Set global options

Convenience function for setting global options such as program paths.

## Usage

``` r
set_global_opts(...)
```

## Details

The following options can be set:

- Required program paths (`samtools`, `minimap2`, `seqtool`, `vsearch`)

- `tmp_dir`: Path to a custom temporary directory for short-lived data
  and not overly large data. On Linux, the default is to place these
  files in a subdirectory of `$XDG_RUNTIME_DIR` (usually
  `/run/user/<userid>`, tmpfs).

Under the hood, this simply does
`options(DadaNanoBC.opt = 'value', ...)`.

The same options can be set as environment variables:
`DadaNanoBC_opt=value` (in R: `Sys.setenv(DadaNanoBC_opt = 'value')`.
`set_global_opt` always has priority over environment variables.

## Examples

``` r
set_global_opts(samtools = '/path/to/samtools', tmp_dir = '/tmp/nanopore-barcoding')
```
