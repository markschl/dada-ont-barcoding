# Download a taxonomic database and save it in UTAX format

Download a taxonomic database and save it in UTAX format

## Usage

``` r
load_taxdb(
  db_url,
  db_type = c("unite", "qiime_fasta", "qiime_qza"),
  db_file,
  db_tax_url = NULL,
  ...
)
```

## Arguments

- db_url:

  Character vector with one or more file URLs; sequences from all URLs
  are combined into one final database.

- db_type:

  Database type (character) (see details)

- db_file:

  Output file (will be GZIPped FASTA, `.fasta.gz`)

- db_tax_url:

  (optional) URL(s) pointing to file(s) containing the taxonomy for the
  sequences in the files at `db_url` (currently implemented for
  *qiime_qza* format)

## Details

Currently, three database types are implemented (more may follow):

- *"unite"*: UNITE ITS database (https://unite.ut.ee/repository.php):
  choose a dataset (e.g. [all
  Eukaryotes](https://doi.plutof.ut.ee/doi/10.15156/BIO/3301243)), click
  on the file at *Downloads* and copy the actual file URL once the file
  starts downloading in the browser. The SH sequence set with a *dynamic
  clustering threshold* is chosen from the archive and imported.

- *"qiime_fasta"*: FASTA file with lineages in headers in the form
  `<seqid> rank1__name; rank2__name; ...` whereby `rank1` and `rank2`
  are one-letter codes such as `k` (for kingdom), `g` (for genus), etc.
  If only sequences are in the file, specify the URL of a tab-delimited
  file mapping sequence IDs to lineages. One example of database
  providing such files [GTDB](https://gtdb.ecogenomic.org) (e.g.
  [here](https://data.gtdb.aau.ecogenomic.org/releases/release226/226.0/genomic_files_reps))

- *"qiime_qza"*: [QIIME2
  QZA](https://amplicon-docs.qiime2.org/en/stable/explanations/archives.html)
  format, which is essentially a ZIP file containing a
  *qiime_fasta*-formatted FASTA file, while another QZA file should
  contain the tab-delimited taxonomy (provide with `db_tax_url`).

This function may be expanded and made more flexible in the future
