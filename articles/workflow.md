# Workflow steps

## Pipeline steps

1.  Search the primers and short sample indexes (located up-/downstream
    of the primer sequences) to group them by sample (de-multiplex);
    remove low-quality reads
    ([`do_trim_demux()`](https://markschl.github.io/DadaNanoBC/reference/do_trim_demux.md))
2.  Infer the barcode sequences using a sequential procedure based on
    [DADA2](https://benjjneb.github.io/dada2) and/or fixed-threshold
    clustering (with
    [DECIPHER](https://www2.decipher.codes/Clusterize.html)); compare
    them with already known sequences if present
    ([`infer_barcode()`](https://markschl.github.io/DadaNanoBC/reference/infer_barcode.md)
    /
    [`do_infer_all_barcodes()`](https://markschl.github.io/DadaNanoBC/reference/do_infer_all_barcodes.md))
3.  Auto-assign the taxonomy and compare with the name given by
    morphological identification (if present) to validate and check for
    contamination
    ([`do_assign_compare_taxonomy()`](https://markschl.github.io/DadaNanoBC/reference/do_taxonomy.md))
4.  Export summary table, which can be further manually curated
    ([`create_excel_report()`](https://markschl.github.io/DadaNanoBC/reference/create_excel_report.md))

## Primer search

It is assumed that sample-specific tags are attached to primers on both
sides (dual indexing):

TAG1FPRIMERBARCODESEQUENCERPRIMERTAG2

**Combining multiple amplicons** is possible. Primers are searched in
the order that amplicons appear in the sample sheet. Therefore, if there
are *nested amplicons*, the shorter one should be placed *at the end* in
the sheet.

## Clustering

The basic procedure for the clustering is implemented by the
[`infer_barcode()`](https://markschl.github.io/DadaNanoBC/reference/infer_barcode.md)
function. The procedure starts with DADA2 clustering, followed by
consensus building. Additional steps may follow depending on whether the
consensus is unambiguous (each alignment column supported by enough
identical bases), or not. The workflow is illustrated in this flowchart.
Relevant configuration options are shown in *grey*.

![](img/cluster-flowchart.png)

### Steps/strategies

#### DADA2 denoising

[DADA2](https://benjjneb.github.io/dada2) is a popular program for
inferring sequence variants (ASV) in amplicon sequence data. The program
requires that there are least a few error-free reads present. If there
is not enough duplication, fixed-threshold clustering (see
[below](#fixed-threshold-clustering)) is applied instead.

![](img/cluster-dada.png)

After obtaining the DADA2 ASV(s), the next step is to infer a [consensus
sequence](#how-is-the-consensus-sequence-determined-) for the ASV, which
should usually be identical with the ASV sequence. Ambiguous bases in
the consensus sequence may indicate unresolved sequence variation, which
is further addressed by attempting [haplotype
splitting](#haplotype-splitting) or by more sensitive denoising (see
[flowchart](#flowchart)).

#### Haplotype splitting

For DADA2 ASVs with an ambiguous consensus, it is checked if there are
**two dominant sequences** (supported by a few identical sequences)
serving as references for splitting the sequencing reads into two
similarly abundant sub-clusters with (ideally) *unambiguous* sequences.
The splitting is done if the total number of ambiguous bases in the two
resulting consensus sequences is smaller than the number of ambiguous
bases in the original consensus of all combined reads (see also
non-public function
[`try_split_haplotypes()`](https://markschl.github.io/DadaNanoBC/reference/try_split_haplotypes.md)).

This simple procedure is often successful and avoids having to re-run
DADA2 with higher sensitivity, or it may alo resolve InDel variation
that DADA2 did not discover.

![](img/split-haplo.png)

#### Fixed-threshold clustering

Fixed-threshold clustering is applied to low-coverage samples with
insufficient sequence duplication (see `dada_min_identical` and
`dada_min_n0` settings in
[`infer_barcode()`](https://markschl.github.io/DadaNanoBC/reference/infer_barcode.md).

By default, single-linkage is applied with a 97% identity threshold
(with [DECIPHER](https://www2.decipher.codes/Clusterize.html)). This can
result in an ambiguous consensus sequence, which can still be manually
inspected.

![](img/cluster-fixed.png)

### How is the consensus sequence determined?

After DADA2 denoising, haplotype splitting and fixed-threshold
clustering, the raw (quality-filtered) reads are aligned to the
representative sequence (ASV or dominant unique sequence) with
[Minimap2](https://lh3.github.io/minimap2).

The consensus sequence is then obtained with [samtools
consensus](https://www.htslib.org/doc/samtools-consensus.html) using a
simple **frequency-based** approach. No elaborate neural network-based
inference (see e.g. [Medaka](https://github.com/nanoporetech/medaka)) is
currently applied.

By default, each base in the consensus sequence needs to be supported by
at least 65% of the reads (weighted by quality scores) to be
*unambiguous*. Otherwise, an [ambiguous
base](https://en.wikipedia.org/wiki/Nucleic_acid_notation) is shown,
which may indicate unresolved sequence variation.

The correctness of the consensus is further ensured by re-mapping the
reads against the consensus sequence. If necessary, this step is
repeated until *samtools consensus* returns the same sequence in two
successive rounds.

The *final reported sequence* is always the consensus sequence (see also
[`vignette("curation")`](https://markschl.github.io/DadaNanoBC/articles/curation.md)).
If the barcode sample has a sufficient read depth and sequencing
quality, it should usually be unambiguous and correct. For low-depth
samples with ambiguities, alignments might be inspected manually.

#### Homopolymers

The length of **long homopolymer** repeats is often not very consistent
between different sequencing reads, which can result in `N` bases in the
consensus. In this case, a simple frequency-based approach is used to
“fix” this situation: the homopolymer length is taken from the *most
abundant sequence* (see
[`fix_homopolymers()`](https://markschl.github.io/DadaNanoBC/reference/fix_homopolymers.md)).

### Taxa grouping

Sequence variation *within* taxa (e.g. haplotypes) is distinguished from
*between-taxa* variation with simple fixed-threshold single-linkage
clustering done with
[DECIPHER](https://www2.decipher.codes/Clusterize.html) (default: 97%
threshold). This usually works well at least with the ITS barcode.

> The most problematic case would be if there is a closely related
> **contaminant** taxon in the DNA extract, which may not be well
> distinghished from the target taxon.

## Taxonomic assignments

Taxonomic names are automatically assigned with the the [SINTAX
algorithm](https://doi.org/10.1101/074161) using a reference database
(see
[`do_assign_compare_taxonomy()`](https://markschl.github.io/DadaNanoBC/reference/do_taxonomy.md)).

Second, taxonomic names associated with the specimens are compared to
the sequence-based taxonomic assignments. From this, a **taxonomic
overlap** is calculated, which is the *proportion of matching ranks in
the taxonomic lineages* (as far as names are defined at these levels).
The lineages are retrieved from the [GBIF
backbone](https://doi.org/10.15468/39omei) taxonomy.

Inconsistent taxonomic labels may indicate errors during the preparation
or identification (or possibly the sequence-based identification), or
the presence of contamination. Potential **dominant contaminants** are
automatically flagged (see
[`do_assign_compare_taxonomy()`](https://markschl.github.io/DadaNanoBC/reference/do_taxonomy.md),
*details* section).
