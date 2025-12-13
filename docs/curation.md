---
hide:
  - toc
---

# Sequence curation

The results of the automated analysis are in *ONT_analysis/<run name\>/output/report.xlsx*. There is one row per barcode, the columns are explained by comments in the document.

The **issues** list should be investigated more closely, although not all of them are critical. The list can be filtered.

![report overview](img/report-overview.png)


## Final sequence set

The *seq1*, [*seq2*, etc.] columns in the **Curation** section contain the final sequences. In case of problems, they can be edited. Changed sequences appear in <span style="color:blue">blue</span>

> FASTA sequences can be exported by a simple copy&paste from the *FA1* [*FA2*, etc.] columns to a text file.

## Common problems

### Inconsistent taxonomy

The **tax-mimatch** issue indicates a mismatch between the taxonomic name in the *morphospecies* entry and the sequence-based identification. Possible reasons:

- mis-identified specimen
- something got mixed up when preparing specimens/extracting/amplifying
- the automatic taxonomy assignment of barcode sequences did not work correctly, or the taxonomic lineages obtained from GBIF are not correct
- contamination not automatically recognized ([how it works](workflow.md#taxonomic-assignments))

It is advisable to check the *details* worksheet: Copy the sample barcode combination (e.g. *bc018-bc183*) and filter the *details* table to list *all* sequences from this sample. The *group* column indicates the taxon cluster (the top being the one shown in the main table).

> If the top taxon group is not the correct one, copy from the "sequence" field(s) to "seq1" ("seq2", ...) in the main table (paste values).

> The **contamination** issue indicates that contamination is very likely present and was automatically down-ranked. It may also be good to investigate these samples more precisely in the *details* tab.

### Ambiguities

If the alignment consensus has ambiguous positions (*ambig-consensus* issue), these may be inspected by clicking on *→ data* to open a folder with alignment files. Use a sequence viewer of your choice, e.g. Geneious, CLC Workbench or UGENE.

The *...seq_comparison.fasta* file may be inspected first. Inconsistencies should be readily visible (example alignment view in Geneious).

![sequence comparison](img/seq-comparison.png)

The BAM file (and the corresponding reference consensus sequences) offer additional insights:

![aln-detail](img/aln-detail.png)

In this case, a long homopolymer stretch results in an unclear sequence.

> Consider copying the "group1_seq1_dominant_seq" (= most abundant sequence in sample) to the *seq1* column in the *Curation* section (without gaps).

## List of issues

-   **ambig-consensus:** At least one of the top taxon consensus
    sequences have ambiguities, suggesting either unresolved sequence
    variation, or sequencing errors e.g. in long homopolmyer regions. To
    investigate, click on "→ data"
    to open a folder with alignment files.

-   **[consensus-diffs]:** At least one of the top taxon consensus
    sequences does not match the dominant unique sequence. This is
    not necessarily an issue (the consensus is always reported),
    but it can happen with only few and/or noisy sequences.

-   **tax-mismatch:** Inconsistency found between the name in the
    *morphospecies* column and the auto-assigned taxon (see also "matching
    ranks"). Strong inconsistencies might also indicate a preparation
    error or contamination, although contamination is usually
    automatically recognized

-   **contamination:** The most abundant taxon is suspected to be a
    contamination and was therefore ignored (still listed in the
    "details" sheet).

-   **known-seq-diffs:** There is at least one mismatch between the top
    taxon and the provided "known" sequence. To investigate, click on "→
    data" to open a folder containing ...seq_comparison.fasta

-   **known-seq-contamination:** The provided "known" sequence appears
    to be a contamination itself

-   **many-variants:** Strong sequence variability for the top taxon (\>
    4 abundant variants above the frequency threshold). Maybe multiple
    individuals are
    present in the DNA mix?

-   **low-coverage:**
    \<20 Nanopore sequences support the top taxon; Low-coverage samples
    may have more ambiguities, and there is a higher risk for errors
    (possibly investigate the BAM alignments)

-  **< 3 times overabundant**:
     The top taxon is less than 3x more abundant than the second taxon.
     It may be worth to check for contamination.

-  **Related taxa in mix**:
     The taxonomy assignment suggests that there are several related taxa
     present in the same sample. It may be worth to check for contamination
     to be sure that the top taxon is actually the one you are looking for.

-  **Re-mapping to consensus gives another consensus!**:
     During the clustering workflow, consensus sequences were not identical
     before and after a re-mapping of reads. Consider inspecting as if this
     was an ambiguous sequence.
     *This should be better handled in the clustering procedure*

  