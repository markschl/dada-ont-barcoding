
install_sequtil_deps <- function() {
  p <- installed.packages()[,c('Package')]
  install.packages(setdiff('BiocManager', p))
  p.bioc = setdiff(c('Biostrings', 'DECIPHER'), p)
  if (length(p.bioc) > 0)
    BiocManager::install(p.bioc)
}

n_ambigs <- function(x) {
  nchar(gsub('[ACGT]', '', as.character(x)))  
}

limit_sample <- function(x, limit) {
  if (length(x) > limit) sample(x) else x
}


#' Get edit distance for best-matching strand
match_both_orient <- function(seq1, seq2, cores = 1) {
  both_orient <- list(seq1, as.character(Biostrings::reverseComplement(Biostrings::DNAStringSet(seq1))))
  diffs <- simplify2array(lapply(both_orient, function(s1) {
    pairwise_align(s1, seq2, count_end_gaps=FALSE, cores = cores)[,'diffs']
  }), except=NA)
  apply(diffs, 1, min)
}

write_dna <- function(seqs, seq_file) {
  stopifnot(!is.null(names(seqs)))
  Biostrings::writeXStringSet(Biostrings::DNAStringSet(seqs), seq_file)
}

read_dna <- function(seq_file) {
  Biostrings::readDNAStringSet(seq_file)
}

#' Does a pairwise alignment between pairs of sequences and reports
#' the number of mismatches and pattern/subject gaps.
#' Sets all these values to Inf if the overlap length relative to 'seq1' 
#' is lower than 'min_overlap'
pairwise_align <- function(seq1,
                          seq2,
                          count_end_gaps = TRUE,
                          min_overlap = 0.5,
                          cores = 1,
                          simplify = TRUE) {
  seq1 <- Biostrings::DNAStringSet(seq1)
  seq2 <- Biostrings::DNAStringSet(seq2)
  aln <- DECIPHER::AlignPairs(
    seq1, seq2,
    # penalties chosen to avoid alignments with very short overlap
    # in case of totally non-matching sequences
    perfectMatch=4, misMatch=-4, gapOpening=-4, gapExtension=-1,
    # helps avoiding mis-alignments
    bandWidth=length(seq1[[1]]),
    processors=cores, verbose=F
  )
  if (!count_end_gaps) {
    # filter end gaps
    # TODO: better way?
    aln$PatternGapLength = lapply(1:nrow(aln), function(i) {
      aln$PatternGapLength[[i]][aln$PatternGapPosition[[i]] > 1 & aln$PatternGapPosition[[i]] <= aln$PatternEnd[i]]
    })
    aln$SubjectGapLength = lapply(1:nrow(aln), function(i) {
      aln$SubjectGapLength[[i]][aln$SubjectGapPosition[[i]] > 1 & aln$SubjectGapPosition[[i]] <= aln$SubjectEnd[i]]
    })
  }
  out <- lapply(1:nrow(aln), function(i) {
    a <- aln[i,]
    seq1_coverage <- (a$Matches + a$Mismatches) / length(seq1[[i]])
    out <- if (seq1_coverage >= min_overlap) {
      c(mismatches = a$Mismatches[1], 
        seq1_gaps = sum(a$PatternGapLength[[1]]),
        seq2_gaps = sum(a$SubjectGapLength[[1]]))
    } else {
      c(mismatches=Inf, seq1_gaps=Inf, seq2_gaps=Inf)
    }
    out['gaps'] = out['seq1_gaps'] + out['seq2_gaps']
    out['diffs'] = out['mismatches'] + out['gaps']
    out
  })
  out <- simplify2array(out, except=if (simplify) c(0L, 1L) else NA)
  if (!is.null(dim(out))) {
    out <- t(out)
  }
  out
}


make_fasta <- function(seqs, nm=names(seqs)) {
  stopifnot(!is.null(nm))
  sprintf('>%s\n%s\n', nm, seqs)
}

