#!/usr/bin/env bash

set -euo pipefail

#### Set options ####

outdir=
overwrite=false
primer_mismatch_rate=0.2
bc_max_mismatches=1
concat_search_err_rate=0.2
min_len=50
compr=zst
st=st
ncores=1

usage="
De-multiplexing of dual-indexed (Nanopore) reads

Requires https://github.com/markschl/seqtool (path can be provided with -s)

$(basename "$0") [-h] [options] reads.fastq[.gz] seq_prefix, reads.fastq[.gz|bz2|lz4|zstd], ...

    -h                show this help text.
    -o                the output directory containing the demultiplexed and 
                      discarded reads [default = <input>_trim]
    -l                the minimum read length [default = $min_len]
    -p                the maximum mismatch rate allowed in primers [default = $primer_mismatch_rate]
    -b                the maximum number of mismatches in barcodes for reads to be
                      included in <outdir>/trimmed.fastq[.gz|...] [default = $bc_max_mismatches]
    -u                the maximum mismatch rate for the primer search in trimmed reads
                      (to remove concatenated products) [default = $concat_search_err_rate]
    -c                the compression format for the output files [default = $compr].
                      Available are: 'none', 'gz', 'bz2', 'lz4', 'zstd'
    -w                Overwrite the output files, forcing re-calculation of all steps [default = false]
    -t                number of threads to use for seqtool [default = $ncores]
                      note that the effective number is slightly higher, since 
                      multiple commnands are piped
    -s                path to the seqtool binary [default = $st]
    seq_prefix        Path prefix to primer and sample index files in FASTA format
                      (<prefix>fwd_primers.fasta, <prefix>fwd_primers.fasta,
                      <prefix>fwd_idx.fasta, <prefix>rev_idx.fasta)
    reads.fastq, ...  One or more FASTQ file(s), optionally compressed

Output directory structure:
outdir/
  trimmed.fastq[.gz|...]             demultiplexed reads
  no_primer.<fwd/rev>.fastq[.gz|...] reads without forward/reverse primers
  no_index.<fwd/rev>.fastq[.gz|...]  reads without forward/reverse sample indexes
  too_short.fastq[.gz|...]           reads shorter than $min_len
  concatenated.fastq[.gz|...]        putative concatenated artifacts (primers found in trimmed reads)
  trim_counts.tsv                    Statistics: read file, count
  length_stats.tsv                   Statistics: read file, length, count
  primer_positions.tsv               Statistics: primer position, count
  qual_stats.tsv                     Statistics: read quality, primer mismatches, barcode mismatches
                                     (without concatenated and too short reads)
"

while getopts "ho:l:p:b:u:c:wt:s:" opt; do
    case "$opt" in
    h)
        echo "$usage" >&2
        exit 0
        ;;
    o) outdir=$OPTARG ;;
    l) min_len=$OPTARG ;;
    p) primer_mismatch_rate=$OPTARG ;;
    b) bc_max_mismatches=$OPTARG ;;
    u) concat_search_err_rate=$OPTARG ;;
    c) compr=$OPTARG ;;
    w) overwrite=true ;;
    t) ncores=$OPTARG ;;
    s) st=$OPTARG ;;
    *)
        echo "Invalid option: -$OPTARG" >&2
        exit 1
        ;;
    esac
done

shift $((OPTIND - 1))

if [ "$#" -lt 2 ]; then
    echo "Incorrect number of arguments (see -h)" >&2
    exit 1
fi

seq_prefix="$1" && shift
fastq_input="$@"

fwd_primers="$seq_prefix"fwd_primers.fasta
rev_primers="$seq_prefix"rev_primers.fasta
fwd_bc="$seq_prefix"fwd_idx.fasta
rev_bc="$seq_prefix"rev_idx.fasta

if [[ -z "$outdir" ]]; then
    outdir="${fastq_input%.*}_trim"
fi

fq=fastq
if [[ $compr != "none" ]]; then
    fq=$fq.$compr
fi
trimmed_out="$outdir/trimmed.$fq"


#### Setup ####

# TODO: balance parallel calculations
# adjust ncores to avoid using too many resources
ncores=$(($ncores >= 2 ? $ncores - 1 : $ncores))

# max. barcode length per direction (usually all should have the same length)
f_bc_len=$($st . --to-tsv seqlen "$fwd_bc" | sort -rnu | head -n1)
r_bc_len=$($st . --to-tsv seqlen "$rev_bc" | sort -rnu | head -n1)

# final possible range within which the barcode can occur is
# the barcode length + the maximum allowed edit distance
f_bc_len=$(($f_bc_len + $bc_max_mismatches))
r_bc_len=$(($r_bc_len + $bc_max_mismatches))

# create reverse complement of primers and barcodes

rev_bc_rc="$outdir/_rev_bc_rc.fasta"
fwd_primers_rc="$outdir/_fwd_primers_rc.fasta"
rev_primers_rc="$outdir/_rev_primers_rc.fasta"

mkdir -p "$outdir"
$st revcomp "$rev_bc" -o "$rev_bc_rc"
$st revcomp "$fwd_primers" -o "$fwd_primers_rc"
$st revcomp "$rev_primers" -o "$rev_primers_rc"

#### Do the searching ####

if [[ ! $overwrite && -f "$trimmed_out" ]]; then
    echo "Output file already exists, skipping" >&2
    exit 0
fi

# find the primers in the forward orientation and trim,
# but still include up- and downstream containing the barcodes

echo "Searching for primers in forward orientation..." >&2

$st set -i {seq_num} $fastq_input | # create short headers to save space
    $st del -d --fq |
    $st filter --fq "seqlen >= $min_len" --dropped "$outdir/too_short.$fq" |
    $st find "file:$fwd_primers" \
        -R $primer_mismatch_rate \
        --in-order \
        -t $ncores --seqtype dna --fq \
        -f --dropped "$outdir/_no_fwd_primer.fastq" \
        -a f={match_start} -a fl={match_len} -a fm={match_diffs} |
    $st find "file:$rev_primers_rc" \
        -R $primer_mismatch_rate \
        --in-order \
        -t $ncores --seqtype dna --fq \
        -f --dropped "$outdir/no_primer.rev.$fq" \
        -a r={match_neg_end} -a rl={match_len} -a rm={match_diffs} |
    $st trim --fq "{num(attr('f')) - $f_bc_len}:{num(attr('r')) + $r_bc_len}" \
        >"$outdir/primers_trimmed.fastq"

echo "Searching for primers in reverse orientation..." >&2

# -> primers_trimmed.fastq always contains reads in the forward orientation
$st find "file:$rev_primers" "$outdir/_no_fwd_primer.fastq" \
    -R $primer_mismatch_rate \
    --in-order \
    -t $ncores --seqtype dna --fq \
    -f --dropped "$outdir/_no_rev_primer_rc.$fq" \
    -a r={match_start} -a rl={match_len} -a rm={match_diffs} |
    $st find "file:$fwd_primers_rc" \
        -R $primer_mismatch_rate \
        --in-order \
        -t $ncores --seqtype dna --fq \
        -f --dropped "$outdir/no_primer.fwd.$fq" \
        -a f={match_neg_end} -a fl={match_len} -a fm={match_diffs} |
    $st trim --fq "{num(attr('r')) - $r_bc_len}:{num(attr('f')) + $f_bc_len}" |
    $st revcomp --fq \
        >>"$outdir/primers_trimmed.fastq"

rm "$outdir/_no_fwd_primer.fastq"
cat "$outdir/_no_rev_primer_rc.$fq" >> "$outdir/no_primer.rev.$fq"
rm "$outdir/_no_rev_primer_rc.$fq"

echo "Searching for sample indexes..." >&2

all_primers="$outdir/_all_primers.fasta"
cat "$fwd_primers" "$rev_primers" "$rev_primers_rc" "$fwd_primers_rc" >"$all_primers"

$st find "file:$fwd_bc" "$outdir/primers_trimmed.fastq" \
    --rng ":$f_bc_len" --anchor-end 0 \
    -D $bc_max_mismatches \
    -t $ncores --seqtype dna \
    -f --dropped "$outdir/no_index.fwd.$fq" \
    -a fi={pattern_name} \
    -a fi.m={match_diffs} -a fi2.m='{match_diffs(1,2)}' |
    $st find "file:$rev_bc_rc" \
        --rng "-$r_bc_len:" --anchor-start 0 \
        -D $bc_max_mismatches \
        -t $ncores --seqtype dna --fq \
        -f --dropped "$outdir/no_index.rev.$fq" \
        -a ri={pattern_name} \
        -a ri.m={match_diffs} -a ri2.m='{match_diffs(1,2)}' |
    $st trim -e --fq "{$f_bc_len + num(attr('fl'))}:{-(num(attr('rl')) + $r_bc_len)}" |
    $st find -e "file:$all_primers" \
        -R $concat_search_err_rate \
        --dropped "$outdir/concatenated.$fq" \
        -t $ncores --seqtype dna --fq |
    $st filter --fq 'opt_attr("fi2.m") === undefined && opt_attr("ri2.m") === undefined' \
        --dropped "$outdir/ambig_index.$fq" |
    $st filter --fq "seqlen >= $min_len" --dropped "$outdir/_too_short2.$fq" \
        -o "$trimmed_out"

cat "$outdir/_too_short2.$fq" >> "$outdir/too_short.$fq"
rm "$outdir/_too_short2.$fq"

#### Statistics ####

echo "Assembling read statistics..." >&2

echo -e "file\tcount" >"$outdir"/trim_counts.tsv
$st count \
  -k filename \
  "$outdir"/no_primer.fwd.$fq \
  "$outdir"/no_primer.rev.$fq \
  "$outdir"/no_index.fwd.$fq \
  "$outdir"/no_index.rev.$fq \
  "$outdir"/concatenated.$fq \
  "$outdir/too_short.$fq" \
  "$outdir/ambig_index.$fq" \
  "$trimmed_out" >>"$outdir"/trim_counts.tsv

header="exp_err\tf_primer_mis\tr_primer_mis\tf_idx_mis\tr_idx_mis\tf_idx2_mis\tr_idx2_mis\tcount"
echo -e "$header" >"$outdir"/qual_stats.tsv
$st count \
  -k 'bin(exp_err, 1)' \
  -k 'opt_attr(fm)' -k 'opt_attr(rm)' \
  -k "opt_attr('fi.m')" -k "opt_attr('ri.m')" \
  -k "opt_attr('fi2.m')" -k "opt_attr('ri2.m')" \
  "$outdir"/no_primer.fwd.$fq \
  "$outdir"/no_primer.rev.$fq \
  "$outdir"/no_index.fwd.$fq \
  "$outdir"/no_index.rev.$fq \
  "$outdir/ambig_index.$fq" \
  "$trimmed_out" \
  >>"$outdir/qual_stats.tsv"

# TODO: 2 entries for same position can exist
echo -e "dir\tpos\tcount" >"$outdir/primer_pos.tsv"
for dir in fwd rev; do
    d=${dir:0:1}
    $st count "$trimmed_out" -k "attr($d)" | 
    awk -F'\t' -v OFS='\t' '{if ($1 > 0) {print "fwd",$1,$2} else {print "rev",-$1,$2}}' \
    >> "$outdir/primer_pos.tsv"
done

echo -e "file\tlength\tcount" > "$outdir/length_stats.tsv"
$st count \
  "$trimmed_out" \
  "$outdir/concatenated.$fq" \
  -k filename -k seqlen \
  >> "$outdir/length_stats.tsv"

# clean up
rm -f "$outdir/primers_trimmed.fastq" \
  "$fwd_primers" "$fwd_primers_rc" \
  "$rev_primers" "$rev_primers_rc" \
  "$all_primers" \
  "$fwd_bc" "$rev_bc" "$rev_bc_rc"
