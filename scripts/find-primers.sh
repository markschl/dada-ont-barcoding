#!/usr/bin/env bash

set -euo pipefail

#### Set options ####

outdir=
overwrite=false
primer_mismatch_rate=0.2
bc_max_mismatches=1
artifact_search_err_rate=0.2
compr=zst
st=st
ncores=1

usage="
De-multiplexing of dual-indexed (Nanopore) reads

Requires https://github.com/markschl/seqtool (path can be provided with -s)

$(basename "$0") [-h] [options] reads.fastq[.gz] seq_prefix, reads.fastq[.gz|bz2|lz4|zstd], ...

    -h                show this help text.
    -o                the output directory containing the demultiplexed and 
                      discarded reads [default = <input>_demux]
    -p                the maximum mismatch rate allowed in primers [default = $primer_mismatch_rate]
    -b                the maximum number of mismatches in barcodes for reads to be
                      included in <outdir>/demux.fastq.$compr [default = $bc_max_mismatches]
    -u                the maximum mismatch rate for the primer search in trimmed reads
                      (to remove concatenated artifacts) [default = $artifact_search_err_rate]
    -f                the compression format for the output files [default = $compr].
                      Available are: 'none', 'gz', 'bz2', 'lz4', 'zstd'
    -w                Overwrite the output files, forcing re-calculation of all steps [default = false]
    -t                number of threads to use for seqtool [default = $ncores]
                      note that the effective number is slightly higher, since 
                      multiple commnands are piped
    -s                path to the seqtool binary [default = $st]
    seq_prefix        Path prefix to primer and sample index files in FASTA format
                      (<prefix>fwd_primers.fasta, <prefix>fwd_primers.fasta,
                      <prefix>fwd_idx.fasta, <prefix>rev_idx.fasta)
    reads.fastq, ...  One or more FASTQ file(s)

Output directory structure:
outdir/
  demux.fastq.$compr                demultiplexed reads
  no_primer.<fwd/rev>.fastq.$compr  reads without forward/reverse primers
  no_index.<fwd/rev>.fastq.$compr   reads without forward/reverse sample indexes
  fusion_artifact.fastq.$compr      possible concatenated artifacts (primers found in trimmed reads)
  read_counts.tsv                   read counts for all the listed FASTQ files
  primer_positions.tsv              primer position statistics
  length_stats.tsv                  TSV file: read file, length, count
  demux_stats.tsv                   read counts for different combinations of
  read_stats.tsv                    read quality / primer mismatches / barcode mismatches
"

# TODO: balance parallel calculations

while getopts "ho:p:b:u:wt:s:" opt; do
    case "$opt" in
    h)
        echo "$usage" >&2
        exit 0
        ;;
    o) outdir=$OPTARG ;;
    p) primer_mismatch_rate=$OPTARG ;;
    b) bc_max_mismatches=$OPTARG ;;
    u) artifact_search_err_rate=$OPTARG ;;
    f) compr=$OPTARG ;;
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
    outdir="${fastq_input%.*}_demux"
fi

fq=fastq
if [[ $compr != "none" ]]; then
    fq=$fq.$compr
fi
demux_out="$outdir/demux.$fq"


#### Setup ####

# adjust ncores to avoid using too many resources
# TODO: not very elegant
ncores=$(($ncores >= 2 ? $ncores - 1 : $ncores))

# barcode length

f_bc_len=$($st . --to-tsv seqlen "$fwd_bc" | uniq)
r_bc_len=$($st . --to-tsv seqlen "$rev_bc" | uniq)

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

if [[ ! $overwrite && -f "$demux_out" ]]; then
    echo "Output file already exists, skipping" >&2
    exit 0
fi

# find the primers in the forward orientation and trim,
# but still include up- and downstream containing the barcodes

echo "Searching for primers in forward orientation..." >&2

$st set -i {seq_num} $fastq_input | # create short headers to save space
    $st del -d --fq |
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
        >$outdir/pre_trim.fastq

echo "Searching for primers in reverse orientation..." >&2

$st find "file:$rev_primers" "$outdir/_no_fwd_primer.fastq" \
    -R $primer_mismatch_rate \
    --in-order \
    -t $ncores --seqtype dna --fq \
    -f --dropped "$outdir/no_primer.rev.rc.$fq" \
    -a r={match_start} -a rl={match_len} -a rm={match_diffs} |
    $st find "file:$fwd_primers_rc" \
        -R $primer_mismatch_rate \
        --in-order \
        -t $ncores --seqtype dna --fq \
        -f --dropped "$outdir/no_primer.fwd.rc.$fq" \
        -a f={match_neg_end} -a fl={match_len} -a fm={match_diffs} |
    $st trim --fq "{num(attr('r')) - $r_bc_len}:{num(attr('f')) + $f_bc_len}" |
    $st revcomp --fq \
        >>"$outdir/pre_trim.fastq"

rm "$outdir/_no_fwd_primer.fastq"

echo "Searching for sample indexes..." >&2

all_primers="$outdir/_all_primers.fasta"
cat "$fwd_primers" "$rev_primers" "$rev_primers_rc" "$fwd_primers_rc" >"$all_primers"

$st find "file:$fwd_bc" "$outdir/pre_trim.fastq" \
    --rng ":$f_bc_len" --anchor-end 0 \
    -D $bc_max_mismatches \
    -t $ncores --seqtype dna \
    -f --dropped "$outdir/no_index.fwd.$fq" \
    -a fbc={pattern_name} \
    -a fbc.m={match_diffs} -a fbc2.m='{match_diffs(1,2)}' | #  -a fbc.start='{match_start}' -a fbc.len='{match_len}'
    $st find "file:$rev_bc_rc" \
        --rng "-$r_bc_len:" --anchor-start 0 \
        -D $bc_max_mismatches \
        -t $ncores --seqtype dna --fq \
        -f --dropped "$outdir/no_index.rev.$fq" \
        -a rbc={pattern_name} \
        -a rbc.m={match_diffs} -a rbc2.m='{match_diffs(1,2)}' | #  -a rbc.end='{match_neg_end}' -a rbc.len='{match_len}'
    $st trim -e --fq "{$f_bc_len + num(attr('fl'))}:{-(num(attr('rl')) + $r_bc_len)}" |
    $st find -e "file:$all_primers" \
        -R $artifact_search_err_rate \
        --dropped "$outdir/fusion_artifact.$fq" \
        -t $ncores --seqtype dna --fq |
    $st filter --fq 'opt_attr("fbc2.m") === undefined && opt_attr("rbc2.m") === undefined' \
        --dropped $outdir/ambiguous_barcode.$fq \
        -o "$demux_out"

#cutadapt -j $cores -e 0 --no-indels -g '^file:fwd_bc.fasta' --suffix ' f={name}' trim.fastq |
#  cutadapt -j $cores -e 0 --no-indels -a 'file$:rev_bc_rc.fasta' --suffix ' r={name}' - |
#  $st split --fq -po 'demux/{attr(f)}-{attr(r)}.fastq.gz'

echo "Assembling read statistics..." >&2

echo -e "file\tcount" >"$outdir"/demux_stats.tsv
$st count \
  -k filestem \
  "$outdir"/no_primer.*.$fq \
  "$outdir"/no_index.*.$fq \
  "$outdir"/fusion_artifact.$fq \
  "$demux_out" >>"$outdir"/demux_stats.tsv

header="exp_err\tf_primer_mis\tr_primer_mis\tf_bc_mis\tr_bc_mis\tf_bc2_mis\tr_bc2_mis\tcount"
echo -e "$header" >"$outdir"/read_stats.tsv
$st count \
  -k 'bin(exp_err, 1)' \
  -k 'opt_attr(fm)' -k 'opt_attr(rm)' \
  -k "opt_attr('fbc.m')" -k "opt_attr('rbc.m')" \
  -k "opt_attr('fbc2.m')" -k "opt_attr('rbc2.m')" \
  "$outdir"/no_primer.*.$fq \
  "$outdir"/no_index.*.$fq \
  "$demux_out" \
  >>"$outdir/read_stats.tsv"

for dir in f r; do
    echo -e "pos\tcount" >"$outdir/$dir"_primer_pos.tsv
    $st count "$demux_out" -k "attr($dir)" >>"$outdir/$dir"_primer_pos.tsv
done

echo -e "file\tlength\tcount" > "$outdir/length_stats.tsv"
$st count \
  "$demux_out" \
  "$outdir/fusion_artifact.$fq" \
  -k filestem -k seqlen \
  >> "$outdir/length_stats.tsv"

# clean up
rm -f "$outdir/pre_trim.fastq" "$outdir"/_*.fasta
