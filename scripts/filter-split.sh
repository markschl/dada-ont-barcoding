#!/usr/bin/env bash
# Do quality filtering and split by sample barcode (assuming header attributes added in find-primers.sh)

set -euo pipefail

exp_err=$1
min_len=$2
out="$3"
seqtool="$4"
shift 4
demux_fq="$@"

rm -f "$out"/*--xxxx.fastq.zst "$out"/*.fastq.gz

$seqtool filter "exp_err <= $exp_err && seqlen >= $min_len" $demux_fq |
  $seqtool split --fq -po "$out/{attr(fi)}--xxxx.fastq.zst"

for f in "$out"/*-xxxx.fastq.zst; do
  $seqtool split "$f" -po "$out/{attr(fi)}--{attr(ri)}.fastq.gz"
done

rm -f "$out"/*--xxxx.fastq.zst
