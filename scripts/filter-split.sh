#!/usr/bin/env bash
# Do quality filtering and split by sample barcode (assuming header attributes added in find-primers.sh)

set -euo pipefail

exp_err=$1
min_len=$2
out="$3"
seqtool="$4"
shift 4
demux_fq="$@"

$seqtool filter "exp_err <= $exp_err && seqlen >= $min_len" $demux_fq |
  $seqtool split --fq -po "$out/{attr(fbc)}-{attr(rbc)}.fastq.gz"
