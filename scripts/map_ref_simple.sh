#!/usr/bin/env bash
# Do Minimap2 mapping against a reference
# (same as in ref_consensus.sh)

set -euo pipefail

reads="$1"
ref="$2"
bam_out="$3"
cores="$4"
minimap2="$5"
samtools="$6"

$minimap2 -ax map-ont -t $cores "$ref" "$reads" 2>/dev/null |
  $samtools view --no-PG -@ $cores -T "$ref" -F 2308 -Su - |   # -F 2308 = without UNMAP,SECONDARY,SUPPLEMENTARY
  $samtools sort --no-PG -@ $cores -o "$bam_out" 2>/dev/null

$samtools index "$bam_out"
rm -f "$ref".fai

$samtools view "$bam_out"
