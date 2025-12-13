#!/usr/bin/env bash
# Do Minimap2 mapping against a reference, then obtain the consensus sequence
# using samtools consensus with the 'simple' base frequency counting method

set -euo pipefail

reads="$1"
ref="$2"
out="$3"
fast="$4"
cores="$5"
minimap2="$6"
samtools="$7"
shift 7

comp_level=6
if [ $fast = "true" ]; then
  comp_level=1
fi

# map to reference
$minimap2 -ax map-ont -t $cores "$ref" "$reads" 2>/dev/null |
  $samtools view --no-PG -@ $cores -T "$ref" -F 2308 -Su - |   # -F 2308 = without UNMAP,SECONDARY,SUPPLEMENTARY
  $samtools sort --no-PG -l $comp_level -@ $cores -o "$out".bam 2>/dev/null

$samtools index "$out".bam
rm -f "$ref".fai

# Wbtain the consensus sequence with ambiguities (-A)
# without any heterzygosity adjustments (-H 2)
# (+ other arguments passed to the script)
$samtools consensus \
  -A \
  -l 0 \
  -H 2 \
  -@ $cores \
  "$@" \
  -o "$out"_consensus.fasta \
  "$out".bam

# output read counts
samtools idxstats "$out".bam
