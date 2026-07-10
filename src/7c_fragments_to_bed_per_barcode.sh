#!/usr/bin/env bash
# Split a 10x-style _fragments.tsv.gz into per-barcode BED files.
#
# Usage: fragments_to_bed_per_barcode.sh <fragments.tsv.gz> <prefix> <bed_output_dir>
#   fragments.tsv.gz  Input fragments file (gzipped). Comment rows starting
#                     with '#' are skipped. First 4 columns must be
#                     chr, start, end, barcode.
#   prefix            Output filename prefix. One BED file is written per
#                     barcode at <bed_output_dir>/<prefix>_<barcode>.bed
#   bed_output_dir    Output directory (created if missing).

set -euo pipefail

if [[ $# -ne 3 ]]; then
    echo "Usage: $0 <fragments.tsv.gz> <prefix> <bed_output_dir>" >&2
    exit 1
fi

frag_file="$1"
prefix="$2"
outdir="$3"

mkdir -p "$outdir"

# Pipeline:
#   1. zcat the fragments
#   2. drop comment rows, keep first 4 cols
#   3. sort by barcode so we only ever have one output file open at a time
#   4. awk writes (chr, start, end) into prefix_<barcode>.bed, closing the
#      previous file each time the barcode changes
zcat "$frag_file" \
    | awk 'BEGIN{OFS="\t"} !/^#/ {print $1, $2, $3, $4}' \
    | sort -k4,4 \
    | awk -v prefix="$prefix" -v outdir="$outdir" 'BEGIN{OFS="\t"} {
        if ($4 != cur_bc) {
            if (cur_out != "") close(cur_out)
            cur_bc = $4
            cur_out = outdir "/" prefix "_" $4 ".bed"
        }
        print $1, $2, $3 > cur_out
    }'

echo "Wrote per-barcode BED files to $outdir" >&2

