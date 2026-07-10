#!/bin/bash
set -uo pipefailZ

# Path to motif file
motif="/home/bx2ur/code/renin_atac/data/H14CORE_meme_format.meme"

# Directory containing input files
input_dir="/project/shefflab/brickyard/results_pipeline/gomez_atac/results_pipeline/differential/geniml_univ_cc/fasta"

# Output directory
output_dir="/project/shefflab/brickyard/results_pipeline/gomez_atac/results_pipeline/differential/geniml_univ_cc/meme_res"

# Number of subsets to process in parallel. Override per run, e.g.:
#   JOBS=8 ./6_motif_enrichment.sh
JOBS="${JOBS:-4}"

# Set FORCE=1 to re-run subsets that already completed successfully.
FORCE="${FORCE:-0}"

run_sea() {
    local input_file="$1"
    local filename subset output_subdir done_marker
    filename=$(basename -- "$input_file")
    subset="${filename%.*}"
    output_subdir="$output_dir/$subset"
    done_marker="$output_subdir/.done"

    # Resume logic: a .done sentinel is written only after sea exits 0,
    # so any subset that was killed mid-run will not have one and will re-run.
    if [[ "$FORCE" != "1" && -f "$done_marker" ]]; then
        echo "[SKIP] $subset (already complete)"
        return 0
    fi

    mkdir -p "$output_subdir"
    echo "[RUN ] $subset"

    # -oc overwrites any partial output from a previous interrupted run.
    if sea --p "$input_file" --n "$input_dir/ENCODE_ATAC_mm10_consensusPeaks.narrowPeak" \
           --m "$motif" \
           -oc "$output_subdir/"; then
        : > "$done_marker"
        echo "[DONE] $subset"
    else
        echo "[FAIL] $subset" >&2
        return 1
    fi
}
export -f run_sea
export motif input_dir output_dir FORCE

# Collect inputs, skipping control.fa.
mapfile -t inputs < <(find "$input_dir" -maxdepth 1 -name '*.fa' ! -name 'control.fa' | sort)

if [[ ${#inputs[@]} -eq 0 ]]; then
    echo "No input fasta files found in $input_dir" >&2
    exit 1
fi

echo "Processing ${#inputs[@]} subsets with $JOBS parallel job(s)..."
printf '%s\n' "${inputs[@]}" | xargs -n 1 -P "$JOBS" -I {} bash -c 'run_sea "$@"' _ {}

