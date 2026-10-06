#!/bin/bash
set -uo pipefail

# MEME suite install. `sea` is not on PATH by default on this cluster.
MEME_BIN="${MEME_BIN:-/home/bx2ur/meme/bin}"
export PATH="$MEME_BIN:$PATH"

# Path to motif file
motif="/home/bx2ur/code/renin_atac/data/H14CORE_meme_format.meme"

# Directory containing input files
input_dir="/project/shefflab/brickyard/results_pipeline/gomez_atac/results_pipeline/differential/geniml_univ_cc_073126/fasta"

# Output directory
output_dir="/project/shefflab/brickyard/results_pipeline/gomez_atac/results_pipeline/differential/geniml_univ_cc_073126/meme_res"

# Control (negative) sequences: FASTA of the non-renin ENCODE consensus peaks.
# Kept in input_dir under its .narrowPeak name on purpose -- the *.fa glob below
# must not pick it up as a positive set.
control="$input_dir/universe_minus_diff_up.fa"

# Number of subsets to process in parallel. Override per run, e.g.:
#   JOBS=8 ./6_motif_enrichment.sh
JOBS="${JOBS:-4}"

# Set FORCE=1 to re-run subsets that already completed successfully.
FORCE="${FORCE:-0}"

# sequences.tsv reaches 6.8 GB on some subsets and nothing downstream reads it
# (6_motif_enrichment.Rmd loads only sea.tsv), so NOSEQS=1 is available purely to
# save disk. It is NOT the fix for the memory blowup on small positive sets --
# measured, sequences.tsv is *smallest* exactly where memory is worst, because
# its size tracks the number of significant motifs, which grows with positive-set
# size. Default stays 0 so output matches the 18 subsets already completed.
NOSEQS="${NOSEQS:-0}"
sea_opts=""
[[ "$NOSEQS" == "1" ]] && sea_opts="--noseqs"

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
    # $sea_opts is deliberately unquoted -- it is a (possibly empty) option list.
    if sea --p "$input_file" --n "$control" \
           --m "$motif" $sea_opts \
           -oc "$output_subdir/"; then
        : > "$done_marker"
        echo "[DONE] $subset"
    else
        echo "[FAIL] $subset" >&2
        return 1
    fi
}
export -f run_sea
export motif input_dir output_dir control FORCE sea_opts

# Preflight: fail loudly up front instead of once per subset.
command -v sea >/dev/null || { echo "sea not found (MEME_BIN=$MEME_BIN)" >&2; exit 1; }
[[ -s "$motif"   ]] || { echo "Motif file not found: $motif" >&2; exit 1; }
[[ -s "$control" ]] || { echo "Control fasta not found: $control" >&2; exit 1; }

# Collect inputs. With no arguments, run every subset (the normal whole-run mode).
# With arguments, run only those fasta files -- used by 6_motif_enrichment_rerun.sbatch
# to give a single subset a whole SLURM task's memory.
#
# Excluded from the glob: control.fa, and the universe -- the universe is a
# background set, not a region set to test for enrichment, so scoring it against
# the ENCODE control is not a contrast we want (and it is one of the two ~260 MB
# inputs). It stays in input_dir because other steps read it from there.
if [[ $# -gt 0 ]]; then
    inputs=("$@")
    for f in "${inputs[@]}"; do
        [[ -s "$f" ]] || { echo "Input fasta not found: $f" >&2; exit 1; }
    done
else
    mapfile -t inputs < <(find "$input_dir" -maxdepth 1 -name '*.fa' \
                               ! -name 'control.fa' ! -name 'universe_*.fa' | sort)
fi

if [[ ${#inputs[@]} -eq 0 ]]; then
    echo "No input fasta files found in $input_dir" >&2
    exit 1
fi

echo "Processing ${#inputs[@]} subsets with $JOBS parallel job(s)..."
printf '%s\n' "${inputs[@]}" | xargs -n 1 -P "$JOBS" -I {} bash -c 'run_sea "$@"' _ {}

