#!/bin/bash

# Path to motif file
motif="/home/bx2ur/code/renin_atac/data/HOCOMOCOv11_full/HOCOMOCOv11_full_MOUSE_mono_meme_format.meme"

# Directory containing input files
input_dir="/project/shefflab/brickyard/results_pipeline/gomez_atac/results_pipeline/differential/gitk_univ/fasta"

# Output directory
output_dir="/project/shefflab/brickyard/results_pipeline/gomez_atac/results_pipeline/differential/gitk_univ/meme_res"

# Loop through input files in the directory
for input_file in "$input_dir"/*.fa; do
    # Skip control.fa
    if [ "$(basename "$input_file")" = "control.fa" ]; then  # Corrected condition
        echo "$input_file"
        continue
    fi
    
    # Extract filename without extension
    filename=$(basename -- "$input_file")
    output_subdir="$output_dir/${filename%.*}"
    
    # Make output subdirectory
    mkdir -p "$output_subdir"
    
    # Run sea command
    sea --p "$input_file" --n "$input_dir/control.fa" \
        --m "$motif" \
        -oc "$output_subdir/"
done
