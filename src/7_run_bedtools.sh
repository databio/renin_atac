# Input directory containing the files
input_directory="/project/shefflab/processed/gomez_atac/results_pipeline/differential/gitk_univ/pseudobulk_data/cell_barcodes/"

# Output directory for intersection files
output_directory="/project/shefflab/processed/gomez_atac/results_pipeline/differential/gitk_univ/pseudobulk_data/cell_barcodes/peaks/"

# Universe bed file directory
universe_bed_dir="/project/gomezlab/processed"

# Create the output directory if it doesn't exist
mkdir -p "$output_directory"

# Loop through files starting with E12 or P30
for file in "$input_directory"/*.bed; do
    # Check if file exists
    if [ -e "$file" ]; then
        # Extract file name without path
        file_name=$(basename "$file")
        
        # Extract prefix of the file name (e.g., "E12" or "P30")
        prefix=$(echo "$file_name" | cut -d'_' -f1)

        # Construct the bed file path based on the prefix
        bed_file="$universe_bed_dir/${prefix}_scATAC-seq/outs/peaks.bed"

        # Output file path for each intersection
        output_file="$output_directory/${file_name%.bed}.bed"

        # Run bedtools intersect and remove duplicates
        bedtools intersect -a "$bed_file" -b "$file" -wa | sort -u > "$output_file"

        echo "Intersection for $file_name is saved in $output_file"
    fi
done