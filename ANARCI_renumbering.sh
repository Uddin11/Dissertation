#!/bin/bash

# Paths to directories
BOUND_DIR="/mnt/clusters/grayling/data/c22111854/ResearchProject/Data/SabDab_all_structures/all_structures/fasta/bound"
UNBOUND_DIR="/mnt/clusters/grayling/data/c22111854/ResearchProject/Data/SabDab_all_structures/all_structures/fasta/unbound"
BOUND_OUTPUT_DIR="/mnt/clusters/grayling/data/c22111854/ResearchProject/Data/SabDab_all_structures/all_structures/fasta/bound_renumbered"
UNBOUND_OUTPUT_DIR="/mnt/clusters/grayling/data/c22111854/ResearchProject/Data/SabDab_all_structures/all_structures/fasta/unbound_renumbered"

# Create output directories if they don't exist
mkdir -p "$BOUND_OUTPUT_DIR"
mkdir -p "$UNBOUND_OUTPUT_DIR"

# Process all FASTA files in the bound directory
for file in "$BOUND_DIR"/*.fasta; do
    filename=$(basename "$file")
    output_file="$BOUND_OUTPUT_DIR/${filename%.fasta}.renumbered"
    echo "Processing $file..."
    ANARCI -i "$file" -o "$output_file" --scheme imgt
done

# Process all FASTA files in the unbound directory
for file in "$UNBOUND_DIR"/*.fasta; do
    filename=$(basename "$file")
    output_file="$UNBOUND_OUTPUT_DIR/${filename%.fasta}.renumbered"
    echo "Processing $file..."
    ANARCI -i "$file" -o "$output_file" --scheme imgt
done

echo "Processing completed."
