#!/bin/bash

# Check for correct number of arguments
if [ "$#" -ne 3 ]; then
    echo "Usage: $0 mapping.txt input.gff output.gff"
    exit 1
fi

mapping_file=$1
input_gff=$2
output_gff=$3

# Create a temporary file for storing modified GFF content
temp_file=$(mktemp)

echo "Mapping file: $mapping_file"
echo "Input GFF file: $input_gff"
echo "Output GFF file: $output_gff"
echo "Temporary file: $temp_file"

# Function to parse mapping file and generate sed commands
generate_sed_commands() {
    local mapping_file=$1
    while IFS=$'\t' read -r protein_id gene_id rest; do
        echo "s/ID=${protein_id}.exon/ID=exon-predicted-${gene_id}-/g"
        echo "s/ID=cds.${protein_id}/ID=cds-predicted-${gene_id}/g"
    done < "$mapping_file"
}

# Generate sed commands from mapping file
generate_sed_commands "$mapping_file" > sed_commands.txt
sed_commands_file="sed_commands.txt"

echo "Generated sed commands stored in: $sed_commands_file"

# Apply sed commands to modify GFF file
sed -f "$sed_commands_file" "$input_gff" > "$temp_file"

# Output modified GFF to the final output file
cat "$temp_file" > "$output_gff"

# Clean up temporary file
rm "$temp_file"
rm "$sed_commands_file"

echo "Modified GFF file saved as $output_gff"
