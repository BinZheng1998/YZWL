#!/bin/bash

# Check number of arguments
if [ "$#" -ne 2 ]; then
    echo "Usage: $0 input.vcf.gz input_sample.txt"
    exit 1
fi

input_vcf="$1"
pop_samples_file="$2"

# Define paths to tools
java="/home/bzheng/software/jdk-22.0.1/bin/java"
beagle_jar="/home/bzheng/software/beagle.27May24.118.jar"
vcftools="vcftools"
bgzip="bgzip"
tabix="tabix"
bcftools="bcftools"

# Check if input files exist
if [ ! -f "$input_vcf" ]; then
    echo "Error: Input VCF file $input_vcf does not exist"
    exit 1
fi
if [ ! -f "$pop_samples_file" ]; then
    echo "Error: Sample file $pop_samples_file does not exist"
    exit 1
fi

# Check if required tools are available
for tool in "$java" "$vcftools" "$bgzip" "$tabix" "$bcftools" parallel; do
    if ! command -v "$tool" &> /dev/null; then
        echo "Error: $tool not found. Please ensure it is installed and in your PATH"
        exit 1
    fi
done
if [ ! -f "$beagle_jar" ]; then
    echo "Error: Beagle JAR file $beagle_jar does not exist"
    exit 1
fi

# Function to process a single population
process_population() {
    population="$1"
    samples="$2"
    echo "Processing population: $population with samples: $samples"

    # Create sample file
    sample_file="./${population}_samples.txt"
    echo "$samples" | tr ',' '\n' > "$sample_file"
    if [ ! -s "$sample_file" ]; then
        echo "Error: Failed to create sample file $sample_file for $population"
        exit 1
    fi

    # Subset VCF with vcftools
    output_vcf="${population}_subset"
    $vcftools --gzvcf "$input_vcf" --keep "$sample_file" --recode --out "$output_vcf" 2> "${output_vcf}.log"
    if [ $? -ne 0 ]; then
        echo "Error: vcftools failed for $population. Check ${output_vcf}.log"
        exit 1
    fi

    # Compress VCF
    $bgzip "${output_vcf}.recode.vcf"
    if [ $? -ne 0 ] || [ ! -f "${output_vcf}.recode.vcf.gz" ]; then
        echo "Error: bgzip failed for ${output_vcf}.recode.vcf"
        exit 1
    fi

    # Index VCF
    $tabix -p vcf "${output_vcf}.recode.vcf.gz"
    if [ $? -ne 0 ] || [ ! -f "${output_vcf}.recode.vcf.gz.tbi" ]; then
        echo "Error: tabix failed for ${output_vcf}.recode.vcf.gz"
        exit 1
    fi

    # Run Beagle imputation, redirecting logs to a separate file
    beagle_output="./${population}_beagle"
    $java -Xmx64g -jar "$beagle_jar" gt="${output_vcf}.recode.vcf.gz" impute=true nthreads=6 out="$beagle_output" >> "${beagle_output}.log" 2>&1
    if [ $? -ne 0 ] || [ ! -f "${beagle_output}.vcf.gz" ]; then
        echo "Error: Beagle failed for $population. Check ${beagle_output}.log"
        exit 1
    fi

    # Index Beagle output
    $tabix -p vcf "${beagle_output}.vcf.gz"
    if [ $? -ne 0 ] || [ ! -f "${beagle_output}.vcf.gz.tbi" ]; then
        echo "Error: tabix failed for ${beagle_output}.vcf.gz"
        exit 1
    fi

    # Output absolute path of Beagle VCF
    realpath "${beagle_output}.vcf.gz"

    # Clean up intermediate files
    rm -f "$sample_file" "${output_vcf}.recode.vcf.gz" "${output_vcf}.recode.vcf.gz.tbi"
    echo "Beagle imputation and indexing completed for $population. Output: ${beagle_output}.vcf.gz"
}

# Export function and variables for parallel
export -f process_population
export input_vcf java beagle_jar vcftools bgzip tabix

# Clear any existing vcf_files.txt
> vcf_files.txt

# Run parallel processing, capturing VCF paths and Beagle logs separately
parallel --halt now,fail=1 -j 12 --colsep ' ' process_population {1} {2} :::: "$pop_samples_file" | tee vcf_paths.txt

# Check if any VCF files were generated
if [ ! -s vcf_paths.txt ]; then
    echo "Error: No VCF files were generated. Check $pop_samples_file format and *.log files"
    exit 1
fi

# Collect absolute paths of all *_beagle.vcf.gz files into vcf_files.txt
find "$(pwd)" -type f -name "*_beagle.vcf.gz" > vcf_files.txt

# Check if vcf_files.txt is empty
if [ ! -s vcf_files.txt ]; then
    echo "Error: No Beagle VCF files found in current directory. Check Beagle processing"
    exit 1
fi

# Concatenate VCF files
merged_vcf="./merged_beagle.vcf.gz"
echo "Concatenating VCF files:"
cat vcf_files.txt

$bcftools merge -Oz -o "$merged_vcf" $(cat vcf_files.txt) 2> "concat.log"
if [ $? -ne 0 ] || [ ! -f "$merged_vcf" ]; then
    echo "Error: bcftools concat failed. Check concat.log"
    exit 1
fi

# Index merged VCF
$tabix -p vcf "$merged_vcf"
if [ $? -ne 0 ] || [ ! -f "${merged_vcf}.tbi" ]; then
    echo "Error: tabix failed for $merged_vcf"
    exit 1
fi

# Clean up temporary files (keep vcf_files.txt for reference)
rm -f vcf_paths.txt

echo "All populations processed, concatenated, and indexed into $merged_vcf"
echo "Beagle logs saved in individual *.log files (e.g., Pop1_beagle.log)"
echo "Processing completed."
