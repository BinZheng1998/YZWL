#!/bin/bash

# Check for input arguments
if [ "$#" -lt 1 ]; then
  echo "Usage: $0 <input_1.fq.gz> [input_2.fq.gz]"
  exit 1
fi

# Input files
READ1=$1
READ2=$2

# Extract sample name
if [ -z "$READ2" ]; then
  SAMPLE_NAME=$(basename "$READ1" | sed 's/.fastq.gz//; s/.fq.gz//')
  PAIRING="single"
else
  SAMPLE_NAME=$(basename "$READ1" | sed 's/_1.fastq.gz//; s/_1.fq.gz//')
  PAIRING="paired"
fi

# Output directory for this sample
SAMPLE_OUTDIR="/home/bzheng/project/06_rnaseq/03_result/$SAMPLE_NAME/fusion_res"

# Create output directory
mkdir -p "$SAMPLE_OUTDIR"

# Reference genome and annotation
GENOME_DIR="/home/bzheng/project/06_rnaseq/01_ref/fusion_ref/STAR_index"
GTF_FILE="/home/bzheng/project/06_rnaseq/01_ref/rnaseq_ref/20250528_chicken_final_noUTR.gtf"
FA_FILE="/home/bzheng/project/06_rnaseq/01_ref/rnaseq_ref/huxu.fa"
BLACKLIST="/home/bzheng/project/06_rnaseq/02_script/blacklist.tsv"
# Step 1: Quality control and trimming with fastp
if [ "$PAIRING" == "paired" ]; then
  fastp -i "$READ1" -I "$READ2" -q 20 -n 15 -u 50 -l 50 -e 20 --thread 4 -o "$SAMPLE_OUTDIR/${SAMPLE_NAME}_trimmed_1.fq.gz" -O "$SAMPLE_OUTDIR/${SAMPLE_NAME}_trimmed_2.fq.gz"
else
  fastp -i "$READ1" -q 20 -n 15 -u 50 -l 50 -e 20 --thread 4 -o "$SAMPLE_OUTDIR/${SAMPLE_NAME}_trimmed.fq.gz"
fi

# Step 2: Alignment with STAR
if [ "$PAIRING" == "paired" ]; then
  STAR --runThreadN 6 \
       --genomeDir "$GENOME_DIR" --genomeLoad NoSharedMemory \
       --readFilesIn "$SAMPLE_OUTDIR/${SAMPLE_NAME}_trimmed_1.fq.gz" "$SAMPLE_OUTDIR/${SAMPLE_NAME}_trimmed_2.fq.gz" --readFilesCommand zcat \
       --outStd BAM_Unsorted --outSAMtype BAM Unsorted --outSAMunmapped Within --outBAMcompression 0 \
       --outFilterMultimapNmax 50 --peOverlapNbasesMin 10 --alignSplicedMateMapLminOverLmate 0.5 --alignSJstitchMismatchNmax 5 -1 5 5 \
       --chimSegmentMin 10 --chimOutType WithinBAM HardClip --chimJunctionOverhangMin 10 --chimScoreDropMax 30 \
       --chimScoreJunctionNonGTAG 0 --chimScoreSeparation 1 --chimSegmentReadGapMax 3 --chimMultimapNmax 50 |
       ~/software/arriba_v2.4.0/arriba -x /dev/stdin \
       -o "$SAMPLE_OUTDIR/${SAMPLE_NAME}_fusions.tsv" \
       -O "$SAMPLE_OUTDIR/${SAMPLE_NAME}_fusions_discarded.tsv" \
       -a "$FA_FILE" \
       -g "$GTF_FILE" \
       -b "$BLACKLIST"
else
  STAR --runThreadN 6 \
       --genomeDir "$GENOME_DIR" --genomeLoad NoSharedMemory \
       --readFilesIn "$SAMPLE_OUTDIR/${SAMPLE_NAME}_trimmed.fq.gz" --readFilesCommand zcat \
       --outStd BAM_Unsorted --outSAMtype BAM Unsorted --outSAMunmapped Within --outBAMcompression 0 \
       --outFilterMultimapNmax 50 --alignSplicedMateMapLminOverLmate 0.5 --alignSJstitchMismatchNmax 5 -1 5 5 \
       --chimSegmentMin 10 --chimOutType WithinBAM HardClip --chimJunctionOverhangMin 10 --chimScoreDropMax 30 \
       --chimScoreJunctionNonGTAG 0 --chimScoreSeparation 1 --chimSegmentReadGapMax 3 --chimMultimapNmax 50 |
       ~/software/arriba_v2.4.0/arriba -x /dev/stdin \
       -o "$SAMPLE_OUTDIR/${SAMPLE_NAME}_fusions.tsv" \
       -O "$SAMPLE_OUTDIR/${SAMPLE_NAME}_fusions_discarded.tsv" \
       -a "$FA_FILE" \
       -g "$GTF_FILE" \
       -b "$BLACKLIST"
fi


# Clean up
rm "$SAMPLE_OUTDIR/${SAMPLE_NAME}_trimmed.fq.gz" \
        "$SAMPLE_OUTDIR/${SAMPLE_NAME}_trimmed_1.fq.gz" \
        "$SAMPLE_OUTDIR/${SAMPLE_NAME}_trimmed_2.fq.gz" \
        "$SAMPLE_OUTDIR/${SAMPLE_NAME}.bam"

echo "Pipeline completed successfully for sample $SAMPLE_NAME."
