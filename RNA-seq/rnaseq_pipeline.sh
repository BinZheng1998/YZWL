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
SAMPLE_OUTDIR="/home/bzheng/project/04_rnaseq/03_result/$SAMPLE_NAME/rnaseq_res"

# Create output directory
mkdir -p $SAMPLE_OUTDIR

# Reference genome and annotation
GENOME_INDEX="/home/bzheng/project/04_rnaseq/01_ref/rnaseq_ref/huxu_hisat2_index"
GTF_FILE="/home/bzheng/project/04_rnaseq/01_ref/rnaseq_ref/20250528_chicken_final_noUTR.gtf"

# Step 1: Quality control and trimming with fastp
if [ "$PAIRING" == "paired" ]; then
  fastp \ 
          -i $READ1 -I $READ2 \ 
          -q 20 -n 15 -u 50 -l 50 -e 20 \ 
          --thread 8 \ 
          -o $SAMPLE_OUTDIR/${SAMPLE_NAME}_trimmed_1.fq.gz \ 
          -O $SAMPLE_OUTDIR/${SAMPLE_NAME}_trimmed_2.fq.gz
else
  fastp \ 
          -i $READ1 \ 
          -q 20 -n 15 -u 50 -l 50 -e 20 \ 
          --thread 8 \ 
          -o $SAMPLE_OUTDIR/${SAMPLE_NAME}_trimmed.fq.gz
fi

# Step 2: Alignment with HISAT2 and converting to BAM with SAMtools
if [ "$PAIRING" == "paired" ]; then
  hisat2 \ 
          -x $GENOME_INDEX \ 
          -p 8 \ 
          -1 $SAMPLE_OUTDIR/${SAMPLE_NAME}_trimmed_1.fq.gz \ 
          -2 $SAMPLE_OUTDIR/${SAMPLE_NAME}_trimmed_2.fq.gz \ 
          --dta | samtools view -bS - | samtools sort \ 
          -o $SAMPLE_OUTDIR/${SAMPLE_NAME}_aligned_sorted.bam
else
  hisat2 \ 
          -x $GENOME_INDEX \ 
          -p 8 \ 
          -U $SAMPLE_OUTDIR/${SAMPLE_NAME}_trimmed.fq.gz \ 
          --dta | samtools view -bS - | samtools sort \ 
          -o $SAMPLE_OUTDIR/${SAMPLE_NAME}_aligned_sorted.bam
fi

# Index the BAM file
samtools index $SAMPLE_OUTDIR/${SAMPLE_NAME}_aligned_sorted.bam

# Step 3: Assembly and quantification with StringTie
stringtie \ 
        $SAMPLE_OUTDIR/${SAMPLE_NAME}_aligned_sorted.bam \ 
        -G $GTF_FILE \ 
        -o $SAMPLE_OUTDIR/${SAMPLE_NAME}_transcripts.gtf \ 
        -p 8 -B -e -A $SAMPLE_OUTDIR/${SAMPLE_NAME}_gene_abund.tab

# Step 4: Read counting with featureCounts
if [ "$PAIRING" == "paired" ]; then
  featureCounts \ 
          -a $GTF_FILE \ 
          -o $SAMPLE_OUTDIR/${SAMPLE_NAME}_counts.txt \ 
          -T 8 -p -t exon -g gene_id $SAMPLE_OUTDIR/${SAMPLE_NAME}_aligned_sorted.bam
else
  featureCounts \ 
          -a $GTF_FILE \ 
          -o $SAMPLE_OUTDIR/${SAMPLE_NAME}_counts.txt \ 
          -T 8 -t exon -g gene_id $SAMPLE_OUTDIR/${SAMPLE_NAME}_aligned_sorted.bam
fi

rm $SAMPLE_OUTDIR/${SAMPLE_NAME}_trimmed.fq.gz \ 
        $SAMPLE_OUTDIR/${SAMPLE_NAME}_trimmed_1.fq.gz \ 
        $SAMPLE_OUTDIR/${SAMPLE_NAME}_trimmed_2.fq.gz \ 
        $SAMPLE_OUTDIR/${SAMPLE_NAME}_aligned_sorted.bam*
echo "Pipeline completed successfully for sample $SAMPLE_NAME."
