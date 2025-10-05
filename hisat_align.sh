#!/bin/bash

# === Paths ===
FASTQ_DIR="/home/ddibina_93/RNA_seq/fastq"
GENOME_INDEX="/home/ddibina_93/RNA_seq/grch38/genome"   # path to HISAT2 index (without .1.ht2 extension)
LOGFILE="/home/ddibina_93/RNA_seq/hisat2_alignment_log.txt"

# Clear or create logfile
> $LOGFILE

# List of FASTQ files (update if you renamed after concatenation)
FILES=(
    "LNCaP_EV_Norm_rep1.fastq.gz"
    "LNCaP_EV_Norm_rep2.fastq.gz"
    "LNCaP_EV_Hyp_rep1.fastq.gz"
    "LNCaP_EV_Hyp_rep2.fastq.gz"
)

# === Loop through each file ===
for f in "${FILES[@]}"; do
    SAMPLE_NAME=$(basename "$f" .fastq.gz)
    echo "Processing $SAMPLE_NAME at $(date)" | tee -a $LOGFILE
    START_TIME=$(date +%s)

    # Run HISAT2 alignment and Samtools sorting/indexing
    hisat2 -q -x $GENOME_INDEX -U $FASTQ_DIR/$f | \
    samtools sort -o ${SAMPLE_NAME}.bam

    samtools index ${SAMPLE_NAME}.bam

    END_TIME=$(date +%s)
    ELAPSED=$((END_TIME - START_TIME))
    echo "Finished $SAMPLE_NAME in $ELAPSED seconds at $(date)" | tee -a $LOGFILE
    echo "--------------------------------------" | tee -a $LOGFILE

done

echo "All files processed successfully at $(date)" | tee -a $LOGFILE
