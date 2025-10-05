#!/bin/bash

# Go to directory with BAM files
cd /home/ddibina_93/RNA_seq

# Make output folder
mkdir -p /home/ddibina_93/RNA_seq/quants

# Loop over all BAM files
for bam in *.bam; do
    start=$(date +%s)  # start time

    echo "Processing $bam ..."
    featureCounts -s 0 \
        -a /home/ddibina_93/RNA_seq/Homo_sapiens.GRCh38.114.gtf \
        -o /home/ddibina_93/RNA_seq/quants/${bam%.bam}_featurecounts.txt \
        "$bam"

    end=$(date +%s)  # end time
    runtime=$(( (end - start) / 60 ))  # in minutes

    echo "✅ Completed $bam in $runtime minutes."
    echo "------------------------------------"
done
