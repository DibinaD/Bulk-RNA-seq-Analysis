# Bulk RNA-seq Analysis Pipeline 

### Author: Dibina  
### Dataset: GSE106305  
### Reference Repo: [erilu/bulk-rnaseq-analysis](https://github.com/erilu/bulk-rnaseq-analysis)  

---

# Overview

Bulk RNA sequencing (RNA-seq) is a high-throughput technique that captures and quantifies the entire transcriptome — the complete set of RNA molecules expressed in a sample. In bulk RNA-seq analysis, RNA from a group of cells or tissue is extracted, sequenced, and computationally analyzed to estimate gene expression levels. The workflow includes several key steps: quality control of raw sequencing reads, adapter and quality trimming, alignment to a reference genome, quantification of gene-level reads, and differential expression analysis using statistical models (like DESeq2). The results reveal which genes are significantly up- or down-regulated between conditions (e.g., healthy vs. disease), offering deep insights into cellular responses, signaling pathways, and regulatory mechanisms.

This repository walks you through a complete bulk RNA-seq analysis pipeline, starting from raw FASTQ files and ending with biologically meaningful interpretation.

It’s designed for beginners and intermediate bioinformaticians who want to understand every decision made along the way, not just run commands blindly.

---

# Tools Used

| Step | Tool | Why This Tool Was Used | Significance |
|------|------|-------------------------|--------------------------------|
| **1. Data Retrieval** | **SRA Toolkit (prefetch, fastq-dump)** | To fetch and convert SRA files from NCBI Sequence Read Archive to FASTQ format | Standard and reliable for downloading public datasets |
| **2. QC - Pre-alignment** | **FastQC + MultiQC** | To assess sequencing quality (per-base quality, GC content, adapter contamination) | Detects low-quality reads before alignment |
| **3. Trimming** | **Trimmomatic** | To remove adapter sequences and low-quality bases from reads | Improves alignment accuracy and reduces mapping errors |
| **4. Alignment** | **HISAT2** | Splice-aware aligner for eukaryotic genomes | Chosen over STAR for better memory efficiency and speed on small datasets |
| **5. BAM Processing** | **Samtools** | Converts SAM to BAM, sorts, and indexes alignments | Necessary for viewing in IGV and counting with featureCounts |
| **6. Alignment QC** | **Qualimap** | Evaluates alignment quality and coverage | Confirms mapping success before quantification |
| **7. Quantification** | **featureCounts (Subread)** | Counts reads mapped to genomic features using GTF annotation | Produces gene-level counts per sample |
| **8. Count Matrix Creation** | **Python** | Aggregates all per-sample counts into one matrix | Prepares data for DESeq2 input |
| **9. DEG Analysis & Visualization** | **DESeq2 R package** | Performs differential expression analysis and visualizes gene-level patterns | Identifies significantly up/downregulated genes and reveals biological insights |

---

# About dataset

Our goal is to find differentially expressed genes in response to hypoxia for the LNCaP and PC3 cell lines. Therefore, we will select the
control samples for both cell lines (Empty\_Vector for LNCaP and siCtrl for PC3) in conditions of either normoxia or hypoxia. The specific
samples we need to download are outlined in the table below:

| Sample Name                                   | GSM Identifier | SRA Identifier (SRX) | SRA Runs (SRR, download these)                     |
|-----------------------------------------------|----------------|----------------------|----------------------------------------------------|
| LNCaP\_RNA-Seq\_Empty\_Vector\_Normoxia\_rep1 | GSM3145509     | SRX4096735           | SRR7179504, SRR7179505, SRR7179506, and SRR7179507 |
| LNCaP\_RNA-Seq\_Empty\_Vector\_Normoxia\_rep2 | GSM3145510     | SRX4096736           | SRR7179508, SRR7179509, SRR7179510, and SRR7179511 |
| LNCaP\_RNA-Seq\_Empty\_Vector\_Hypoxia\_rep1  | GSM3145513     | SRX4096739           | SRR7179520, SRR7179521, SRR7179522, and SRR7179523 |
| LNCaP\_RNA-Seq\_Empty\_Vector\_Hypoxia\_rep2  | GSM3145514     | SRX4096740           | SRR7179524, SRR7179525, SRR7179526, and SRR7179527 |
| PC3\_RNA-Seq\_siCtrl\_Normoxia\_rep1          | GSM3145517     | SRX4096743           | SRR7179536                                         |
| PC3\_RNA-Seq\_siCtrl\_Normoxia\_rep2          | GSM3145518     | SRX4096744           | SRR7179537                                         |
| PC3\_RNA-Seq\_siCtrl\_Hypoxia\_rep1           | GSM3145521     | SRX4096747           | SRR7179540                                         |
| PC3\_RNA-Seq\_siCtrl\_Hypoxia\_rep2           | GSM3145522     | SRX4096748           | SRR7179541                                         |

---

# Data Download

    python fastq_download.py
  
Run this script from the folder where you want to save the FASTQ files. After execution, all raw FASTQ files will be downloaded into a fastq/ directory, ready for downstream analysis.

Obtaining raw sequencing reads is the first essential step in any RNA-seq analysis. These FASTQ files contain the actual reads produced by the sequencer, which will later be processed, aligned, and quantified. Keeping them organized in a dedicated folder ensures reproducibility and easy access for downstream steps.

After downloading, verify that all expected samples and runs are present. Missing or incomplete downloads may lead to incomplete datasets and affect later analysis.

---

# Concatenating FASTQ Files

    cat SRR7179504_pass.fastq.gz SRR7179505_pass.fastq.gz SRR7179506_pass.fastq.gz SRR7179507_pass.fastq.gz  > LNCAP_Normoxia_S1.fastq.gz
    cat SRR7179508_pass.fastq.gz SRR7179509_pass.fastq.gz SRR7179510_pass.fastq.gz SRR7179511_pass.fastq.gz  > LNCAP_Normoxia_S2.fastq.gz
    cat SRR7179520_pass.fastq.gz SRR7179521_pass.fastq.gz SRR7179522_pass.fastq.gz SRR7179523_pass.fastq.gz  > LNCAP_Hypoxia_S1.fastq.gz
    cat SRR7179524_pass.fastq.gz SRR7179525_pass.fastq.gz SRR7179526_pass.fastq.gz SRR7179527_pass.fastq.gz  > LNCAP_Hypoxia_S2.fastq.gz

    mv SRR7179536_pass.fastq.gz PC3_Normoxia_S1.fastq.gz
    mv SRR7179537_pass.fastq.gz PC3_Normoxia_S2.fastq.gz
    mv SRR7179540_pass.fastq.gz PC3_Hypoxia_S1.fastq.gz
    mv SRR7179541_pass.fastq.gz PC3_Hypoxia_S2.fastq.gz


For samples with multiple SRA runs (e.g., LNCaP), the resulting FASTQ files are concatenated into a single file per sample using cat. Single-run samples (e.g., PC3) are simply renamed with mv. After concatenation and renaming, only the final 8 FASTQ files remain, ready for alignment.

Consolidating multiple runs ensures accurate representation of each sample and prevents errors during alignment. Proper naming also reduces confusion and makes the pipeline easier to follow.

Confirm the number of reads before and after concatenation to ensure no data is lost.
Remove individual SRA-run files if no longer needed (rm SRR*) to save storage.
After this step, the folder should contain exactly 8 final FASTQ files (4 LNCaP + 4 PC3), ready for alignment and quality control.

---

# Quality Control of the raw files

## FASTQC

    fastqc fastq/*.fastq.gz -o fastqc_results/ --threads 8 
    
FASTQC performs quality control (QC) checks on raw sequencing reads. It generates per-sample reports that summarize metrics such as per-base sequence quality, GC content, sequence duplication levels, adapter content, and overrepresented sequences.

QC ensures that the sequencing data is reliable before any downstream analysis. Poor-quality reads, contamination, or sequencing biases can distort alignment and gene quantification, leading to false biological interpretations.

Key parameters to inspect:

Per-base sequence quality: Should ideally be above Q30 across most of the read. Low-quality tails may indicate the need for trimming.
Per-sequence GC content: Should reflect expected human transcriptome (~40–60% GC). Deviations may suggest contamination or bias.
Adapter content: High adapter presence indicates leftover sequencing adapters that must be trimmed.
Sequence duplication levels: Excessive duplication can indicate PCR bias or over-amplification.

Decision-making:

If the reads have consistently high quality, low adapter contamination, and expected GC content, you can proceed directly to alignment.
If there are low-quality tails, high adapter content, or other anomalies, proceed with trimming to remove low-quality bases and adapters.

## MultiQC

    multiqc fastqc_results/ -o multiqc_report/ 

Aggregates all FASTQC quality control results from fastqc_results/ into a single comprehensive report in multiqc_report/ for easy visualization and assessment of sequencing quality.

Why it’s important:

For large datasets, inspecting individual FASTQC reports is cumbersome. MultiQC allows you to quickly spot global trends, such as consistently low-quality reads in certain samples or systematic adapter contamination.

Decision-making:

Look for samples that deviate from the expected quality metrics.
Samples with extreme issues may need re-sequencing or additional trimming.
Overall trends can guide your trimming strategy (e.g., whether to use aggressive quality trimming or adapter removal).

## Trimming

    java -jar Trimmomatic-0.39/trimmomatic-0.39.jar SE -threads 4 fastq/LNCAP_Hypoxia_S1.fastq.gz fastq/LNCAP_Hypoxia_S1_trimmed.fastq TRAILING:10 -phred33 

Trimmomatic to trim low-quality bases from FASTQ files (e.g., LNCAP_Hypoxia_S1.fastq.gz) with TRAILING:10. The output is a cleaned FASTQ file (*_trimmed.fastq) ready for downstream analysis. Run FastQC on trimmed files to compare quality improvements.

Why it’s important:

Trimming improves alignment accuracy, reduces false positive counts, and ensures that downstream quantification reflects true biological signal rather than technical artifacts.

Key parameters to focus on:

TRAILING / LEADING quality thresholds: Remove low-quality bases from ends.
MINLEN: Discards reads shorter than a minimum length (commonly 20–30 bp) to avoid mapping errors.
Adapter removal: Ensures adapters do not map to the genome falsely.

Decision-making:

After trimming, rerun FASTQC to confirm improvements in per-base quality and removal of adapters.
If significant improvement is seen (higher average quality, reduced adapter content), proceed to alignment.
If low-quality bases persist or reads are too short, adjust trimming parameters (higher quality threshold, longer MINLEN) or consider additional preprocessing.

Possible outcomes in human transcriptome data:
Most reads pass QC → proceed directly to alignment.
Low-quality tails or adapters detected → trim reads.
Extremely poor quality → consider re-sequencing.

---

# Alignment & Quantification

## Reference Genome Setup

    wget https://genome-idx.s3.amazonaws.com/hisat/grch38_genome.tar.gz
    tar -xvzf grch38_genome.tar.gz
    wget ftp://ftp.ensembl.org/pub/release-114/gtf/homo_sapiens/Homo_sapiens.GRCh38.114.gtf.gz
    gunzip Homo_sapiens.GRCh38.114.gtf.gz

We start by downloading the GRCh38 human genome and its corresponding annotation file (GTF format). The genome file is extracted using tar, and the GTF file is unzipped. This reference genome and annotation are essential for aligning RNA-seq reads to the correct genomic coordinates and for counting reads per gene accurately. Without a properly indexed genome and annotation, downstream alignment and quantification would be impossible or inaccurate.

## Read Alignment

We align our cleaned FASTQ files to the reference genome using the provided **hisat2alignment.sh** script. This step produces BAM files, which are binary files storing aligned reads. Accurate alignment is crucial for correct quantification and downstream differential expression analysis. The alignment may take several hours depending on the dataset size and computational resources.

## Quality Control of Alignments

    ./qualimap_v2.3/2.4/qualimap rnaseq -bam yourpath/sample.bam -gtf Homo_sapiens.GRCh38.114.gtf –outdir yourpath/outputfolder –java-mem-size=8G

After alignment, we run Qualimap RNA-seq to assess the quality of BAM files. This tool provides metrics such as mapping rate, coverage uniformity, and gene body coverage, helping identify any technical issues or biases in the sequencing data. High-quality alignments are necessary to ensure confidence in read counts and subsequent analyses.

## Read Quantification

Using featureCounts, we count the number of reads mapping to each gene from the BAM files. **./featurecounts.sh** this script generates a per-sample raw read count file, which is the foundational data for RNA-seq analysis. Accurate counting depends on high-quality alignments and correct annotation, as miscounts can lead to false interpretations in differential expression studies.

## Count Matrix Generation

Finally, we combine individual sample counts into a single count matrix using the provided Python script (countsmatrix_wholedata.py). The count matrix is a table where rows represent genes and columns represent samples, forming the required input for downstream analyses like DESeq2 in R. This consolidation allows for easier handling, normalization, and statistical testing across all samples.

---

# Differential Expression Analysis (DESeq2)

This step focuses on identifying genes that are differentially expressed between experimental conditions (e.g., normoxia vs. hypoxia, LNCaP vs. PC3) using DESeq2, a robust R package for bulk RNA-seq analysis.

## Creating the DESeq2 Object

A DESeq2 object is constructed using the raw count matrix and a sample annotation table. This object integrates counts, experimental design, and metadata, forming the foundation for all downstream normalization, statistical testing, and visualization.

Why it’s important: Creating a structured DESeq2 object ensures that gene counts and sample information are consistently linked, which is critical for reproducibility and accurate statistical inference.

## Annotation file

Gene annotations are retrieved from genome.gtf file using Create_AnnotationFile.R script, adding meaningful metadata such as gene symbols, biotypes, and geneids.
Annotated genes allow for biologically interpretable results, proper labeling in visualizations, and easier downstream functional analysis.

## Exploratory Data Analysis (Sample Variability)

Before performing differential expression, sample-level variability is assessed to detect outliers or batch effects. Common visualizations include:
Distance Plot: Measures pairwise distances between samples to identify clustering patterns.
Variable Genes Heatmap: Highlights genes with the highest variance across samples, showing patterns of biological or technical variability.
PCA Plot: Principal Component Analysis summarizes sample relationships in a low-dimensional space, revealing clusters corresponding to experimental conditions.

These checks ensure that biological differences dominate over technical variation and help detect potential outlier samples.

## Subsetting DESeq2 Objects

Subsetting the DESeq2 object allows targeted analysis of specific comparisons (e.g., only LNCaP samples or specific experimental groups).

Focused comparisons reduce complexity, improve statistical power, and make downstream visualization clearer.

## Extracting Differential Expression Results

DESeq2 calculates log2 fold changes, p-values, and adjusted p-values for each gene, identifying those significantly up- or down-regulated between conditions.

This step generates the primary data used for biological interpretation, ranking genes based on significance and magnitude of change.

## Visualizations for DE Results

To interpret and communicate results, multiple visualizations are generated:

PlotCounts (upgraded): Shows normalized expression of specific genes across conditions to validate differential expression.
Differential Gene Heatmap: Displays expression patterns of top DE genes across all samples.
Volcano Plot: Combines log2 fold change and significance to quickly highlight the most relevant DE genes.
LogFoldChange Comparison Plot: Compares fold changes across multiple conditions or subsets to detect shared or unique DE genes.
Gene Set Enrichment Analysis (GSEA): Identifies pathways or functional categories enriched in up- or down-regulated genes, providing biological context to the DE results.

Visualization and functional analysis make results interpretable, support hypothesis generation, and facilitate reporting for publications or clinical interpretation.

---
 
# Conclusions

This repository demonstrates a complete bulk RNA-sequencing analysis pipeline, starting from downloading and processing raw sequencing files to visualizing differentially expressed genes and pathways.

From analysis:

Cell line differences: LNCaP and PC3 cells show substantial differences in their global gene expression profiles.
Shared hypoxia response: Both cell lines upregulate a common set of genes under hypoxic conditions.
Metabolic adaptation: Hypoxia upregulates glycolysis-related genes and downregulates oxidative phosphorylation genes, indicating a metabolic switch.
Differential expression results: Comprehensive DE gene lists for each cell line are generated, which can be further explored for novel genes associated with hypoxia response.

I hope this pipeline helps you confidently explore and interpret your bulk RNA-sequencing data.

---





