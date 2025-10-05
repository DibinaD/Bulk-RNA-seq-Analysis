library(rtracklayer)

# Import your GTF file
gtf <- rtracklayer::import("Homo_sapiens.GRCh38.114.gtf")

# Keep only gene-level entries
genes <- subset(gtf, type == "gene")

# Build annotation table
annotation <- data.frame(
  ensembl_gene_id = gsub("\\..*", "", genes$gene_id),   # strip version numbers
  gene_name       = genes$gene_name,
  chr             = as.character(seqnames(genes)),
  start           = start(genes),
  end             = end(genes),
  strand          = strand(genes),
  gene_biotype    = genes$gene_biotype
)

# Save to CSV
write.csv(annotation, "GRCh38annotation.csv", row.names = FALSE)