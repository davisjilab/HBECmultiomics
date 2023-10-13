# load required R libraries

library(tximportData)
library(tximport)
library(DESeq2)
library(dplyr)
library(tidyr)
library(pheatmap)

# Load in sample data from RSEM/Bowtie2

setwd("/hbec_rnaseq/alignments")
dir <- "/hbec_rnaseq/alignments"
samples <- read.table("sample_info.txt", header = TRUE)
files <- file.path(dir, paste0(samples$run, ".genes.results.gz"))
names(files) <- c("WTSal1", "WTSal2", "WTSal3", "WTSal4", "KDSal1", "KDSal2", "KDSal3", "KDSal4", "WTHDM1", "WTHDM2", "KDHDM1", "KDHDM2")
txi.rsem <- tximport(files, type = "rsem", txIn = FALSE, txOut = FALSE)
head(txi.rsem$counts)

sampleTable <- data.frame(condition = factor(c("WTSal", "WTSal", "WTSal", "WTSal", "KDSal", "KDSal", "KDSal", "KDSal", "WTHDM", "WTHDM", "KDHDM", "KDHDM")))
rownames(sampleTable) <- colnames(txi.rsem$counts)
txi.rsem$length[txi.rsem$length == 0] <- 1
dds <- DESeqDataSetFromTximport(txi.rsem, sampleTable, ~condition)

# Perform differential expression analysis

dds <- estimateSizeFactors(dds)
dds <- estimateDispersions(dds)
dds <- nbinomWaldTest(dds, maxit=500)

# Perform pairwise comparisons
 
deseq2.res_WTHDMvsWTSal <- lfcShrink(dds, contrast=c("condition","WTHDM","WTSal"))
deseq2.res_KDSalvsWTSal <- lfcShrink(dds, contrast=c("condition","KDSal","WTSal"))
deseq2.res_KDHDMvsWTHDM <- lfcShrink(dds, contrast=c("condition","KDHDM","WTHDM"))

# Order results by adjusted p-value

deseq2.res_WTHDMvsWTSalordered <- deseq2.res_WTHDMvsWTSal[order(deseq2.res_WTHDMvsWTSal$padj),]
deseq2.res_KDSalvsWTSalordered <- deseq2.res_KDSalvsWTSal[order(deseq2.res_KDSalvsWTSal$padj),]
deseq2.res_KDHDMvsWTHDMordered <- deseq2.res_KDHDMvsHDM[order(deseq2.res_KDHDMvsWTHDM$padj),]

# Keep only significant results in a new table

deseq2.res_WTHDMvsWTSalordered.Sig <- subset(deseq2.res_WTHDMvsWTSalordered, padj <= 0.05 & abs(log2FoldChange) >= log2(1.2))
deseq2.res_KDSalvsWTSalordered.Sig <- subset(deseq2.res_KDSalvsWTSalordered, padj <= 0.05 & abs(log2FoldChange) >= log2(1.2))
deseq2.res_KDHDMvsWTHDMordered.Sig <- subset(deseq2.res_KDHDMvsWTHDMordered, padj <= 0.05 & abs(log2FoldChange) >= log2(1.2))

# Write all results into a csv

write.csv(as.data.frame(deseq2.res_WTHDMvsWTSalordered), 
          file="WTHDM_vs_WTSal_genes_DESeq2.csv")
write.csv(as.data.frame(deseq2.res_KDSalvsWTSalordered), 
          file="KDSal_vs_WTSal_genes_DESeq2.csv")
write.csv(as.data.frame(deseq2.res_KDHDMvsWTHDMordered), 
          file="KDHDM_vs_WTHDM_genes_DESeq2.csv")

# Write significant results into a separate csv

write.csv(as.data.frame(deseq2.res_WTHDMvsWTSalordered.Sig), 
          file="WTHDM_vs_WTSal_genes_DESeq2_significant.csv")
write.csv(as.data.frame(deseq2.res_KDSalvsWTSalordered.Sig), 
          file="KDSal_vs_WTSal_genes_DESeq2_significant.csv")
write.csv(as.data.frame(deseq2.res_KDHDMvsWTHDMordered.Sig), 
          file="KDHDM_vs_WTHDM_genes_DESeq2_significant.csv")

# Perform PCA

vsd <- vst(dds)
pdf("PCA_allsamples_condition.pdf")
plotPCA(vsd, "condition")
dev.off()

# Prepare data for heatmap

samplenames <- c("WTSal1", "WTSal2", "WTSal3", "WTSal4", "KDSal1", "KDSal2", "KDSal3", "KDSal4", "WTHDM1", "WTHDM2", "KDHDM1", "KDHDM2")
sample_type <- ("WTSal", "WTSal", "WTSal", "WTSal", "KDSal", "KDSal", "KDSal", "KDSal", "WTHDM", "WTHDM", "KDHDM", "KDHDM")
meta <- data.frame(cbind(sample_type))
rownames(meta) <- samplenames

vsd_mat <- assay(vsd)
vsd_mat_variable_genes <- vsd_mat[apply(vsd_mat, MARGIN = 1, FUN = function(x) sd(x) != 0),]

# Make heatmap

pdf("heatmap_rnaseq_allsamples.pdf")
pheatmap(vsd_mat_variable_genes, annotation_col = meta, show_rownames = FALSE, show_colnames = FALSE, 
	color = rev(RColorBrewer::brewer.pal(11, name = "RdBu")),
	border_color = "grey", fontsize = 10, 
	scale = "row",
	main = "RNA-seq clustering")
dev.off()
