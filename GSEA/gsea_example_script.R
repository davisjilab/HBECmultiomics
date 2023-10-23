library(GSEA)

## ATAC peaks vs. Expression correlation
## .rnk file (in this case) includes all genes ranked by expression p-value, first column is gene name, second column is log2foldchange
## .cls file includes sample number information
##  .gmt file (in this case) includes all differentially accessible genes

GSEA("KDSalvsWTSal_ranked_pval_expression.rnk",
     "KDSalvsWTSal_clsfile_GSEA_correlate_RNAseq_ATACseq.cls", gs.db = "db_file_gsea_ATAC_all_MANorm2.gmt", 
     gs.size.threshold.max = 100000, random.seed = 42, gsea.type = "preranked", reshuffling.type = "gene.labels")
