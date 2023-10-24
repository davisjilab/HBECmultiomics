# This is an example of how we used bedtools to find overlaps between two sets of coordinates

# ATAC overlap with WGBS (KDSal vs WTSal)

# This line looks for coordinate overlaps between two sets of coordinates

bedtools intersect -a KDSalvsWTSal_ATAC.bed -b KDSalvsWTSal_WGBS.bed -wa -wb > KDSalvsWTSal_ATAC_WGBS_alloverlaps.txt

# This next line calculates significance of the overlap between these two sets of coordinates

bedtools fisher -a KDSalvsWTSal_ATAC.bed -b KDSalvsWTSal_WGBS.bed -g hg38.genome > bedtools_fisher_KDSalvsWTSal.txt
