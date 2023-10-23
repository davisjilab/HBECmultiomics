library(MAnorm2)
library(readxl)
library(tidyverse)
library(dplyr)

setwd("differential-analysis/MAnorm2")

hbec_profile <- read_xlsx("hbec_profile_bins.xlsx",sheet = "hbec_profile_bins")

head(hbec_profile)
colnames(hbec_profile)
################################################################################
# Construct a bioCond for each group of samples.
################################################################################s
# conds_raw_reads <- list(WTSal = bioCond(log(hbec_profile[4:5] + 0.5, base = 2), hbec_profile[12:13], name = "WTSal"),
#               WTHDM = bioCond(log(hbec_profile[6:7] + 0.5, base = 2), hbec_profile[14:15], name = "WTHDM"),
#               KDHDM = bioCond(log(hbec_profile[8:9] + 0.5, base = 2), hbec_profile[16:17], name = "KDHDM"),
#               KDSal = bioCond(log(hbec_profile[10:11] + 0.5, base = 2), hbec_profile[18:19], name = "KDSal"))

################################################################################
# MA plots before normalization.
################################################################################
# pdf("MA_WTHDM_vs_WTSal_before_norm.pdf", height=6, width=6)
# MAplot(conds_raw_reads[[1]], conds_raw_reads[[2]], ylim = c(-12, 12), main = "WTHDM vs. WTSal Before normalization")
# abline(h = 0, lwd = 2, lty = 5)
# dev.off()
# 
# 
# pdf("MA_KDSal_vs_WTSal_before_norm.pdf", height=6, width=6)
# MAplot(conds_raw_reads[[1]], conds_raw_reads[[4]], ylim = c(-12, 12), main = "KDSal vs. WTSal Before normalization")
# abline(h = 0, lwd = 2, lty = 5)
# dev.off()
# 
# pdf("MA_KDHDM_vs_WTHDM_before_norm.pdf", height=6, width=6)
# MAplot(conds_raw_reads[[3]], conds_raw_reads[[2]], ylim = c(-12, 12), main = "KDHDM vs. WTHDM Before normalization")
# abline(h = 0, lwd = 2, lty = 5)
# dev.off()
# 
################################################################################
# Perform within-group normalization.
################################################################################
norm <- normalize(hbec_profile, count = 4:5, occupancy = 12:13)
norm <- normalize(norm, count = 6:7, occupancy = 14:15)
norm <- normalize(norm, count = 8:9, occupancy = 16:17)
norm <- normalize(norm, count = 10:11, occupancy = 18:19)

head(norm)
################################################################################
# Construct a bioCond for each group of samples.
################################################################################
conds <- list(WTSal = bioCond(norm[4:5], norm[12:13], name = "WTSal"),
              WTHDM = bioCond(norm[6:7], norm[14:15], name = "WTHDM"),
              KDHDM=bioCond(norm[8:9], norm[16:17], name = "KDHDM"),
              KDSal=bioCond(norm[10:11], norm[18:19], name = "KDSal"))
conds
conds <- normBioCond(conds, baseline=NULL)

names(attributes(norm))
attributes(conds)
################################################################################
# MA plots After normalization.
################################################################################
# pdf("MA_WTHDM_vs_WTSal_after_norm.pdf", height=6, width=6)
# MAplot(conds[[1]], conds[[2]], ylim = c(-12, 12), main = "WTHDM vs. WTSal")
# abline(h = 0, lwd = 2, lty = 5)
# dev.off()
# 
# pdf("MA_KDSal_vs_WTSal_after_norm.pdf", height=6, width=6)
# MAplot(conds[[1]], conds[[4]], ylim = c(-12, 12), main = "KDSal vs. WTSal")
# abline(h = 0, lwd = 2, lty = 5)
# dev.off()
# 
# pdf("MA_KDHDM_vs_WTHDM_after_norm.pdf", height=6, width=6)
# MAplot(conds[[3]], conds[[2]], ylim = c(-12, 12), main = "KDHDM vs. WTHDM")
# abline(h = 0, lwd = 2, lty = 5)
# dev.off()
# 
################################################################################
# Fit an MVC.
# The "parametric" method sometimes requires the users to explicitly specify
# initial coefficients. Try setting init.coef = c(0.1, 10) in these cases.
conds <- fitMeanVarCurve(conds, method = "parametric", occupy.only = TRUE,init.coef = c(0.1, 10))


################################################################################
# Perform differential tests.
################################################################################
#WTHDM vs. WTSal
################################################################################
res <- diffTest(conds[[1]], conds[[2]])
head(res)

#write results to file
out <- as.data.frame(res)
hbec_profile_df <- as.data.frame(hbec_profile[1:3])
out_with_coords <- merge(hbec_profile_df, out, by.x = 0, by.y = 0, all.x = TRUE, all.y = TRUE,sort = FALSE)
out_with_coords$Row.names <- NULL
write.table(out_with_coords, file="WTHDM_vs_WTSal.tsv", sep="\t", quote=F, row.names=F)

WTHDM_enriched <- filter(out_with_coords,  padj < 0.05 & Mval > 1)
write.table(WTHDM_enriched, file="WTHDM_enriched_WTHDM_vs_WTSal.tsv", sep="\t", quote=F, row.names=F)
WTSal_enriched <-  filter(out_with_coords,  padj < 0.05 & Mval < -1)
write.table(WTSal_enriched, file="WTSal_enriched_WTHDM_vs_WTSal.tsv", sep="\t", quote=F, row.names=F)
common_peaks <- filter(out_with_coords, abs(Mval) < 1 )
write.table(common_peaks, file="common_WTHDM_vs_WTSal.tsv", sep="\t", quote=F, row.names=F)


# visualize test results
pdf("MA_WTHDM_vs_WTSal_diff.pdf", height=6, width=6)
MAplot(res, padj = 0.05, main = "WTHDM vs. WTSal")
abline(h = 0, lwd = 2, lty = 5, col = "green3")
dev.off()

################################################################################
#KDSal vs. WTSal
################################################################################
res <- diffTest(conds[[1]], conds[[4]])
head(res)

#write results to file
out <- as.data.frame(res)
hbec_profile_df <- as.data.frame(hbec_profile[1:3])
out_with_coords <- merge(hbec_profile_df, out, by.x = 0, by.y = 0, all.x = TRUE, all.y = TRUE,sort = FALSE)
out_with_coords$Row.names <- NULL
write.table(out_with_coords, file="KDSal_vs_WTSal.tsv", sep="\t", quote=F, row.names=F)

KDSal_enriched <- filter(out_with_coords,  padj < 0.05 & Mval > 1)
write.table(KDSal_enriched, file="KDSal_enriched_KDSal_vs_WTSal.tsv", sep="\t", quote=F, row.names=F)
WTSal_enriched <-  filter(out_with_coords,  padj < 0.05 & Mval < -1)
write.table(WTSal_enriched, file="WTSal_enriched_KDSal_vs_WTSal.tsv", sep="\t", quote=F, row.names=F)
common_peaks <- filter(out_with_coords, abs(Mval) < 1 )
write.table(common_peaks, file="common_KDSal_vs_WTSal.tsv", sep="\t", quote=F, row.names=F)


# visualize test results
pdf("MA_KDSal_vs_WTSal_diff.pdf", height=6, width=6)
MAplot(res, padj = 0.05, main = "KDSal vs. WTSal")
abline(h = 0, lwd = 2, lty = 5, col = "green3")
dev.off()

################################################################################
#KDHDM vs. WTHDM
################################################################################
res <- diffTest(conds[[2]], conds[[3]])
head(res)

#write results to file
out <- as.data.frame(res)
hbec_profile_df <- as.data.frame(hbec_profile[1:3])
out_with_coords <- merge(hbec_profile_df, out, by.x = 0, by.y = 0, all.x = TRUE, all.y = TRUE,sort = FALSE)
out_with_coords$Row.names <- NULL
write.table(out_with_coords, file="KDHDM_vs_hdm.tsv", sep="\t", quote=F, row.names=F)

KDHDM_enriched <- filter(out_with_coords,  padj < 0.05 & Mval > 1)
write.table(KDHDM_enriched, file="KDHDM_enriched_KDHDM_vs_WTHDM.tsv", sep="\t", quote=F, row.names=F)
hdm_enriched <-  filter(out_with_coords,  padj < 0.05 & Mval < -1)
write.table(WTHDM_enriched, file="WTHDM_enriched_KDHDM_vs_hdm.tsv", sep="\t", quote=F, row.names=F)
common_peaks <- filter(out_with_coords, abs(Mval) < 1 )
write.table(common_peaks, file="common_KDHDM_vs_WTHDM.tsv", sep="\t", quote=F, row.names=F)


# visualize test results
pdf("MA_KDHDM_vs_WTHDM_diff.pdf", height=6, width=6)
MAplot(res, padj = 0.05, main = "KDHDM vs. WTHDM")
abline(h = 0, lwd = 2, lty = 5, col = "green3")
dev.off()

