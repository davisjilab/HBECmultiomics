library(DMRichR)

## Run DMRichR separately for each comparison
## Edit sample_info.xlsx to include samples in the comparison of interest (e.g. KDSal vs. WTSal, WTHDM vs. WTSal, or KDHDM vs. WTHDM) 

DMRichR::DM.R(genome = "hg38", testCovariate = "SampleType"
	      coverage = 1, perGroup = 1, minCpGs = 5, maxPerms = 10, cutoff = 0.05)
