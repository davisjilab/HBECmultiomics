cd "/directory/containing/fastq/files/"

# "--quality 20" trims low-quality bases from the ends of reads in addition to adapter sequences
#  "--stringency 6" defines the amount of overlap with adapter sequence required to trim a sequence.
# The default value for "stringency" is 1, meaning that even a single base pair matching the adapter at the end of a read will be trimmed
# "--gzip" compresses the output file, saving storage space
# "--length 50" discards reads that get shorter than 50 base pairs due to trimming (default is 20)
# "--paired" is used when you have paired-end reads
# "--clip_R1 10" trims 10 bases off of the 5' end of read one. This is commonly necessary because the beginning of the reads are usually poor quality
# "--clip_R2 10" trims 10 bases off of the 5' end of read two. This is commonly necessary because the beginning of the reads are usually poor quality
# "--fastqc_args" will cause it to run fastqc once it is done trimming, with the arguments "--noextract" and "--nogroup" which keeps the output compressed and lets you see QC data for each individual base position rather than groups of positions.

trim_galore --quality 20 --stringency 6 --gzip --length 50 --paired --clip_R1 10 --clip_R2 10 --fastqc_args "--noextract --nogroup" read1_here.fastq.gz read2_here.fastq.gz