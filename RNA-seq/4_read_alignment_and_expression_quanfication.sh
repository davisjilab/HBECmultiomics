# fastq files need to be unzipped for processing, but we will zip them again at the end to save storage space

gunzip sample1.read1.fastq.gz
gunzip sample1.read2.fastq.gz

# calculate expression using RSEM and bowtie2
# "-p 4" indicates you want to use 4 CPUs/cores ("multithreading") which will speed up the process
# "--bowtie2" indicates that you want to use bowtie2 for alignment
# "--estimate-rspd" enables RSEM to learn from data how the reads are distributed across a transcript. The learned statistics can help us assess if any positional biases are shown in the data
# "--paired-end" indicates that the data is paired-end sequencing data (take this out if that is not true)
# "--append-names" will add gene names to the output files
# "--output-genome-bam" will output the bam file as well as the expression data. Bam files contain information about how each read aligned to the reference, so this could be useful later

rsem-calculate-expression -p 4 --bowtie2  --estimate-rspd --paired-end --append-names --output-genome-bam \
sample1.read1.fastq sample1.read2.fastq /path/to/prepared/reference/from/previous/script /path/to/output/directory

# Compress the fastq files again to save space

gzip sample1.read1.fastq
gzip sample1.read2.fastq