# Prepare reference for bowtie2
# "--gtf" indicates that the annotation file is a gtf
# "--bowtie2" indicates that we are preparing the reference for processing by bowtie2
# Adding "\" allows you to continue writing the script on a separate line. If you just added an enter key without the "\", then the script would not work

rsem-prepare-reference --gtf reference_annotation_hg38.gtf --bowtie2 \
reference_assembly_hg38.fasta \
/path/where/you/want/the/output