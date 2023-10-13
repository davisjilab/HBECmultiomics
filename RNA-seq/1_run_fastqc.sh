cd "/directory/containing/fastq/files"

# "nogroup" forces fastqc to output the results on a per base basis, rather than grouping some base positions together
# "noextract" keeps the output compressed, which saves storage space on the cluster

fastqc * --nogroup --noextract