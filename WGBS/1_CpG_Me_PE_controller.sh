# Command line arguments set genome and array variables
# Provide a task_samples.txt file of sample ids (no file extensions) with one per a line in working directory and a raw_sequences folder with paired fastq files (.fq.gz)

genome=$1
sample=`sed "${SLURM_ARRAY_TASK_ID}q;d" task_samples.txt`

###################
# Run Information #
###################

start=`date +%s`

hostname

########
# Trim #
########

jid1=$(sbatch \
--job-name=Trim \
--output=Trim_${sample}_%j.out \
--error=Trim_${sample}_%j.err \
--partition=${partition} \
--ntasks=15 \
--mem=12000 \
--time=0-03:00:00 \
CpG_Me_PE_switch.sh \
trim \
${genome} \
| cut -d " " -f 4)

####################
# Screen and Align #
####################

# Set threads for fastqscreen in config file
# Each multicore needs 3 cores and 5 GB RAM per a core for directional libraries

jid2=$(sbatch \
--job-name=Align \
--output=Align_${sample}_%j.out \
--error=Align_${sample}_%j.err \
--partition=${partition} \
--dependency=afterok:${jid1} \
--ntasks=18 \
--mem=64000 \
--time=3-00:00:00 \
CpG_Me_PE_switch.sh \
align \
${genome} \
| cut -d " " -f 4)

#########################
# Remove PCR Duplicates #
#########################

jid3=$(sbatch \
--job-name=Dedup \
--output=Dedup_${sample}_%j.out \
--error=Dedup_${sample}_%j.err \
--partition=${partition} \
--dependency=afterok:${jid2} \
--ntasks=1 \
--mem=30000 \
--time=2-00:00:00 \
CpG_Me_PE_switch.sh \
deduplicate \
${genome} \
| cut -d " " -f 4) 

#######################
# Insert Size Metrics #
#######################

jid4=$(sbatch \
--job-name=Insert \
--output=Insert_${sample}_%j.out \
--error=Insert_${sample}_%j.err \
--partition=${partition} \
--dependency=afterok:${jid3} \
--ntasks=2 \
--mem=10000 \
--time=0-04:00:00 \
CpG_Me_PE_switch.sh \
insert \
| cut -d " " -f 4) 

#######################
# Nucleotide Coverage #
#######################

jid5=$(sbatch \
--job-name=Coverage \
--output=Coverage_${sample}_%j.out \
--error=Coverage_${sample}_%j.err \
--partition=${partition} \
--dependency=afterok:${jid3} \
--ntasks=1 \
--mem=4000 \
--time=2-00:00:00 \
CpG_Me_PE_switch.sh \
coverage \
${genome} \
| cut -d " " -f 4) 

#######################
# Extract Methylation #
#######################

# Each multicore needs 3 cores, 2GB overhead on buffer --split_by_chromosome \
jid6=$(sbatch \
--job-name=Extract \
--output=Extract_${sample}_%j.out \
--error=Extract_${sample}_%j.err \
--partition=${partition} \
--dependency=afterok:${jid3} \
--ntasks=18 \
--mem-per-cpu=2000 \
--time=2-00:00:00 \
CpG_Me_PE_switch.sh \
extract \
${genome} \
| cut -d " " -f 4)

###########################
# Cytosine and QC reports #
###########################

sbatch \
--job-name=Report \
--output=Report_${sample}_%j.out \
--error=Report_${sample}_%j.err \
--partition=${partition} \
--dependency=afterok:${jid5}:${jid6} \
--ntasks=3 \
--mem-per-cpu=2000 \
--time=2-00:00:00 \
CpG_Me_PE_switch.sh \
cytosineReport \
${genome}

###################
# Run Information #
###################

end=`date +%s`
runtime=$((end-start))
echo ${runtime}

squeue -u $USER -o "%.8A %.4C %.10m %.20E"