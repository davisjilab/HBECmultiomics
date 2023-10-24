#!/bin/bash

###################
# Run Information #
###################

start=`date +%s`

hostname

THREADS=${SLURM_NTASKS}
MEM=$((${SLURM_MEM_PER_CPU}/1024))

echo "Allocated threads: ${THREADS}"
echo "Allocated memory:  ${MEM}"

#########
# Tidy  #
#########

mkdir slurm_logs
mv {*.out,*.err} ./slurm_logs

###########
# MultiQC #
###########

call="multiqc \
. \
--ignore slurm_logs/ \
--ignore raw_sequences/ \
--config multiqc_config_PE.yaml"

echo ${call}
eval ${call}

###########
# Bismark #
###########

call="bismark2summary \
"$(find `.` -name '*_pe.bam' -print | tr '\n' ' ')""

echo ${call}
eval ${call}

#########
# Tidy  #
#########

# Remove non-deduplicated BAM files
if [ -f "bismark_summary_report.html" ]
then
    find . -type f -name "*_pe.bam" -exec rm -f {} \;
fi

# Copy cytosine reports to central directory
mkdir cytosine_reports
find \
. \
-name '*cov.gz.CpG_report.txt.gz' \
-type f \
-not -path "./cytosine_reports" \
-print0 | \
xargs -0 \
cp -t \
"./cytosine_reports" 

# Copy merged cytosine reports to central directory
mkdir cytosine_reports_merged
find \
. \
-name '*merged_CpG_evidence.cov.gz' \
-type f \
-not -path "./cytosine_reports_merged" \
-print0 | \
xargs -0 \
cp -t \
"./cytosine_reports_merged" 

###################
# Run Information #
###################

end=`date +%s`
runtime=$((end-start))
echo ${runtime}
