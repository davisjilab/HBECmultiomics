#!/bin/bash

# on HPC, make sure that Caper's conf ~/.caper/default.conf is correctly configured to work with your HPC
# the following command will submit Caper as a leader job to SLURM with Singularity
caper hpc submit atac.wdl -i "${INPUT_JSON}" --singularity --leader-job-name ANY_GOOD_LEADER_JOB_NAME
