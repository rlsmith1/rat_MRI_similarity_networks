#!/bin/bash
#SBATCH --mem=20g
#SBATCH --array=0-301
#SBATCH -e /data/CamRat/WHS_preprocessing/scripts/logs/run-06b/%a.err
#SBATCH -o /data/CamRat/WHS_preprocessing/scripts/logs/run-06b/%a.out
#SBATCH -t 00:05:00
#SBATCH --mail-type=ALL

# set base directory
base_dir=/data/CamRat/WHS_preprocessing

# load R module
module load R/4.3.0

# identify subject and session
sub_ses_table=($(cat ${base_dir}/sub_ses.txt))
export row=${sub_ses_table[${SLURM_ARRAY_TASK_ID}]}

IFS=',' read sub ses <<< ${row}

# set as environment variables
export ARG1="$sub"
export ARG2="$ses"

# Run R script
Rscript ${base_dir}/scripts/06b.subset-csv-for-MIND.R "$sub" "$ses"
