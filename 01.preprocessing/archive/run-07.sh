#!/bin/bash
#SBATCH --mem=50g
#SBATCH --array=0-301
#SBATCH -e /data/CamRat/WHS_preprocessing/scripts/logs/run-07/%a.err
#SBATCH -o /data/CamRat/WHS_preprocessing/scripts/logs/run-07/%a.out
#SBATCH -t 01:30:00
#SBATCH --mail-type=ALL

# set base directory
base_dir=/data/CamRat/WHS_preprocessing

# load python module
module load python/3.10

# identify subject & session
sub_ses_table=($(cat ${base_dir}/sub_ses.txt))
export row=${sub_ses_table[${SLURM_ARRAY_TASK_ID}]}

IFS=',' read sub ses <<< ${row}

# set as environment variables
export ARG1="$sub"
export ARG2="$ses"

# Run R script and send terminal output to logs folder
python ${base_dir}/scripts/07.calc_MIND.py "$sub" "$ses"
