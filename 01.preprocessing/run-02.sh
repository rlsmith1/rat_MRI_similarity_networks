#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --array=0-903
#SBATCH -e /data/smithral/cam_psych_rat/WHS_preprocessing/scripts/logs/%a.err
#SBATCH -o /data/smithral/cam_psych_rat/WHS_preprocessing/scripts/logs/%a.out
#SBATCH -t 00:20:00
#SBATCH --mail-type=ALL

base_dir=/data/CamRat/WHS_preprocessing/
script_dir=${base_dir}/scripts

scan_table=($(cat ${base_dir}/scan_table_for_registration.txt))
export row=${scan_table[${SLURM_ARRAY_TASK_ID}]}

IFS=',' read sub ses cont init_scale <<< ${row}

echo ${sub}
echo ${ses}
echo ${cont}

# Run alignment script using values from table
bash ${script_dir}/02.align-to-atlas.sh ${sub} ${ses} ${cont}
