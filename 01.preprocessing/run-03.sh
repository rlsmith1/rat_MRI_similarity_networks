#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=100g
#SBATCH --array=0-17
#SBATCH -e /data/CamRat/WHS_preprocessing/scripts/logs/run-03/%a_qc3.err
#SBATCH -o /data/CamRat/WHS_preprocessing/scripts/logs/run-03/%a_qc3.out
#SBATCH -t 12:00:00
#SBATCH --mail-type=ALL

base_dir=/data/CamRat/WHS_preprocessing
script_dir=${base_dir}/scripts

# indicate QC round (leave blank '' if round 1)
qc_round=4

scan_table=($(cat ${base_dir}/scan_table_for_registration${qc_round}.txt))
export row=${scan_table[${SLURM_ARRAY_TASK_ID}]}

IFS=',' read sub ses cont init_scale <<< ${row}

echo ${sub}
echo ${ses}
echo ${cont}
echo ${init_scale}

# start with just MTw!
if [ ${cont} == 'MTw' ] 
then
	# assign directory for @animal_warper outputs
	out_dir=${base_dir}/data/derivatives/${sub}/${ses}/anat/03.register_${cont}

	# erase what we had before, if it exists (turn on for QC reruns
	#if [[ -d ${out_dir} ]]; then
        #	rm -r ${out_dir}
	#fi

	# run animal warper on scans using specified settings	
	bash ${script_dir}/03.register-to-atlas.sh ${sub} ${ses} ${cont} ${init_scale}
else
	echo  "only running MTw scans right now!"
fi
