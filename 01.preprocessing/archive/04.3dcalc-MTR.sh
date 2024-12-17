#!/bin/bash
#SBATCH --mem=10g
#SBATCH --array=0-301
#SBATCH -e /data/CamRat/WHS_preprocessing/scripts/logs/04-3dcalc/%a.err
#SBATCH -o /data/CamRat/WHS_preprocessing/scripts/logs/04-3dcalc/%a.out
#SBATCH -t 00:02:00
#SBATCH --mail-type=ALL

# load AFNI module
module load afni/current-py3

# set paths to directories
base_dir=/data/CamRat/WHS_preprocessing
data_dir=${base_dir}/data/derivatives

# read subject and session based on array ID
sub_ses_table=($(cat ${base_dir}/sub_ses.txt))
export row=${sub_ses_table[${SLURM_ARRAY_TASK_ID}]}

IFS=',' read sub ses <<< ${row}

# Identify PDw and MTw files for this subject and session
mtw=${data_dir}/${sub}/${ses}/anat/${sub}_${ses}_MTw_deobliqued_RIP_shft.nii
pdw=${data_dir}/${sub}/${ses}/anat/${sub}_${ses}_PDw_deobliqued_RIP_shft.nii

# Check that files exist
if [[ -f ${mtw} ]] && [[ -f ${pdw} ]]; then
	echo "Calculating MTR for "${sub}" "${ses}
else
	echo "Cannot find MTw and/or PDw for "${sub}" "${ses}
fi

# run 3dcalc using MTR equation (PDw - MTw)/PDw
output=${data_dir}/${sub}/${ses}/anat/${sub}_${ses}_MTR.nii
3dcalc -a ${pdw} -b ${mtw} -expr '(a - b)/a' -prefix ${output}

# orient to align with atlas in native space
#mset=${data_dir}/${sub}/${ses}/anat/03.register_MTw/WHS_atlas_in_${sub}_${ses}_MTw.nii.gz
#mset_orient=$(3dinfo -orient ${mset})
#mtr_orient=$(3dinfo -orient ${output})
#if [[ ${mset_orient} != ${mtr_orient} ]]; then
#	3dcopy ${output} ${data_dir}/${sub}/${ses}/anat/${sub}_${ses}_MTR_${mset_orient}.nii
#	3drefit -orient ${mset_orient} ${data_dir}/${sub}/${ses}/anat/${sub}_${ses}_MTR_${mset_orient}.nii
#else
#	echo "mset and MTR file align"
#fi

echo "finished! check output files"
