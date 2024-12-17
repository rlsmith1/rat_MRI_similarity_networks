#!/bin/bash
#SBATCH --mem=5g
#SBATCH --array=0-301
#SBATCH -e /data/CamRat/WHS_preprocessing/scripts/logs/05-3dROIstats/%a.err
#SBATCH -o /data/CamRat/WHS_preprocessing/scripts/logs/05-3dROIstats/%a.out
#SBATCH -t 00:01:00
#SBATCH --mail-type=ALL

# point path to updated AFNI binaries
export PATH=$PATH:${base_dir}/AFNIbin/

# set paths to directories
base_dir=/rds/project/rds-QhAmE41a88w
data_dir=${base_dir}/data/derivatives
out_dir=${base_dir}/outputs/3dROIstats/PDw

# read subject and session based on array ID
sub_ses_table=($(cat ${base_dir}/sub_ses.txt))
export row=${sub_ses_table[${SLURM_ARRAY_TASK_ID}]}

IFS=',' read sub ses <<< ${row}

scan_dir=${data_dir}/${sub}/${ses}/anat

# Identify MTR file for this subject and session
mtr=${scan_dir}/${sub}_${ses}_MTR.nii

# Identify atlas in native MTw space
mset=${scan_dir}/03.register_MTw/WHS_atlas_in_${sub}_${ses}_MTw.nii.gz

# calculate median MTR for each ROI
if [[ -f ${mtr} ]] && [[ -f ${mset} ]]; then

	echo "Calculating median ROI MTR for "${sub}" "${ses}
	3dROIstats -nzmedian -nzvolume -nomeanout -mask ${mset} ${mtr} > ${out_dir}/${sub}_${ses}.txt

else
	echo "Can't find MTR image or mset for "${sub}" "${ses}". Check registration directory"
fi
