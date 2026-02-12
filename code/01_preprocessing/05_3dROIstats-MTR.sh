#!/bin/bash
#SBATCH --mem=5g
#SBATCH --array=0-301
#SBATCH -e /data/CamRat/WHS_preprocessing/scripts/logs/05-3dROIstats/%a.err
#SBATCH -o /data/CamRat/WHS_preprocessing/scripts/logs/05-3dROIstats/%a.out
#SBATCH -t 00:01:00
#SBATCH --mail-type=ALL

# load AFNI module
module load afni/current-py3

# set paths to directories
base_dir=/data/CamRat/WHS_preprocessing
data_dir=${base_dir}/data/derivatives
out_dir=${base_dir}/outputs/3dROIstats_MTR_scaled

# create output directory if it doesn't exist
if [ ! -d "$out_dir" ]; then
    mkdir -p "$out_dir"
    echo "Directory created: $out_dir"
else
    echo "Directory already exists: $out_dir"
fi

# read subject and session based on array ID
sub_ses_table=($(cat ${base_dir}/sub_ses.txt))
export row=${sub_ses_table[${SLURM_ARRAY_TASK_ID}]}

IFS=',' read sub ses <<< ${row}

scan_dir=${data_dir}/${sub}/${ses}/anat

# Identify MTR file for this subject and session
mtr=${scan_dir}/${sub}_${ses}_MTR_scaled.nii

# Identify atlas in native MTw space
mset=${scan_dir}/03.register_MTw/WHS_atlas_in_${sub}_${ses}_MTw.nii.gz

# calculate median MTR for each ROI
if [[ -f ${mtr} ]] && [[ -f ${mset} ]]; then

	echo "Calculating median ROI MTR for "${sub}" "${ses}
	3dROIstats -nzmedian -nzvolume -nomeanout -mask ${mset} ${mtr} > ${out_dir}/${sub}_${ses}.txt

else
	echo "Can't find MTR image or mset for "${sub}" "${ses}". Check registration directory"
fi
