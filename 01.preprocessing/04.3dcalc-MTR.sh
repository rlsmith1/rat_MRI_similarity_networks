#!/bin/bash
#SBATCH --mem=10g
#SBATCH --array=0-301
#SBATCH -e /data/CamRat/WHS_preprocessing/scripts/logs/04-3dcalc/%a.err
#SBATCH -o /data/CamRat/WHS_preprocessing/scripts/logs/04-3dcalc/%a.out
#SBATCH -t 00:02:00
#SBATCH --mail-type=ALL

## Load AFNI module
module load afni/current-py3

## Set paths to directories
base_dir=/data/CamRat/WHS_preprocessing
data_dir=${base_dir}/data/derivatives

## Read subject and session based on array ID
sub_ses_table=($(cat ${base_dir}/sub_ses.txt))
export row=${sub_ses_table[${SLURM_ARRAY_TASK_ID}]}

IFS=',' read sub ses <<< ${row}

## Identify PDw and MTw files for this subject and session
mtw=${data_dir}/${sub}/${ses}/anat/${sub}_${ses}_MTw_deobliqued_RIP_shft.nii
pdw=${data_dir}/${sub}/${ses}/anat/${sub}_${ses}_PDw_deobliqued_RIP_shft.nii

## Check that files exist
if [[ -f ${mtw} ]] && [[ -f ${pdw} ]]; then
	echo "Calculating MTR for "${sub}" "${ses}
else
	echo "Cannot find MTw and/or PDw for "${sub}" "${ses}
fi

## Read acqps_RG.txt to pull recieved gain (RG) values for each scan
rg_values=${base_dir}/acqps_RG.txt

## Filter table for current subject & session
sub_abbr=$(echo "$sub" | awk -F'-' '{print $2}')
ses_abbr=$(echo "$ses" | awk -F'-' '{gsub(/[^0-9]/,"",$2); printf "%d", $2}')
extract_RG_value() {
    local contrast="$1"
    awk -F'\t' -v subject="$sub_abbr" -v timepoint="$ses_abbr" -v contrast="$contrast" '$1 == subject && $2 == timepoint && $4 == contrast {print $6}' ${rg_values}
}

## Extract RG value for each contrast
RG_MTw=$(extract_RG_value "MTw")
RG_PDw=$(extract_RG_value "PDw")

## Calculate scaling factor based on RG values for each scan
# If for image i=1,2,3 for MT, PD, T1 and j is the subject such that RGij is the RG for image i for the jth subject, then
# IMij = IMij / (RGij/max( RG1j, RG2j, RG3j ))

# define a bash function to find the maximum value:
find_max_value() {
    local values=("$@")
    sorted_values=($(printf "%s\n" "${values[@]}" | sort -n))
    max_value=${sorted_values[-1]}
    echo "$max_value"
}
max_RG_val=$(find_max_value "$RG_MTw" "$RG_PDw")

# identify scaling factor (SF) for each scan
sf_MTw=$(echo "$RG_MTw / $max_RG_val" | bc -l)
sf_PDw=$(echo "$RG_PDw / $max_RG_val" | bc -l)

## Rescale MTw and PDw images using their respective scaling factors
mtw_scaled=${data_dir}/${sub}/${ses}/anat/${sub}_${ses}_MTw_deobliqued_RIP_shft_scaled.nii
3dcalc -a ${mtw} -expr "a/$sf_MTw" -prefix ${mtw_scaled}
pdw_scaled=${data_dir}/${sub}/${ses}/anat/${sub}_${ses}_PDw_deobliqued_RIP_shft_scaled.nii
3dcalc -a ${pdw} -expr "a/$sf_PDw" -prefix ${pdw_scaled}

## Run 3dcalc using MTR equation (PDw - MTw)/PDw on the scaled images
output=${data_dir}/${sub}/${ses}/anat/${sub}_${ses}_MTR_scaled.nii
3dcalc -a ${pdw_scaled} -b ${mtw_scaled} -expr '(a - b)/a' -prefix ${output}

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
