#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH -t 00:05:00
#SBATCH --mail-type=ALL

# point path to updated AFNI
#export PATH=$PATH:/data/smithral/cam_psych_rat/AFNIbin/

module load afni/current-py3

# set paths
base_dir=/data/CamRat/WHS_preprocessing
data_dir=${base_dir}/data
derived_dir=${base_dir}/data/derivatives
template=${base_dir}/atlas/WHS_SD_rat_T2star_v1.01.nii

mkdir -p ${derived_dir}

#subject [sub], session [ses] and contrast [cont] selected in parent shell script (run-02.sh)
sub=$1
ses=$2
cont=$3

# set infile (native scan)
infile=${data_dir}/${sub}/${ses}/anat/${sub}_${ses}_${cont}.nii

# create corresponding directories in derived_dir
out_dir=${derived_dir}/${sub}/${ses}/anat
mkdir -p ${out_dir}
	
echo "reorienting scan "${infile}
	
# deoblique if scan is oblique
is_oblique=$( 3dinfo -is_oblique ${infile} )
if [[ ${is_oblique} ]]; then
	3dWarp -oblique2card -prefix ${out_dir}/${sub}_${ses}_${cont}_deobliqued.nii ${infile}
else
        cp ${infile} ${out_dir}/${sub}_${ses}_${cont}_deobliqued.nii
fi

# rotate scan to align orientation with WHS template
3dresample -orient RAI -prefix ${out_dir}/${sub}_${ses}_${cont}_deobliqued_RAI.nii -input ${out_dir}/${sub}_${ses}_${cont}_deobliqued.nii
3dcopy ${out_dir}/${sub}_${ses}_${cont}_deobliqued_RAI.nii ${out_dir}/${sub}_${ses}_${cont}_deobliqued_RIP.nii
3drefit -orient RIP ${out_dir}/${sub}_${ses}_${cont}_deobliqued_RIP.nii

# Center of mass align
@Align_Centers -base ${template} -dset ${out_dir}/${sub}_${ses}_${cont}_deobliqued_RIP.nii
