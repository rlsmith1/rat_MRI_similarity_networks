#!/bin/bash

# point path to updated AFNI
#export PATH=$PATH:/rds/project/rds-cAQcxgoLHoA/AFNIbin/

module load afni/current-py3

# set paths to atlas, template, and data
base_dir=/data/CamRat/WHS_preprocessing
template=${base_dir}/atlas/WHS_SD_rat_T2star_v1.01.nii
atlas=${base_dir}/atlas/WHS_SD_rat_atlas_v4.nii

data_dir=${base_dir}/data/derivatives

# session, contrast, and subject are defined in the 'parent' shell script (run-03.sh)
sub=$1
ses=$2
cont=$3
init_scale=$4

subj_dir=${data_dir}/${sub}/${ses}/anat

# assign output directory
out_dir=${subj_dir}/03.register_${cont}

# erase what we had before, if it exists (turn on for QC reruns)
if [[ -d ${out_dir} ]]; then
       	rm -r ${out_dir}
fi

# set infile based on these parameters
infile=${subj_dir}/${sub}_${ses}_${cont}_deobliqued_RIP_shft.nii

# run @animal_warper
if [[ -f ${infile} ]]; then
       	echo "Warping file "${infile}
       	@animal_warper \
               	-input ${infile} \
               	-base  ${template} \
               	-atlas  ${atlas} \
               	-outdir ${out_dir} \
		-cost lpa+ZZ \
		-init_scale ${init_scale} \
		-input_abbrev ${sub}_${ses}_${cont} \
		-base_abbrev WHS_template \
		-atlas_abbrevs WHS_atlas
else
       	echo "Can't find "${infile}
fi

