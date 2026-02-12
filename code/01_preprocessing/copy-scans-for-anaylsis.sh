#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --array=0-47
#SBATCH --cpus-per-task=1
#SBATCH -t 00:01:00
#SBATCH --mail-type=ALL

base_dir=/rds/project/rds-cAQcxgoLHoA/afni_preprocessing
data_dir=${base_dir}/data/derived

study='EDA' #'EDA' 'JWD'
ses='ses-PND300' # 'ses-PND020' 'ses-PND035' 'ses-PND063' 'ses-PND300'
cont='T1w' # 'T1w' 'PDw' 'MTR'

outdir=${base_dir}/outputs/${cont}_scans
mkdir ${outdir}

study_outdir=${outdir}/${study}
mkdir ${study_outdir}

ses_outdir=${study_outdir}/${ses}
mkdir ${ses_outdir}

mset_outdir=${ses_outdir}/mset
infile_outdir=${ses_outdir}/infile
mask_outdir=${ses_outdir}/mask

mkdir ${mset_outdir}
mkdir ${infile_outdir}
mkdir ${mask_outdir}

subject_list=($(cat ${base_dir}/${study}_subject_list.txt))
export sub=${subject_list[${SLURM_ARRAY_TASK_ID}]}

echo ${sub}
echo ${ses}
echo ${cont}

scan_dir=${data_dir}/${sub}/${ses}/anat/02.warp_struct_${cont}
mset=${scan_dir}/SIGMA_Anatomical_Brain_Atlas_in_${sub}_${ses}_${cont}_deobliqued_RIP_shft.nii.gz
infile=${scan_dir}/${sub}_${ses}_${cont}_deobliqued_RIP_shft.nii.gz
mask=${scan_dir}/${sub}_${ses}_${cont}_deobliqued_RIP_shft_mask.nii.gz

cp ${mset} ${mset_outdir}/SIGMA_Anatomical_Brain_Atlas_in_${sub}_${ses}_${cont}_deobliqued_RIP_shft.nii.gz
cp ${infile} ${infile_outdir}/${sub}_${ses}_${cont}_deobliqued_RIP_shft.nii.gz
cp ${mask} ${mask_outdir}/${sub}_${ses}_${cont}_deobliqued_RIP_shft_mask.nii.gz
