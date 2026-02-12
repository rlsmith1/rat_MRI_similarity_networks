#!/bin/bash

# Set paths to directories where data are
base_dir=/data/CamRat/WHS_preprocessing
derived_dir=${base_dir}/data/derivatives

# Indicate round of preprocessing (i.e., number of times you have run the pipeline)
qc_round=5

# Set output directory (where you want QC images to go)
out_dir=${base_dir}/outputs/registration_QC/QC${qc_round}
mkdir ${out_dir}

# Read scan table for list of scans
scan_table=$(cat ${base_dir}/scan_table_for_registration${qc_round}.txt)

# Loop through scans
for row in ${scan_table}
do

	# Pull subject, session, and contrast info from table
	IFS=',' read sub ses cont init_scale <<< ${row}

	echo "Copying "${sub}" "${ses}" "${cont}

	# Locate QC directory
	qc_dir=${derived_dir}/${sub}/${ses}/anat/03.register_${cont}/QC

	# Identify QC images of interest (in this case, qc_03)
	images=${qc_dir}/qc_03.input+wrpd_ATL_WHS_atlas.${sub}_${ses}_${cont}.*.png

	# Rename QC images with sub, ses, and cont, and copy to output directory
	for i in ${images}
	do 
		plane=$(echo $(basename $i) | cut -d '.' -f 4)
		cp $i ${out_dir}/${sub}_${ses}_${cont}_${plane}.png
	done

done

echo "finished"
