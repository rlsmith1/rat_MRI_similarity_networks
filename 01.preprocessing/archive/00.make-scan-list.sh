#!/bin/bash

###################################################################

## Generates .txt file table of scans for preprocessing pipeline ##

###################################################################

# Set original (Steve) data directory
orig_dir=/rds/project/rds-83zLpPieDV8/q10014/nii

# Set output directory
out_dir=/rds/project/rds-QhAmE41a88w

# Set study ID prefix
study='B31007\|B3526'

# Get list of subjects
subjects=$(ls ${orig_dir} | grep ${study})

# Create txt where column 1 is subject ID and column 2 is scan date
for sub in ${subjects}
do 
	for ses in $(ls ${orig_dir}/$sub)
	do 
		echo "$(echo "$sub" | sed -r 's/[,_]+/-/g'),$ses"
	done
done > ${out_dir}/scan_list.txt
