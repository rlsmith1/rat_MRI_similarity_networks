#!/bin/bash

### ***************** NOTE: FOR CAMBRIDGE HPC ONLY ************************

# set base directory
base_dir=/rds/project/rds-cAQcxgoLHoA/afni_preprocessing/WHS_preprocessing

# load R module
module load R/4.2.2

# Run R script and send terminal output to logs folder
Rscript ${base_dir}/scripts/01.data-to-bids.R > ${base_dir}/scripts/logs/01.out 
