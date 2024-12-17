#!/usr/bin/env Rscript

# libraries ---------------------------------------------------------------

library(tidyverse)
library(oro.nifti)
library(janitor)

# set directories ---------------------------------------------------------

base_dir <- "/data/CamRat/WHS_preprocessing/"
data_dir <- paste0(base_dir, "data/derivatives/")
out_dir <- paste0(base_dir, "outputs/")


# read sub_ses table ------------------------------------------------------


df_sub_ses <- read_csv(paste0(base_dir, "sub_ses.txt"), col_names = FALSE)


# loop through subjects and sessions, read mask, and calculate TBV --------

df_tbv <- tibble()
for (i in 1:nrow(df_sub_ses)) {

	sub <- df_sub_ses[i,] %>% pull(1)
	ses <- df_sub_ses[i,] %>% pull(2)

	print(paste0("Calcuting TBV for ", sub, " ", ses))

	# read mask
	mask_path <- paste0(data_dir, sub, "/", ses, "/anat/03.register_MTw/", 
			    sub, "_", ses, "_MTw_mask.nii.gz"
			    )
	
	if (!file.exists(mask_path)) {
		
		print("Cannot find mask  :(")
		tbv <- NA

	} else {
	
		mask <- readNIfTI(mask_path)
        	mask_data <- mask@.Data

        	# count non-zero voxels in mask
        	non_zero_vox  <- mask_data[mask_data != 0] %>% length

        	# calculate TBV (scan resolution is 0.16 mm^3), note: this includes Spinal Cord!!
        	tbv <- non_zero_vox * (0.16^3)

	}

	# return tibble
	df_tmp <- tibble(
			 subject = sub,
			 session = ses,
			 tbv = tbv
			 )

	# combine with other scan data
	df_tbv <- df_tbv %>% bind_rows(df_tmp)

}

# save tibble
save(df_tbv, file = paste0(out_dir, "27Oct2023_df_tbv.RDS"))
print("finished!")
