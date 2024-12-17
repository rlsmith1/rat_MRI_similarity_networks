#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
sub <- args[1]
ses <- args[2]

print(paste0("Subsetting MIND input csv for ", sub, " ", ses))

# libraries ---------------------------------------------------------------

# INSTALL PACKAGES IF IT NOT INSTALLED
#install.packages(setdiff(c("tidyverse", "janitor"), rownames(installed.packages())))

library(tidyverse)
library(janitor)

# set directories ---------------------------------------------------------

base_dir <- "/data/CamRat/WHS_preprocessing/"
in_dir <- paste0(base_dir, "outputs/MIND_files/MTR/all_ROIs/input/")
out_dir <- paste0(base_dir, "outputs/MIND_files/MTR/GM_only_noBrainstemCerebellum/input/")
#dir.create(out_dir, recursive = TRUE)

atlas_dir <- paste0(base_dir, "atlas/")

# read data ---------------------------------------------------

# atlas hierarchy
df_hierarchy <- read_csv(paste0(atlas_dir, "WHS_hierarchy.csv"))

# sub-ses input file with all ROIs
df_all_rois <- read_csv(paste0(in_dir, sub, "_", ses, ".csv"))

# Subset input file for regions of interest and export to new directory for MIND --------

gray_matter_rois <- df_hierarchy %>% 
	filter(level1 == "Gray matter") %>% 
	pull(level8) %>%
       	unique()

df_subset_rois <- df_all_rois %>% 
	filter(Label %in% gray_matter_rois & 
	       !(str_detect(Label, regex('cerebel', ignore_case = TRUE))) &
	       !(str_detect(Label, regex('brainstem', ignore_case = TRUE)))
       )

# export
write.csv(df_subset_rois, 
	  file = paste0(out_dir, sub, "_", ses, ".csv"), 
	  row.names = FALSE)

print("Finished! Please check output :)")
