#!/usr/bin/env Rscript

print(paste0("Subsetting MIND input csv for WHS T2star"))

# libraries ---------------------------------------------------------------

# INSTALL PACKAGES IF IT NOT INSTALLED
#install.packages(setdiff(c("tidyverse", "janitor"), rownames(installed.packages())))

library(tidyverse)
library(janitor)

# set directories ---------------------------------------------------------

base_dir <- "/data/CamRat/WHS_preprocessing/"
in_dir <- paste0(base_dir, "outputs/MIND_files/WHS_T2star/all_ROIs/input/")
out_dir <- paste0(base_dir, "outputs/MIND_files/WHS_T2star/CortexOB_only/input/")
dir.create(out_dir, recursive = TRUE)

atlas_dir <- paste0(base_dir, "atlas/")

# read data ---------------------------------------------------

# atlas hierarchy
df_hierarchy <- read_csv(paste0(atlas_dir, "WHS_hierarchy.csv"))

# WHS T2star input file with all ROIs
df_all_rois <- read_csv(paste0(in_dir, "WHS_SD_rat_T2star_v1.01.csv"))

# Subset input file for regions of interest and export to new directory for MIND --------

# FOR CORTEX (& OB) ONLY:
keep_rois <- df_hierarchy %>%
        filter(level3 == "Laminated pallium") %>%
        pull(level8)

df_subset_rois <- df_all_rois %>% filter(Label %in% keep_rois)

# export
write.csv(df_subset_rois, 
	  file = paste0(out_dir, "WHS_SD_rat_T2star_v1.01.csv"), 
	  row.names = FALSE)

print("Finished! Please check output :)")
