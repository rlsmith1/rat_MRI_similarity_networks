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
library(readxl)

# set directories ---------------------------------------------------------

base_dir <- "/data/CamRat/WHS_preprocessing/"
in_dir <- paste0(base_dir, "outputs/MIND_files/MTR_scaled/all_ROIs/input/")
out_dir <- paste0(base_dir, "outputs/MIND_files/MTR_scaled/cerebral_cortex/input/")
#dir.create(out_dir, recursive = TRUE)

atlas_dir <- paste0(base_dir, "atlas/")

# read data ---------------------------------------------------

# atlas hierarchy
df_hierarchy <- read_csv(paste0(atlas_dir, "WHS_hierarchy.csv"))
#df_hierarchy <- read_xlsx(paste0(atlas_dir, "WHS_hierarchy_RLS.xlsx"))

# sub-ses input file with all ROIs
df_all_rois <- read_csv(paste0(in_dir, sub, "_", ses, ".csv"))

# Subset input file for regions of interest and export to new directory for MIND --------

## USING WAXHOLM ATLAS
#keep_rois <- df_hierarchy %>% 
#    filter(level1 == "Gray matter" &
#             level2 != "Rhombencephalon" &
#             level3 != "Brainstem, unspecified" &
#	     #level4 != "Olfactory bulb" & # did not register well
#	     #level8 != "Glomerular layer of the olfactory bulb" & # did not register well
#	     level8 != "Nucleus of the stria medullaris" # did not register in every scan
#    ) %>% 
#    pull(level8)

# FOR CORTEX (& OB) ONLY:
keep_rois <- df_hierarchy %>%
	#filter(level3 == "Laminated pallium") %>% # includes olfactory bulb ROIs
	filter(level4 == "Cerebral cortex") %>% # excludes olfactory bulb ROIS
	pull(level8)

df_subset_rois <- df_all_rois %>% filter(Label %in% keep_rois)

# export
write.csv(df_subset_rois, 
	  file = paste0(out_dir, sub, "_", ses, ".csv"), 
	  row.names = FALSE)

print("Finished! Please check output :)")
