#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
sub <- args[1]
ses <- args[2]

print(paste0("Generating MIND input csv for ", sub, " ", ses))

# libraries ---------------------------------------------------------------

# INSTALL PACKAGES IF IT NOT INSTALLED
#install.packages(setdiff(c("tidyverse", "RNifti", "oro.nifti", "janitor"), rownames(installed.packages())))

library(tidyverse)
#library(RNifti)
library(oro.nifti)
library(janitor)

# set directories ---------------------------------------------------------

base_dir <- "/data/CamRat/WHS_preprocessing/"
out_dir <- paste0(base_dir, "outputs/MIND_files/MTR_scaled/all_ROIs/input/")
#dir.create(out_dir, recursive = TRUE)

# read atlas label table -------------------------------------------------

df_labels <- read.table(paste0(base_dir, "atlas/WHS_SD_rat_atlas_v4.label")) %>% 
  as_tibble() %>% 
  row_to_names(1)

ROIs <- df_labels %>% filter(IDX != 0) %>% pull(IDX)

# Load Nifti data ---------------------------------------------------------

scan_dir <- paste0(base_dir, "/data/derivatives/", sub, "/", ses, "/anat/")

# Read MTR in native space (calculated in 05)
print("reading MTR infile")
infile_path <- paste0(scan_dir, sub, "_", ses, "_MTR_scaled.nii")
if (!file.exists(infile_path)) {stop("Cannot find MTR infile")}
infile <- readNIfTI(infile_path)
infile_data <- infile@.Data

# flip axis to align with mset
infile_data <- aperm(infile_data, c(1, 3, 2))

# Read atlas in MTw native space (registered in 03)
print("reading mset")
mset_path <- paste0(scan_dir, "03.register_MTw/WHS_atlas_in_", sub, "_", ses, "_MTw.nii.gz")
if (!file.exists(mset_path)) {stop("Cannot find mset")}
mset <- readNIfTI(mset_path)
mset_data <- mset@.Data

# Find distribution of MTR values for each ROI ------------------------------

if (all(dim(infile_data) == dim(mset_data))) { # check that atlas aligns with native MTw scan

  df_mind <- ROIs %>%
    map_dfr(
      ~tibble(IDX = .x,
              MTR = infile_data[mset_data == .x]) # pull all MTR values for ROI index
    ) %>% filter(MTR != 0)

} else {
  print("Native scan and atlas in native space do not align")
}

# Format df_mind and export as .csv ------------------------------------------

# Remove regions with 2 or less voxels (not a distribution)
null_regions <- df_mind %>% count(IDX) %>% filter(n <= 2) %>% pull(IDX) %>% unique()

# add atlas labels
df_mind <- df_mind %>%
  left_join(df_labels, by = join_by(IDX)) %>%
  filter(!(IDX %in% null_regions)) %>%
  dplyr::select(LABEL, MTR) %>%
  dplyr::rename("Label" = "LABEL")

write.csv(df_mind, file = paste0(out_dir, sub, "_", ses, ".csv"), row.names = FALSE)

print("finished!")

