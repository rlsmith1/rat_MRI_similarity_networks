
#==============================================================================#
# Generate atlas for ggplot from WHS atlas nifti images 
# The right and left hemisphers are generated separately using AFNI 3dcalc
# Since the Waxholm atlas itself doesn't separate them
#==============================================================================#


# libraries ---------------------------------------------------------------

library(tidyverse)
library(sf)
library(RNifti)
library(janitor)
library(concaveman)
library(viridis)
library(lidR)
library(readxl)


# Read data ---------------------------------------------------------------

base_dir <- "~/Documents/PhD/projects/CamRat/CamRat/"
objects_dir <- paste0(base_dir, "outputs/objects")

atlas_dir <- "~/Documents/PhD/projects/CamRat/atlases/WHS_SD_rat_atlas_v4_pack/"

# atlas ROI data
LH_path <- paste0(atlas_dir, "left_mask.nii")
LH_data <- asNifti(LH_path)

RH_path <- paste0(atlas_dir, "right_mask.nii")
RH_data <- asNifti(RH_path)

# brain mask data
mask_path <- paste0(atlas_dir, "WHS_SD_v2_brainmask_bin.nii")
mask_data <- asNifti(mask_path)

# Atlas labels
df_labels <- read.table(paste0(atlas_dir, "WHS_SD_rat_atlas_v4.label")) %>% 
  dplyr::rename_all(~c("IDX", "-R-",  "-G-", "-B-", "-A--",  "VIS", "MSH", "region_of_interest")) %>% 
  as_tibble()
  
# Atlas hierarchy
df_system_hierarchy <- read_xlsx(path = paste0(base_dir, "data/WHS_analysis/WHS_hierarchy_RLS.xlsx")) %>% 
  dplyr::select(system2, system, region_of_interest)




# Write function to convert Nifti data to tibble --------------------------


nifti2tibble <- function(atlas_data, mask_data = mask_data) {
  
  # Set any voxels that are not part of the brain mask to 0
  atlas_data_masked <- atlas_data
  atlas_data_masked[mask_data == 0] <- NA
  
  # Get the dimensions of the 3D matrix
  dim_x <- dim(atlas_data_masked)[1]
  dim_y <- dim(atlas_data_masked)[2]
  dim_z <- dim(atlas_data_masked)[3]
  
  # Create a tibble with all combinations of x, y, z coordinates
  # Combine the coordinates and matrix values into a data frame
  df_atlas_data <- expand.grid(x = 1:dim_x, y = 1:dim_y, z = 1:dim_z) %>% 
    mutate(IDX = as.vector(atlas_data_masked)) %>% 
    filter(!is.na(IDX)) # remove voxels that were not included in the mask
  
  # Return atlas object
  return(df_atlas_data)
  
}


# Format atlas image data -------------------------------------------------

## Run nifti2tibble on both left and right hemisphere data
df_lh_data <- nifti2tibble(LH_data, mask_data)
df_rh_data <- nifti2tibble(RH_data, mask_data)

# combine left and right hemisphere data
df_whs <- df_lh_data %>% 
  mutate(hemisphere = "left", .before = 1) %>% 
  bind_rows(
    df_rh_data %>% 
      mutate(hemisphere = "right", .before = 1)
  ) %>% 
  
  # add atlas label information
  left_join(df_labels %>% dplyr::select(IDX, region_of_interest)) %>% 
  
  # remove uninteresting regions (ugly to plot)
  distinct() %>% 
  na.omit() %>% 
  filter(
    !(region_of_interest %in% c("Clear Label", "Periventricular gray", "Central canal")) & 
      !str_detect(region_of_interest, regex("spinal", ignore_case = TRUE))
  )


# Identify 'noise' in point cloud to remove as outliers ----------------------------------------

## Convert all ROI point clouds to LAS objects
df_clouds <- df_whs %>% 
  mutate(LABEL = paste0(hemisphere, "_", region_of_interest)) %>% 
  dplyr::select(-c(IDX, hemisphere)) %>% 
  mutate_at(vars(c(x, y, z)), ~ as.numeric(.x)) %>% 
  group_by(LABEL) %>% 
  nest() %>%
  
  mutate(
    cloud = map(
      .x = data,
      .f = ~ LAS(.x)
    )
  )

## Classify noisy points
df_clouds <- df_clouds %>% 
  mutate(
    denoised_cloud = map(
      .x = cloud,
      .f = ~ classify_noise(.x, algorithm = ivf(res = 4, n = 25))@data %>% 
        as_tibble %>% 
        clean_names %>% 
        mutate(classification = factor(classification))
    )
  )


# Convert to 2d for plotting ----------------------------------------------

## Convert to 2d for plotting, add jitter to points
df_clouds_2d <- df_clouds %>% 
  mutate(
    
    axial_2d = map(
      .x = denoised_cloud,
      .f = ~ .x %>% 
        filter(classification == 1) %>% 
        dplyr::select("x", "y") %>% 
        dplyr::rename_all( ~ c("x", "y")) %>% 
        mutate(x = jitter(x), y = jitter(y))
    ),
    coronal_2d = map(
      .x = denoised_cloud,
      .f = ~ .x %>% 
        filter(classification == 1) %>% 
        dplyr::select("x", "z") %>% 
        dplyr::rename_all( ~ c("x", "y")) %>% 
        mutate(x = jitter(x), y = jitter(y))
    ),
    sagittal_2d = map(
      .x = denoised_cloud,
      .f = ~ .x %>% 
        filter(classification == 1) %>% 
        dplyr::select("y", "z") %>% 
        dplyr::rename_all( ~ c("x", "y")) %>% 
        mutate(x = jitter(x), y = jitter(y))
    )
    
  )


# Create polygons for atlas -----------------------------------------------

## Generate hulls from 2d point clouds
doParallel::registerDoParallel()
df_hulls <- df_clouds_2d %>% 
  mutate(
    axial_hull = map(
      .x = axial_2d,
      .f = ~ concaveman(.x, concavity = 3)
    ),
    coronal_hull = map(
      .x = coronal_2d,
      .f = ~ concaveman(.x, concavity = 3)
    ),
    sagittal_hull = map(
      .x = sagittal_2d,
      .f = ~ concaveman(.x, concavity = 3)
    )
  )

## Plot
df_hulls %>% 
  filter(!str_detect(LABEL, "Ventricular system")) %>% 
  unnest(cols = c(coronal_hull)) %>% 
  
  ggplot() +
  geom_polygon(aes(x = x, y = y, fill = LABEL), color = "black", alpha = 0.75) +
  #facet_wrap(vars(LABEL), scales = "free") +
  theme_minimal() +
  labs(title = "Smooth Hull Around 3D Point Cloud (Projected to 2D)") +
  theme(legend.position = "none")

## Pivot longer and save
df_whs_atlas <- df_hulls %>% 
  dplyr::select(-c(data, cloud, denoised_cloud)) %>% 
  pivot_longer(2:ncol(.)) %>% 
  separate(name, into = c("plane", "data")) %>% 
  filter(data == "hull") %>% 
  dplyr::select(-data) %>% 
  separate(LABEL, into = c("hemisphere", "region_of_interest"), sep = "_") %>% 
  mutate(region_of_interest = str_trim(region_of_interest))


## Save!
save(df_whs_atlas, file = paste0(objects_dir, "16Dec2024_WHS_atlas_for_plotting.RDS"))
