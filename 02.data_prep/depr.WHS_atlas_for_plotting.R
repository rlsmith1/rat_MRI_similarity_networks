
################################################################################

## Generate atlas for ggplot from WHS atlas nifti image
## Also calculate centroid from point cloud

################################################################################



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
atlas_dir <- "~/Documents/PhD/projects/CamRat/atlases/WHS_SD_rat_atlas_v4_pack/"

# atlas ROI data
atlas_path <- paste0(atlas_dir, "WHS_SD_rat_atlas_v4.nii")
atlas_data <- asNifti(atlas_path)

# brain mask data
mask_path <- paste0(atlas_dir, "WHS_SD_v2_brainmask_bin.nii")
mask_data <- asNifti(mask_path)

# Atlas labels
df_labels <- read.table(paste0(atlas_dir, "WHS_SD_rat_atlas_v4.label")) %>% 
  dplyr::rename_all(~c("IDX", "-R-",  "-G-", "-B-", "-A--",  "VIS", "MSH", "LABEL")) %>% 
  as_tibble()

# Atlas hierarchy
df_system_hierarchy <- read_xlsx(path = paste0(base_dir, "data/WHS_analysis/WHS_hierarchy_RLS.xlsx")) %>% 
  dplyr::select(system2, system, region_of_interest)


# write labels csv
#write.csv(df_labels, file = paste0(base_dir, "data/WHS_analysis/WHS_labels.csv"), row.names = FALSE)


# Format atlas image data -------------------------------------------------

## Exclude any voxels that are not in brain mask (set to NA)
atlas_data_masked <- atlas_data
atlas_data_masked[mask_data == 0] <- NA

## Convert 3D matrix into 4-column table, with x, y, z coordinates of each voxel and the respective value

# Get the dimensions of the 3D matrix
dim_x <- dim(atlas_data_masked)[1]
dim_y <- dim(atlas_data_masked)[2]
dim_z <- dim(atlas_data_masked)[3]

# Create a tibble with all combinations of x, y, z coordinates
# Combine the coordinates and matrix values into a data frame

df_atlas_data <- expand.grid(x = 1:dim_x, y = 1:dim_y, z = 1:dim_z) %>% 
  mutate(IDX = as.vector(atlas_data_masked)) %>% 
  filter(!is.na(IDX)) # remove voxels that were not included in the mask

## Combine image data with atlas labels by IDX
df_whs <- df_labels %>% 
  dplyr::select(IDX, LABEL) %>% 
  left_join(df_atlas_data, by = join_by(IDX)) %>% 
  distinct() %>% 
  na.omit() %>% 
  filter(!(LABEL %in% c("Clear Label", "Periventricular gray", "Central canal")) & !str_detect(LABEL, regex("spinal", ignore_case = TRUE)))


# Identify the midline (to split by hemisphere) ------------------------------

# ## Plot
# df_whs %>% 
#   #filter(str_detect(LABEL, regex("olfactory", ignore_case = TRUE))) %>% 
#   #filter(LABEL == "Intermediodorsal thalamic nucleus") %>% 
#   # filter(x > 235 & x < 250 &
#   #          y > 500 & y < 600 &
#   #          z > 250 & z < 350
#   # ) %>% 
#   #filter(LABEL == "Cornu ammonis 1") %>% # & y < 600) %>% 
#   
#   ggplot(aes(x = x, y = y, color = z)) +
#   geom_point() +
#   geom_vline(aes(xintercept = 245.5)) +
#   scale_color_viridis(option = "A") +
#   scale_x_continuous(breaks = seq(225, 275, 5)) +
#   theme(legend.position = "none")

# arbitrarily select midline
df_whs_hemi <- df_whs %>% 
  mutate(hemisphere = ifelse(x <= 245, "left", "right")) %>% 
  mutate(LABEL = paste0(hemisphere, "_", LABEL))
df_whs_hemi


# Identify 'noise' in point cloud to remove as outliers ----------------------------------------

## Convert all ROI point clouds to LAS objects
df_clouds <- df_whs_hemi %>% 
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


# Identify centroids -------------------------------------------------------


## Calculate centroid as mean x, y, z coordinate for each ROI
df_centroids <- df_clouds %>% 
  mutate(
    centroid = map(
      .x = denoised_cloud,
      .f = ~ .x %>% 
        filter(classification == 1) %>% 
        dplyr::select(-classification) %>% 
        summarise_all( ~ mean(.x))
    )
  ) %>% 
  unnest(cols = c(centroid)) %>% 
  dplyr::select(LABEL, x, y, z) %>% 
  separate(LABEL, into = c("hemisphere", "region_of_interest"), sep = "_") %>% 
  mutate(region_of_interest = str_trim(region_of_interest))

# # quick plot to check
# df_centroids %>% 
#   dplyr::rename_at(vars(c(x, y, z)), ~ paste0(.x, "_centroid")) %>% 
#   filter(str_detect(LABEL, "Amygdaloid area, unspecified")) %>% 
#   unnest(cols = c(denoised_cloud)) %>% 
#   
#   ggplot(aes(x = x, y = z)) +
#   geom_point(aes(color = y)) +
#   geom_point(aes(x = x_centroid, y = z_centroid)) +
#   scale_color_viridis(option = "A")

## Calculate system centroid as mean x, y, z coordinate for each cortical system
df_system_centroids <- df_clouds %>% 
  
  # Group by cortical system instead of ROI
  separate(LABEL, into = c("hemisphere", "region_of_interest"), sep = "_") %>% 
  left_join(df_system_hierarchy, join_by(region_of_interest)) %>% 
  filter(!is.na(system)) %>% 
  mutate(LABEL = paste0(hemisphere, "_", system)) %>% 
  dplyr::select(LABEL, denoised_cloud) %>% 
  unnest(cols = c(denoised_cloud)) %>% 
  group_by(LABEL) %>% 
  nest() %>% 
  
  # Calculate the centroid of each cortical region
  mutate(
    centroid = map(
      .x = data,
      .f = ~ .x %>% 
        filter(classification == 1) %>% 
        dplyr::select(-classification) %>% 
        summarise_all( ~ mean(.x))
    )
  ) %>% 
  unnest(cols = c(centroid)) %>% 
  dplyr::select(LABEL, x, y, z) %>% 
  separate(LABEL, into = c("hemisphere", "system"), sep = "_") %>% 
  mutate(system = str_trim(system))


## SAVE
save(df_centroids, df_system_centroids,
     file = paste0(base_dir, "objects/WHS_centroids.RData")
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
  unnest(cols = c(sagittal_hull)) %>% 
  
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

save(df_whs_atlas, file = paste0(base_dir, "objects/20June2024_WHS_atlas_for_plotting.RDS"))



# DEPR: using sf ----------------------------------------------------------


df_sagittal_atlas <- df_whs %>%
  st_as_sf(coords = c("y", "z", "x")) %>% # swap x and z coords
  group_by(LABEL) %>%
  summarize(geometry = st_union(geometry)) %>%
  st_concave_hull(ratio = 0.35)

load(paste0(base_dir, "objects/07June2024_WHS_atlas_for_plotting.RDS"))
df_sagittal_atlas %>% 
  filter(st_geometry_type(geometry) %in% c("POLYGON", "MULTIPOLYGON")) %>% 
  
  ggplot() +
  geom_sf(aes(geometry = geometry, group = -1, fill = LABEL), lwd = 0.2) +
  theme_void() +
  theme(legend.position = "none")
