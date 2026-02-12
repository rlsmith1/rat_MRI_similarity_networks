# libraries ---------------------------------------------------------------

library(tidyverse)
library(sf)
library(RNifti)
library(janitor)


# read data ---------------------------------------------------------------

base_dir <- "/data/CamRat/WHS_preprocessing/"
atlas_dir <- paste0(base_dir, "atlas/")

# atlas ROI data
atlas_path <- paste0(atlas_dir, "WHS_SD_rat_atlas_v4.nii")
atlas_data <- asNifti(atlas_path)

# atlas labels
df_labels <- read.table(paste0(atlas_dir, "WHS_SD_rat_atlas_v4.label")) %>% 
  #dplyr::rename_all(~c("IDX", "-R-",  "-G-", "-B-", "-A--",  "VIS", "MSH", "LABEL")) %>% 
  row_to_names(1) %>%	  
  as_tibble() %>%
  mutate(IDX = as.numeric(IDX))

# write labels csv
#write.csv(df_labels, file = paste0(base_dir, "data/WHS_analysis/WHS_labels.csv"), row.names = FALSE)



# format atlas image data -------------------------------------------------

# CONVERT 3D MATRIX INTO 4 COLUMN DATAFRAME, WHERE VALUE IS ROI IDX
#doParallel::registerDoParallel()

print("Converting atlas nifti file to tibble")

l_atlas_data <- 1:dim(atlas_data)[1] %>% 
  map( ~ atlas_data[.x, , ] %>% 
         as_tibble %>% 
         mutate(x = .x, .before = 1)
  )

df_atlas_data <- 1:length(l_atlas_data) %>% 
  map( ~  l_atlas_data[[.x]] %>%     
         mutate(y = row_number(), .before = 2) %>% 
         pivot_longer(3:ncol(.), names_to = "z", values_to = "IDX") %>% 
         mutate(z = str_remove(z, "V") %>% as.numeric()) %>% 
         filter(IDX != 0)
  ) %>% 
  bind_rows()

# COMBINE IMAGE DATA WITH ATLAS LABELS BY IDX

print("Adding atlas labels to coordinate idx tibble")

df_whs <- df_labels %>% 
  dplyr::select(IDX, LABEL) %>% 
  left_join(df_atlas_data, by = join_by(IDX)) %>% 
  distinct() %>% 
  filter(!is.na(x))


# Create polygon atlas for plotting -----------------------------------------------------

# FIND BORDERS FOR ALL REGIONS ACROSS ALL PLANES AND CREATE CONVEX HULL AS 3D OUTLINE

# coronal

print("generating coronal hull")

y_slices <- df_whs %>% filter(!is.na(y)) %>% arrange(y) %>% pull(y) %>% unique
df_coronal_atlas <- y_slices %>% 
  map_df( ~ df_whs %>%
            filter(y == .x) %>% 
            st_as_sf(coords = c("x", "z")) %>% 
            group_by(LABEL) %>%
            summarize(geometry = st_union(geometry)) %>%
            st_convex_hull() %>% 
            dplyr::mutate(side = "coronal",
                          slice = .x,
                          .before = 1
            )
  )

# sagittal

print("generating sagittal hull")

x_slices <- df_whs %>% filter(!is.na(x)) %>% arrange(x) %>% pull(x) %>% unique
df_sagittal_atlas <- x_slices %>% 
  map_df( ~ df_whs %>%
            filter(x == .x) %>% 
            st_as_sf(coords = c("y", "z")) %>% 
            group_by(LABEL) %>%
            summarize(geometry = st_union(geometry)) %>%
            st_convex_hull() %>% 
            dplyr::mutate(side = "sagittal",
                          slice = .x,
                          .before = 1
            )
  )

# axial

print("generating axial hull")

z_slices <- df_whs %>% filter(!is.na(z)) %>% arrange(z) %>% pull(z) %>% unique
df_axial_atlas <- z_slices %>% 
  map_df( ~ df_whs %>%
            filter(z == .x) %>% 
            st_as_sf(coords = c("x", "y")) %>% 
            group_by(LABEL) %>%
            summarize(geometry = st_union(geometry)) %>%
            st_convex_hull() %>% 
            dplyr::mutate(side = "axial",
                          slice = .x,
                          .before = 1
            )
  )

# COMBINE
df_whs_atlas <- df_axial_atlas %>% 
  bind_rows(df_sagittal_atlas) %>% 
  bind_rows(df_coronal_atlas) %>% 
  na.omit()

# SAVE
save(df_whs_atlas, file = paste0(atlas_dir, "23Oct2023_WHS_atlas_for_plotting.RDS"))

