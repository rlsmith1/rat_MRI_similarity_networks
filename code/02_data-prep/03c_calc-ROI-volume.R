
#==============================================================================#
# Read 3dROIstats output files and pull ROI volumes for each subject
#==============================================================================#


# setup -------------------------------------------------------------------

## LIBRARIES
library(tidyverse)

## SET DIRECTORIES
base_dir <- "~/Documents/PhD/projects/CamRat/CamRat/"
atlas_dir <- "~/Documents/PhD/projects/CamRat/atlases/WHS_SD_rat_atlas_v4_pack/"
data_dir <- paste0(base_dir, "data/WHS_analysis/3dROIstats_MTR/")
objects_dir <- paste0(base_dir, "outputs/objects/")


## LOAD ATLAS DATA

# atlas labels
df_labels <- read.table(paste0(atlas_dir, "WHS_SD_rat_atlas_v4.label")) %>% 
  dplyr::rename_all(~c("IDX", "-R-",  "-G-", "-B-", "-A--",  "VIS", "MSH", "LABEL")) %>% 
  as_tibble() %>% 
  mutate(LABEL = str_trim(LABEL))

# WHS hierachy
df_hierarchy <- read_xlsx(paste0(base_dir, "data/WHS_analysis/WHS_hierarchy.xlsx"))


# read 3dROIstats files --------------------------------------------------------------

files <- list.files(data_dir)
df_vol <- map_dfr(
  .x = files,
  .f = ~ read_table(paste0(data_dir, .x)) %>% 
    dplyr::select(contains("Volume")) %>%
    pivot_longer(1:ncol(.), names_to = "IDX", values_to = "volume") %>%
    mutate(IDX = str_remove(IDX, "Volume_") %>% as.numeric) %>%
    left_join(df_labels %>% dplyr::select(IDX, LABEL), by = join_by(IDX)) %>%
    mutate(subject = .x %>%
             str_remove("_.*") %>%
             str_remove("sub-"),
           timepoint = .x %>%
             str_remove(".*_") %>%
             str_remove(".txt") %>%
             str_remove("ses-PND") %>%
             as.numeric %>%
             factor(levels = c(20, 35, 63, 300)),
           .before = 1
    )
  
)


# format & save -----------------------------------------------------------

cortical_rois <- df_hierarchy %>% filter(level4 == "Cerebral cortex") %>% pull(level8)

# subset for cortical ROIs
df_vol_cortex <- df_vol %>% 
  filter(LABEL %in% cortical_rois) %>% 
  dplyr::rename("region_of_interest" = "LABEL") %>% 
  dplyr::select(subject, timepoint, region_of_interest, volume)

## SAVE
save(df_vol_cortex, file = paste0(objects_dir, "03June2024_df_vol_cortex.RDS"))

