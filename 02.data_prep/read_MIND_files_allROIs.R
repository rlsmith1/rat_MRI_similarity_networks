
################################################################################

## Read MIND output files (for all ROIs, not just cortex) 
## combine into df with metadata

################################################################################


# setup ---------------------------------------------------------------

## LIBRARIES
library(tidyverse)

## SET DIRECTORIES
base_dir <- "~/Documents/PhD/projects/CamRat/CamRat/"
data_dir <- paste0(base_dir, "data/WHS_analysis/MIND_files/MTR_scaled/all_ROIs/output/")

## LOAD DATA
load(paste0(base_dir, "objects/df_metadata.RDS")) # df_metadata (generated in combine_metadata.R)


# Read & format MIND nets -------------------------------------------------


files <- list.files(data_dir)
df_mind_cortex_allROIs <- map_dfr(
  .x = files,
  .f = ~ read_csv(paste0(data_dir, .x)) %>% 
    mutate(R1 = colnames(.), .before = 1) %>% 
    pivot_longer(2:ncol(.), names_to = "R2", values_to = "weight") %>% 
    mutate(subject = .x %>% 
             str_remove("_.*") %>% 
             str_remove("sub-"),
           timepoint = .x %>% 
             str_remove(".*_") %>% 
             str_remove(".csv") %>% 
             str_remove("ses-PND") %>% 
             as.numeric %>% 
             factor(levels = c(20, 35, 63, 300)),
           .before = 1
    )
) %>%
  
  mutate(R1 = str_trim(R1), R2 = str_trim(R2)) %>% 
  
  # define edges
  mutate(edge = paste0(pmin(R1, R2), sep = " - ", pmax(R1, R2))) %>% 
  
  # add metadata, including only subjects/timepoints included in df_data_format
  inner_join(df_metadata %>% 
               ungroup %>% 
               dplyr::select(subject, timepoint, age, sex, group, tbv, scan_date) %>% 
               distinct(),
             by = join_by(subject, timepoint)
  )

## SAVE
save(df_mind_cortex_allROIs, file = paste0(base_dir, "objects/11Mar2025_df_mind_cortex_allROIs.RDS")) # df_mind_cortex


