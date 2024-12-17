
################################################################################

## Read MIND output files and combine into df (with metadata)

################################################################################


### SET UP ###

## libraries
library(tidyverse)

## set directories
base_dir <- "~/Documents/PhD/projects/CamRat/CamRat/"
data_dir <- paste0(base_dir, "data/WHS_analysis/MIND_files/MTR_scaled/cerebral_cortex_withOB/output/")

## load data
load(paste0(base_dir, "objects/df_metadata.RDS")) # df_metadata (generated in combine_metadata.R)


### READ AND FORMAT MIND NETS ###

files <- list.files(data_dir)
df_mind_cortex_withOB <- map_dfr(
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
  ) # n = 57 ROIs

### SAVE ###
save(df_mind_cortex_withOB, file = paste0(base_dir, "objects/11July2024_df_mind_cortex_withOB.RDS")) # df_mind_cortex_withOB

