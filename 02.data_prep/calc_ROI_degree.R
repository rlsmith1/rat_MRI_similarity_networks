
################################################################################

## Read MIND data and calculate ROI strength for each subject

################################################################################


# setup -------------------------------------------------------------------

## LIBRARIES
library(tidyverse)

## SET DIRECTORIES
base_dir <- "~/Documents/PhD/projects/CamRat/CamRat/"

## LOAD DATA
load(paste0(base_dir, "objects/03June2024_df_mind_cortex.RDS")) # df_mind_cortex (generated in read_MIND_files.R)


# calculate node strengths ------------------------------------------------


df_degree_cortex <- df_mind_cortex %>% 
  group_by(subject, timepoint) %>% 
  dplyr::select(subject, timepoint, R1, R2, weight) %>% 
  nest() %>% 
  
  mutate(
    strength = map(
      .x = data,
      .f = ~ .x %>% 
        pivot_wider(id_cols = R1, names_from = R2, values_from = weight) %>% 
        dplyr::select(-R1) %>% 
        colSums %>% 
        enframe %>% 
        dplyr::rename_all( ~ c("region_of_interest", "degree")) %>% 
        arrange(region_of_interest)
    )
  ) %>% 
  unnest(cols = c(strength)) %>% 
  dplyr::select(-data)

# SAVE
save(df_strength_cortex, file = paste0(base_dir, "objects/03June2024_df_degree_cortex.RDS")) # df_strength

