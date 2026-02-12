
#==============================================================================#
# Combine calculated ROI values with metadata for downstream analyses
#==============================================================================#


# setup -------------------------------------------------------------------

## LIBRARIES
library(tidyverse)

## SET DIRECTORIES
base_dir <- "~/Documents/PhD/projects/CamRat/CamRat/"
objects_dir <- paste0(base_dir, "outputs/objects/")


## LOAD DATA
load(paste0(objects_dir, "df_metadata.RDS")) # df_metadata (generated in combine_metadata.R)
load(paste0(objects_dir, "03June2024_df_vol_cortex.RDS")) # df_vol_cortex (generated in calc_ROI_volume.R)
load(paste0(objects_dir, "03June2024_df_mtr_cortex.RDS")) # df_mtr_cortex (generated in calc_ROI_MTR.R)
load(paste0(objects_dir, "03June2024_df_degree_cortex.RDS")) # df_strength_cortex (generated in calc_ROI_strength.R)


## SCAN QUALITY INFO

# Note: These are hand-generated files

# did scan had gadolinium ?
df_gad <- read_xlsx(paste0(base_dir, "data/gadolinium_use_records.xlsx")) %>%
  dplyr::select(subject, scan_date, gad_by_record)

# did scan register well?
df_qc <- read_xlsx(paste0(base_dir, "registration_QC/qc4_review.xlsx")) %>% 
  filter(registration_success == 0) %>% 
  mutate(subject = str_remove(subject, "sub-"),
         timepoint = str_remove(session, "ses-PND") %>% 
           as.numeric %>% 
           factor(levels = c(20, 35, 63, 300)),
         .before = 1
  )



# Combine & format data ---------------------------------------------------


df_data_cortex <- df_vol_cortex %>% 
  left_join(df_mtr_cortex, join_by(subject, timepoint, region_of_interest)) %>% 
  left_join(df_degree_cortex, join_by(subject, timepoint, region_of_interest)) %>% 
  
  # pivot longer
  pivot_longer(c("volume", "MTR", "degree"), names_to = "feature", values_to = "value") %>% 
  
  # add metadata
  left_join(df_metadata, by = join_by(subject, timepoint)) %>% 
  
  # reorder columns
  dplyr::select(study, subject, timepoint, sex, group, age, tbv, scan_date,
                region_of_interest, feature, value) %>% 
  
  # remove scans that had gadolinium
  anti_join(
    df_gad %>% filter(gad_by_record == 1),
    by = join_by(subject, scan_date)
  ) %>% 
  
  # remove scans that failed registration
  anti_join(df_qc, by = join_by(subject, timepoint)) %>% 
  distinct()


## SAVE for downstream analyses 
save(df_data_cortex, file = paste0(objects_dir, "03June2024_combined_data_cortex.RDS")) 
  

