
#==============================================================================#
# Combine TBV and subject info as metadata
#==============================================================================#

# libraries ---------------------------------------------------------------

library(tidyverse)
library(readxl)

# data --------------------------------------------------------------------

base_dir <- "~/Documents/PhD/projects/CamRat/CamRat/"
objects_dir <- paste0(base_dir, "outputs/objects/")
data_dir <- paste0(base_dir, "data/")

load(paste0(objects_dir, "27Oct2023_df_tbv.RDS")) # df_tbv (generated in calc-tbv-from-masks.R on cluster)
load(paste0(objects_dir, "ratdb_scans.Rdata")) # df_scans (metadata from collaborator)

# subject metadata
df_eda_subject_info <- read_csv(paste0(data_dir, "subject_info.csv")) %>% 
  na.omit() %>% 
  clean_names %>% 
  dplyr::rename("sex" = "gender") %>% 
  mutate(subject = id %>% 
           as.character %>% 
           ifelse(nchar(.) < 2, paste0(0, .), .) %>% 
           paste0("sub-EDAA", .),
         .before = 1
  ) %>% 
  dplyr::select(-id)



# combine and export metadata -----------------------------------------------------------------------

df_metadata <- df_tbv %>% # start with TBV
  
  # add group and sex info
  left_join(df_eda_subject_info, by = join_by(subject)) %>% 
  mutate(sex = ifelse(str_detect(subject, "JWD"), "Male", sex),
         group = ifelse(str_detect(subject, "JWD"), "normative", group),
         study = ifelse(str_detect(subject, "JWD"), "MRC", "GSK"),
         timepoint = str_remove(session, "ses-PND") %>% as.numeric %>% as.factor,
         subject = str_remove(subject, "sub-"),
         .before = 1) %>% 
  
  # add ages
  left_join(df_scans %>% 
              mutate(pattern = ifelse(str_detect(ids, "EDA"), ".*_", ".*-"),
                     subject = str_remove(ids, pattern),
                     age = dages) %>% 
              dplyr::rename("scan_date" = "dates") %>% 
              dplyr::select(subject, age, scan_date) %>% 
              mutate(timepoint = cut(age, 
                                     breaks = c(0, 25, 35, 70, 310), 
                                     labels = c("20", "35", "63", "300")) 
              ),
            by = join_by(subject, timepoint)
  ) %>% 
  distinct()


## SAVE
save(df_metadata, file = paste0(objects_dir, "df_metadata.RDS"))

