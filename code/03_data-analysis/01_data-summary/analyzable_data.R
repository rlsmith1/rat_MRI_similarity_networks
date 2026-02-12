
#==============================================================================#
# Summarize study data
#==============================================================================#

## Set up
base_dir <- "~/Documents/PhD/projects/CamRat/CamRat/"
source(paste0(base_dir, "code/03_data-analysis/setup.R"))


# Write table to summarize data -------------------------------------------

df_data_cortex %>% 
  filter(!is.na(age)) %>% 
  filter(!(timepoint %in% c(20, 300) & study == "GSK")) %>% 
  dplyr::select(study, subject, timepoint, group) %>% 
  distinct() %>% 
  count(study, timepoint, group) %>% 
  
  mutate(cohort = case_when(
    study == "GSK" & group == "control" ~ "Experimental stress, Control",
    study == "GSK" & group == "MS" ~ "Experimental stress, RMS",
    study == "MRC" ~ "Normative development"
  )
  ) %>% 
  dplyr::select(cohort, timepoint, n) %>% 
  arrange(str_detect(cohort, "Experimental"), timepoint) %>%
  
  write_xlsx(path = paste0(tables_dir, "Table1.analyzable_samples.xlsx"))
