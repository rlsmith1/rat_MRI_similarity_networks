
###
### Generate scan table for next round of registration from QC results
###


# load libraries ----------------------------------------------------------

library(tidyverse)
library(readxl)
library(xlsx)


# read data ---------------------------------------------------------------

qc_round <- 4
base_dir <- "~/Documents/PhD/projects/CamRat/CamRat/registration_QC/"
df_qc <- read_xlsx(paste0(base_dir, "qc", qc_round, "_review.xlsx"))


# df_qc %>% 
#   filter(contrast == "MTw") %>% 
#   write.xlsx(file = paste0(base_dir, "qc", qc_round, "_review.xlsx"))


# format ------------------------------------------------------------------

df_qc %>% 
  filter(registration_success != 1) %>% 
  filter(!is.na(init_scale)) %>% 
  
  # if scan didn't register intially, rerun using original init_scale
  mutate(
    init_scale2 = ifelse(registration_success == "NA", init_scale, init_scale2)
  ) %>% 
  
  dplyr::select(-c(init_scale, registration_success, notes)) %>% 
  dplyr::rename("init_scale" = "init_scale2") %>% 
  
  # add scan that was missing
  # bind_rows(
  #   expand_grid(
  #     sub = "sub-19",
  #     ses = "ses-T1",
  #     cont = c("MTR", "PDw", "T1w"),
  #     init_scale = 0.95
  #   )
  # ) %>% 
  # arrange(sub, ses, cont) %>% 
  # filter(!is.na(init_scale)) %>% 
  
  # write table 
  write.table(paste0(base_dir, "scan_table_for_registration", (qc_round + 1), ".txt"),
              sep = ",",  col.names = FALSE, row.names = FALSE, quote = FALSE)
