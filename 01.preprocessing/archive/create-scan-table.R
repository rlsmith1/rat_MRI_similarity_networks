
#### ***NOTE: TO BE RUN LOCALLY** ###

########################################################################

## Generate scan table for next round of registration from QC results ##

########################################################################

# load libraries ----------------------------------------------------------

library(tidyverse)
library(readxl)
library(readxl)

# read data ---------------------------------------------------------------

qc_round <- 5
base_dir <- "~/Documents/PhD/projects/CamRat/preprocessing/Olivia/registration_QC/"
df_qc <- read_xlsx(paste0(base_dir, "qc", qc_round, "_review.xlsx"))


# df_qc %>% 
#   filter(contrast == "MTw") %>% 
#   write.xlsx(file = paste0(base_dir, "qc", qc_round, "_review.xlsx"))

# format ------------------------------------------------------------------

df_qc %>% 
  filter(registration_success != 1) %>% 

  # reset inti_scale as init_scale2 for next registration round
  dplyr::select(-c(init_scale, registration_success, notes)) %>% 
  dplyr::rename("init_scale" = "init_scale2") %>% 
  

# write table 
  write.table(paste0(base_dir, "scan_table_for_registration", (qc_round + 1), ".txt"),
              sep = ",",  col.names = FALSE, row.names = FALSE, quote = FALSE)
