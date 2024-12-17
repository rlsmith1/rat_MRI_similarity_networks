### ***************** NOTE: FOR CAMBRIDGE HPC ONLY ******************

#####################################################################

## Convert scan format as it is outputted from the scanner to BIDS ##
## Copies data over to specified directory			   ##

#####################################################################

# INSTALL TIDYVERSE IF IT IS NOT INSTALLED
install.packages(setdiff(c("tidyverse"), rownames(installed.packages())))

# LOAD LIBRARIES
library(tidyverse)

# SET DIRECTORIES
orig_dir <- "/rds/project/rds-83zLpPieDV8/q10014/nii"
base_dir <- "/rds/project/rds-QhAmE41a88w"
data_dir <- paste0(base_dir, "/data")

# READ IN SCAN LIST AS TIBBLE
df_scan_list <- read_csv(paste0(base_dir, "/scan_list.txt"))

# change column names to id & date
df_scan_list <- df_scan_list %>% 
  bind_rows(tibble(colnames(.)[1], 
                   colnames(.)[2]) %>% 
              dplyr::rename_all(~colnames(df_scan_list)), 
            .) %>% 
  dplyr::rename_all(~c("id", "date"))

# ASSIGN SUBJECT AND SESSION NUMBER
df_scan_list <- df_scan_list %>%

	# add subject ID by removing study ID
	mutate(sub = paste0("sub-", id)) %>%

	# add session number
	group_by(id) %>% 
	mutate(ses = paste0("ses-", row_number())) %>%
	ungroup()

# LOOP OVER SUBJECTS AND SESSIONS (each line in df_scan_list) TO COPY SCANS TO OUR DIRECTORY IN BIDS FORMAT
for (i in 1:nrow(df_scan_list)) {
  
  # identify subject & session number
  id <- df_scan_list[i, ] %>% pull(id)
  date <- df_scan_list[i, ] %>% pull(date)
  sub <- df_scan_list[i, ] %>% pull(sub)
  ses <- df_scan_list[i, ] %>% pull(ses)
  print(paste0('Copying: ', sub, '; ' , ses))
  
  # set up BIDS folder structure in our data dir
  current_data_dir <- paste0(data_dir, '/', sub, '/', ses, '/anat/')
  dir.create(current_data_dir, recursive = T)  
  #dir.create(paste0(out_dir, sub, '/', ses, '/func/'), recursive = T) # no functional scans for ex-vivo
  
  # copy over MT, T1, and PD files
  
    # MT
    print("Copying MT scan")
    mt_old <- paste0(orig_dir, "/", id, "/", date, '/Series_000_/mt.nii')
    mt_new <- paste0(current_data_dir, "/", sub, "_", ses, "_MTw.nii")
    file.copy(mt_old, mt_new) %>% try()
    
    # T1
    print("Copying T1 scan")
    t1_old <- paste0(orig_dir, "/", id, "/", date, '/Series_000_/t1.nii')
    t1_new <- paste0(current_data_dir, "/", sub, "_", ses, "_T1w.nii")
    file.copy(t1_old, t1_new) %>% try()
    
    # PD
    print("Copying PD scan")
    pd_old <- paste0(orig_dir, "/", id, "/", date, '/Series_000_/pd.nii')
    pd_new <- paste0(current_data_dir, "/", sub, "_", ses, "_PDw.nii")
    file.copy(pd_old, pd_new) %>% try()
  
}

print("finished!")

# CREATE SCAN TABLE FOR DATA PREP AND REGISTRATION
df_scan_table <- df_scan_list %>% 

	# only register most recent scan (assumes earlier scans are corrupt)
	mutate(ses = factor(ses, levels = c("ses-1", "ses-2", "ses-3"))) %>%
	group_by(sub) %>%
	slice_max(n = 1, order_by = ses) %>% 
	ungroup %>%

	# add contrasts & init_scale for AFNI pipeline
	expand_grid(contrast = c("MTw", "PDw", "T1w")) %>% 
	mutate(init_scale = 0.9) %>% 
	dplyr::select(-id, -date)

write.table(df_scan_table, paste0(base_dir, "/scan_table_for_registration.txt"),
            sep = ",",  col.names = FALSE, row.names = FALSE, quote = FALSE)

print("end script")

