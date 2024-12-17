###
### Convert scan format as it is outputted from the scanner to BIDS
###

### ***************** NOTE: FOR CAMBRIDGE HPC ONLY ************************

# INSTALL TIDYVERSE IF IT IS NOT INSTALLED
install.packages(setdiff(c("tidyverse"), rownames(installed.packages())))

# LOAD LIBRARIES
library(tidyverse)

# SET DIRECTORIES
orig_dir <- "/rds/project/rds-83zLpPieDV8/q10014/nii"
base_dir <- "/rds/project/rds-cAQcxgoLHoA/afni_preprocessing/WHS_preprocessing"
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

# LOAD INFORMATION ABOUT SCANS FROM STEVE'S RATDB
load(paste0(base_dir, "/ratdb_scans.Rdata")) # df_scans

# format to get subject ID, scan date, and age/timepoint info
df_scans <- df_scans %>%
	dplyr::select(ids, dates, dages) %>%
	dplyr::rename_all(~c("id", "date", "age")) %>% 
	filter(str_detect(id, "JWD|EDA")) %>%
	mutate(timepoint = cut(age, 
				breaks = c(0, 30, 50, 150, 400), 
				labels = c("PND020", "PND035", "PND063", "PND300")
				)
		 ) %>%
	distinct() %>%
	mutate(id = str_replace_all(id,  ",_", "-"))

# ASSIGN SUBJECT AND SESSION NUMBER
df_scan_list <- df_scan_list %>%
	mutate(date = str_remove(date, "_:") %>% as.Date(format = "%Y%m%d")) %>%
	left_join(df_scans, by = c("id", "date")) %>%
	mutate(age = ifelse(is.na(age), 62, age),
		timepoint = ifelse(is.na(timepoint), "PND063", as.character(timepoint))
	) %>%
	mutate(sub = str_remove(id, ".*-") %>% paste0("sub-", .),
		ses = paste0("ses-", timepoint)
	)

# LOOP OVER SUBJECTS AND SESSIONS (each line in df_scans) TO COPY SCANS TO OUR DIRECTORY IN BIDS FORMAT
for (i in 1:nrow(df_scan_list)) {
  
  # identify subject & session number
  id <- df_scan_list[i, ] %>% pull(id) %>% str_replace("-EDAA", ",_EDAA")
  date <- df_scan_list[i, ] %>% pull(date) %>% str_remove_all("-") %>% paste0(., "_")
  sub <- df_scan_list[i, ] %>% pull(sub)
  ses <- df_scan_list[i, ] %>% pull(ses)
  print(paste0('Copying: ', sub, '; ' , ses))
  
  # set up BIDS folder structure in our data dir
  current_data_dir <- paste0(data_dir, "/", sub, '/', ses, '/anat/')
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
  expand_grid(contrast = c("MTw", "PDw", "T1w")) %>% 
  mutate(init_scale = case_when(
	
	timepoint == "PND020" ~ 0.7,
	timepoint == "PND035" ~ 0.75,
	timepoint == "PND063" ~ 0.8,
	timepoint == "PND300" ~ 0.85

		)
	) %>% 
  dplyr::select(-id, -date, -age, -timepoint)

write.table(df_scan_table, paste0(base_dir, "/scan_table_for_registration.txt"),
            sep = ",",  col.names = FALSE, row.names = FALSE, quote = FALSE)

print("end script")

