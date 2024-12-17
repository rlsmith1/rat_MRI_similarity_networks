
######################################################################################

## Calculate the volume of each ROI across subjects using the WHS standard template ##

######################################################################################

# libraries ---------------------------------------------------------------

library(tidyverse)
library(janitor)
library(readxl)
library(writexl)
library(ggrepel)

# data --------------------------------------------------------------------

base_dir <- "~/Documents/PhD/projects/CamRat/preprocessing/Olivia/"
atlas_dir <- paste0(base_dir, "atlas/WHS_SD_rat_atlas_v4_pack/")
data_dir <- paste0(base_dir, "data/3dROIstats/")

# atlas labels
df_labels <- read.table(paste0(atlas_dir, "WHS_SD_rat_atlas_v4.label")) %>% 
    dplyr::rename_all(~c("IDX", "-R-",  "-G-", "-B-", "-A--",  "VIS", "MSH", "LABEL")) %>% 
    as_tibble() %>% 
    mutate(LABEL = str_trim(LABEL))

# list of scans registered
df_scan_list <- read_csv(paste0(base_dir, "scan_list.txt"), col_names = FALSE) %>% 
    dplyr::rename_all( ~ c("subject_id", "scan_date")) %>% 
    mutate(scan_date = as.Date(scan_date, format = "%Y%m%d_", "%Y-%m-%d")) %>% 
    separate(subject_id, into = c("study", "id"), sep = "-") %>% 
    
    # add session number
    group_by(study, id) %>% 
    mutate(ses = paste0("ses-", row_number())) %>%
    ungroup()


# read files --------------------------------------------------------------

contrasts <- c("PDw", "MTw", "T1w")
df_volumes <- tibble()
for (cont in contrasts) {
    
    print("Reading volumes from ", cont, " scans")
    
    # List all files in contrast data directory
    files <- list.files(paste0(data_dir, cont))
    
    # Read each file and format to align IDX with atlas LABEL
    df_tmp <- map_dfr(
        .x = files,
        .f = ~ read_table(paste0(paste0(data_dir, cont, "/"), .x)) %>% 
            dplyr::select(contains("Volume")) %>%
            pivot_longer(1:ncol(.), names_to = "IDX", values_to = "volume") %>%
            mutate(IDX = str_remove(IDX, "Volume_") %>% as.numeric) %>%
            left_join(df_labels %>% dplyr::select(IDX, LABEL), by = join_by(IDX)) %>%
            mutate(sub_ses = .x %>% str_remove(".txt"), contrast = cont, .before = 1) %>% 
            dplyr::select(sub_ses, contrast, IDX, LABEL, volume)
    )
    
    # Combine with volume data from other contrasts
    df_volumes <- df_volumes %>% bind_rows(df_tmp)
    
}

df_volumes %>% 
    count(sub_ses) # 108 * 3 contrasts each = 324 total scans


# format and and combine volume data with scan info --------------------------------------

df_volumes_formatted <- df_volumes %>% 
    separate(sub_ses, into = c("subject_id", "ses"), sep = "_") %>% 
    mutate(subject_id = str_remove(subject_id, "sub-")) %>% 
    separate(subject_id, into = c("study", "id"), sep = "-") %>% 
    
    # add scan date
    left_join(df_scan_list) %>% 
    dplyr::select(study, id, ses, scan_date, contrast, IDX, LABEL, volume)

df_volumes_formatted %>% filter(study == "B3526")
df_scan_list %>% filter(study == "B3526")

# filter out scans that failed registration! ------------------------------

# Last updated: 21 March 2024
qc_round <- 4
df_qc <- read_xlsx(paste0(base_dir, "registration_QC/qc", qc_round, "_review.xlsx")) %>% 
    mutate(subject = str_remove(subject, "sub-")) %>% 
    separate(subject, into = c("study", "id"), sep = "-") %>% 
    dplyr::rename("ses" = "session")
df_failed_scans <- df_qc %>% 
    filter(registration_success == 0) %>% 
    dplyr::select(study, id, ses, contrast)

# Filter out failed scans from formatted volume data
df_volumes_filt <- df_volumes_formatted %>% 
    anti_join(df_failed_scans)

df_volumes %>% count(sub_ses, contrast) %>% print(n = nrow(.))
df_volumes_formatted %>% count(study, id, ses, contrast) %>% print(n = nrow(.))
df_volumes_filt %>% count(study, id, ses, contrast) %>% print(n = nrow(.))

# check for outliers with quick exploratory plots -------------------------

# Histograms
df_volumes_filt %>% 
    filter(contrast == "PDw" & study == "B3526") %>% 
    ggplot(aes(x = volume)) +
    geom_histogram() +
    facet_wrap(vars(LABEL), scales = "free")

# PCA
pca_res <- df_volumes_filt %>% 
    filter(contrast == "PDw" & study == "B3526") %>% 
    pivot_wider(id_cols = c(id), names_from = LABEL, values_from = volume) %>% 
    clean_names %>% 
    
    # impute missing values with median
    mutate(across(where(is.numeric), 
                  ~ replace(., is.na(.), median(., na.rm = TRUE)))
    ) %>% 
    # run PCA
    column_to_rownames("id") %>%
    prcomp(center = TRUE, scale = TRUE)

# plot subjects
df_pca_subjects <- pca_res$x %>% 
    as.data.frame() %>% 
    rownames_to_column("id") %>% 
    as_tibble() %>% 
    left_join(df_volumes_filt %>% 
                  filter(study == "B3526") %>% 
                  dplyr::select(id, scan_date) %>% 
                  distinct) %>% 
    dplyr::select(id, scan_date, everything()) %>% 
    filter(!is.na(scan_date))

df_pca_subjects %>% 
    ggplot(aes(x = PC1, y = PC2)) +
    geom_point(aes(color = scan_date)) +
    geom_text_repel(aes(label = id))

cor.test(df_pca_subjects$PC1, as.numeric(df_pca_subjects$scan_date))
cor.test(df_pca_subjects$PC2, as.numeric(df_pca_subjects$scan_date))

# plot ROIs
pca_res$rotation %>% 
    as.data.frame() %>% 
    rownames_to_column("LABEL") %>% 
    as_tibble() %>% 

    ggplot(aes(x = PC1, y = PC2)) +
    geom_point() +
    geom_text_repel(aes(label = LABEL))



# Export to spreadsheets to send ------------------------------------------


# B3526
df_volumes_filt %>% 
    filter(study == "B3526" & contrast == "PDw") %>% 
    pivot_wider(id_cols = c(id, scan_date, ses), names_from = LABEL, values_from = volume) %>% 
    
    write_xlsx(paste0(base_dir, "outputs/B3526_PDw_volumes.xlsx"))

# B3526
df_volumes_filt %>% 
    filter(study == "B31007" & contrast == "PDw") %>% 
    pivot_wider(id_cols = c(id, scan_date, ses), names_from = LABEL, values_from = volume) %>% 
    
    write_xlsx(paste0(base_dir, "outputs/B31007_PDw_volumes.xlsx"))




