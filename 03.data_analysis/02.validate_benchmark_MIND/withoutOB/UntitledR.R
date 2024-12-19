base_dir <- "~/Documents/PhD/projects/CamRat/CamRat/"
data_dir <- paste0(base_dir, "data/AMBA/")

df_merfish_expr <- read.csv(paste0(data_dir, "2024Feb28_MERFISH_exp_mat.csv")) %>% as_tibble()

read_csv(paste0(data_dir, "MERFISH_cell_megadata_with_parcellation_annotation.csv"))
