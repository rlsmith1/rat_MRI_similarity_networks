
################################################################################

## Read MIND output files and combine into df (with metadata)

################################################################################


# setup ---------------------------------------------------------------

## LIBRARIES
library(tidyverse)

## SET DIRECTORIES
base_dir <- "~/Documents/PhD/projects/CamRat/CamRat/"
data_dir <- paste0(base_dir, "data/WHS_analysis/MIND_files/MTR_scaled/cerebral_cortex/output/")

## LOAD DATA
load(paste0(base_dir, "objects/df_metadata.RDS")) # df_metadata (generated in combine_metadata.R)


# Read & format MIND nets -------------------------------------------------


files <- list.files(data_dir)
df_mind_cortex <- map_dfr(
  .x = files,
  .f = ~ read_csv(paste0(data_dir, .x)) %>% 
    mutate(R1 = colnames(.), .before = 1) %>% 
    pivot_longer(2:ncol(.), names_to = "R2", values_to = "weight") %>% 
    mutate(subject = .x %>% 
             str_remove("_.*") %>% 
             str_remove("sub-"),
           timepoint = .x %>% 
             str_remove(".*_") %>% 
             str_remove(".csv") %>% 
             str_remove("ses-PND") %>% 
             as.numeric %>% 
             factor(levels = c(20, 35, 63, 300)),
           .before = 1
    )
) %>%
  
  mutate(R1 = str_trim(R1), R2 = str_trim(R2)) %>% 
  
  # define edges
  mutate(edge = paste0(pmin(R1, R2), sep = " - ", pmax(R1, R2))) %>% 
  
  # add metadata, including only subjects/timepoints included in df_data_format
  inner_join(df_metadata %>% 
               ungroup %>% 
               dplyr::select(subject, timepoint, age, sex, group, tbv, scan_date) %>% 
               distinct(),
             by = join_by(subject, timepoint)
  )

## SAVE
save(df_mind_cortex, file = paste0(base_dir, "objects/03June2024_df_mind_cortex.RDS")) # df_mind_cortex




#### an aside (to delete later) -----------------------------------------------


## Correlate with old MIND network (with unscaled MTR data)
df_mind_cortex_scaled <- df_mind_cortex
load(paste0(base_dir, "objects/archive/29Apr2024_df_mind_cortex.RDS"))

df_compare_mind_nets <- df_mind_cortex %>% 
  dplyr::select(subject, timepoint, R1, R2, edge, weight) %>% 
  filter(!str_detect(edge, regex("olfactory", ignore_case = TRUE))) %>% 
  left_join(
    df_mind_cortex_scaled %>% 
      dplyr::select(subject, timepoint, R1, R2, edge, weight) %>% 
      dplyr::rename("weight_scaled" = "weight")
  ) %>% 
  mutate(study = ifelse(str_detect(subject, "EDA"), "GSK", "MRC"), .before = 1) %>% 
  distinct()

df_compare_mind_nets_median <- df_compare_mind_nets %>% 
  group_by(study, timepoint, R1, R2, edge) %>% 
  summarise_if(is.numeric, ~ median(.x, na.rm = TRUE))

# Pearson & cophenetic correlations
f_calc_coph_cor <- function(df_net1, df_net2) {
  
  ## Generate adjacency matrices for both networks
  net1 <- df_net1 %>% 
    arrange(R1, R2) %>% 
    pivot_wider(id_cols = R1, names_from = R2, values_from = weight) %>% 
    column_to_rownames("R1")
  net2 <- df_net2 %>% 
    arrange(R1, R2) %>% 
    pivot_wider(id_cols = R1, names_from = R2, values_from = weight) %>% 
    column_to_rownames("R1")
  
  # Convert adjacency matrices to distance matrices
  # Perform hierarchical clustering
  # Compute cophenetic distance matrices
  net1_coph <- hclust(as.dist(1 - net1), method = "average") %>% cophenetic
  net2_coph <- hclust(as.dist(1 - net2), method = "average") %>% cophenetic
  
  ## Calculate cophenetic correlation coefficient between two networks
  coph_cor <- cor.test(net1_coph, net2_coph)
  
  ## Return cophenetic correlation test results
  return(coph_cor)
  
}

df_compare_mind_nets_median %>% 
  group_by(study, timepoint) %>% 
  nest() %>%
  mutate(
    pearsons_r = map(
      .x = data,
      .f = ~ cor.test(.x$weight, .x$weight_scaled)$estimate
    )
  ) %>% 
  unnest(cols = c(pearsons_r))

