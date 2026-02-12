
#==============================================================================#
# Calculate the median PND 63 network in GSK control male animals 
# for comparison with normative cohort
#==============================================================================#

## Set up
base_dir <- "~/Documents/PhD/projects/CamRat/CamRat/"
source(paste0(base_dir, "code/03_data-analysis/setup.R"))


# Define the control PND 63 connectome in GSK cohort ----------------------


## Define the network without the OB
## Include PND 63 control subjects only

df_experimental_mind <- df_mind_cortex %>%
  filter(timepoint == 63 & str_detect(subject, "EDA") & group == "control") %>%
  group_by(R1, R2) %>%
  summarise(median_weight = median(weight, na.rm = TRUE),
            sd_weight = sd(weight, na.rm = TRUE),
            z_score = median_weight/sd_weight
  ) %>% 
  arrange(-median_weight) %>% 
  
  # add system-level info
  left_join(df_system_hierarchy %>% dplyr::rename_all( ~ c("S1", "R1"))) %>% 
  left_join(df_system_hierarchy %>% dplyr::rename_all( ~ c("S2", "R2"))) %>% 
  dplyr::select(starts_with("S", ignore.case = FALSE), everything()) %>% 
  ungroup

## Align & combine with MRC normative network
df_norm_exp_mind <- df_experimental_mind %>% 
  dplyr::select(-sd_weight, -z_score) %>% 
  dplyr::rename("exp_weight" = "median_weight") %>% 
  
  left_join(
    df_normative_mind %>% 
      dplyr::select(-sd_weight, -z_score) %>% 
      dplyr::rename("norm_weight" = "median_weight")
  ) %>% 
  
  # remove self-connections (skews network correlations)
  filter(R1 != R2)

## Edge correlation P-value
exp_norm_cor <- cor.test(df_norm_exp_mind$norm_weight, df_norm_exp_mind$exp_weight)
t_val <- exp_norm_cor$statistic
df_val <- exp_norm_cor$parameter
pt(t_val, df = df_val, lower.tail = FALSE)



# Calculate median strength for PND 63 control GSK subjects --------------------------


## Identify control male subjects included in GSK MIND network
gsk_control_subjs <- df_mind_cortex %>%
  filter(timepoint == 63 & str_detect(subject, "EDA") & group == "control") %>% 
  pull(subject) %>% unique

## Take median across normative individuals to define control strength in the experimental cohort
df_experimental_strength <- df_strength_cortex %>% 
  filter(subject %in% gsk_control_subjs & timepoint == 63) %>% 
  group_by(system, region_of_interest) %>% 
  summarise(median_value = median(strength),
            sd_value = sd(strength),
            z_score = median_value/sd_value
  ) %>% 
  arrange(-median_value) %>% 
  ungroup

## Align and combine with MRC normative MIND strength values
df_norm_exp_strength <- df_experimental_strength %>% 
  dplyr::select(-sd_value, -z_score) %>% 
  dplyr::rename("exp_strength" = "median_value") %>% 
  
  left_join(
    df_normative_strength %>% 
      dplyr::select(-sd_value, -z_score) %>% 
      dplyr::rename("norm_strength" = "median_value")
  )

## Strength correlation P-value
exp_norm_cor <- cor.test(df_norm_exp_strength$exp_strength, df_norm_exp_strength$norm_strength)
t_val <- exp_norm_cor$statistic
df_val <- exp_norm_cor$parameter
pt(t_val, df = df_val, lower.tail = FALSE)




# Correlate each network with each other network ----------------


## Define sex & group metadata
df_subject_metadata <- df_data_cortex %>% 
  dplyr::select(subject, timepoint, study, sex, group) %>% 
  distinct() %>% 
  filter(!(str_detect(subject, "EDAA") & timepoint == 20))

## Create data frame of all subject-timepoint combinations
df_subj_timepoint <- df_mind_cortex %>% 
  dplyr::select(subject, timepoint) %>% 
  distinct() %>% 
  filter(!(str_detect(subject, "EDAA") & timepoint == 20)) # remove GSK PND20 scans (not included in these analyses)
df_combos <- df_subj_timepoint %>% 
  dplyr::rename_all( ~ paste0(.x, 1)) %>% 
  expand_grid(
    df_subj_timepoint %>% 
      dplyr::rename_all( ~ paste0(.x, 2))
  ) %>% 
  filter(paste0(subject1, timepoint1) != paste0(subject2, timepoint2)) %>% 
  
  # Remove duplicates in reverse order
  mutate(id1 = paste0(subject1, " PND", timepoint1),
         id2 = paste0(subject2, " PND", timepoint2)
  ) %>% 
  group_by(id = paste0(pmin(id1, id2), sep = "-", pmax(id1, id2))) %>% 
  slice(1) %>% 
  ungroup %>% 
  dplyr::select(-contains("id"))

## Write function to calculate the cophenetic correlation between two networks 
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

## Calculate correlations
doParallel::registerDoParallel()
df_cor_all_nets <- df_combos %>%
  mutate(
    
    # Pearson correlation
    pearson_cor_test = pmap(
      .l = list(..1 = subject1, ..2 = timepoint1, ..3 = subject2, ..4 = timepoint2),
      .f = ~ cor.test(
        df_mind_cortex %>% filter(subject == ..1 & timepoint == ..2) %>% arrange(edge) %>% pull(weight),
        df_mind_cortex %>% filter(subject == ..3 & timepoint == ..4) %>% arrange(edge) %>% pull(weight)
      )
    ),
    pearson_r = map(
      .x = pearson_cor_test,
      .f = ~ .x$estimate
    ),
    pearson_p_val = map(
      .x = pearson_cor_test,
      .f = ~ .x$p.value
    ),
    
    # Cophenetic correlation (considers network structure)
    cophenetic_cor_test = pmap(
      .l = list(..1 = subject1, ..2 = timepoint1, ..3 = subject2, ..4 = timepoint2),
      .f = ~ f_calc_coph_cor(
        df_mind_cortex %>% filter(subject == ..1 & timepoint == ..2),
        df_mind_cortex %>% filter(subject == ..3 & timepoint == ..4)
      )
    ),
    cophenetic_cor = map(
      .x = cophenetic_cor_test,
      .f = ~ .x$estimate
    ),
    cophenetic_p_val = map(
      .x = cophenetic_cor_test,
      .f = ~ .x$p.value
    )
    
  ) %>% 
  unnest(cols = c(pearson_r, pearson_p_val, cophenetic_cor, cophenetic_p_val)) %>% 
  
  # add back metadata
  left_join(df_subject_metadata %>% dplyr::rename_all( ~ paste0(.x, 1))) %>%
  left_join(df_subject_metadata %>% dplyr::rename_all( ~ paste0(.x, 2))) %>% 
  arrange(-cophenetic_cor) %>% 
  
  # Remove subject-timepoint combos that did not pass QC filtering
  filter(!is.na(study1) & !is.na(study2) & !(study1 == "GSK" & timepoint1 == 20) & !(study2 == "GSK" & timepoint2 == 20)) %>% 
  
  # Change names/organization for plotting
  mutate(group1 = ifelse(group1 == "MS", "RMS", group1)) %>% 
  mutate(group2 = ifelse(group2 == "MS", "RMS", group2))
  



# Save objects for plotting ---------------------------------------

save(
  df_norm_exp_mind, df_norm_exp_strength, df_cor_all_nets,
  file = paste0(objects_dir, "07June2024_GSK_MIND_comparison.RData")
)

