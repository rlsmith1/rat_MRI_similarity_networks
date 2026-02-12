
##############################################################################################################

## Take another look at behavioral data from the RMS study

##############################################################################################################

base_dir <- "~/Documents/PhD/projects/CamRat/CamRat/"
source(paste0(base_dir, "code/03.data_analysis/setup.R"))
tables_dir <- paste0(base_dir, "outputs/tables/")
analysis_objects_dir <- paste0(base_dir, "outputs/objects/") 



# Load & format behavior data ------------------------------------------------------


## Read in behavior data
df_behavior <- read_xlsx(paste0(base_dir, "data/GSK_E1_all_behaviour__collapsed_wide.xlsx"), sheet = 1) %>% 
  mutate_at(c(2:ncol(.)), ~ as.numeric(.x))
df_behavior_key <- read_xlsx(paste0(base_dir, "data/GSK_E1_all_behaviour__collapsed_wide.xlsx"), sheet = 2) %>% 
  clean_names

## Combine with metadata on subjects
df_subj_id <- df_data_cortex %>% 
  filter(str_detect(subject, "EDA")) %>% 
  dplyr::select(subject, sex, group) %>% 
  distinct() %>% 
  mutate(ID = str_remove(subject, "EDAA") %>% as.numeric, .before = 1) %>% 
  mutate(group = ifelse(group == "MS", "RMS", "Control"))

df_behavior <- df_behavior %>% 
  pivot_longer(2:ncol(.), names_to = "ethan_code", values_to = "behavior_score") %>% 
  left_join(df_subj_id) %>% 
  
  # add behavior key for plotting
  left_join(df_behavior_key) %>% 
  dplyr::select(subject, sex, group, task_group, subtask, behavior_score) %>% 
  
  # remove all but one of the FRPR PR (PR4, PR8, and PR16 are all the same)
  filter(!str_detect(subtask, "PR4|PR16"))

## Rename group colors to match
names(group_cols) <- c("Control", "RMS")


# Test for group differences in behavior scores ----------------------------------------------


## Run models
df_behavior_caseControl <- df_behavior %>% 
  mutate(task = paste0(task_group, " ", subtask)) %>% 
  group_by(task) %>% 
  mutate(behavior_score_scaled = scale(behavior_score)[,1]) %>% 
  filter(!is.na(group)) %>% 
  nest() %>% 
  mutate(
    lm = map(
      .x = data,
      .f = ~ lm(behavior_score_scaled ~ group + sex, data = .x)
    ),
    lm_res = map(lm, ~ tidy(.x) %>% filter(term == "groupRMS") %>% clean_names)
  ) %>% 
  unnest(cols = c(lm_res)) %>% 
  arrange(-abs(statistic)) %>% 
  ungroup %>% 
  mutate(p_adj = p.adjust(p_value, method = "fdr"))

## Plot
df_behavior_caseControl %>% 
  filter(p_value < 0.05) %>% 
  unnest(cols = c(data)) %>% 
  
  ggplot(aes(x = group, y = behavior_score_scaled)) +
  geom_point(aes(color = group), size = 0.75, position = position_jitter(width = 0.2)) +
  geom_boxplot(aes(fill = group, color = group), alpha = 0.1, outlier.shape = NA) +
  geom_text(aes(label = paste("t = ", round(statistic, 2))), x = 2, y = 4, check_overlap = TRUE) +
  facet_wrap(vars(task)) +
  scale_color_manual(values = group_cols, guide = "none") +
  scale_fill_manual(values = group_cols, guide = "none") +
  labs(x = NULL, y = "Behavior score (scaled)")
  



# Correlate all behavior scores with regional features --------------------


df_behavior_feature_cor <- df_behavior %>% 
  left_join(df_data_cortex %>% 
              filter(str_detect(subject, "EDAA")) %>% 
              mutate(group = ifelse(group == "MS", "RMS", "Control"))
  ) %>% 
  mutate(task = paste0(task_group, " ", subtask)) %>% 
  group_by(task, region_of_interest, feature, timepoint, group) %>% 
  nest() %>% 
  filter(!is.na(region_of_interest) & timepoint %in% c(63, 300)) %>% 
  
  mutate(
    lm = map(
      .x = data,
      .f = ~ lm(behavior_score ~ value, data = .x)
    ),
    lm_res = map(lm, ~ tidy(.x) %>% filter(term == "value") %>% clean_names %>% dplyr::select(-term))
  ) %>% 
  unnest(cols = c(lm_res))

## Plot
df_behavior_feature_cor %>%
  filter(timepoint == 63 & feature == "degree") %>%
  
  group_by(region_of_interest, feature) %>% 
  mutate(p_adj = p.adjust(p_value, method = "fdr")) %>% 
  arrange(-abs(statistic)) %>% 
  mutate(region_of_interest = factor(region_of_interest, levels = roi_order),
         label = ifelse(p_value < 0.05, "*", ""),
         label = ifelse(p_adj < 0.05, "**", label)
  ) %>% 
  
  ggplot(aes(x = task, y = region_of_interest)) +
  geom_tile(aes(fill = statistic)) +
  geom_text(aes(label = label), vjust = 0.75) +
  facet_wrap(vars(group)) +
  scale_fill_gradientn(colors = RMS_effect_size_scale, limits = c(-4.4, 4.4)) +
  labs(x = NULL, y = NULL,
       title = "Relationship between regional PND 63 degree and behavior score"
  ) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))

## Calculate identity line residuals to see what effect sizes were the most different between the two groups
distance_point_to_line <- function(x, y, m, b) {
  distance <- abs(m * x - y + b) / sqrt(m^2 + 1)
  return(distance)
}
df_behavior_degree_dist <- df_behavior_feature_cor %>%
  filter(timepoint == 63 & feature == "degree") %>% 
  pivot_wider(id_cols = c(task, region_of_interest), names_from = group, values_from = statistic) %>% 
  mutate(distance = distance_point_to_line(Control, RMS, 1, 0)) %>% 
  arrange(-abs(distance))

## Plot correlation between control and RMS effect sizes
df_behavior_degree_dist %>%
  mutate(distance = case_when(
    RMS < 0 & Control < 0 ~ NA_real_,
    RMS > 0 & Control > 0 ~ NA_real_,
    TRUE ~ distance
  )
  ) %>% 

  ggplot(aes(x = Control, y = RMS)) +
  geom_hline(yintercept = 0, color = "gray") +
  geom_vline(xintercept = 0, color = "gray") +
  geom_point(aes(color = distance), size = 0.75) +
  geom_abline(lty = 2) +
  geom_smooth(method = "lm", color = "black") +
  stat_cor() +
  scale_color_gradientn(colors = RMS_effect_size_scale, limits = c(-3.6, 3.6), guide = "none") +
  labs(x = "Control", y = "RMS",
       title = "Correlation of heatmaps"
  )

## Plot comparisons with highest differences in effect size
df_behavior_degree_dist %>% 
  ungroup %>% 
  top_n(20, distance) %>% 
  ggplot(aes(x = distance, y = paste0(task, ", ", region_of_interest) %>% reorder(distance))) +
  geom_col(aes(fill = distance)) +
  scale_fill_gradientn(colors = RMS_effect_size_scale, limits = c(-3.6, 3.6), guide = "none") +
  labs(x = "Distance from effect size to y=x", y = NULL,
       title = "Task-regional degree effect sizes with \ngreatest case-control distance"
       )


## Plot case-control relationships in affected regions
sig_tasks <- df_behavior_caseControl %>% 
  filter(p_value < 0.05) %>% 
  pull(task)
df_behavior_feature_cor %>% 
  filter(timepoint == 63 & feature == "degree" & task %in% sig_tasks) %>% 
  group_by(region_of_interest, feature) %>% 
  mutate(p_adj = p.adjust(p_value, method = "fdr")) %>% 
  arrange(-abs(statistic)) %>% 
  mutate(region_of_interest = factor(region_of_interest, levels = roi_order),
         label = ifelse(p_value < 0.05, "*", ""),
         label = ifelse(p_adj < 0.05, "**", label)
  ) %>% 
  
  ggplot(aes(x = task, y = region_of_interest)) +
  geom_tile(aes(fill = statistic)) +
  geom_text(aes(label = label), vjust = 0.75) +
  facet_wrap(vars(group)) +
  scale_fill_gradientn(colors = RMS_effect_size_scale, limits = c(-4.4, 4.4)) +
  labs(x = NULL, y = NULL,
       title = "Relationship between regional PND 63 degree and behavior score"
  ) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))

## Plot comparisons with highest differences in effect size for significant tasks
df_behavior_degree_dist %>% 
  filter(task %in% sig_tasks) %>% 
  ungroup %>% 
  top_n(20, distance) %>% 
  ggplot(aes(x = distance, y = paste0(task, ", ", region_of_interest) %>% reorder(distance))) +
  geom_col(aes(fill = distance)) +
  scale_fill_gradientn(colors = RMS_effect_size_scale, limits = c(-3.6, 3.6), guide = "none") +
  labs(x = "Distance from effect size to y=x", y = NULL,
       title = "Task-regional degree effect sizes with \ngreatest case-control distance"
  )








  
