
#==============================================================================#
# Calculate developmental slope of systems-level nodal degree slope in early and later life using LME model,
# where the beta coefficient for age is the slope, accounting for random effects by subject & ROI
#==============================================================================#

### SETUP ###

base_dir <- "~/Documents/PhD/projects/CamRat/CamRat/"
source(paste0(base_dir, "code/03_data-analysis/setup.R"))


### RUN MODELS ###

## Identify systems that are comprised of only one ROI
singletons <- df_system_hierarchy %>% 
  count(system, sort = TRUE) %>% 
  filter(n == 1) %>% 
  pull(system)

## Run LM/LME on each feature for each combination of timepoints
df_node_slopes_system <- df_data_cortex %>% 
  filter(study == "MRC") %>% # not running this on GSK cohort because we're ignoring PND 20 scans
  
  # define early development as PND 20 --> 35, and aging as PND 63 --> 230
  mutate(period = ifelse(timepoint %in% c(20, 35), "early", "aging")) %>% 
  
  # add system-level information
  left_join(df_system_hierarchy) %>% 
  group_by(system, feature, period) %>% 
  nest() %>% 
  
  # run models (LME with region_of_interest as random effect if there is more than one region in the system, regular lm if not)
  mutate(
    model_res = ifelse(
      system %in% singletons,
      map(
        .x = data,
        .f = ~ f_lme_model(
          dataframe = .x,
          independentVar = "age",
          dependentVar = "value",
          covariates = "tbv",
          randomEffects = c("subject")
        )$results 
      ),
      map(
        .x = data,
        .f = ~ f_lme_model(
          dataframe = .x,
          independentVar = "age",
          dependentVar = "value",
          covariates = "tbv",
          randomEffects = c("region_of_interest", "subject")
        )$results
      )
    ) 
  ) %>% 
  
  # format results
  unnest(cols = c(model_res)) %>% 
  mutate(period = factor(period, levels = c("early", "aging"))) %>% 
  dplyr::rename("slope" = "estimate") %>% 
  arrange(period, -abs(slope)) %>% 
  mutate(significant = ifelse(abs(t_value) > 2, 1, 0) %>% as.factor) %>% 
  filter(term != "(Intercept)")

df_node_slopes_system %>%
  filter(feature == "degree" & term == "age") %>% 
  pivot_wider(id_cols = system, names_from = period, values_from = t_value) %>% 
  arrange(-early, aging)

### SAVE ###

save(
  df_node_slopes_system,
  file = paste0(objects_dir, "27Sept2024_node_slopes_SYS.RDS")
)
