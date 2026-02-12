
#==============================================================================#
# Calculate developmental slope of systems-level edge weight slope in early and later life using LME model,
# where the beta coefficient for age is the slope, accounting for random effects by subject & ROI-level edge
#==============================================================================#

### SETUP ###

base_dir <- "~/Documents/PhD/projects/CamRat/CamRat/"
source(paste0(base_dir, "code/03_data-analysis/setup.R"))


### RUN MODELS ###

## Run LM/LME on each ROI-level edge for each combination of timepoints
doParallel::registerDoParallel()
df_edge_slopes_roi <- df_mind_cortex %>%
  filter(R1 != R2) %>% 
  
  # define early development as PND 20 --> 35, and aging as PND 63 --> 230
  mutate(period = ifelse(timepoint %in% c(20, 35), "early", "aging")) %>% 
  group_by(edge, period) %>% 
  nest() %>% 
  
  # run models (LME with region_of_interest as random effect if there is more than one region in the system, regular lm if not)
  mutate(
    model_res = 
      map(
        .x = data,
        .f = ~ f_lme_model(
          dataframe = .x,
          independentVar = "age",
          dependentVar = "weight",
          covariates = "tbv",
          randomEffects = c("subject")
        )$results
      )
  ) %>% 
  
  # format results
  unnest(cols = c(model_res)) %>% 
  mutate(period = factor(period, levels = c("early", "aging"))) %>% 
  dplyr::rename("slope" = "estimate") %>% 
  arrange(period, -abs(slope)) %>% 
  mutate(significant = ifelse(abs(t_value) > 2, 1, 0) %>% as.factor) %>% 
  filter(term != "(Intercept)")

df_edge_slopes_roi %>%
  filter(term == "age") %>% 
  pivot_wider(id_cols = edge, names_from = period, values_from = t_value) %>% 
  arrange(-early, aging)


### SAVE ###

save(
  df_edge_slopes_roi,
  file = paste0(objects_dir, "27Sept2024_edge_slopes_ROI.RDS")
)
