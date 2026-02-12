
#==============================================================================#
# Calculate developmental slope of systems-level edge weight slope in early and later life using LME model,
# where the beta coefficient for age is the slope, accounting for random effects by subject & ROI-level edge
#==============================================================================#

### SETUP ###

base_dir <- "~/Documents/PhD/projects/CamRat/CamRat/"
source(paste0(base_dir, "code/03_data-analysis/setup.R"))


### FORMAT DATA ###

## Add system-level info to MIND network
df_mind_cortex_system <- df_mind_cortex %>% 
  filter(R1 != R2 & str_detect(subject, "JWD")) %>% # only considering MRC (normative development) study for this analysis
  
  ## remove identical edges in reverse order
  mutate(edge = paste0(pmin(R1, R2), " - ", pmax(R1, R2))) %>% 
  group_by(subject, timepoint, edge) %>%
  slice(1) %>% 
  
  ## add system-level connection info
  left_join(df_system_hierarchy %>% dplyr::select(-system2), by = join_by("R1" == "region_of_interest")) %>% 
  dplyr::rename("S1" = "system") %>% 
  left_join(df_system_hierarchy %>% dplyr::select(-system2), by = join_by("R2" == "region_of_interest")) %>% 
  dplyr::rename("S2" = "system") %>% 
  mutate(system_edge = paste0(pmin(S1, S2), " - ", pmax(S1, S2)))



### RUN MODELS ###

## Identify system_edges that are comprised of only one ROI-level edge
edge_singletons <- df_mind_cortex_system %>% 
  ungroup() %>% 
  dplyr::select(system_edge, edge) %>% 
  distinct() %>% 
  count(system_edge, sort = TRUE) %>% 
  filter(n == 1) %>% 
  pull(system_edge)

## Run LM/LME on each system-level edge for each combination of timepoints
doParallel::registerDoParallel()
df_edge_slopes_system <- df_mind_cortex_system %>% 

  # define early development as PND 20 --> 35, and aging as PND 63 --> 230
  mutate(period = ifelse(timepoint %in% c(20, 35), "early", "aging")) %>% 
  group_by(system_edge, period) %>% 
  nest() %>% 
  
  # run models (LME with region_of_interest as random effect if there is more than one region in the system, regular lm if not)
  mutate(
    model_res = 
    ifelse(
      system_edge %in% edge_singletons,
      map(
        .x = data,
        .f = ~ f_lme_model(
          dataframe = .x,
          independentVar = "age",
          dependentVar = "weight",
          covariates = "tbv",
          randomEffects = c("subject")
        )$results
      ),
      map(
        .x = data,
        .f = ~ f_lme_model(
          dataframe = .x,
          independentVar = "age",
          dependentVar = "weight",
          covariates = "tbv",
          randomEffects = c("edge", "subject")
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



### STATS ###

## Correlate early development and aging slopes
df_edge_slopes_wide <- df_edge_slopes_system %>%
    filter(term == "age") %>% 
    pivot_wider(id_cols = system_edge, names_from = period, values_from = t_value) %>% 
    arrange(-early, aging)
cor.test(df_edge_slopes_wide$early, df_edge_slopes_wide$aging)



### SAVE ###

save(
  df_edge_slopes_system,
  file = paste0(objects_dir, "27Sept2024_edge_slopes_SYS.RDS")
)
