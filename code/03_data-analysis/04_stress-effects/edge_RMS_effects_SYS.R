#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Calculate case-control differences in systems-level edge weight following RMS using LME model,
# where the beta coefficient for group is the effect size, accounting for random effects by ROI-level edge
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



# Setup -------------------------------------------------------------------

## Set directories & load data
base_dir <- "~/Documents/PhD/projects/CamRat/CamRat/"
source(paste0(base_dir, "code/03.data_analysis/setup.R"))
analysis_objects_dir <- paste0(base_dir, "outputs/objects/") # to save analysis objects for figures

## Write function to calculate z-score (for later)
f_calc_zscore <- function(x) {(x - mean(x))/sd(x)}

## Add system-level info to MIND network
df_mind_cortex_system <- df_mind_cortex %>% 
  filter(R1 != R2 & str_detect(subject, "EDA")) %>% # Use Experimental stress cohort only for this analysis
  
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

## Identify system_edges that are comprised of only one ROI-level edge
edge_singletons <- df_mind_cortex_system %>% 
  ungroup() %>% 
  dplyr::select(system_edge, edge) %>% 
  distinct() %>% 
  count(system_edge, sort = TRUE) %>% 
  filter(n == 1) %>% 
  pull(system_edge)


# Calculate actual case-control effect sizes across system-level edges -----------------


## Run LM/LME on each system-level edge for each combination of timepoints
doParallel::registerDoParallel()
df_edge_stress_system <- df_mind_cortex_system %>%
  filter(timepoint == 63) %>% # only looking at the young adult timepoint following RMS
  group_by(timepoint, system_edge) %>% 
  nest() %>% 
  
  # run models depending on number of regions in system
  mutate(
    model_res = ifelse(
      system_edge %in% edge_singletons,
      map(
        .x = data,
        .f = ~ f_lm_model(
          dataframe = .x,
          independentVar = "group",
          dependentVar = "weight",
          covariates = c("tbv", "sex", "age")
        )$results 
      ),
      map(
        .x = data,
        .f = ~ f_lme_model(
          dataframe = .x,
          independentVar = "group",
          dependentVar = "weight",
          covariates = c("tbv", "sex", "age"),
          randomEffects = c("edge")
        )$results
      )
    ) 
  ) %>% 
  
  # format results
  unnest(cols = c(model_res)) %>% 
  arrange(timepoint, t_value) %>% 
  filter(term != "(Intercept)")

df_edge_stress_system %>%
  filter(term == "groupMS:sexMale") %>% 
  dplyr::select(system_edge, t_value) %>% 
  arrange(t_value)

##
df_edge_stress_system %>%
  filter(term == "groupMS") %>% 
  ungroup %>% 
  filter(system_edge == "Hippocampal formation - Orbitofrontal cortex") %>% 
  unnest(cols = c(data)) %>% 
  
  ggplot(aes(x = group, y = weight)) +
  geom_point(aes(color = sex)) +
  geom_boxplot(aes(color = sex))
  


# Run permutations: mix up group label and re-calculate case-control effect size --------------------


## Method: Mix up group label vector and calculate case-control effect size for each ROI

## Set up permutation
n_perm <- 1000
df_edge_stress_system_permAll <- tibble()
doParallel::registerDoParallel()
set.seed(0927)

## Run permutation in loop
for (i in 1:n_perm) {
  
  if (i %% 5 == 0) {print(paste0("perm ", i))}
  
  # Run model using on permuted group assignments
  df_tmp <- df_mind_cortex_system %>%
    filter(timepoint == 63) %>% # only looking at the young adult timepoint following RMS
    
    ## PERMUTE GROUP HERE ##
    mutate(group = sample(group)) %>% 
    
    # group by system_edge to run each model on
    group_by(timepoint, system_edge) %>% 
    nest() %>% 

    # run models depending on number of regions in system
    mutate(
      model_res = ifelse(
        system_edge %in% edge_singletons,
        map(
          .x = data,
          .f = ~ f_lm_model(
            dataframe = .x,
            independentVar = "group",
            dependentVar = "weight",
            covariates = c("tbv", "sex", "age")
          )$results 
        ),
        map(
          .x = data,
          .f = ~ f_lme_model(
            dataframe = .x,
            independentVar = "group",
            dependentVar = "weight",
            covariates = c("tbv", "sex", "age"),
            randomEffects = c("edge")
          )$results
        )
      ) 
    ) %>% 
    
    # format results
    unnest(cols = c(model_res)) %>% 
    arrange(timepoint, t_value) %>% 
    filter(term != "(Intercept)")
  
  # bind rows with tibble
  df_edge_stress_system_permAll <- df_edge_stress_system_permAll %>% 
    bind_rows(
      df_tmp %>% 
        mutate(perm = i, .before = 1)
    )
  
}


## Calculate the z-score of the *actual* effect size relative to the permuted null distribution
df_edge_stress_system_perm <- df_edge_stress_system %>% 
  ungroup %>% 
  filter(term == "groupMS") %>% 
  dplyr::select(system_edge, t_value) %>% 
  mutate(perm = "actual", .before = 1) %>% 
  bind_rows(
    df_edge_stress_system_permAll %>% 
      filter(term == "groupMS") %>% 
      dplyr::select(perm, system_edge, t_value) %>% 
      mutate(perm = paste0("perm", perm))
  ) %>% 
  mutate(perm_effect_size = f_calc_zscore(t_value)) %>% 
  
  # filter for our actual case-control effect size
  filter(perm == "actual") %>% 
  dplyr::select(-perm) %>% 
  dplyr::rename("actual_effect_size" = "t_value") %>% 
  arrange(-abs(perm_effect_size))



# Save actual and permuted effect sizes -----------------------------------

save(
  df_edge_stress_system, df_edge_stress_system_perm,
  file = paste0(analysis_objects_dir, "27Sept2024_edge_RMS_effects_SYS.RData")
)
