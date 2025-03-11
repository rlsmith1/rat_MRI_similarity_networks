#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Calculate case-control differences in systems-level nodal degree following RMS using LME model,
# where the beta coefficient for group is the effect size, accounting for random effects by ROI
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


# Setup -------------------------------------------------------------------

## Set directories & load data
base_dir <- "~/Documents/PhD/projects/CamRat/CamRat/"
source(paste0(base_dir, "code/03.data_analysis/setup.R"))
analysis_objects_dir <- paste0(base_dir, "outputs/objects/") # to save analysis objects for figures

## Identify systems that are comprised of only one ROI
singletons <- df_system_hierarchy %>% 
  count(system, sort = TRUE) %>% 
  filter(n == 1) %>% 
  pull(system)

## Write function to calculate z-score (for later)
f_calc_zscore <- function(x) {(x - mean(x))/sd(x)}


# Calculate actual case-control effect sizes across systems --------------------------------------------------------


## Run models to identify case-control differences at PND 63 (following RMS)

df_node_stress_system <- df_data_cortex %>%
  filter(study == "GSK" & timepoint == 63) %>% # only looking at the young adult timepoint following RMS
  
  # add system-level information
  left_join(df_system_hierarchy) %>% 
  group_by(timepoint, system, feature) %>% 
  nest() %>% 
  
  # run models depending on number of regions in system
  mutate(
    model_res = ifelse(
      system %in% singletons,
      map(
        .x = data,
        .f = ~ f_lm_model(
          dataframe = .x,
          independentVar = "group",
          dependentVar = "value",
          covariates = c("tbv", "sex", "age")
        )$results 
      ),
      map(
        .x = data,
        .f = ~ f_lme_model(
          dataframe = .x,
          independentVar = "group",
          dependentVar = "value",
          covariates = c("tbv", "sex", "age"),
          randomEffects = c("region_of_interest")
        )$results
      )
    ) 
  ) %>% 
  
  # format results
  unnest(cols = c(model_res)) %>% 
  arrange(timepoint, t_value) %>% 
  filter(term != "(Intercept)")

df_node_stress_system %>% 
  filter(term == "groupMS") %>% 
  pivot_wider(id_cols = system, names_from = feature, values_from = t_value) %>% 
  arrange(degree)



# Run permutations: mix up group label and re-calculate case-control effect size --------------------


## Method: Mix up group label vector and calculate case-control effect size for each ROI

## Set up permutation
n_perm <- 1000
df_node_stress_system_permAll <- tibble()
doParallel::registerDoParallel()
set.seed(0927)

## Run permutation in loop
for (i in 1:n_perm) {
  
  if (i %% 5 == 0) {print(paste0("perm ", i))}
  
  # Run model using on permuted group assignments
  df_tmp <- df_data_cortex %>%
    filter(study == "GSK" & timepoint == 63) %>% # only looking at the young adult timepoint following RMS
    
    ## PERMUTE GROUP HERE ##
    mutate(group = sample(group)) %>% 
    
    # add system-level information
    left_join(df_system_hierarchy) %>% 
    group_by(timepoint, system, feature) %>% 
    nest() %>% 
    
    # run models depending on number of regions in system
    mutate(
      model_res = ifelse(
        system %in% singletons,
        map(
          .x = data,
          .f = ~ f_lm_model(
            dataframe = .x,
            independentVar = "group",
            dependentVar = "value",
            covariates = c("tbv", "sex", "age")
          )$results 
        ),
        map(
          .x = data,
          .f = ~ f_lme_model(
            dataframe = .x,
            independentVar = "group",
            dependentVar = "value",
            covariates = c("tbv", "sex", "age"),
            randomEffects = c("region_of_interest")
          )$results
        )
      ) 
    ) %>% 
    
    # format results
    unnest(cols = c(model_res)) %>% 
    arrange(timepoint, t_value) %>% 
    filter(term != "(Intercept)")
  
  # bind rows with tibble
  df_node_stress_system_permAll <- df_node_stress_system_permAll %>% 
    bind_rows(
      df_tmp %>% 
        mutate(perm = i, .before = 1)
    )
  
}


## Calculate the z-score of the *actual* effect size relative to the permuted null distribution
df_node_stress_system_perm <- df_node_stress_system %>% 
  ungroup %>% 
  filter(term == "groupMS") %>% 
  dplyr::select(feature, system, t_value) %>% 
  mutate(perm = "actual", .before = 1) %>% 
  bind_rows(
    df_node_stress_system_permAll %>% 
      filter(term == "groupMS") %>% 
      dplyr::select(perm, feature, system, t_value) %>% 
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
  df_node_stress_system, df_node_stress_system_perm,
  file = paste0(analysis_objects_dir, "27Sept2024_node_RMS_effects_SYS.RData")
)
