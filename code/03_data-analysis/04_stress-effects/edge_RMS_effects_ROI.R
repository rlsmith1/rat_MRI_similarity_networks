
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Calculate case-control differences in ROI-level edge weight following RMS using linear model,
# where the beta coefficient for group is the effect size
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


# Setup -------------------------------------------------------------------

## Set directories & load data
base_dir <- "~/Documents/PhD/projects/CamRat/CamRat/"
source(paste0(base_dir, "code/03.data_analysis/setup.R"))
analysis_objects_dir <- paste0(base_dir, "outputs/objects/") # to save analysis objects for figures

## Write function to calculate z-score (for later)
f_calc_zscore <- function(x) {(x - mean(x))/sd(x)}



# Calculate actual case-control effect sizes across roi-level edges -----------------


## Run LM/LME on each roi-level edge for each combination of timepoints
doParallel::registerDoParallel()
df_edge_stress_roi <- df_mind_cortex %>%
  filter(R1 != R2 & str_detect(subject, "EDA") & timepoint == 63) %>% # Use Experimental stress cohort only for this analysis; timepoint 63  
  group_by(timepoint, edge) %>% 
  nest() %>% 
  
  # run linear model on each ROI edge
  mutate(
    model_res = map(
      .x = data,
      .f = ~ f_lm_model(
        dataframe = .x,
        independentVar = "group",
        dependentVar = "weight",
        covariates = c("tbv", "sex", "age")
      )$results 
    ) 
  ) %>% 
  
  # format results
  unnest(cols = c(model_res)) %>% 
  arrange(timepoint, t_value) %>% 
  filter(term != "(Intercept)")

df_edge_stress_roi %>%
  filter(term == "groupMS") %>% 
  dplyr::select(edge, t_value) %>% 
  arrange(t_value)


# Run permutations: mix up group label and re-calculate case-control effect size --------------------


## Method: Mix up group label vector and calculate case-control effect size for each ROI

## Set up permutation
n_perm <- 1000
df_edge_stress_roi_permAll <- tibble()
doParallel::registerDoParallel()
set.seed(0927)

## Run permutation in loop
for (i in 1:n_perm) {
  
  if (i %% 5 == 0) {print(paste0("perm ", i))}
  
  # Run model using on permuted group assignments
  df_tmp <- df_mind_cortex %>%
    filter(R1 != R2 & str_detect(subject, "EDA") & timepoint == 63) %>% # Use Experimental stress cohort only for this analysis; timepoint 63  
    
    ## PERMUTE GROUP HERE ##
    mutate(group = sample(group)) %>% 
    
    # group by ROI-level edge to run model on each
    group_by(timepoint, edge) %>% 
    nest() %>% 
    
    # run linear model on each ROI edge
    mutate(
      model_res = map(
        .x = data,
        .f = ~ f_lm_model(
          dataframe = .x,
          independentVar = "group",
          dependentVar = "weight",
          covariates = c("tbv", "sex", "age")
        )$results 
      ) 
    ) %>% 
    
    # format results
    unnest(cols = c(model_res)) %>% 
    arrange(timepoint, t_value) %>% 
    filter(term != "(Intercept)")
  
  # bind rows with tibble
  df_edge_stress_roi_permAll <- df_edge_stress_roi_permAll %>% 
    bind_rows(
      df_tmp %>% 
        mutate(perm = i, .before = 1)
    )
  
}


## Calculate the z-score of the *actual* effect size relative to the permuted null distribution
df_edge_stress_roi_perm <- df_edge_stress_roi %>% 
  ungroup %>% 
  filter(term == "groupMS") %>% 
  dplyr::select(edge, t_value) %>% 
  mutate(perm = "actual", .before = 1) %>% 
  bind_rows(
    df_edge_stress_roi_permAll %>% 
      filter(term == "groupMS") %>% 
      dplyr::select(perm, edge, t_value) %>% 
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
  df_edge_stress_roi, df_edge_stress_roi_perm,
  file = paste0(analysis_objects_dir, "27Sept2024_edge_RMS_effects_ROI.RData")
)
