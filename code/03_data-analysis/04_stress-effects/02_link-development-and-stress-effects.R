#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Are the effects of RMS embedded in normative development of the adolescent rat brain?
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


### SETUP ###

base_dir <- "~/Documents/PhD/projects/CamRat/CamRat/"
source(paste0(base_dir, "code/03.data_analysis/setup.R"))
analysis_objects_dir <- paste0(base_dir, "outputs/objects/")

## Load objects from development and stress analyses
objects <- list.files(analysis_objects_dir)
objects <- objects[!str_detect(objects, "null|archive")]
for (obj in objects){
  print(paste0("loading ", obj, "...."))
  load(paste0(analysis_objects_dir, obj))
}



### COMBINE EARLY DEVELOPMENT AND STRESS EFFECTS ###

## System-level edge
df_edge_slopes_stress_system <- df_edge_slopes_system %>% 
  ungroup %>% 
  filter(term == "age") %>% 
  dplyr::select(period, system_edge, t_value) %>% 
  dplyr::rename("development_slope" = "t_value") %>% 
  
  left_join(
    df_edge_stress_system_perm %>% 
      ungroup %>% 
      dplyr::select(system_edge, actual_effect_size, perm_effect_size) %>% 
      dplyr::rename("RMS_effect_size" = "perm_effect_size")
  )

## System-level node
df_node_slopes_stress_system <- df_node_slopes_system %>% 
  ungroup %>% 
  filter(term == "age") %>% 
  dplyr::select(period, feature, system, t_value) %>% 
  dplyr::rename("development_slope" = "t_value") %>% 
  
  left_join(
    df_node_stress_system_perm %>% 
      ungroup %>% 
      dplyr::select(feature, system, actual_effect_size, perm_effect_size) %>% 
      dplyr::rename("RMS_effect_size" = "perm_effect_size"),
    by = join_by(system, feature)
  )

## ROI-level edge
df_edge_slopes_stress_roi <- df_edge_slopes_roi %>% 
  ungroup %>% 
  filter(term == "age") %>% 
  dplyr::select(period, edge, t_value) %>% 
  dplyr::rename("development_slope" = "t_value") %>% 
  
  left_join(
    df_edge_stress_roi %>% 
      ungroup %>% 
      filter(term == "groupMS") %>% 
      dplyr::select(edge, t_value) %>% 
      dplyr::rename("RMS_effect_size" = "t_value")
  )


### CORRELATE EDGE-LEVEL EFFECTS ###
cor_dev <- cor.test(
  df_edge_slopes_stress_system %>% filter(period == "early") %>% pull(development_slope), 
  df_edge_slopes_stress_system %>% filter(period == "early") %>% pull(RMS_effect_size), 
)$estimate
cor_age <- cor.test(
  df_edge_slopes_stress_system %>% filter(period == "aging") %>% pull(development_slope), 
  df_edge_slopes_stress_system %>% filter(period == "aging") %>% pull(RMS_effect_size), 
)$estimate


### PERMUTATIONS ###
# Permute edge-labeling of RMS effects and correlate with unpermuted dev & aging effects
# Is the observed correlation greater than 95% of permuted correlations?

## Run permutations
nperm <- 10000
set.seed(20241101)
df_stress_dev_cor_perm <- tibble(
  perm = 1:nperm,
  data = map(
    .x = perm,
    .f = ~ df_edge_slopes_system %>% 
      ungroup %>% 
      filter(term == "age") %>% 
      dplyr::select(period, system_edge, t_value) %>% 
      dplyr::rename("development_slope" = "t_value") %>% 
      
      left_join(
        df_edge_stress_system_perm %>% 
          dplyr::select(system_edge, actual_effect_size, perm_effect_size) %>% 
          dplyr::rename("RMS_effect_size" = "actual_effect_size") %>% 
          mutate(RMS_effect_size = sample(RMS_effect_size)), ## Permute here!!
        by = join_by(system_edge)
      ) %>% 
      mutate(perm = .x, .before = 1)
  )
) %>% 
  mutate(
    cor_test_dev = map(
      .x = data,
      .f = ~ cor.test(
        .x %>% filter(period == "early") %>% pull(development_slope), 
        .x %>% filter(period == "early") %>% pull(RMS_effect_size), 
      )
    ),
    cor_dev = map(
      .x = cor_test_dev,
      .f = ~ .x$estimate
    ),
    cor_test_age = map(
      .x = data,
      .f = ~ cor.test(
        .x %>% filter(period == "aging") %>% pull(development_slope), 
        .x %>% filter(period == "aging") %>% pull(RMS_effect_size), 
      )
    ),
    cor_age = map(
      .x = cor_test_age,
      .f = ~ .x$estimate
    )
  ) %>% 
  unnest(cols = c(cor_dev, cor_age))


### SAVE OBJECTS ###
save(
  cor_dev, cor_age,
  df_edge_slopes_stress_system, df_stress_dev_cor_perm,
  file = paste0(analysis_objects_dir, "27Sept2024_link_stress_and_development.RData")
)

