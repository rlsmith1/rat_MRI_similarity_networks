
#==============================================================================#
# Generate null networks considering weight bins for subsequent validation analyses
#==============================================================================#

## NOTE: Run after 00_edge-distance.R, because we use distance bins for permutations to preserve network structure

## Load distance data for permutation bins
base_dir <- "~/Documents/PhD/projects/CamRat/CamRat/"
source(paste0(base_dir, "code/03_data-analysis/setup.R"))
load(paste0(objects_dir, "08July2024_normative_mind_distance_withOB.RDS")) # df_normative_mind_distance_withOB

## Generate 10000 null networks to compare validation results against
set.seed(08072024)
df_null_nets_withOB <- df_normative_mind_distance_withOB %>% 
  dplyr::select(-edge) %>% 
  expand_grid(network = paste0("null", 1:10000)) %>% 
  group_by(network, distance_bin) %>% 
  nest() %>% 
  mutate(
    data = map(
      .x = data,
      .f = ~ .x %>% mutate(weight = sample(median_weight)) %>% dplyr::select(-sd_weight, -z_score, -median_weight)
    )
  ) %>% 
  unnest(cols = c(data)) %>% 
  arrange(network, distance_bin)

## Save for benchmarking scripts
save(df_null_nets_withOB, 
     file = paste0(objects_dir, "08July2024_null_nets_withOB.RDS")
)
