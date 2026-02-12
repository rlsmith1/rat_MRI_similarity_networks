
#==============================================================================#
# Relate MIND edge weight to edge distance between centroids
# Hypothesis: Regions that are closer together will have higher MIND similarity
#==============================================================================#

base_dir <- "~/Documents/PhD/projects/CamRat/CamRat/"
source(paste0(base_dir, "code/03_data-analysis/setup.R"))


## Write function to calculate Euclidean distance between two points
CalculateEuclideanDistance <- function(vect1, vect2) sqrt(sum((vect1 - vect2)^2))

## Calculate distance between each pair of WHS atlas ROI centroids (calculated in data_prep/WHS_ROI_centroids.R)
df_edge_distance_withOB <- df_mind_cortex_withOB %>% 
  ungroup %>% 
  dplyr::select(R1, R2) %>% 
  distinct() %>% 
  
  # calculate Euclidean distance between centroid of each ROI
  mutate(distance = map2(
    .x = R1,
    .y = R2,
    .f = ~ CalculateEuclideanDistance(
      vect1 = c(df_centroids %>% filter(region_of_interest == .x) %>% pull(y),
                df_centroids %>% filter(region_of_interest == .x) %>% pull(z)
      ),
      vect2 = c(df_centroids %>% filter(region_of_interest == .y) %>% pull(y),
                df_centroids %>% filter(region_of_interest == .y) %>% pull(z)
      )
    )
    
  )
  ) %>% 
  unnest(cols = c(distance))

## Combine distance data with normative MIND weights
df_normative_mind_distance_withOB <- df_normative_mind_withOB %>% 
  
  # add distance info
  left_join(df_edge_distance_withOB) %>% 
  
  # remove duplicate edges in reverse order
  mutate(R1 = as.character(R1), R2 = as.character(R2)) %>% 
  group_by(edge = paste0(pmin(R1, R2), sep = " - ", pmax(R1, R2))) %>% 
  slice(1) %>% 
  ungroup %>% 
  
  # remove self connections
  filter(R1 != R2)

## Divide edge distances into 3 evenly sized bins for permutation tests
df_normative_mind_distance_withOB <- df_normative_mind_distance_withOB %>% 
  mutate(distance_bin = cut_number(distance, 3)) %>% 
  mutate(distance_bin = as.numeric(distance_bin),
         distance_bin = 
           case_when(
             distance_bin == 1 ~ "proximal",
             distance_bin == 2 ~ "intermediate",
             distance_bin == 3 ~ "distal"
           ) %>% factor(levels = c("proximal", "intermediate", "distal"))
  )

## SAVE OBJECTS
save(
  df_normative_mind_distance_withOB,
  file = paste0(objects_dir, "08July2024_normative_mind_distance_withOB.RDS")
)
