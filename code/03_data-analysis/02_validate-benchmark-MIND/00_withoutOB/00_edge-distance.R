
#==============================================================================#
# Relate MIND edge weight to edge distance between centroids
# Hypothesis: Regions that are closer together will have higher MIND similarity
#==============================================================================#

base_dir <- "~/Documents/PhD/projects/CamRat/CamRat/"
source(paste0(base_dir, "code/03_data-analysis/setup.R"))


## Write function to calculate Euclidean distance between two points
# d = √(x1 − x0)^2 + (y1 − y0)^2 + (z1 − z0)^2
CalculateEuclideanDistance <- function(vect0, vect1) sqrt(sum((vect1 - vect0)^2))

## Calculate distance between each pair of WHS atlas ROI centroids (calculated in data_prep/WHS_ROI_centroids.R)
df_edge_distance <- df_mind_cortex %>% 
  ungroup %>% 
  dplyr::select(R1, R2) %>% 
  distinct() %>% 
  
  # calculate Euclidean distance between centroid of each ROI
  mutate(distance = map2(
    .x = R1,
    .y = R2,
    .f = ~ CalculateEuclideanDistance(
      vect0 = c(df_roi_centers %>% filter(region_of_interest == .x & hemisphere == "right") %>% pull(x),
                df_roi_centers %>% filter(region_of_interest == .x & hemisphere == "right") %>% pull(y),
                df_roi_centers %>% filter(region_of_interest == .x & hemisphere == "right") %>% pull(z)
      ),
      vect1 = c(df_roi_centers %>% filter(region_of_interest == .y & hemisphere == "right") %>% pull(x),
                df_roi_centers %>% filter(region_of_interest == .y & hemisphere == "right") %>% pull(y),
                df_roi_centers %>% filter(region_of_interest == .y & hemisphere == "right") %>% pull(z)
      )
    )
    
  )
  ) %>% 
  unnest(cols = c(distance))

## Combine distance data with normative MIND weights
df_normative_mind_distance <- df_normative_mind %>% 
  
  # add distance info
  left_join(df_edge_distance) %>% 
  
  # remove duplicate edges in reverse order
  mutate(R1 = as.character(R1), R2 = as.character(R2)) %>% 
  group_by(edge = paste0(pmin(R1, R2), sep = " - ", pmax(R1, R2))) %>% 
  slice(1) %>% 
  ungroup %>% 
  
  # remove self connections
  filter(R1 != R2)
  

## Divide edge distances into 3 evenly sized bins for permutation tests
df_normative_mind_distance <- df_normative_mind_distance %>% 
  mutate(distance_bin = cut_number(distance, 3)) %>% 
  mutate(distance_bin = as.numeric(distance_bin),
         distance_bin = 
           case_when(
             distance_bin == 1 ~ "proximal",
             distance_bin == 2 ~ "intermediate",
             distance_bin == 3 ~ "distal"
           ) %>% factor(levels = c("proximal", "intermediate", "distal"))
  )

## Correlation P-value
edge_distance_cor <- cor.test(df_normative_mind_distance$median_weight, df_normative_mind_distance$distance)
t_val <- edge_distance_cor$statistic
df_val <- edge_distance_cor$parameter

# One-tailed p-value from t distribution
pt(t_val, df = df_val, lower.tail = TRUE)


## SAVE OBJECTS
save(
  df_normative_mind_distance,
  file = paste0(objects_dir, "17Dec2024_normative_mind_distance.RDS")
)

