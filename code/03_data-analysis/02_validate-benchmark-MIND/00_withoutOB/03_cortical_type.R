
#==============================================================================#
# Compare normative rat connectome to cortical cell type
# Hypothesis: Regions that are within the same cortical type will have higher MIND similarity
# source: https://link.springer.com/article/10.1007/s00429-022-02548-0#Tab4 
#==============================================================================#

base_dir <- "~/Documents/PhD/projects/CamRat/CamRat/"
source(paste0(base_dir, "code/03_data-analysis/setup.R"))


# Load & format cortical type data -------------------------------------------------------------------

## Load atlas mappings (Zilles --> WHS)
data_dir <- paste0(base_dir, "data/WHS_analysis/")
df_cortical_type_all <- read_xlsx(paste0(data_dir, "GarciaCabeza2023_all_cortical_types_RLS.xlsx")) %>% 
  clean_names

## Identify regions in the WHS atlas that map to the Zilles atlas
whs_rois <- df_cortical_type_all %>% pull(whs) %>% unique


# Define function to calculate intra-class overlap ------------------------

## Define number of edges to calculate overlap for (what percent of edges?)
n_edges <- nrow(df_mind_cortical_type_all)
max_n <- ceiling(n_edges*0.1)

## At each network density threshold, define the number of "present" edges that are intra (within) vs inter (between) cortex
f_calc_class_overlap <- function(df, n_edges_subset) {
  df %>% 
    arrange(-weight) %>% 
    head(n_edges_subset) %>% 
    count(cortex_comparison) %>% 
    mutate(percent_edge = n/sum(n)) %>% 
    filter(cortex_comparison == "within") %>% 
    dplyr::rename("n_edges_in_class" = "n") %>% 
    mutate(n_edges_total = n_edges_subset,
           density = n_edges_subset/n_edges, .before = 1
    )
}


## Plot
types <- c("archicortex", "agranular", "dysgranular", "paleocortex")
df_mind_cortical_type_all %>% 
  filter(cortical_type1 %in% types & cortical_type2 %in% types) %>% 
  mutate(cortex_edge = case_when(
    cortex_edge == "allo-allo" ~ "Allocortex -\nAllocortex",
    cortex_edge == "allo-meso" ~ "Allocortex -\nMesocortex",
    cortex_edge == "meso-meso" ~ "Mesocortex -\nMesocortex"
  )) %>% 
  
  ggplot(aes(x = cortical_type_edge, y = weight)) +
  geom_point(aes(fill = type_comparison), shape = 21, position = position_jitter(width = 0.2)) +
  geom_boxplot(aes(color = type_comparison), fill = "transparent", outlier.shape = NA) +
  annotate(geom = "text", label = "Within-class \nconnection", color = "maroon", x = 1, y = 0.65, size = 3.5)  +
  scale_color_manual(values = c("within" = "maroon", "between" = "darkgrey"), guide = "none") +
  scale_fill_manual(values = c("within" = "maroon", "between" = "darkgrey"), guide = "none") +
  labs(x = "", y = "Normative MIND edge weight",
       title = "Edge weight distributions by connection type")


# Calculate intra-class proportions across network density thresholds (Normative MIND) ----------------

df_intraclass_overlap_acutal <- map_dfr(
  .x = seq(1, max_n, 1),
  .f = ~ f_calc_class_overlap(df_mind_cortical_type_all, .x) %>% 
    mutate(network = "normative MIND", .before = 1)
)


# Calculate intra-class proportions across network density thresholds (Null nets) ----------------

## Load null network data
load(paste0(analysis_objects_dir, "03June2024_null_nets.RDS")) # df_null_nets (generated in data_analysis/01.validate_MIND/withOB/generate_null_nets.R)

## Add cortical type info to null networks
df_null_nets_corticalType <- df_null_nets %>% 
  
  # add cortical cell type data
  left_join(df_cortical_type_all, by = join_by("R1" == "whs")) %>% 
  dplyr::rename_at(c("cortical_area", "cortex", "cortical_type"), ~ paste0(.x, "1")) %>% 
  left_join(df_cortical_type_all, by = join_by("R2" == "whs")) %>% 
  dplyr::rename_at(c("cortical_area", "cortex", "cortical_type"), ~ paste0(.x, "2")) %>% 
  
  # define if cortical type class is same or different
  mutate(type_comparison = ifelse(cortical_type1 == cortical_type2, "within", "between"),
         cortex_comparison = ifelse(cortex1 == cortex2, "within", "between")
  ) %>% 
  
  # nest for analysis  
  group_by(network) %>% 
  nest()

## Calculate intra-class edge overlap for all null networks
df_combos <- expand_grid(
  network = df_null_nets_corticalType %>% pull(network),
  n_edges_total = seq(1, max_n, 1)
)
doParallel::registerDoParallel()
df_intraclass_overlap_null <- map2_dfr(
  .x = df_combos %>% pull(network),
  .y = df_combos %>% pull(n_edges_total),
  .f = ~ f_calc_class_overlap(df_null_nets_corticalType %>% filter(network == .x) %>% pull(data) %>% .[[1]], .y) %>% 
    mutate(network = .x, .before = 1)
)


# Combine & save intra-class proportion data -----------------------------------------------------

df_intraclass_overlap <- df_intraclass_overlap_acutal %>% 
  bind_rows(df_intraclass_overlap_null) %>% 
  mutate(color = ifelse(str_detect(network, "null"), "gray", "#B8DE29FF"))

save(
  df_mind_cortical_type_all, df_intraclass_overlap,
  file = paste0(objects_dir, "03June2024_cortical_type.RData")
)




