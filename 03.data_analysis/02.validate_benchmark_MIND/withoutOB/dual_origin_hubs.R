
base_dir <- "~/Documents/PhD/projects/CamRat/CamRat/"
source(paste0(base_dir, "scripts/final/data_analysis/setup_noOB.R"))
analysis_objects_dir <- paste0(base_dir, "scripts/final/figures/objects/systemsLevel/") # to save analysis objects for figures
figures_dir <- paste0(base_dir, "scripts/final/figures/dualOriginExplore/")

### Ed's comment on draft:

# Is it coincidental that the "hub" regions - most highly similar with all other areas of cortex - are theoretically the "dual origins" of cortex? 
# It would be interesting to see if piriform and hippocampal hubs have different or similar patterns of similarity across the cortex. 
# In principal, they could both be hubs for one of two reasons, I think: 
#   (i)  because they are both somewhat similar to ALL other cortical areas and 
#   (ii) because they are each highly similar to an anatomically distinct subset of cortical areas. 
# If the reality looked more like (ii) than (i) that would be consistent with dual origin theory, and another angle on biological validation...


# Plot connectome by with 15 systems ----------------------------------------


df_normative_mind %>%  
  
  # set row & column order
  mutate(R1 = factor(R1, levels = roi_order),
         R2 = factor(R2, levels = roi_order)
  ) %>% 
  mutate(median_weight = ifelse(R1 == R2, NA, median_weight)) %>% #pull(R1) %>% unique
  
  # plot
  ggplot(aes(x = R1, y = R2, fill = median_weight)) +
  geom_tile() +
  
  # add lines
  #geom_hline(yintercept = system_lines, color = "white", linewidth = 1) +
  #geom_vline(xintercept = system_lines, color = "white", linewidth = 1) +
  
  # add module colors as side bars
  
  # right
  annotate(ymin = c(-Inf, system_lines),
           ymax = c(system_lines, length(roi_order) + 0.5),
           xmin = length(roi_order) + 1.25, xmax = length(roi_order) + 2.5,
           geom = "rect",
           fill = system_colors) +
  # top
  annotate(xmin = c(0.5, system_lines),
           xmax = c(system_lines, length(roi_order) + 0.5),
           ymin = length(roi_order) + 1.25, ymax = length(roi_order) + 2.75,
           geom = "rect",
           fill = system_colors) +
  
  # # label systems
  # annotate(geom = "text", x = length(roi_order) + 3, y = system_annotations,
  #          label = names(system_annotations), hjust = 0, vjust = 0.25, size = 3
  # ) +
  
  # plot aesthetics
  #coord_cartesian(xlim = c(0, length(roi_order) + 15), clip = "off") +
  labs(x = "", y = "", title = "A | Adult normative cortical connectome") +
  scale_fill_viridis(option = "viridis", na.value = "white") +
  guides(fill = guide_colourbar(title.position = "left", title.hjust = 0.5, title = "edge weight")) +
  theme(axis.text.x = element_blank(),
        axis.text.y = element_text(size = 9),
        axis.ticks.x = element_blank(),
        axis.line = element_blank(), 
        legend.position = "bottom",
        legend.justification = "center",
        legend.title = element_text(size = 8),
        legend.text = element_text(size = 8),
        legend.key.height = unit(0.2, "cm"),
        legend.key.width = unit(1.5, "cm"),
        legend.margin = margin(t = -20)
  )
map(
  .x = c(".png", ".pdf"),
  .f = ~ ggsave(paste0(figures_dir, "MIND_connectome_15sys", .x), width = 9, height = 6)
)


# Plot new systems legend ------------------------------------------------------

#legend_color_names <- str_remove(names(system_colors), " .*")
#legend_colors <- system_colors
#names(legend_colors) <- legend_color_names
system_legend <- (
  system_colors %>% 
    enframe("system", "color") %>% 
    mutate(system = factor(system, levels = rev(system_order))) %>% 
    #mutate(system = str_remove(system, " .*") %>% factor(levels = rev(system_order) %>% str_remove(" .*"))) %>% 
    ggplot(aes(x = system, y = color)) +
    geom_point(aes(color = system)) +
    guides(colour = guide_legend(title = NULL, ncol = 1, override.aes = list(size = 4))) +
    scale_color_manual(values = system_colors) +
    theme(legend.margin = margin(t = 0, r = 0, b = 0, l = 0))
) %>%
  get_legend()
as.ggplot(system_legend)
map(
  .x = c(".png", ".pdf"),
  .f = ~ ggsave(paste0(figures_dir, "system_legend", .x), width = 1.75, height = 3.25)
)


# Plot connectome by median system similarity ----------------------------------------

df_normative_mind %>%  
  group_by(S1, S2) %>% 
  mutate(median_weight = median(median_weight, na.rm = TRUE)) %>% 
  mutate(median_weight = ifelse(
    (S1 == "Frontal association cortex" & S2 == "Frontal association cortex") |
      (S1 == "Temporal association cortex" & S2 == "Temporal association cortex"),
    NA_real_,
    median_weight
  )
  ) %>% 
  
  # set row & column order
  mutate(R1 = factor(R1, levels = roi_order),
         R2 = factor(R2, levels = roi_order)
  ) %>% 

  # plot
  ggplot(aes(x = R1, y = R2, fill = median_weight)) +
  geom_tile() +
  
  # add lines
  geom_hline(yintercept = system_lines, color = "white", linewidth = 1) +
  geom_vline(xintercept = system_lines, color = "white", linewidth = 1) +
  
  # add module colors as side bars
  
  # right
  annotate(ymin = c(-Inf, system_lines),
           ymax = c(system_lines, length(roi_order) + 0.5),
           xmin = length(roi_order) + 1.25, xmax = length(roi_order) + 2.5,
           geom = "rect",
           fill = system_colors) +
  # top
  annotate(xmin = c(0.5, system_lines),
           xmax = c(system_lines, length(roi_order) + 0.5),
           ymin = length(roi_order) + 1.25, ymax = length(roi_order) + 2.75,
           geom = "rect",
           fill = system_colors) +
  
  # # label systems
  # annotate(geom = "text", x = length(roi_order) + 3, y = system_annotations,
  #          label = names(system_annotations), hjust = 0, vjust = 0.25, size = 3
  # ) +
  
  # plot aesthetics
  #coord_cartesian(xlim = c(0, length(roi_order) + 15), clip = "off") +
  labs(x = "", y = "", title = "Median system-system connectivity") +
  scale_fill_viridis(option = "viridis", na.value = "white") +
  guides(fill = guide_colourbar(title.position = "top", title.hjust = 0.5, title = "median edge weight")) +
  theme(axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.line = element_blank(), 
        legend.position = "bottom",
        legend.justification = "center",
        legend.title = element_text(size = 8),
        legend.text = element_text(size = 8),
        legend.key.height = unit(0.2, "cm"),
        legend.key.width = unit(1.5, "cm"),
        legend.margin = margin(t = -20)
  )
map(
  .x = c(".png", ".pdf"),
  .f = ~ ggsave(paste0(figures_dir, "MIND_connectome_15sys_median", .x), width = 6, height = 6)
)



# Do hippocampal and piriform regions have high similarity with specific other subsystems? --------------------------------

## Identify Hippocampal and piriform regions
hip_pir_rois <- df_system_hierarchy %>% 
  filter(system %in% c("Hippocampal formation", "Piriform cortex")) %>% 
  pull(region_of_interest)


## Add system edge info
df_normative_mind_sys <- df_normative_mind %>% 
  
  # remove identical edges in reverse order
  mutate(edge = paste0(pmin(R1, R2), " - ", pmax(R1, R2))) %>% 
  group_by(edge) %>% 
  slice(1) %>% 
  
  # remove self connections
  filter(R1 != R2) %>% 
  
  # add system edge
  dplyr::select(-S1, -S2) %>% 
  left_join(df_system_hierarchy_rls, by = join_by("R1" == "region_of_interest")) %>% 
  dplyr::rename("S1" = "system") %>% 
  left_join(df_system_hierarchy_rls, by = join_by("R2" == "region_of_interest")) %>% 
  dplyr::rename("S2" = "system") %>% 
  mutate(system_edge = ifelse(
    paste0(pmin(S1, S2), " - ", pmax(S1, S2)) != "Hippocampal formation - Piriform cortex",
    paste0(pmin(S1, S2), " - ", pmax(S1, S2)),
    paste0(S1, " - ", S2)
  ) 
  )

## Filter for Hippocampal & Piriform edges
df_normative_mind_sys_HipPir <- df_normative_mind_sys %>% 
  filter(str_detect(system_edge, "Hippocampal formation|Piriform cortex")) %>% 
  mutate(system = ifelse(str_detect(system_edge, "Piriform cortex"), "Piriform cortex", "Hippocampal formation")) %>% 
  mutate(system = ifelse(system_edge == "Piriform cortex - Hippocampal formation", "Hippocampal formation", system )) %>% # special case for both systems included
  mutate(system2 = str_remove(system_edge, system) %>% str_remove(" - "))


## Plot  
df_normative_mind_sys_HipPir %>% 
  ggplot(aes(x = median_weight, y = reorder_within(system2, by = median_weight, within = system))) +
  geom_jitter(aes(color = system2), size = 0.5) +
  geom_boxplot(aes(color = system2, fill = system2), alpha = 0.2, outlier.shape = NA) +
  facet_wrap(vars(system), scales = "free_y") +
  scale_y_reordered() +
  scale_color_manual(values = system_colors, guide = "none") +
  scale_fill_manual(values = system_colors, guide = "none") +
  labs(x = "Edge weight", y = NULL,
       title = "Connectivity of hippocampal and piriform regions \nwith other cortical subsystems")

map(
  .x = c(".png", ".pdf"),
  .f = ~ ggsave(paste0(figures_dir, "Hip_pir_similarity_dist", .x), width = 6.5, height = 4)
)

## Stats
df_normative_mind_sys_HipPir_stats <- df_normative_mind_sys_HipPir %>% 
  group_by(system2) %>% 
  nest() %>% 
  
  # run model on each comparison
  mutate(
    model = map(.x = data, .f = ~ lm(median_weight ~ system, data = .x)),
    model_res = map(
      .x = model,
      .f = ~ .x %>% 
        summary %>% 
        coefficients %>% 
        as.data.frame() %>% 
        rownames_to_column("term") %>% 
        as_tibble() %>% 
        clean_names
    )
  ) %>% 
  unnest(cols = c(model_res)) %>% 
  filter(term != "(Intercept)") %>% 
  dplyr::rename("p_value" = "pr_t") %>% 
  mutate(p_adj = p.adjust(p_value, method = "fdr"))


## Compare directly
df_normative_mind_sys_HipPir_stats %>% 
  mutate(comparison = ifelse(t_value < 0, "Hip > Pir", "Pir > Hip")) %>% 
  unnest(cols = c(data)) %>% 
  
  # plot
  ggplot(aes(x = median_weight, y = reorder_within(system2, abs(t_value), comparison))) +
  geom_point(aes(color = system), position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.75), size = 0.5) +
  geom_boxplot(aes(color = system, fill = system), alpha = 0.2, outlier.shape = NA) +
  facet_wrap(vars(comparison), scales = "free_y") +
  scale_y_reordered() +
  scale_color_manual(values = system_colors, guide = "none") +
  scale_fill_manual(values = system_colors, guide = "none") +
  labs(x = "Edge weight", y = NULL,
       title = "Connectivity of hippocampal and piriform regions \nwith other cortical subsystems")
map(
  .x = c(".png", ".pdf"),
  .f = ~ ggsave(paste0(figures_dir, "Hip_pir_similarity_comparison", .x), width = 6.5, height = 3)
)



# Brain map of t-values ---------------------------------------------------

## Map t-values to ROIs for brain map
df_HipPir_diff <- df_normative_mind_sys_HipPir_stats %>% 
  dplyr::select(system2, t_value) %>% 
  dplyr::rename("system" = "system2") %>% 
  left_join(df_system_hierarchy_rls)

## Identify non-cortical ROIs
na_rois <- df_whs_atlas %>% 
  anti_join(df_HipPir_diff) %>% 
  pull(region_of_interest) %>%
  unique

## Set layer order for brain plot based on normative strength
HipPir_roi_order <- df_HipPir_diff %>% 
  slice_max(order_by = abs(t_value), n = 1) %>% 
  arrange(abs(t_value)) %>% 
  pull(region_of_interest)

## Combine atlas information with ROI weights based on nodal strength slopes
df_HipPir_diff_plot <- df_whs_atlas %>% 
  
  # select plane for plotting
  #filter(plane == "sagittal" & hemisphere == "right") %>% 
  #dplyr::select(-plane) %>% 
  unnest(cols = c(value)) %>% 
  
  # add strength data
  left_join(df_HipPir_diff) %>% 

  # arrange ROIs for plot
  mutate(region_of_interest = factor(region_of_interest, levels = c(na_rois, HipPir_roi_order)))

## Plot
df_HipPir_diff_plot %>% filter(plane == "coronal") %>% 
  ggplot() +
  geom_polygon(aes(x = x, y = y, group = region_of_interest, fill = t_value, color = t_value), 
               alpha = 1, lwd = 0.25) +
  scale_fill_gradient2(low = "blue", mid = "white", high = system_colors[["Piriform cortex"]],
                      na.value = "lightgray") +
  scale_color_gradient2(low = "blue", mid = "white", high = system_colors[["Piriform cortex"]],
                      na.value = "lightgray", guide = "none") +
  guides(fill = guide_colorbar(title = NULL, title.position = "top", title.hjust = 0.5)) +
  #labs(title = "E | Anatomy of developmental degree slopes") +
  theme_void() +
  theme(#plot.title = element_text(face = "bold", size = 12),
    strip.text = element_text(size = 11, face = "bold"),
    strip.clip = "off",
    legend.title = element_text(size = 9),
    legend.text = element_text(size = 8),
    legend.position = "bottom",
    legend.justification = "center",
    legend.key.width = unit(1.5, "cm"),
    legend.key.height = unit(0.1, "cm"),
    plot.margin = margin(l = 5)
  )

map(
  .x = c(".png", ".pdf"),
  .f = ~ ggsave(paste0(figures_dir, "Hip_pir_similarity_brainmaps_cor", .x), width = 2.5, height = 2)
)



# PCA ---------------------------------------------------------------------


## Run PCA
pca_mind <- df_normative_mind %>% 
  pivot_wider(id_cols = R1, names_from = R2, values_from = median_weight) %>% 
  column_to_rownames("R1") %>% 
  prcomp(scale = TRUE, center = TRUE)

## Extract results
df_loadings <- pca_mind$rotation %>% 
  as.data.frame() %>% 
  rownames_to_column("region_of_interest") %>% 
  as_tibble()
df_scores <- pca_mind$x %>% 
  as.data.frame() %>% 
  rownames_to_column("region_of_interest") %>% 
  as_tibble()
summary(pca_mind)$importance %>% 
  as.data.frame %>% 
  rownames_to_column("feature") %>% 
  as_tibble()

## Scatterplot
df_loadings %>%
  left_join(df_system_hierarchy) %>%
  ggplot(aes(x = PC1, y = PC2)) +
  geom_hline(yintercept = 0, color = "gray") +
  geom_vline(xintercept = 0, color = "gray") +
  geom_point(aes(fill = system), shape = 21, size = 3) +
  scale_fill_manual(values = system_colors) +
  labs(title = "Normative MIND PCA loadings")
map(
  .x = c(".png", ".pdf"),
  .f = ~ ggsave(paste0(figures_dir, "PCA_loadings_scatter", .x), width = 5, height = 4)
)

df_scores %>%
  left_join(df_system_hierarchy) %>%
  ggplot(aes(x = PC1, y = PC2)) +
  geom_hline(yintercept = 0, color = "gray") +
  geom_vline(xintercept = 0, color = "gray") +
  geom_point(aes(fill = system), shape = 21, size = 3) +
  scale_fill_manual(values = system_colors) +
  labs(title = "Normative MIND PCA scores")
map(
  .x = c(".png", ".pdf"),
  .f = ~ ggsave(paste0(figures_dir, "PCA_scores_scatter", .x), width = 5, height = 4)
)


# PCA brain plots ---------------------------------------------------------


## Identify non-cortical ROIs
na_rois <- df_whs_atlas %>% 
  anti_join(df_scores) %>% 
  pull(region_of_interest) %>%
  unique

## Set layer order for brain plot based on normative strength
PC1_roi_order <- df_loadings %>% 
  arrange(abs(PC1)) %>% 
  pull(region_of_interest)

## Combine atlas information with PC information
df_PC1_plot <- df_whs_atlas %>% 
  
  # select plane for plotting
  filter(plane == "sagittal" & hemisphere == "right") %>% 
  dplyr::select(-plane) %>% 
  unnest(cols = c(value)) %>% 
  
  # add PCA data
  left_join(df_loadings) %>% 
  
  # arrange ROIs for plot
  mutate(region_of_interest = factor(region_of_interest, levels = c(na_rois, PC1_roi_order)))

## Plot
df_PC1_plot %>% #filter(plane == "coronal") %>% 
  ggplot() +
  geom_polygon(aes(x = x, y = y, group = region_of_interest, fill = PC1, color = PC1), 
               alpha = 1, lwd = 0.25) +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red",
                       na.value = "lightgray") +
  scale_color_gradient2(low = "blue", mid = "white", high = "red",
                        na.value = "lightgray", guide = "none") +
  guides(fill = guide_colorbar(title = NULL, title.position = "top", title.hjust = 0.5)) +
  labs(title = "PC1 loadings") +
  theme_void() +
  theme(#plot.title = element_text(face = "bold", size = 12),
    strip.text = element_text(size = 11, face = "bold"),
    strip.clip = "off",
    legend.title = element_text(size = 9),
    legend.text = element_text(size = 8),
    legend.position = "bottom",
    legend.justification = "center",
    legend.key.width = unit(1.25, "cm"),
    legend.key.height = unit(0.1, "cm"),
    plot.margin = margin(l = 5)
  )

map(
  .x = c(".png", ".pdf"),
  .f = ~ ggsave(paste0(figures_dir, "PC1_loadings_brainmap_sag", .x), width = 3, height = 2)
)



# Correlate PCA with distance from HPC or PIR -----------------------------

df_scores

## Function to calculate distance 
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
      vect0 = c(df_centroids %>% filter(region_of_interest == .x & hemisphere == "right") %>% pull(x),
                df_centroids %>% filter(region_of_interest == .x & hemisphere == "right") %>% pull(y),
                df_centroids %>% filter(region_of_interest == .x & hemisphere == "right") %>% pull(z)
      ),
      vect1 = c(df_centroids %>% filter(region_of_interest == .y & hemisphere == "right") %>% pull(x),
                df_centroids %>% filter(region_of_interest == .y & hemisphere == "right") %>% pull(y),
                df_centroids %>% filter(region_of_interest == .y & hemisphere == "right") %>% pull(z)
      )
    )
    
  )
  ) %>% 
  unnest(cols = c(distance))


df_scores

## Identify HPC and PIR edges
df_edge_distance_HipPir <- df_edge_distance %>% 
  
  # add system info
  left_join(df_system_hierarchy, by = join_by("R1" == "region_of_interest")) %>% 
  dplyr::rename("S1" = "system") %>% 
  left_join(df_system_hierarchy, by = join_by("R2" == "region_of_interest")) %>% 
  dplyr::rename("S2" = "system") %>% 
  mutate(system_edge = ifelse(
    paste0(pmin(S1, S2), " - ", pmax(S1, S2)) != "Hippocampal formation - Piriform cortex",
    paste0(pmin(S1, S2), " - ", pmax(S1, S2)),
    paste0(S1, " - ", S2)
  ) 
  ) %>% 

  # filter for HIP and PIR
  filter(str_detect(system_edge, "Hippocampal formation|Piriform cortex")) %>% 
  mutate(system = ifelse(str_detect(system_edge, "Piriform cortex"), "Piriform cortex", "Hippocampal formation")) %>% 
  mutate(system = ifelse(system_edge == "Piriform cortex - Hippocampal formation", "Hippocampal formation", system )) %>% # special case for both systems included
  mutate(system2 = str_remove(system_edge, system) %>% str_remove(" - "))

# for each ROI, calculate average distance to hippocampal formation and piriform
df_edge_distance_HipPir <- df_edge_distance_HipPir %>% 
  group_by(system, R2) %>% 
  summarize(avg_distance = median(distance)) %>% 
  arrange(system, -avg_distance) %>% 
  dplyr::rename("region_of_interest" = "R2")

### Plot on rat brain ###

## Identify non-cortical ROIs
na_rois <- df_whs_atlas %>% 
  anti_join(df_edge_distance_HipPir) %>% 
  pull(region_of_interest) %>%
  unique

## Set layer order for brain plot based on normative strength
HipPir_roi_order <- df_edge_distance_HipPir %>% 
  
  filter(system == "Piriform cortex") %>% 
  
  arrange(abs(avg_distance)) %>% 
  pull(region_of_interest) %>% 
  unique

## Combine atlas information with PC information
df_HipPir_plot <- df_whs_atlas %>% 
  
  # select plane for plotting
  filter(plane == "sagittal" & hemisphere == "right") %>% 
  dplyr::select(-plane) %>% 
  unnest(cols = c(value)) %>% 
  
  # add PCA data
  left_join(df_edge_distance_HipPir %>% filter(system == "Piriform cortex")) %>% 
  
  # arrange ROIs for plot
  mutate(region_of_interest = factor(region_of_interest, levels = c(na_rois, HipPir_roi_order)))

## Plot
df_HipPir_plot %>% #filter(plane == "coronal") %>% 
  ggplot() +
  geom_polygon(aes(x = x, y = y, group = region_of_interest, fill = avg_distance, color = avg_distance), 
               alpha = 1, lwd = 0.25) +
  scale_fill_viridis(option = "A", na.value = "lightgray") +
  scale_color_viridis(option = "A", na.value = "lightgray", guide = "none") +
  guides(fill = guide_colorbar(title = NULL, title.position = "top", title.hjust = 0.5)) +
  labs(title = "Distance to HIP") +
  theme_void() +
  theme(#plot.title = element_text(face = "bold", size = 12),
    strip.text = element_text(size = 11, face = "bold"),
    strip.clip = "off",
    legend.title = element_text(size = 9),
    legend.text = element_text(size = 8),
    legend.position = "bottom",
    legend.justification = "center",
    legend.key.width = unit(1.25, "cm"),
    legend.key.height = unit(0.1, "cm"),
    plot.margin = margin(l = 5)
  )

## Correlate with PC1 & PC2
df_edge_distance_HipPir %>% 
  left_join(df_scores %>% 
              dplyr::select(region_of_interest, PC1, PC2) %>% 
              pivot_longer(contains("PC"), names_to = "PC", values_to = "score")
  ) %>% 
  
  ggplot(aes(x = avg_distance, y = score)) +
  geom_point() +
  geom_smooth(method = "lm") +
  stat_cor(label.sep = "\n", label.y.npc = "bottom", vjust = 0) +
  facet_grid2(PC ~ system, scales = "free", independent = "all") +
  labs(title = "Distance to origin vs PC score")

map(
  .x = c(".png", ".pdf"),
  .f = ~ ggsave(paste0(figures_dir, "PC_score_vs_dualOrigin_distance", .x), width = 5, height = 4)
)



