
## Setup 
base_dir <- "~/Documents/PhD/projects/CamRat/CamRat/"
source(paste0(base_dir, "code/03.data_analysis/setup.R"))
analysis_objects_dir <- paste0(base_dir, "outputs/objects/")

### GLASS BRAIN NETWORK NODES & EDGES ###

# Generate brain atlas outline
df_atlas_outline <- df_whs_atlas %>% 
  filter(!str_detect(region_of_interest, "corticofugal")) %>% 
  dplyr::select(-hemisphere, -region_of_interest) %>% 
  unnest(cols = c(value)) %>% 
  group_by(plane) %>% 
  nest() %>% 
  mutate(
    hull = map(
      .x = data,
      .f = ~ concaveman(as.matrix(.x), concavity = 2.15) %>% 
        as.data.frame %>% 
        as_tibble() %>% 
        dplyr::rename_all(~ c("x", "y"))
    )
  ) %>% 
  dplyr::select(-data)
  
## Write function to plot glass brain with one or more system
f_plot_glass_brain <- function(system, 
                               plane_to_plot = "sagittal",
                               outline_width = 1.5,
                               node_size = 1,
                               edge_size = 0.75
) {
  
  # Identify system to plot  
  current_system <- system
  
  # Identify regions within that system
  current_regions <- df_system_hierarchy %>% 
    filter(system %in% current_system) %>% 
    pull(region_of_interest)
  
  # Identify system color for plotting
  if (length(current_system) > 1) {
    brain_color <- "gray"
    brain_fill <- "#d3d3d333"
  } else {
    brain_color <- system_colors[current_system]
    brain_fill <- str_replace(system_colors[current_system], "FF$", "33") # set alpha to 20% in hex code
    
  }
  
  # Define nodes
  df_nodes <- df_roi_centers %>% 
    filter(region_of_interest %in% current_regions & hemisphere == "right") %>% 
    left_join(df_system_hierarchy, by = join_by(region_of_interest)) %>% 
    left_join(enframe(system_colors, name = "system", value = "color"), by = join_by(system))
  
  if (plane_to_plot == "axial") {
    df_nodes <- df_nodes %>% dplyr::select(hemisphere, region_of_interest, color, x, y)
  } else if (plane_to_plot == "coronal") {
    df_nodes <- df_nodes %>% dplyr::select(hemisphere, region_of_interest, color, x, z) %>% dplyr::rename("y" = "z")
  } else if (plane_to_plot == "sagittal") {
    df_nodes <- df_nodes %>% dplyr::select(region_of_interest, color, y, z) %>% dplyr::rename("y" = "z", "x" = "y")
  }
  
  # Define edges (all connections between regions within system)
  df_edges <- df_normative_mind %>% 
    
    # filter for regions in current system
    filter(S1 %in% current_system & S2 %in% current_system) %>% 
    filter(S1 == S2) %>% 
    
    # add colors
    left_join(enframe(system_colors, name = "S1", value = "color"), by = join_by(S1)) %>% 
    
    # remove duplicates in reverse order
    group_by(grp = paste0(pmin(R1, R2), sep = " - ", pmax(R1, R2))) %>% 
    slice(1) %>% 
    ungroup %>% dplyr::select(-grp) %>% 
    
    # add node position information
    left_join(df_nodes, by = c("R1" = "region_of_interest")) %>% 
    dplyr::rename_at(vars(x, y), ~ paste0(.x, 1)) %>% 
    left_join(df_nodes, by = c("R2" = "region_of_interest")) %>% 
    dplyr::rename_at(vars(x, y), ~ paste0(.x, 2)) %>% 
    filter(R1 != R2)
  
  ## PLOT ##
  p <- ggplot() +
    
    # brain outline
    geom_polygon(data = df_atlas_outline %>% filter(plane == plane_to_plot) %>% unnest(cols = c(hull)),
                 aes(x = x, y = y), 
                 fill = brain_fill, color = brain_color
    ) +
    
    # network
    geom_segment(data = df_edges,
                 mapping = aes(x = x1, xend = x2, y = y1, yend = y2, color = I(color)),
                 alpha = 0.75, linewidth = edge_size
    ) +
    geom_point(data = df_nodes,
               aes(x = x, y = y, fill = I(color)),
               shape = 21, color = "black", size = node_size
    ) +
    
    # plot aesthetics
    scale_alpha_continuous(range = c(0, 0.75), limits = c(0.6, 1)) +
    scale_linewidth_continuous(range = c(0, 0.75), limits = c(0.6, 1)) +
    #labs(title = str_remove(current_system, " .*")) +
    theme_void() +
    theme(legend.position = "none",
          plot.title = element_text(hjust = 0.5, size = 10)
    )
  return(p)
  
}

# Example
map(.x = rev(system_order), 
    .f = ~ f_plot_glass_brain(.x, 
                              outline_width = 1.5,
                              node_size = 1,
                              edge_size = 0.25,
                              plane_to_plot = "sagittal")
) %>% 
  wrap_plots(ncol = 3)

## SAVE for Fig 2 ##
save(
  df_atlas_outline, f_plot_glass_brain,
  file = paste0(analysis_objects_dir, "27Sept2024_glass_brain_plot.RData")
)

