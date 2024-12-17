
#################################################################

## Function to plot network circle plots to show effect sizes ##

#################################################################

base_dir <- "~/Documents/PhD/projects/CamRat/CamRat/"
source(paste0(base_dir, "scripts/03.data_analysis/setup.R"))
analysis_objects_dir <- paste0(base_dir, "outputs/objects/") # to save analysis objects for figures

df_edge_slopes_system

# Node positions ---------------------------------------------------------------

## Define nodes for circle plot
nodes <- system_order
n <- length(nodes)

## Create a set of points of length n, so that the points are evenly spaced around a circle with radius = 1 and center = (0, 0)
inner_radius <- 0.1
circumference <- 2*pi*inner_radius
arc_length <- circumference / n # determine the arc length between each pair of points
theta <- -arc_length / inner_radius # the angle associated with each arc length (in radians)

## Initiate first point at (0, 1)
df_node_positions <- tibble(
  x = 0,
  y = 1
)

## Loop through to determine the position of the rest of the nodes
for (i in 2:n) {
  
  # set last point as (x0, y0)
  row <- df_node_positions[nrow(df_node_positions),]
  x0 <- row$x
  y0 <- row$y
  
  # calculate next point based on this position
  x1 <- x0*cos(theta) - y0*sin(theta)
  y1 <- x0*sin(theta) + y0*cos(theta)
  
  # combine with previous points
  df_node_positions <- df_node_positions %>% 
    bind_rows(
      tibble(
        x = x1,
        y = y1
      )
    )
  
}

## Assign each cortical ROI to a point in the circle
df_node_positions <- df_node_positions %>% 
  mutate(system = nodes)


# Node labels -------------------------------------------------------------

## Write function to define point features
f_calc_radius <- function(x, y) { sqrt(x^2 + y^2) }
f_calc_cos <- function(x, y) { x / sqrt(x^2 + y^2) }
f_calc_sin <- function(x, y) { y / sqrt(x^2 + y^2) }

## Define node label positions as an extension of the circle radius at the same angle as the associated node
df_node_angles <- df_node_positions %>% 

  # calculate node radius (should be inner_radius = 1)
  mutate(node_radius = f_calc_radius(x, y)) %>% 
  
  # calculate sine and cosine of nodes using x & y coordinates of each
  mutate(sin_angle = f_calc_sin(x, y),
         cos_angle = f_calc_cos(x, y)
  )
  

# Define edge positions ------------------------------------------------------------


df_edge_positions <- df_mind_cortex %>% 
  
  # select edges used in analysis
  dplyr::select(R1, R2) %>%
  distinct() %>% 
  
  # add system-level connection info
  left_join(df_system_hierarchy %>% dplyr::select(-system2), by = join_by("R1" == "region_of_interest")) %>% 
  dplyr::rename("S1" = "system") %>% 
  left_join(df_system_hierarchy %>% dplyr::select(-system2), by = join_by("R2" == "region_of_interest")) %>% 
  dplyr::rename("S2" = "system") %>% 
  dplyr::select(S1, S2) %>% 
  distinct() %>% 

  # add node positions for R1 and R2 
  left_join(df_node_positions %>% 
              dplyr::select(system, x, y) %>% 
              dplyr::rename_all( ~ c("S1", "x1", "y1"))
  ) %>% 
  left_join(df_node_positions %>% 
              dplyr::select(system, x, y) %>% 
              dplyr::rename_all( ~ c("S2", "x2", "y2"))
  )


# Plotting function -----------------------------------------------------------

f_circle_plot <- function(
    df_node_weights, # tibble containing nodes (systems), their associated weights (in a column labeled node_effect_size), and their desired color (node_color)
    df_edge_weights, # tibble containing edges and their associated weights (in a column labeled edge_effect_size)
    edge_significance_threshold = 3.3, # effect size threshold for edges to plot (to avoid over-saturation of plot)
    node.size = 1.5, # size of nodes on circle
    node_labels = nodes_to_label, # a vector of regions of interest to label
    sf_label = 1.25, # ratio of label radius to node radius
    sf_segment = 1.25, # ratio of segment radius to node radius
    color_scale = slope_effect_size_scale, # color scale for nodes & edges
    label.text.size = 3, # size of text for node labels
    label.width = 20, # str_wrap width for node labels
    plot.title = "Development slope", # plot title
    include.caption = c(TRUE, FALSE), # include the effect size caption
    plot.margins = margin(t = 0, r = 0, b = 0, l = 0), # margins for the plot (to remove white space around circle but still include node labels)
    legend.title = "Effect size", # legend title
    legend.orientation = c("vertical", "horizontal", "none"), # should the legend be vertical or horizontal?
    legend.position.x = 0.1, # legend position on the x-axis
    legend.position.y = 0.2 # legend position on the y-axis
) {
  
  ### DEFINE NODE WEIGHTS ###
  df_nodes <- df_node_weights %>% 
    
    # add points here (assign each cortical system to a point on the circle, in order)
    left_join(df_node_positions)
  
  ### DEFINE EDGE WEIGHTS ###
  df_edges <- df_edge_weights %>%
    
    # filter for edges that reach significance_threshold
    filter(abs(edge_effect_size) > edge_significance_threshold) %>% 
    
    # separate edge to join with node positions
    separate(edge, into = c("S1", "S2"), sep = " - ") %>% 
    filter(S1 != S2) %>% 
    
    # add edge position information
    left_join(df_edge_positions) %>% 
    
    # arrange to plot highest absolute edge weights last
    arrange(-abs(edge_effect_size))
  
  ### DEFINE PLOT LABELS ###
  df_node_labels <- df_node_angles %>% 
    
    # filter for select nodes that we want to label
    filter(system %in% node_labels) %>% 
    
    # wrap text to condense for plot
    # mutate(region_of_interest = str_replace(region_of_interest, "Primary somatosensory area", "S1")) %>% 
    # mutate(region_of_interest = str_replace(region_of_interest, "Cornu ammonis ", "CA")) %>% 
    mutate(system = str_wrap(system, width = label.width)) %>% 
    
    # calculate position for labels and connecting segments based on desired radius length
    mutate(label_x = sf_label * (node_radius * cos_angle),
           segment_x = sf_segment * (node_radius * cos_angle),
           
           label_y = sf_label * (node_radius * sin_angle),
           segment_y = sf_segment * (node_radius * sin_angle)
    )
  
  ### PLOT AESTHETICS ###
  abs_scale_min <- df_edges %>% 
    filter(abs(edge_effect_size) > edge_significance_threshold) %>% 
    pull(edge_effect_size) %>% 
    abs() %>% 
    min %>% 
    floor
  scale_max <- df_edges %>% 
    pull(edge_effect_size) %>% 
    abs %>% 
    max %>% 
    ceiling
  
  ### PLOT LABELS ###
  if (include.caption == TRUE) {
    plot.caption = paste0("abs(edge effect size) > ", edge_significance_threshold)
  } else (plot.caption = NULL)
  
  ### PLOT ###
  p <- ggplot() +
    
    # edges
    geom_curve(data = df_edges,
               mapping = aes(x = x1, y = y1, xend = x2, yend = y2,
                             color = edge_effect_size, linewidth = abs(edge_effect_size), alpha = abs(edge_effect_size)
               ),
               #curvature = -1, angle = 90
                curvature = 0.2, angle = 90
    ) +

    # nodes
    geom_point(data = df_nodes, aes(x = x, y = y, fill = node_effect_size, color = I(node_color)), 
               shape = 21, size = node.size) +
    
    # labels
    geom_segment(data = df_node_labels, mapping = aes(x = segment_x, xend = x, y = segment_y, yend = y)) +
    geom_node_text(data = df_node_labels, mapping = aes(x = label_x, y = label_y, label = system),
                   color = "black", size = label.text.size, 
    ) +
    
    # plot aesthetics
    scale_fill_gradientn(colors = color_scale, limits = c(-scale_max, scale_max), guide = "none") +
    scale_color_gradientn(colors = color_scale, limits = c(-scale_max, scale_max)) +
    scale_linewidth_continuous(range = c(0.1, 2), limits = c(abs_scale_min, scale_max), guide = "none") +
    scale_alpha_continuous(range = c(0.1, 1), limits = c(abs_scale_min, scale_max), guide = "none") +
    
    # plot labels
    labs(title = plot.title,
         caption = plot.caption
    ) +
    
    # plot_theme
    coord_cartesian(xlim = c(-1.65, 1.65), ylim = c(-1.65, 1.65), clip = "off") +
    theme_void() +
    theme(plot.title = element_text(face = "bold", size = 10, margin = margin(t = 0, r = 0, b = -50, l = 0)),
          axis.title = element_blank(),
          plot.margin = plot.margins
    )
  
  ### PLOT LEGEND ###
  if (legend.orientation == "vertical") {
    p <- p + 
      guides(color = guide_colorbar(title = legend.title, title.position = "left", title.hjust = 0.5)) +
      theme( 
        legend.position = c(legend.position.x, legend.position.y),
        legend.direction = "vertical",
        legend.title = element_text(angle = 90),
        legend.key.width = unit(0.125, "cm"),
        legend.key.height = unit(0.65, "cm")
      )
  } else if (legend.orientation == "horizontal") {
    p <- p + 
      guides(color = guide_colorbar(title = legend.title, title.position = "top", title.hjust = 0.5)) +
      theme(
        legend.position = c(legend.position.x, legend.position.y),
        legend.direction = "horizontal",
        legend.key.width = unit(0.65, "cm"),
        legend.key.height = unit(0.125, "cm") 
      )
  } else if (legend.orientation == "none") {
    p <- p + 
      guides(color = guide_colorbar(title = legend.title, title.position = "top", title.hjust = 0.5)) +
      theme( 
        legend.position = "none",
      )
  }
  
  ### RETURN PLOT ###
  return(p)
  
}


# Save objects & plotting function to generate in figures script ----------

# save(
#   df_node_positions, df_node_angles, df_edge_positions, df_arcs, f_circle_plot,
#   file = paste0(analysis_objects_dir, "20May2024_circle_plot.RData")
# )

