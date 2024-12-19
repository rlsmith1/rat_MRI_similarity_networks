
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Generate ggplot rendering of rat brain flatmap
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


# Setup -------------------------------------------------------------------

# Libraries
library(tidyverse)
library(xml2)
library(rvest)
library(sf)

# Identify the SVG file path
base_dir <- "~/Documents/PhD/projects/CamRat/CamRat/"
flatmap_svg <- paste0(base_dir, "depr.scripts/final/flatmap.svg")
file.exists(flatmap_svg)


# Read and parse the XML ----------------------------------------------------------------------

# Load and parse the SVG file
svg_xml <- read_xml(flatmap_svg)

# Print the first few lines to see the structure
print(xml_children(svg_xml))

# Check for namespaces in the SVG
namespaces <- xml_ns(svg_xml)
print(namespaces) # If a default namespace like "http://www.w3.org/2000/svg" is present, you’ll need to modify your query to include this namespace in the search.

# Extract all <path> elements
paths <- xml_find_all(svg_xml, ".//d1:path | .//d1:ellipse", namespaces)

# Create tibble from elements
svg_df <- data.frame(
  name = xml_name(paths),
  id = xml_attr(paths, "id"),
  d = xml_attr(paths, "d"),
  class = xml_attr(paths, "class"),
  stringsAsFactors = FALSE
) %>% as_tibble()

# Use regular expression to extract class names and fill colors
style_vector <- xml_find_all(svg_xml, ".//d1:style", namespaces)
hex_vector <- regmatches(style_vector, gregexpr("\\.st[0-9]+\\{fill:#[0-9A-Fa-f]+;", style_vector))[[1]]
df_hex <- tibble(hex_vector) %>% 
  mutate(class = sub("\\{.*", "", hex_vector) %>% str_remove("."),
         hex = sub(".*fill:(#[0-9A-Fa-f]+);", "\\1", hex_vector)
  ) %>% 
  dplyr::select(-hex_vector)

# Combine hex codes with IDs & d paths
df_svg_hex <- df_hex %>% 
  left_join(svg_df) 



# Read polygons -----------------------------------------------------------


# Read WKT renderings of SVG paths
wkt_collection = paste(readLines("/Users/smithral/Documents/PhD/projects/CamRat/CamRat/depr.scripts/final/flatmap_wkt.txt"), collapse = "")

# Convert to geometry collection
geometry_collection_sf <- st_as_sfc(wkt_collection)

# Extract the GEOMETRYCOLLECTION (assuming there's only one feature)
geometry <- st_geometry(geometry_collection_sf)[[1]]

# Loop through each geometry in the collection and store in a list
geometry_list <- list()
for (i in seq_along(geometry)) {
  individual_geometry <- geometry[[i]]  # Access the i-th geometry
  geometry_list[[i]] <- individual_geometry  # Store in the list
}

# Create a tibble with geometry ID and individual geometries
df_geometry <- tibble(
  geometry_id = seq_along(geometry_list),  # index for each geometry
  geometry = geometry_list                 # Individual geometries
) %>% st_as_sf # convert to geometry

# Convert geometries to polygons and add region of interest ID and hex code
df_svg_geometry <- df_geometry %>% 
  st_cast("POLYGON") %>% 
  bind_cols(df_svg_hex %>% 
              arrange(name) %>% # the geometry WKT has ellipses first
              dplyr::select(name, id, hex) %>% 
              filter(!is.na(id))
  )

# Plot!
df_svg_geometry %>% 
  ggplot(mapping = aes(geometry = geometry)) +
  geom_sf(aes(fill = I(hex)))



# Add full ROI names ------------------------------------------------------


df_swanson_abbr <- read_xlsx(paste0(base_dir, "data/swanson_abbreviations.xlsx"))
df_whs_to_swanson <- read_xlsx(paste0(base_dir, "data/WHS_analysis/WHS_to_Swanson_cortex.xlsx")) %>% 
  dplyr::rename("roi_abbreviation" = "Swanson_abbr")

# Add BM4 ROI labels (full names)
df_flatmap <- df_svg_geometry %>% 
  mutate(roi_abbreviation = str_remove(id, "_.*"),
         hemisphere = ifelse(str_detect(id, "Right"), "Right", "Left")
  ) %>% 
  left_join(df_swanson_abbr, by = join_by(roi_abbreviation)) %>%
  dplyr::select(region_of_interest, roi_abbreviation, hemisphere, id, hex, geometry) 

# Add WHS labels
df_flatmap <- df_flatmap %>% 
  dplyr::rename("Swanson_name" = "region_of_interest") %>% 
  left_join(df_whs_to_swanson) %>% 
  dplyr::rename("region_of_interest" = "WHS_name") %>% 
  dplyr::select(hemisphere, Swanson_name, region_of_interest, id, hex, geometry)


# Save for plotting! ------------------------------------------------------

save(df_flatmap,
     file = paste0(base_dir, "objects/flatmap_df.RData")) # df_flatmap



# Make some plots! --------------------------------------------------------

figures_dir <- paste0(base_dir, "outputs/figures/try.flatmaps/")

## Plot normative degree
df_flatmap %>% 
  dplyr::rename("Swanson_name" = "region_of_interest") %>% 
  left_join(df_whs_to_swanson) %>% 
  dplyr::rename("region_of_interest" = "WHS_name") %>% 
  expand_grid(timepoint = c(20, 35, 63, 300) %>% factor(levels = c(20, 35, 63, 300))) %>% 
  left_join(
    df_data_cortex %>% 
      filter(study == "MRC" & feature == "degree") %>% 
      group_by(region_of_interest, timepoint) %>% 
      summarise(median_degree = median(value, na.rm = TRUE)),
    by = join_by(region_of_interest, timepoint)
  ) %>% 
  
  ggplot(mapping = aes(geometry = geometry)) +
  geom_sf(aes(fill = median_degree, color = median_degree)) +
  facet_wrap(vars(timepoint), nrow = 1) +
  scale_fill_viridis(na.value = "white") +
  scale_color_viridis(na.value = "darkgray") +
  guides(fill = guide_colorbar(title = "degree", title.position = "top", title.hjust = 0.5),
         color = guide_colorbar(title = "degree", title.position = "top", title.hjust = 0.5)) +
  labs(title = "Normative degree") +
  theme_void() +
  theme(strip.text = element_blank(),
        legend.position = "right",
        legend.key.width = unit(0.1, "cm"),
        legend.key.height = unit(0.65, "cm")
  )
map(
  .x = c(".png", ".pdf"),
  .f = ~ ggsave(paste0(figures_dir, "normative_degree_development", .x), width = 8, height = 2.5)
)

## Plot MTR
df_flatmap %>% 
  dplyr::rename("Swanson_name" = "region_of_interest") %>% 
  left_join(df_whs_to_swanson) %>% 
  dplyr::rename("region_of_interest" = "WHS_name") %>% 
  expand_grid(timepoint = c(20, 35, 63, 300) %>% factor(levels = c(20, 35, 63, 300))) %>% 
  left_join(
    df_data_cortex %>% 
      filter(study == "MRC" & feature == "MTR") %>% 
      group_by(region_of_interest, timepoint) %>% 
      summarise(median_mtr = median(value, na.rm = TRUE)),
    by = join_by(region_of_interest, timepoint)
  ) %>% 
  
  ggplot(mapping = aes(geometry = geometry)) +
  geom_sf(aes(fill = median_mtr, color = median_mtr)) +
  facet_wrap(vars(timepoint), nrow = 1) +
  scale_fill_viridis(na.value = "white") +
  scale_color_viridis(na.value = "darkgray") +
  guides(fill = guide_colorbar(title = "MTR", title.position = "top", title.hjust = 0.5),
         color = guide_colorbar(title = "MTR", title.position = "top", title.hjust = 0.5)) +
  labs(title = "Normative MTR") +
  theme_void() +
  theme(strip.text = element_blank(),
        legend.position = "right",
        legend.key.width = unit(0.1, "cm"),
        legend.key.height = unit(0.65, "cm")
  )
map(
  .x = c(".png", ".pdf"),
  .f = ~ ggsave(paste0(figures_dir, "normative_MTR_development", .x), width = 8, height = 2.5)
)

## Plot volume
df_flatmap %>% 
  dplyr::rename("Swanson_name" = "region_of_interest") %>% 
  left_join(df_whs_to_swanson) %>% 
  dplyr::rename("region_of_interest" = "WHS_name") %>% 
  left_join(
    df_data_cortex %>%
      filter(study == "MRC" & timepoint == 63 & feature == "volume") %>% 
      group_by(region_of_interest) %>% 
      summarise(median_volume = median(value, na.rm = TRUE))
  ) %>% 
  
  ggplot(mapping = aes(geometry = geometry)) +
  geom_sf(aes(fill = median_volume, color = median_volume)) +
  scale_fill_viridis(na.value = "white") +
  scale_color_viridis(na.value = "darkgray") +
  guides(fill = guide_colorbar(title = "volume", title.position = "top", title.hjust = 0.5),
         color = guide_colorbar(title = "volume", title.position = "top", title.hjust = 0.5)) +
  labs(title = "Normative volume") +
  theme_void() +
  theme(legend.position = c(0.75, 0.85),
        legend.direction = "horizontal"
  )
map(
  .x = c(".png", ".pdf"),
  .f = ~ ggsave(paste0(figures_dir, "normative_volume", .x), width = 5, height = 3)
)

## Plot PC2
df_flatmap %>% 
  dplyr::rename("Swanson_name" = "region_of_interest") %>% 
  left_join(df_whs_to_swanson) %>% 
  dplyr::rename("region_of_interest" = "WHS_name") %>% 
  left_join(df_scores %>% dplyr::select(region_of_interest, PC2)) %>%
  
  ggplot(mapping = aes(geometry = geometry)) +
  geom_sf(aes(fill = PC2, color = PC2)) +
  scale_fill_gradientn(colors = rev(brewer.piyg(100)), limits = c(-2.78, 2.78), na.value = "white") +
  scale_color_gradientn(colors = rev(brewer.piyg(100)), limits = c(-2.78, 2.78), na.value = "gray") +
  guides(fill = guide_colorbar(title.position = "top", title.hjust = 0.5),
         color = guide_colorbar(title.position = "top", title.hjust = 0.5)) +
  labs(title = "PC2 scores") +
  theme_void() +
  theme(legend.position = c(0.75, 0.85),
        legend.direction = "horizontal"
  )
map(
  .x = c(".png", ".pdf"),
  .f = ~ ggsave(paste0(figures_dir, "PC2", .x), width = 5, height = 3)
)

## Plot PC1
df_flatmap %>% 
  dplyr::rename("Swanson_name" = "region_of_interest") %>% 
  left_join(df_whs_to_swanson) %>% 
  dplyr::rename("region_of_interest" = "WHS_name") %>% 
  left_join(df_scores %>% dplyr::select(region_of_interest, PC1)) %>%
  
  ggplot(mapping = aes(geometry = geometry)) +
  geom_sf(aes(fill = PC1, color = PC1)) +
  scale_fill_gradientn(colors = rev(brewer.piyg(100)), limits = c(-6.65, 6.65), na.value = "white") +
  scale_color_gradientn(colors = rev(brewer.piyg(100)), limits = c(-6.65, 6.65), na.value = "gray") +
  guides(fill = guide_colorbar(title.position = "top", title.hjust = 0.5),
         color = guide_colorbar(title.position = "top", title.hjust = 0.5)) +
  labs(title = "PC1 scores") +
  theme_void() +
  theme(legend.position = c(0.75, 0.85),
        legend.direction = "horizontal"
  )
map(
  .x = c(".png", ".pdf"),
  .f = ~ ggsave(paste0(figures_dir, "PC1", .x), width = 5, height = 3)
)

## Plot cortical systems
df_flatmap %>% 
  dplyr::rename("Swanson_name" = "region_of_interest") %>% 
  left_join(df_whs_to_swanson) %>% 
  dplyr::rename("region_of_interest" = "WHS_name") %>% 
  left_join(df_system_hierarchy) %>% 
  
  ggplot(mapping = aes(geometry = geometry)) +
  geom_sf(aes(fill = system), color = "black") +
  scale_fill_manual(values = system_colors, na.value = "white") +
  guides(fill = guide_legend(title = NULL)) +
  labs(title = "Cortical system") +
  theme_void() +
  theme(legend.position = "bottom",
        legend.key.size = unit(0.25, "cm"),
        legend.text = element_text(size = 5)
  )
map(
  .x = c(".png", ".pdf"),
  .f = ~ ggsave(paste0(figures_dir, "cortical_system", .x), width = 5, height = 3)
)

## Spatial distribution of development stuff
development_objects_dir <- paste0(base_dir, "outputs/objects/Fig3/") # to save analysis objects for figures
objects <- list.files(development_objects_dir)
objects <- objects[!str_detect(objects, "null|archive")]
for (obj in objects){
  print(paste0("loading ", obj, "...."))
  load(paste0(development_objects_dir, obj))
}

df_flatmap %>% 
  dplyr::rename("Swanson_name" = "region_of_interest") %>% 
  left_join(df_whs_to_swanson) %>% 
  dplyr::rename("region_of_interest" = "WHS_name") %>% 
  left_join(df_node_slopes_roi %>% 
              filter(feature == "degree" & term == "age") %>% 
              mutate(hemisphere = ifelse(period == "aging", "Left", "Right")) %>% 
              dplyr::select(hemisphere, region_of_interest, t_value)
  ) %>% 
  
  ggplot(mapping = aes(geometry = geometry)) +
  geom_sf(aes(fill = t_value, color = t_value)) +
  annotate(geom = "text", label = "Age-related change in degree", x = 100, y = -460, hjust = 0.19, size = 4.35) +
  scale_fill_gradientn(colors = slope_effect_size_scale, na.value = "white", limits = c(-6.4, 6.4)) +
  scale_color_gradientn(colors = slope_effect_size_scale, na.value = "gray", limits = c(-6.4, 6.4)) +
  guides(fill = guide_colorbar(title = "Developmental slope", title.position = "top", title.hjust = 0.5),
         color = guide_colorbar(title = "Developmental slope", title.position = "top", title.hjust = 0.5)) +
  labs(title = "Early development change in degree") +
  coord_sf(clip = "off") +
  theme_void() +
  theme(legend.position = c(0.75, 0.85),
        legend.direction = "horizontal"
  )
map(
  .x = c(".png", ".pdf"),
  .f = ~ ggsave(paste0(figures_dir, "development_slopes", .x), width = 5, height = 3)
)


## Spatial distribution of ELS stuff
stress_objects_dir <- paste0(base_dir, "outputs/objects/Fig4/") # to save analysis objects for figures
objects <- list.files(stress_objects_dir)
objects <- objects[!str_detect(objects, "null|archive")]
for (obj in objects){
  print(paste0("loading ", obj, "...."))
  load(paste0(stress_objects_dir, obj))
}

df_flatmap %>% 
  dplyr::rename("Swanson_name" = "region_of_interest") %>% 
  left_join(df_whs_to_swanson) %>% 
  dplyr::rename("region_of_interest" = "WHS_name") %>% 
  left_join(df_node_stress_roi %>% 
              filter(feature == "degree" & term == "groupMS") %>% 
              dplyr::select(region_of_interest, t_value)
  ) %>% 
  
  ggplot(mapping = aes(geometry = geometry)) +
  geom_sf(aes(fill = t_value, color = t_value)) +
  scale_fill_gradientn(colors = RMS_effect_size_scale, na.value = "white", limits = c(-3.15, 3.15)) +
  scale_color_gradientn(colors = RMS_effect_size_scale, na.value = "gray", limits = c(-3.15, 3.15)) +
  guides(fill = guide_colorbar(title = "Effect size", title.position = "top", title.hjust = 0.5),
         color = guide_colorbar(title = "Effect size", title.position = "top", title.hjust = 0.5)) +
  labs(title = "RMS effects in young adulthood") +
  theme_void() +
  theme(legend.position = c(0.75, 0.85),
        legend.direction = "horizontal"
  )
map(
  .x = c(".png", ".pdf"),
  .f = ~ ggsave(paste0(figures_dir, "stress_effects", .x), width = 5, height = 3)
)



