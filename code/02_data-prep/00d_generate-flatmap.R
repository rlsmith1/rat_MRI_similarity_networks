
#==============================================================================#
# Generate ggplot rendering of rat brain flatmap
#==============================================================================#


# Setup -------------------------------------------------------------------

# Libraries
library(tidyverse)
library(xml2)
library(rvest)
library(sf)

# Identify the SVG file path
base_dir <- "~/Documents/PhD/projects/CamRat/CamRat/"
objects_dir <- paste0(base_dir, "outputs/objects")

flatmap_svg <- paste0(base_dir, "data/flatmap.svg")


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
wkt_collection = paste(readLines(flatmap_svg), collapse = "")

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
     file = paste0(objects_dir, "flatmap_df.RData")) # df_flatmap


