
#==============================================================================#
# Format the Waxholm atlas hierarchy so it functions as a tibble
#==============================================================================#

## See [Kleven et al *Nature Methods* (2023)](https://www.nature.com/articles/s41592-023-02034-3) for atlas download info


# libraries ---------------------------------------------------------------

library(tidyverse)
library(XML)


# read data ---------------------------------------------------------------

base_dir <- "~/Documents/PhD/projects/CamRat/CamRat/"
atlas_dir <- "~/Documents/PhD/projects/CamRat/atlases/MBAT_WHS_SD_rat_atlas_v4_pack/Data/"

xml_data <- xmlParse(paste0(atlas_dir, "WHS_SD_rat_atlas_v4_labels.ilf"))


# parse XML data ----------------------------------------------------------


# Create a function to recursively extract label information
extract_labels <- function(node) {
  label_data <- xmlAttrs(node)
  children <- xmlChildren(node)
  
  if (length(children) > 0) {
    child_data <- lapply(children, extract_labels)
    return(list(label_data = label_data, children = child_data))
  } else {
    return(list(label_data = label_data, children = NULL))
  }
}

# Find the root element (typically <milf>)
root <- xmlRoot(xml_data)

# Locate the starting point for label information (modify as needed)
start_node <- root[[2]][[1]]  # Adjust this based on the actual structure

# Extract the label hierarchy
label_hierarchy <- extract_labels(start_node)

# Define a function to flatten the hierarchy
flatten_hierarchy <- function(hierarchy, parent_label = "") {
  label_data <- hierarchy$label_data
  children <- hierarchy$children
  
  if (length(children) > 0) {
    child_data <- lapply(children, function(child) {
      flatten_hierarchy(child, parent_label = label_data[["name"]])
    })
    return(do.call(rbind, c(list(data.frame(label_data, parent = parent_label, stringsAsFactors = FALSE)), child_data)))
  }else {
    return(data.frame(label_data, parent = parent_label, stringsAsFactors = FALSE))
  }
}

# Flatten the label hierarchy
flat_data <- flatten_hierarchy(label_hierarchy) %>% 
  rownames_to_column("term") %>% 
  as_tibble() %>% 
  
  # remove "Brain" label (we know it's a brain)
  filter(str_detect(term, "label"))


# get abbreviation, color, and id for each ROI name --------------------------


df_atlas_info <- flat_data %>% 
  mutate(level = paste0("level", str_count(term, "label"))) %>% 
  mutate(term = str_remove_all(term, "label.")) %>% 
  mutate(number = str_extract(term, "\\d+"),
         number = ifelse(is.na(number), 0, number)
  ) %>% 
  mutate(term = str_remove_all(term, "[0-9]")) %>% 
  filter(!str_detect(term, "description"))


# combine information in flattened hierarchy to generate a wide csv -----------


df_tmp <- df_atlas_info %>% 
  filter(term == "name", level == "level1") %>% 
  pivot_wider(id_cols = c(number, parent), names_from = level, values_from = label_data) %>% 
  dplyr::select(-parent) %>% 
  
  left_join(
    df_atlas_info %>% 
      filter(term == "name", level == "level2") %>% 
      pivot_wider(id_cols = c(number, parent), names_from = level, values_from = label_data) %>% 
      dplyr::rename("level1" = "parent"),
    by = join_by(level1)
  ) %>% 
  
  left_join(
    df_atlas_info %>% 
      filter(term == "name", level == "level3") %>% 
      pivot_wider(id_cols = c(number, parent), names_from = level, values_from = label_data) %>% 
      dplyr::rename("level2" = "parent"),
    by = join_by(level2)
  ) %>% 
  
  left_join(
    df_atlas_info %>% 
      filter(term == "name", level == "level4") %>% 
      pivot_wider(id_cols = c(number, parent), names_from = level, values_from = label_data) %>% 
      dplyr::rename("level3" = "parent"),
    by = join_by(level3)
  ) %>% 
  
  left_join(
    df_atlas_info %>% 
      filter(term == "name", level == "level5") %>% 
      pivot_wider(id_cols = c(number, parent), names_from = level, values_from = label_data) %>% 
      dplyr::rename("level4" = "parent"),
    by = join_by(level4)
  ) %>% 
  
  left_join(
    df_atlas_info %>% 
      filter(term == "name", level == "level6") %>% 
      pivot_wider(id_cols = c(number, parent), names_from = level, values_from = label_data) %>% 
      dplyr::rename("level5" = "parent"),
    by = join_by(level5)
  ) %>% 

  left_join(
    df_atlas_info %>% 
      filter(term == "name", level == "level7") %>% 
      pivot_wider(id_cols = c(number, parent), names_from = level, values_from = label_data) %>% 
      dplyr::rename("level6" = "parent"),
    by = join_by(level6)
  ) %>% 
  
  left_join(
    df_atlas_info %>% 
      filter(term == "name", level == "level8") %>% 
      pivot_wider(id_cols = c(number, parent), names_from = level, values_from = label_data) %>% 
      dplyr::rename("level7" = "parent"),
    by = join_by(level7)
  ) %>% 
  
  dplyr::select(-contains("number"))

# fill NAs with last value to the left

df_hierarchy <- tibble()
for (i in 1:nrow(df_tmp)) {
  
  print(i)
  row <- df_tmp[i,]

  if (rowSums(is.na(row)) != 0) {
    
    na_cols <- which(is.na(row))
    first_na_idx <- na_cols[1]
    n_reps <- length(na_cols)
    fill_value <- row[first_na_idx - 1] %>% deframe
    
    row <- row %>% mutate(across(everything() , ~ifelse(is.na(.), fill_value, .)))
    
  }
  
  df_hierarchy <- df_hierarchy %>% bind_rows(row)
  
}



# export files as CSVs ----------------------------------------------------


df_atlas_info %>% write.csv(file = paste0(base_dir, "data/WHS_analysis/WHS_atlas_metadata.csv"), row.names = FALSE)
df_hierarchy %>% write.csv(file = paste0(base_dir, "data/WHS_analysis/WHS_hierarchy.csv"), row.names = FALSE)



