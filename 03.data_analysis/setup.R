
################################################################################

## Load and format data objects and set theme for cortex-only analyses
## DOES NOT include olfactory regions

################################################################################


# libraries ---------------------------------------------------------------

library(tidyverse)
library(tidytext)
library(janitor)
library(readxl)
library(writexl)
library(nlme)
library(lme4)
library(sf)
library(patchwork)
library(viridis)
library(ggpubr)
library(ggrepel)
library(ggExtra)
library(ggraph)
library(ggh4x)
library(ggborderline)
library(GGally)
library(ggforce)
library(ggplotify)
library(Polychrome)
library(paletteer)
library(pals)
library(igraph)
library(tidygraph)
library(cowplot)
library(corrr)
library(tidymodels)
library(formattable)
library(dunn.test)
library(RColorBrewer)
library(pals)
library(RNifti)
library(concaveman)


# load data ---------------------------------------------------------------

## Set directories
base_dir <- "~/Documents/PhD/projects/CamRat/CamRat/"

## MIND data
load(paste0(base_dir, "objects/03June2024_df_mind_cortex.RDS")) # df_mind_cortex (generated in read_MIND_files.R)
load(paste0(base_dir, "objects/03June2024_combined_data_cortex.RDS")) # df_data_cortex (generated in combine_vol_MTR_strength.R)
load(paste0(base_dir, "objects/11July2024_df_mind_cortex_withOB.RDS")) # df_mind_cortex_withOB (generated in read_MIND_files_withOB.R) -- use for mouse gene expression comparison only!!

## Atlas
load(paste0(base_dir, "objects/16Dec2024_WHS_atlas_for_plotting.RDS")) # df_whs_atlas (generated in WHS_atlas_for_plotting.R)
load(paste0(base_dir, "objects/WHS_centers.RDS")) # df_roi_centers (generated in WHS_atlas_for_plotting.R)
load(paste0(base_dir, "objects/flatmap_df.RData")) # df_flatmap

#df_hierarchy <- read_excel(paste0(base_dir, "data/WHS_analysis/WHS_hierarchy.xlsx"))

## Functions
source(paste0(base_dir, "functions/linear_models.R"))


# format data -------------------------------------------------------------

## Identify cortical ROIs in MIND analysis
cortical_ROIs <- df_data_cortex %>% pull(region_of_interest) %>% unique
length(cortical_ROIs) # n = 53

## Define system-level groupings for cortical ROIs
# df_system_hierarchy <- df_hierarchy %>% 
#   dplyr::select(level5, level8) %>% 
#   dplyr::rename_all( ~ c("system", "region_of_interest")) %>% 
#   filter(region_of_interest %in% cortical_ROIs)
df_system_hierarchy <- read_xlsx(path = paste0(base_dir, "data/WHS_analysis/WHS_hierarchy_RLS.xlsx")) %>% 
  dplyr::select(system2, system, region_of_interest)


# get WHS atlas abbreviations ---------------------------------------------

df_abbreviations <- read_csv(paste0(base_dir, "data/WHS_analysis/WHS_atlas_metadata.csv")) %>% 
  dplyr::select(term, label_data) %>% 
  filter(term %in% c("abbreviation", "name")) %>% 
  mutate(region = map(1:(608/2), ~rep(.x, 2)) %>% unlist, 
         .before = 1
  ) %>% 
  pivot_wider(id_cols = region, names_from = term, values_from = label_data) %>% 
  filter(name %in% cortical_ROIs) %>% 
  dplyr::select(name, abbreviation) %>% 
  dplyr::rename("region_of_interest" = "name")


# plot theme --------------------------------------------------------------

theme_set(
  theme_cowplot() +
    theme(plot.title = element_text(size = 11),
          strip.text = element_text(size = 10, color = "black"),
          axis.title = element_text(size = 10),
          axis.text = element_text(size = 9),
          legend.title = element_text(size = 9),
          legend.text = element_text(size = 8),
          legend.margin = margin(t = 0, r = 0, b = 0, l = 0),  # Remove legend margin
          plot.margin = margin(t = 0, r = 5, b = 0, l = 0)
    )
)

## Figure height & width (for saving as png/pdf)
fig_width = 2.5
fig_height = 2

# for small plots
small_fig_width = 1.625
small_fig_height = 1.3

## Set slope colors
slope_effect_size_scale <- brewer.brbg(100)

## RMS effect size colors
RMS_effect_size_scale <- rev(brewer.rdbu(100))

## Group colors
#group_cols <- c("MS" = "#4781B4", "control" = "#FFAF01") #DEPR
group_cols <- c("control" = "#15508E", "MS" = "#9C1027") # NEW


# Calculate normative MIND & set ROI order --------------


## Define the normative network without the OB
df_normative_mind <- df_mind_cortex %>%
  filter(timepoint == 63 & str_detect(subject, "JWD")) %>%
  group_by(R1, R2) %>%
  summarise(median_weight = median(weight, na.rm = TRUE),
            sd_weight = sd(weight, na.rm = TRUE),
            z_score = median_weight/sd_weight
  ) %>% 
  arrange(-median_weight) %>% 
  
  # add system-level info
  left_join(df_system_hierarchy %>% dplyr::select(-system2) %>% dplyr::rename_all( ~ c("S1", "R1"))) %>% 
  left_join(df_system_hierarchy %>% dplyr::select(-system2) %>% dplyr::rename_all( ~ c("S2", "R2"))) %>% 
  dplyr::select(starts_with("S", ignore.case = FALSE), everything()) %>% 
  ungroup


## Calculate degree from the network (with no OB ROIs)
df_degree_cortex <- df_mind_cortex %>% 
  group_by(subject, timepoint) %>% 
  dplyr::select(subject, timepoint, R1, R2, weight) %>% 
  nest() %>% 
  mutate(
    degree = map(
      .x = data,
      .f = ~ .x %>% 
        pivot_wider(id_cols = R1, names_from = R2, values_from = weight) %>% 
        dplyr::select(-R1) %>% 
        colSums %>% 
        enframe %>% 
        dplyr::rename_all( ~ c("region_of_interest", "degree")) %>% 
        arrange(region_of_interest)
    )
  ) %>% 
  unnest(cols = c(degree)) %>% 
  dplyr::select(-data) %>% 
  
  # add system-level info
  left_join(df_system_hierarchy) %>% 
  dplyr::select(subject, timepoint, system2, system, everything()) %>% 
  ungroup

## Take median across normative individuals to define normative degree
df_normative_degree <- df_degree_cortex %>% 
  filter(str_detect(subject, "JWD") & timepoint == 63) %>% 
  group_by(system2, system, region_of_interest) %>% 
  summarise(median_degree = median(degree),
            sd_degree = sd(degree),
            z_score = median_degree/sd_degree
  ) %>% 
  arrange(-median_degree) %>% 
  ungroup

hubs <- df_normative_degree %>% top_n(10, median_degree) %>% pull(region_of_interest)



# Define normative network for MIND calculation with OB -------------------

df_normative_mind_withOB <- df_mind_cortex_withOB %>%
  filter(timepoint == 63 & str_detect(subject, "JWD")) %>%
  group_by(R1, R2) %>%
  summarise(median_weight = median(weight, na.rm = TRUE),
            sd_weight = sd(weight, na.rm = TRUE),
            z_score = median_weight/sd_weight
  ) %>% 
  arrange(-median_weight) %>% 
  
  # add system-level info
  left_join(df_system_hierarchy %>% dplyr::select(-system2) %>% dplyr::rename_all( ~ c("S1", "R1"))) %>% 
  left_join(df_system_hierarchy %>% dplyr::select(-system2) %>% dplyr::rename_all( ~ c("S2", "R2"))) %>% 
  dplyr::select(starts_with("S", ignore.case = FALSE), everything()) %>% 
  ungroup


# For plotting heatmaps: define ROI order & system-level line placement -----------------------------------------------------------


## Define system order based on broader system grouping & median degree (lowest to highest)
system_order <- df_normative_degree %>%
  group_by(system2) %>%
  mutate(system2_degree = median(median_degree)) %>%
  group_by(system) %>%
  mutate(system_degree = median(median_degree)) %>%
  arrange(system2_degree, system_degree) %>%
  pull(system) %>% 
  unique

## Define ROI order based on system strength and within-system strength
roi_order <- df_normative_degree %>%
  mutate(system = factor(system, levels = system_order)) %>%
  arrange(system, median_degree) %>%
  pull(region_of_interest) %>% 
  unique

## System colors for plotting
system_colors <- paletteer::paletteer_d("ggthemes::Tableau_20")[1:length(system_order)] %>% c
names(system_colors) <- system_order

## Create lines & colors for plotting heatmap with system
system_lines <- df_normative_degree %>% 
  mutate(region_of_interest = factor(region_of_interest, levels = roi_order)) %>% 
  arrange(region_of_interest) %>% 
  mutate(system = factor(system, levels = unique(.$system))) %>% 
  count(system) %>% 
  mutate(cumsum = cumsum(n)) %>% 
  head(nrow(.) - 1) %>% 
  pull(cumsum) + 0.5

## Specify where system-level annotations should be positioned on the heatmap
system_annotations <- df_normative_degree %>% 
  mutate(region_of_interest = factor(region_of_interest, levels = roi_order)) %>% 
  arrange(region_of_interest) %>% 
  mutate(system = factor(system, levels = unique(.$system))) %>% 
  count(system) %>% 
  mutate(cumsum = cumsum(n)) %>% 
  mutate(position = lag(cumsum) + n/2,
         position = ifelse(is.na(position), n/2, position)
  ) %>% 
  mutate(system = str_remove(system, " .*")) %>% 
  dplyr::select(system, position) %>% 
  deframe

