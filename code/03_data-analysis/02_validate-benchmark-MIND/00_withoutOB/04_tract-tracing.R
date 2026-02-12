
#==============================================================================#
# Compare normative rat connectome to axonal connection tract-tracing data
# Hypothesis: Regions with more similar tract-tracing connectivity profiles will have higher MIND similarity
# source: https://www.pnas.org/doi/10.1073/pnas.2017733117
#==============================================================================#

# sheet 3 (male subsystem bins) key:

# 0 Same origin & termination	identity: Intra-region connections are not included in inter-region analysis.
# 0 No data	no article: Following an extensive search, no article in the primary literature was found with data relevant to the connection and suitable for network analysis.
# 0 Unclear: It is unclear from a report if the connection exists or not. Typically, this is applied only when no clear data is available for a connection, and from the available data it is reasonable to infer that the connection is likely weak (at most), very weak, or absent.
# 0 Absent: A report provided evidence to indicate that a connection is absent (does not exist).
# 2 Axons-of-passage: A report indicated the existence of axons-of-passage, from which it is inferred that some (weak) connectivity is possible.
# 1 Very weak: The lowest possible reported value for a present connection is very weak.
# 2	Weak: The reported value for the present connection is weak.
# 3	Weak to moderate: The reported value for the present connection is weak to moderate.
# 4	Present: The connection is reported to exist (it is present), but its strength is indeterminate.
# 4	Moderate: The reported value for the present connection is moderate.
# 5	Moderate to strong: The reported value for the present connection moderate to strong.
# 6	Strong: The reported value for the present connection is strong.
# 7 Very strong: The highest possible reported value for a present connection is very strong.

base_dir <- "~/Documents/PhD/projects/CamRat/CamRat/"
source(paste0(base_dir, "code/03_data-analysis/setup.R"))


# Load Swanson-Sporns tract-tracing data -------------------------

## Take both right & left sides (can bilaterally average, the WHS atlas does not discriminate between hemispheres)
df_tmp <- read_xlsx(paste0(data_dir, "swanson_sporns_2020.xlsx"), sheet = 3) %>%
  dplyr::select(!starts_with("...")) %>% 
  clean_names() %>% 
  
  # specify ipsi- and contralateral connections (1 = right, 2 = left)
  dplyr::rename_all( ~ case_when(
    str_detect(.x, "x1") ~ paste0("right_", str_remove(.x, ".*_")),
    str_detect(.x, "x2") ~ paste0("left_", str_remove(.x, ".*_")),
    TRUE ~ .x
  )
  ) %>% 
  mutate_at(1, ~ case_when(.x == 1 ~ "right",
                           .x == 2 ~ "left",
                           TRUE ~ .x
  )
  ) %>% 
  
  # remove top 3 rows (unnecessary information for now)
  tail(nrow(.) - 3) %>%
  
  # add ROI abbreviation to each bilateral specification (rows)
  dplyr::rename_at(1, ~ "from") %>% 
  mutate(from = paste0(from, sep = "_", to_side)) %>% 
  dplyr::select(-to_side) %>% 
  filter(from != "NA_Abbr.")

## Format tract-tracing data frame
df_swanson <- df_tmp %>% 
  
  # add ROI abbreviation to each bilateral specification (columns)
  dplyr::rename_at(2:ncol(.), ~ c(df_tmp %>% pull(from))) %>% 
  pivot_longer(2:ncol(.), names_to = "to", values_to = "tract_weight") %>% 
  mutate(tract_weight = as.numeric(tract_weight),
         log10_tract_weight = log10(tract_weight + 0.00001)
  )


## Load WHS to Swanson atlas mappings (cortical only)
df_whs_to_swanson <- read_xlsx(paste0(data_dir, "WHS_analysis/WHS_to_Swanson_cortex.xlsx")) # n = 68 (Swanson + WHS combined)

## Identify regions defined in WHS and Swanson atlases for this analysis
#WHS_cortical_ROIs <- df_whs_to_swanson %>% pull(WHS_name) %>% unique # n = 57
Swanson_cortical_ROIs <- df_whs_to_swanson %>% 
  filter(WHS_name %in% cortical_ROIs) %>% 
  pull(Swanson_abbr) %>% 
  unique # n = 52 (no OB)


# Format Swanson-Sporns tract-tracing data --------------------------------

## In this analysis, we don't care about directionality of connections (input vs output) or laterality (intra- vs inter-hemispheric)
## Thus, we can average across directions and across hemispheres

df_swanson_avg <- df_swanson %>% 
  mutate_at(vars(c(from, to)), ~ str_remove(.x, ".*_")) %>% 
  dplyr::rename_at(vars(c(from, to)), ~ c("swanson_R1", "swanson_R2")) %>% 
  
  # include only regions in cortical analysis
  filter(swanson_R1 %in% Swanson_cortical_ROIs & swanson_R2 %in% Swanson_cortical_ROIs) %>% 
  
  # Add 'broad' Swanson labels and WHS atlas labels
  left_join(df_whs_to_swanson %>% 
              dplyr::rename_all( ~ paste0(.x, 1)), 
            by = join_by("swanson_R1" == "Swanson_abbr1")
  ) %>% 
  left_join(df_whs_to_swanson %>% 
              dplyr::rename_all( ~ paste0(.x, 2)), 
            by = join_by("swanson_R2" == "Swanson_abbr2")
  ) %>% 
  
  # Take median across hemispheres, directions, and 'broad' regions - manually annotated where the Swanson atlas has more subdivisions than WHS
  group_by(edge = paste0(pmin(broad_Swanson_abbr1, broad_Swanson_abbr2), sep = "-", pmax(broad_Swanson_abbr1, broad_Swanson_abbr2))) %>% 
  mutate(
    tract_weight = median(tract_weight),
    log10_tract_weight = median(log10_tract_weight)
  ) %>% 
  dplyr::select(edge, broad_Swanson_name1, broad_Swanson_name2, WHS_name1, WHS_name2, tract_weight, log10_tract_weight) %>% 
  distinct() %>% 
  ungroup


# Combine tract-tracing with MIND data ------------------------------------

## Combine Swanson's data with MIND weights
df_swanson_mind <- df_swanson_avg %>% 
  
  # add MIND weights
  left_join(#df_normative_mind %>% dplyr::select(-sd_weight, -z_score), # n = 57 cortical ROIs
    df_mind_cortex %>% filter(str_detect(subject, "JWD")) %>% group_by(timepoint,R1, R2) %>% summarise(median_weight = median(weight)),
    by = join_by("WHS_name1" == "R1", "WHS_name2" == "R2")
  ) %>% 
  dplyr::rename("mind_weight" = "median_weight") %>% 
  
  # remove all self-connections (these are 0 by default in MIND, but not tract-tracing)
  #filter(broad_Swanson_name1 != broad_Swanson_name2) %>% 
  mutate_if(is.numeric, ~ ifelse(broad_Swanson_name1 == broad_Swanson_name2, NA_real_, .x))


## Atlas groupings - group by Swanson atlas ('broad' groupings remove anatomical subdivisions that don't exist in WHS), 
## average across any sub-divisions in the other
df_swanson_mind_SwansonLabs <- df_swanson_mind %>% 
  group_by(timepoint, edge, broad_Swanson_name1, broad_Swanson_name2) %>% 
  summarise_if(is.numeric, ~ median(.x)) %>% 
  ungroup %>% 
  
  # scale weights for correlation & plotting
  mutate(mind_weight_scaled = scale(mind_weight)[,1],
         log10_tract_weight_scaled = scale(log10_tract_weight)[,1]
  ) 



# Arrange Swanson labels by system for plotting ---------------------------

df_hierarchy_Swanson <- df_system_hierarchy %>% 
  left_join(df_whs_to_swanson %>% 
              dplyr::select(WHS_name, broad_Swanson_name),
            by = join_by("region_of_interest" == "WHS_name")
  )

Swanson_roi_order <- df_hierarchy_Swanson %>% 
  mutate(system = factor(system, levels = system_order)) %>% 
  arrange(system) %>% 
  filter(!is.na(broad_Swanson_name)) %>% 
  pull(broad_Swanson_name) %>% 
  unique



# Write networks to tables for comparison ---------------------------------

## Set self-connections to NA because in MIND these are automatically 0

# write MIND with Swanson labels
df_swanson_mind_SwansonLabs %>% 
  arrange(broad_Swanson_name1, broad_Swanson_name2) %>% 
  pivot_wider(id_cols = broad_Swanson_name1, names_from = broad_Swanson_name2, values_from = mind_weight) %>% 
  write_xlsx(paste0(tables_dir, "normative_MIND_SwansonLabs.xlsx"))

# write tract-tracing with Swanson labels
df_swanson_mind_SwansonLabs %>% 
  arrange(broad_Swanson_name1, broad_Swanson_name2) %>% 
  pivot_wider(id_cols = broad_Swanson_name1, names_from = broad_Swanson_name2, values_from = tract_weight) %>% 
  write_xlsx(paste0(tables_dir, "tract_tracing_SwansonLabs.xlsx"))


### break #


# Plot heatmaps without binarizing? ---------------------------------------

## Edge correlation
df_swanson_mind_SwansonLabs %>% 
  dplyr::select(edge, log10_tract_weight_scaled, mind_weight_scaled) %>% 
  distinct() %>% 
  correlate(method = "pearson")


# Count absent/present connections ----------------------------------------

## Classify Swanson connection as either present (not 0) or absent (== 0 )
df_swanson_mind_SwansonLabs_present <- df_swanson_mind_SwansonLabs %>% 
  mutate(swanson_present = ifelse(tract_weight != 0, "present", "absent") %>% as.factor)

## Count the number of present connections in the Swanson data
present_v_absent <- df_swanson_mind_SwansonLabs_present %>% 
  dplyr::select(-contains("broad")) %>% 
  distinct() %>% 
  count(swanson_present) %>% 
  deframe


# 1. Binarized networks -----------------------------------------------------

## Bin MIND network so there are equivalent proportions of present vs absent connections, 
## as compared to the tract-tracing data (0 = absent, anything else = present)

df_swanson_mind_SwansonLabs_present_binned <- df_swanson_mind_SwansonLabs_present %>% 
  arrange(-mind_weight) %>% 
  
  # identify top n edges in MIND network as "present"
  mutate(mind_present = mind_weight %in% head(sort(mind_weight, decreasing = TRUE), present_v_absent["present"]*2),
         mind_present = ifelse(mind_present == TRUE, "present", "absent") %>% as.factor()
  )


# 2. Connection probability --------------------------------------------------

df_swanson_mind_SwansonLabs_present_connectProb <- df_swanson_mind_SwansonLabs_present %>% 
  
  # Split MIND weights into even bins split by 0.1
  mutate(mind_bin = cut(mind_weight, breaks = seq(0, 1, 0.1), include.lowest = TRUE)) %>% 
  
  # count the number of Swanson connections that are present vs absent in each MIND bin
  count(mind_bin, swanson_present) %>%
  
  # format for plotting
  pivot_wider(id_cols = mind_bin, names_from = swanson_present, values_from = n) %>% 
  mutate(relative_frequency = (present / (present + absent)) * sum(present_v_absent),
         label = (present / (present + absent)) %>% round(2)) %>% 
  pivot_longer(c(absent, present), names_to = "swanson_bin", values_to = "frequency")


# 3. Generate Jaccard similarity matrix from exact atlas mappings, and correlate Jaccard with MIND weight for each edge ---------------------

## Define function to calculate jaccard index
jaccard <- function(a, b) {
  intersection = length(intersect(a, b))
  union = length(union(a, b))
  #union = length(a) + length(b) - intersection
  return (intersection/union)
}

## Calculate the Jaccard index (similarity based on axonal connections with other ROIs) for each pair of ROIs 
df_jaccard <- df_swanson_mind_SwansonLabs %>% 
  arrange(broad_Swanson_name1, broad_Swanson_name2) %>% 
  
  # calculate Jaccard by timepoint
  group_by(timepoint) %>% 
  nest() %>% 
  
  # format data so each cell contains the ROI and its connection value
  mutate(
    data = map(
      .x = data,
      .f = ~ .x %>% 
        pivot_wider(id_cols = broad_Swanson_name1, names_from = broad_Swanson_name2, values_from = tract_weight) %>% 
        mutate_at(2:ncol(.), ~ paste0(broad_Swanson_name1, "_", .x)) %>% 
        dplyr::select(-broad_Swanson_name1)
    ),
    jaccard = map(   # calculate Jaccard index for each pair of columns (ROIs)
      .x = data,
      .f = ~ .x %>% 
        colpair_map(jaccard) %>% 
        pivot_longer(2:ncol(.), names_to = "broad_Swanson_name2", values_to = "jaccard") %>% 
        dplyr::rename("broad_Swanson_name1" = "term") %>% 
        filter(broad_Swanson_name1 != broad_Swanson_name2)
    )
  ) %>% 
  unnest(cols = c(jaccard))
  


## Combine tract-tracing Jaccard index with MIND weight
df_jaccard_mind <- df_jaccard %>% 
  
  # Add WHS labels
  left_join(df_whs_to_swanson %>% 
              dplyr::select(contains(c("WHS", "broad"))) %>% 
              distinct() %>% 
              dplyr::rename_all( ~ paste0(.x, 1))
  ) %>% 
  left_join(df_whs_to_swanson %>% 
              dplyr::select(contains(c("WHS", "broad"))) %>% 
              distinct() %>% 
              dplyr::rename_all( ~ paste0(.x, 2))
  ) %>% 
  
  # add MIND weights
  left_join(
    #df_normative_mind, 
    df_mind_cortex %>% filter(str_detect(subject, "JWD")) %>% group_by(timepoint,R1, R2) %>% summarise(median_weight = median(weight)),
    by = join_by("WHS_name1" == "R1", "WHS_name2" == "R2", timepoint)
  ) %>% 
  dplyr::rename("mind_weight" = "median_weight") %>% 
  
  # group by Swanson and combine
  group_by(timepoint, broad_Swanson_name1, broad_Swanson_name2) %>% 
  summarise(jaccard = median(jaccard),
            mind_weight = median(mind_weight)
  ) %>% 
  ungroup
  
  # remove duplicate edges in reverse order
  #group_by(grp = paste0(pmin(broad_Swanson_name1, broad_Swanson_name2), sep = "-", pmax(broad_Swanson_name1, broad_Swanson_name2))) %>% 
  #slice(1) %>% 
  #ungroup %>% dplyr::select(-grp)

## Calculate p-value of spearman correlation
tract_tracing_rho <- cor.test(df_jaccard_mind$jaccard, df_jaccard_mind$mind_weight, method = "spearman")
t_val <- tract_tracing_rho$statistic
df_val <- tract_tracing_rho$parameter
pt(t_val, df = df_val, lower.tail = TRUE)



# plot by timepoint -------------------------------------------------------

df_jaccard_mind %>% 
  ggplot(aes(x = jaccard, y = mind_weight)) +
  geom_point() +
  geom_smooth(method = "lm") +
  stat_cor() +
  facet_wrap(vars(timepoint))



# 3b. Correlate actual tract-tracing Jaccard with null nets ----------------------------------------

## Load null network data
load(paste0(objects_dir, "03June2024_null_nets.RDS")) # df_null_nets (generated in data_analysis/01.validate_MIND/generate_null_nets.R)

## Add Swanson atlas information to null networks
df_null_nets_tractTracing <- df_null_nets %>% 
  
  # add Swanson atlas labels
  left_join(df_whs_to_swanson %>% 
              dplyr::select(contains(c("WHS", "broad"))) %>% 
              distinct() %>% 
              dplyr::rename_all( ~ paste0(.x, 1)),
            by = join_by("R1" == "WHS_name1")
  ) %>% 
  left_join(df_whs_to_swanson %>% 
              dplyr::select(contains(c("WHS", "broad"))) %>% 
              distinct() %>% 
              dplyr::rename_all( ~ paste0(.x, 2)),
            by = join_by("R2" == "WHS_name2")
  )

## Calculate correlation between actual tract-tracing Jaccard index and permuted MIND weights
df_jaccard_mind_null <- df_null_nets_tractTracing %>% 
  ungroup %>% 
  dplyr::select(network, broad_Swanson_name1, broad_Swanson_name2, weight) %>% 
  left_join(df_jaccard) %>% 
  
  # calculate Spearman's rho for each null network
  group_by(network) %>% 
  nest() %>% 
  mutate(
    rho = map(
      .x = data,
      .f = ~ cor.test(.x$weight, .x$jaccard, method = "spearman")$estimate
    )
  ) %>% 
  unnest(cols = c(rho)) %>% 
  
  # add actual (normative) MIND-Jaccard correlations (without OB)
  bind_rows(
    df_jaccard_mind %>% 
      nest() %>% 
      mutate(
        rho = map(
          .x = data,
          .f = ~ cor.test(.x$mind_weight, .x$jaccard, method = "spearman")$estimate
        )
      ) %>% 
      unnest(cols = c(rho)) %>% 
      mutate(network = "normative MIND"), .
  ) %>% 
  ungroup %>% 
  mutate(zscore = (rho - mean(rho))/sd(rho))


## Calculate permutation test p-value
n_perms <- df_jaccard_mind_null %>% 
    filter(network != "normative MIND") %>% 
    mutate(network = str_remove(network, "null") %>% as.numeric) %>% 
    pull(network) %>% max
rho_normative <- df_jaccard_mind_null %>%
    filter(network == "normative MIND") %>%
    pull(rho)
rho_nulls <- df_jaccard_mind_null %>%
    filter(network != "normative MIND") %>%
    pull(rho)
1 - sum(rho_normative > rho_nulls)/n_perms # One-tailed p-value: proportion of null phi values >= normative phi



# Save analysis objects for plotting --------------------------------------

save(
  Swanson_roi_order, df_swanson_mind_SwansonLabs, 
  df_swanson_mind_SwansonLabs_present_binned, 
  df_swanson_mind_SwansonLabs_present_connectProb, 
  df_jaccard_mind, df_jaccard_mind_null,
  file = paste0(objects_dir, "03June2024_tract_tracing.RData")
)

