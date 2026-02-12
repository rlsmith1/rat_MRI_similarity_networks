
#==============================================================================#
# Calculate the rich club coefficient of the normative connectome compared to 10000 null nets
#==============================================================================#

## Set up
base_dir <- "~/Documents/PhD/projects/CamRat/CamRat/"
source(paste0(base_dir, "code/03_data-analysis/setup.R"))


### HOW TO CALCULATE RICH CLUB

# 1. Identify prominent nodes in network (in this case, define as hubs from our PND 63 normative network)
# 2. Calculate the total sum of weights attached to ties among the prominent nodes
# 3. Calculate the total sum of same number of connections, but the strongest ones, in the network

###


# Calculate rich club coefficient of normative network --------------------

n_connections <- length(hubs)*length(hubs) - length(hubs)

df_rich_club <- df_normative_mind %>% 
  nest() %>% 

  mutate(
    hub_sum = map(
      .x = data,
      .f = ~ .x %>% 
        filter(R1 %in% hubs & R2 %in% hubs & R1 != R2) %>% 
        summarise(hub_sum = sum(median_weight, na.rm = TRUE)) %>% 
        pull(hub_sum)
    ),
    strong_sum = map(
      .x = data,
      .f = ~ .x %>% 
        filter(R1 != R2) %>% 
        slice_max(order_by = median_weight, n = n_connections) %>% 
        summarise(strong_sum = sum(median_weight, na.rm = TRUE))
    )
  ) %>% 
  unnest(cols = c(hub_sum, strong_sum)) %>% 
  mutate(phi = hub_sum/strong_sum)



# Calculate rich club of null networks ------------------------------------


## Load null network data
load(paste0(objects_dir, "03June2024_null_nets.RDS")) # df_null_nets (generated in data_analysis/01.validate_MIND/withOB/generate_null_nets.R)

## Calculate the rich club coefficient of each null network
df_rich_club_nulls <- df_null_nets %>%
  group_by(network) %>% 
  nest() %>% 
  
  mutate(
    hub_sum = map(
      .x = data,
      .f = ~ .x %>% 
        filter(R1 %in% hubs & R2 %in% hubs & R1 != R2) %>% 
        summarise(hub_sum = sum(weight, na.rm = TRUE)) %>% 
        pull(hub_sum)
    ),
    strong_sum = map(
      .x = data,
      .f = ~ .x %>% 
        filter(R1 != R2) %>% 
        slice_max(order_by = weight, n = n_connections) %>% 
        summarise(strong_sum = sum(weight, na.rm = TRUE))
    )
  ) %>% 
  unnest(cols = c(hub_sum, strong_sum)) %>% 
  mutate(phi = hub_sum/strong_sum)



# Combine actual with nulls & save! ---------------------------------------

## Combine normative connectome rich club with null networks
df_rich_club_all <- df_rich_club %>% 
  mutate(network = "normative", .before = 1) %>% 
  bind_rows(df_rich_club_nulls) %>% 
  dplyr::select(-data) %>% 
  mutate(z_score = (phi - mean(phi))/sd(phi)
         )

## Calculate permutation test p-value
n_perms <- df_rich_club_all %>% 
    filter(network != "normative") %>% 
    mutate(network = str_remove(network, "null") %>% as.numeric) %>% 
    pull(network) %>% max
phi_normative <- df_rich_club_all %>%
    filter(network == "normative") %>%
    pull(phi)
phi_nulls <- df_rich_club_all %>%
    filter(network != "normative") %>%
    pull(phi)
1 - sum(phi_normative > phi_nulls)/n_perms # One-tailed p-value: proportion of null phi values >= normative phi


## Save to analysis dir for plotting
save(df_rich_club_all, file = paste0(objects_dir, "01102024_rich_club_distributions.RDS"))
