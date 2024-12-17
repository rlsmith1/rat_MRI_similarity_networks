
## Load libraries
library(tidyverse)

## Set paths
base_dir <- "~/Documents/PhD/projects/CamRat/CamRat/"
load(paste0(base_dir, "objects/ratdb_scans.Rdata")) # df_scans

## Write txt file of scans with file path to collect RG from
df_scans %>%
  dplyr::select(ids, dages, type, acqps) %>% 
  filter(str_detect(ids, "JWD|EDAA") & type %in% c("m", "p", "t")) %>% 
  mutate(acqps = str_remove(acqps, ".*scannerfiles/"), # remove path, extract only filename
         subject = str_remove(ids, ".*-") %>% str_remove(".*_"),
         contrast = case_when(
           type == "m" ~ "MTw",
           type == "p" ~ "PDw",
           type == "t" ~ "T1w"
         )
  ) %>%
  # add timepoint
  mutate(timepoint = .bincode(dages, c(0, 25, 35, 70, 310))) %>% 
  mutate(timepoint = case_when(
    
    timepoint == 1 ~ 20,
    timepoint == 2 ~ 35,
    timepoint == 3 ~ 63,
    timepoint == 4 ~ 300
    
  ) %>% factor(levels = c(20, 35, 63, 300))) %>% 
  dplyr::rename("age" = "dages") %>% 
  dplyr::select(subject, timepoint, age, contrast, acqps) %>% 
  arrange(subject, timepoint, contrast) #%>% 

  write.table(file = paste0(base_dir, "scripts/acqps_file_paths.txt"), sep = "\t", row.names = FALSE, quote = FALSE)

### Run extract_receiver_gains.sh on Cambridge HPC ###

## Read results
df_rg <- read_table(paste0(base_dir, "scripts/acqps_RG.txt")) %>% 
  
  # add timepoint for plotting
  mutate(timepoint = .bincode(dages, c(0, 25, 35, 70, 310))) %>% 
  mutate(timepoint = case_when(
    
    timepoint == 1 ~ 20,
    timepoint == 2 ~ 35,
    timepoint == 3 ~ 63,
    timepoint == 4 ~ 300
    
  ) %>% factor(levels = c(20, 35, 63, 300))) %>%
  filter(!is.na(timepoint)) %>% 
  mutate(study = ifelse(str_detect(ids, "EDA"), "GSK", "MRC"))

## Plot receiver gains

unique_rgs <- df_rg %>% pull(RG) %>% unique
df_rg %>% pull(RG) %>% table


# overall
df_rg %>% 
  ggplot(aes(x = RG)) + 
  geom_histogram()

# by timepoint
df_rg %>% 
  ggplot(aes(x = RG)) +
  geom_histogram() +
  facet_wrap(vars(timepoint))

# by scan contrast
df_rg %>% 
  ggplot(aes(x = RG)) +
  geom_histogram() +
  facet_wrap(vars(type))

# by study
df_rg %>% 
  ggplot(aes(x = RG)) +
  geom_histogram() +
  facet_wrap(vars(study))

# all
df_rg %>% 
  ggplot(aes(x = RG)) +
  geom_histogram(aes(fill = type), position = "dodge") +
  facet_grid2(type ~ timepoint, scales = "free", independent = "all") +
  scale_x_continuous(breaks = unique_rgs)

## Do receiver gain values differ across subject-scans at the same timepoint but different contrasts?
df_rg_compare <- df_rg %>% 
  group_by(ids, timepoint, type) %>% 
  mutate(row = row_number()) %>% 
  dplyr::rename("age" = "dages") %>% 
  mutate(type = 
           case_when(
             type == "m" ~ "MT",
             type == "p" ~ "PD",
             type == "t" ~ "T1"
           ) %>% 
           paste0("_RG")
  ) %>% 

  pivot_wider(id_cols = c(study, ids, age, timepoint, row), names_from = type, values_from = RG) %>% 
  dplyr::select(-row) %>% 
  mutate(same_RG_across_scans = MT_RG == PD_RG & MT_RG == T1_RG)
  
# export table
df_rg_compare %>% write_csv(paste0(base_dir, "scripts/acqs_RG_comparison.csv"))

df_rg_compare %>% ungroup %>% filter(!(str_detect(ids, "EDA") & timepoint == 20)) %>% count(study, timepoint, same_RG_across_scans) %>% arrange(study, timepoint)

