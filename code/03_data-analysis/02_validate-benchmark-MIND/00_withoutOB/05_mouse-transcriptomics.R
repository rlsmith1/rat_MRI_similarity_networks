
#==============================================================================#
# Compare normative MIND network with region-based single-cell gene expression
# profiles in the mouse brain using Allen Mouse Brain Atlas data
# Hypothesis: Regions with more similar gene expression profiles will have higher MIND similarity
# source: https://www.nature.com/articles/s41586-023-06812-z
#==============================================================================#

base_dir <- "~/Documents/PhD/projects/CamRat/CamRat/"
source(paste0(base_dir, "code/03_data-analysis/setup.R"))


# Load AMBA data ----------------------------------------------------------

## Point to AMBA analysis data directory
data_dir <- paste0(base_dir, "data/AMBA/")

## AMBA atlas abbreviations
df_amba_acronyms <- read_xlsx(paste0(data_dir, "AMBA_parcellation_definitions.xlsx"))

## Load spatial transcriptomic expression data generated on cluster (rename structure with MIND alignments)
df_merfish_expr <- read.csv(paste0(data_dir, "2024Feb28_MERFISH_exp_mat.csv")) %>% as_tibble()

## AMBA to WHS atlas mappings
df_whs_to_aba <- read_xlsx(paste0(base_dir, "data/WHS_analysis/BMrat_to_ABAmouse.xlsx"), sheet = 1) %>% 
  dplyr::select(-parcellation_division) %>% 
  filter(!is.na(parcellation_structure)) %>% # remove regions not present in ABA data
  distinct() %>%
  
  # include cortical regions and OB regions for now (can remove later, this is for supplement)
  filter(WHS %in% cortical_ROIs | WHS %in% olfactory_ROIs)


## NOTE: group by ABA except amygdala; endopiriform nucleus, pretectal region; Secondary visual area, lateral part;
## Secondary visual area, medial part, Septal region; Subparafascicular nucleus; Ventral pallidum;
## Ventral posterior nucleus of the thalamus, parvicellular part; ventral striatal region, unspecified;
## Retrosplenial dysgranular area


# Format expression data --------------------------------------------------

## Identify genes that are not expression across any regions
gene_sums <- df_merfish_expr %>% 
  dplyr::select(starts_with("ENSMUSG")) %>% 
  colSums
no_expr_genes <- gene_sums[gene_sums == 0] %>% names
length(no_expr_genes) # 63 genes have 0 expression across any regions

## Remove genes that are not expressed at all
df_merfish_expr <- df_merfish_expr %>% dplyr::select(-all_of(no_expr_genes))

## Test for normality in the expression distribution of each gene (to determine if Pearson vs Spearman correlation is appropriate)
df_merfish_expr_shapiro <- df_merfish_expr %>% 
  dplyr::select(starts_with("ENSMUSG")) %>% 
  map_dfr(~shapiro.test(.x)$p.value) %>% 
  pivot_longer(1:ncol(.), names_to = "gene_id", values_to = "p_value") %>% 
  mutate(p_adj = p.adjust(p_value, method = "fdr"))
df_merfish_expr_shapiro %>% filter(p_adj > 0.05) # no genes are normally expressed, use spearman correlation


# Calculate expression similarity network using only ROIs in MIND (cortical) ---------------------------------------

## Filter expression data for regions that are included in MIND
df_merfish_expr_mind <- df_whs_to_aba %>% 
  left_join(
    df_merfish_expr %>% 
      dplyr::select(-parcellation_structure),
    by = join_by(parcellation_substructure)
  ) %>% 
  dplyr::select(parcellation_index, parcellation_division, parcellation_structure, parcellation_substructure, everything()) %>% 
  filter(!is.na(parcellation_substructure)) # remove Medial geniculate body, marginal zone (no match in WHS)

## Generate gene expression similarity network

# calculate expression similarity matrix as Spearman correlation of expression profiles across genes
merfish_cor <- df_merfish_expr_mind %>% 
  group_by(parcellation_structure) %>% 
  summarise_at(vars(contains("ENSMUSG")), ~ median(.x, na.rm = TRUE)) %>% 
  column_to_rownames("parcellation_structure") %>% 
  t() %>% 
  correlate(method = "spearman") %>% 
  column_to_rownames("term") %>% as.matrix

# convert to df
df_merfish_cor <- merfish_cor %>% 
  as.data.frame %>% 
  rownames_to_column("R1") %>% 
  as_tibble() %>% 
  pivot_longer(2:ncol(.), names_to = "R2", values_to = "value")

# define ROI order for plot
# merfish_roi_order <- merfish_cor[hclust(dist(merfish_cor))$order,] %>% rownames
# length(merfish_roi_order) # n = 43 (without OB regions)

df_hierarchy_merfish <- df_system_hierarchy %>% 
    left_join(df_whs_to_aba %>% dplyr::select(WHS, parcellation_structure),
              by = join_by("region_of_interest" == "WHS")
    )
merfish_roi_order <- df_hierarchy_merfish %>% 
    mutate(system = factor(system, levels = system_order)) %>% 
    arrange(system) %>% 
    filter(!is.na(parcellation_structure)) %>% 
    pull(parcellation_structure) %>% 
    unique
merfish_roi_order_fullName <- df_amba_acronyms %>% 
    mutate(acronym = factor(acronym, levels = merfish_roi_order)) %>%
    filter(!is.na(acronym)) %>% 
    arrange(acronym) %>% 
    pull(label)



# Combine expression similarity network with MIND weights ------------------

df_mind_merfish <- df_normative_mind %>% 
  
  # add ABA labels
  left_join(df_whs_to_aba %>% 
              dplyr::select(-parcellation_substructure) %>% 
              dplyr::rename_at(vars(contains("parcellation")), ~ paste0(.x, 1)), 
            by = join_by("R1" == "WHS")
  ) %>% 
  left_join(df_whs_to_aba %>% 
              dplyr::select(-parcellation_substructure) %>% 
              dplyr::rename_at(vars(contains("parcellation")), ~ paste0(.x, 2)), 
            by = join_by("R2" == "WHS")
  ) %>% 
  
  # remove non-existent and self connections
  filter(!is.na(parcellation_structure1) & !is.na(parcellation_structure2) & parcellation_structure1 != parcellation_structure2) %>% 
  dplyr::select(-sd_weight, -z_score) %>% 
  dplyr::rename("mind_weight" = "median_weight") %>% 
  
  # add gene expression similarity data
  left_join(df_merfish_cor, 
            by = join_by("parcellation_structure1" == "R1",
                         "parcellation_structure2" == "R2"
            )
  ) %>% 
  dplyr::rename("gene_cor" = "value") %>% 
  distinct() %>% 
  
  # take median across parcellation_structure levels (from AMBA)
  group_by(parcellation_structure1, parcellation_structure2) %>% 
  summarise_if(is.numeric, ~ median(.x, na.rm. = TRUE)) %>% 
  ungroup %>% 
  mutate(gene_cor = ifelse(parcellation_structure1 == parcellation_structure2, 0, gene_cor),
         parcellation_structure1 = factor(parcellation_structure1, levels = merfish_roi_order),
         parcellation_structure2 = factor(parcellation_structure2, levels = merfish_roi_order)
  )


# Validate MIND-gene expression correlation against null networks --------------------------

## Load distance-permuted null networks
load(paste0(objects_dir, "03June2024_null_nets.RDS")) # df_null_nets (generated in data_analysis/01.validate_MIND/withoutOB/generate_null_nets.R)

## Map atlas labels in null networks and combine with actually gene expression similarity values (by edge)
df_null_nets_mouseTranscriptome <- df_null_nets %>% 
  
  # add ABA labels
  left_join(df_whs_to_aba %>% 
              dplyr::select(-parcellation_substructure) %>% 
              dplyr::rename_at(vars(contains("parcellation")), ~ paste0(.x, 1)), 
            by = join_by("R1" == "WHS")
  ) %>% 
  left_join(df_whs_to_aba %>% 
              dplyr::select(-parcellation_substructure) %>% 
              dplyr::rename_at(vars(contains("parcellation")), ~ paste0(.x, 2)), 
            by = join_by("R2" == "WHS")
  ) %>% 
  filter(!is.na(parcellation_structure1) & !is.na(parcellation_structure2)) %>% 
  distinct() %>% 
  
  # add gene expression correlation values
  ungroup %>% 
  dplyr::select(network, parcellation_structure1, parcellation_structure2, weight) %>% 
  dplyr::rename_at(2:3, ~ c("R1", "R2")) %>% 
  left_join(df_merfish_cor)

## Calculate correlation of each null network with the mouse transcriptomic similarity vector
df_null_nets_cor <- df_null_nets_mouseTranscriptome %>% 
  group_by(network) %>% 
  nest() %>% 
  mutate(spearmans_rho = map(
    .x = data,
    .f = ~ cor.test(.x %>% filter(R1 != R2) %>% filter(R1 != R2) %>% pull(weight), 
                    .x %>% filter(R1 != R2) %>% filter(R1 != R2) %>% pull(value), 
                    method = "spearman") %>% 
      .$estimate
  )
  ) %>% 
  unnest(cols = c(spearmans_rho))

## Combine with actual correlation data
spearman_rho <- cor.test(
  df_mind_merfish %>% filter(parcellation_structure1 != parcellation_structure2) %>% filter(mind_weight != 0) %>% pull(mind_weight),
  df_mind_merfish %>% filter(parcellation_structure1 != parcellation_structure2) %>% filter(mind_weight != 0) %>% pull(gene_cor),
  method = "spearman"
) %>% .$estimate
df_null_nets_cor_mind <- df_null_nets_cor %>% 
  dplyr::select(-data) %>% 
  bind_rows(
    tibble(
      network = "normative MIND",
      spearmans_rho = spearman_rho
    ), .
  ) %>% 
  mutate(zscore = (spearmans_rho - mean(spearmans_rho))/sd(spearmans_rho))


# Save analysis objects for plotting --------------------------------------

save(
    df_amba_acronyms, df_whs_to_aba, merfish_roi_order, 
    merfish_roi_order_fullName, df_mind_merfish, df_null_nets_cor_mind,
    file = paste0(objects_dir, "08July2024_gene_expression.RData")
)
