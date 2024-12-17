
library(pheatmap)
df_normative_mind %>% 
  mutate(median_weight = ifelse(R1 == R2, NA_real_, median_weight)) %>% 
  pivot_wider(id_cols = R1, names_from = R2, values_from = median_weight) %>% 
  column_to_rownames("R1") %>% 
  pheatmap()

# hclust ------------------------------------------------------------------

## Convert to matrix
normative_mind <- df_normative_mind %>% 
  pivot_wider(id_cols = R1, names_from = R2, values_from = median_weight) %>% 
  column_to_rownames("R1")

## Run hclust
hclust_mind <- normative_mind %>% 
  as.dist() %>% 
  hclust()

## Define ROI order
roi_order <- normative_mind[hclust_mind$order,] %>% rownames

## Define n modules
cutree(hclust_mind, k = 2) %>% 
  enframe(name = "region_of_interest", value = "module") %>% 
  group_by(module) %>% 
  mutate(row = row_number(), module = paste0("M", module)) %>% 
  pivot_wider(id_cols = row, names_from = module, values_from = region_of_interest) %>% 
  dplyr::select(-row) %>% 
  mutate_all( ~ ifelse(is.na(.x), "", .x))
  

# Plot module solution ----------------------------------------------------


df_normative_mind %>%
  mutate(R1 = factor(R1, levels = roi_order),
         R2 = factor(R2, levels = roi_order)) %>%
  mutate(median_weight = ifelse(R1 == R2, NA_real_, median_weight)) %>% 
  
  # plot
  ggplot(aes(x = R1, y = R2, fill = median_weight)) +
  geom_tile() +
  
  # # add module lines
  # geom_vline(xintercept = module_lines, color = "white", linewidth = 0.5) +
  # geom_hline(yintercept = module_lines, color = "white", linewidth = 0.5) +
  # 
  # # # add module colors as side bars
  # annotate(ymin = c(-Inf, module_lines),
  #          ymax = c(module_lines, Inf),
  #          xmin = length(roi_order) + 0.5, xmax = length(roi_order) + 3,
  #          geom = "rect",
  #          fill = module_colors) +
  # annotate(xmin = c(-Inf, module_lines),
  #          xmax = c(module_lines, Inf),
  #          ymin = length(roi_order) + 0.5, ymax = length(roi_order) + 3,
  #          geom = "rect",
  #          fill = module_colors) +
  
  # plot aesthetics
  labs(x = "", y = "", title = "A | Normative adult (PND 63) connectome") +
  scale_fill_viridis(na.value = "white") +
  guides(fill = guide_colourbar(title.position = "top", title.hjust = 0.5, title = "edge weight")) +
  theme(axis.text = element_blank(),
        axis.ticks = element_blank(),
        legend.position = "bottom",
        legend.key.height = unit(0.4, "cm"),
        legend.key.width = unit(2.3, "cm")
  )  
