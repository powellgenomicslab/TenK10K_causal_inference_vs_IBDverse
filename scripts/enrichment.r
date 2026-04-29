library(paletteer)
library(tidyverse)
library(data.table)
library(patchwork)
library(ggrepel)
library(scales)

input.dir <- "pi1_mr_novelty_results"
df_prop_blood <- fread(paste0(input.dir, "/prop_blood_ibdverse.csv"))
tenk_not_novel_res <- fread(paste0(input.dir, "/tenk_not_novel_res.tsv"))
tenk_not_novel_rel_res <- fread(paste0(input.dir, "/tenk_not_novel_res_relative.tsv")) |> 
    filter(is.finite(OR)) |> 
    left_join(df_prop_blood, by = c("category" = "predicted_category"))

df_annot <- fread(paste0(input.dir, "/IBDverse_annotation_map.tsv")) 

df_annot[, cell_num := seq_along(label_machine), by = Level]

# Are IBDverse cat genes enriched for novelty vs mr_genes?

tenk_not_novel_res = tenk_not_novel_res %>% 
    filter(is.finite(OR)) # MESENCHYMAL has OR of Inf

lm_fit <- lm(log(OR) ~ log10(prop_blood), data = tenk_not_novel_res)
slope <- signif(coef(lm_fit)[[2]], 3)
pval <- signif(summary(lm_fit)$coefficients["log10(prop_blood)", "Pr(>|t|)"], 3)
# Print
print(paste0("Slope: ", slope, ". p=", pval))

# Plot this
tenk_not_novel_res$annotation_type = "Major population"
ybreaks_or = c(0.25, 0.5, 0.75, 1, 1.33, 2, 4)
ggplot(tenk_not_novel_res, aes(x = log10(prop_blood), y = log(OR))) + 
  geom_point(aes(fill = log(OR)), size = 2, shape = 21, stroke = 0.7) +
  geom_smooth(method = "lm", se = FALSE, color = "black", linetype = "dashed", linewidth = 0.5) + 
  geom_hline(yintercept = 0, linetype = "dashed", color = "lightgrey", linewidth = 0.5) +
  geom_text_repel(aes(label = category), size = 4, max.overlaps = Inf, point.padding = unit(5, "lines")) +
  # scale_size_manual(values = annot.class.sizes, name = 'Annotation granularity') +
  # scale_fill_manual(values = umap.category.palette, name = 'Major population') +
  labs(
    x = "log10 fraction of cells derived from  blood samples in IBDverse",
    y = "Odds ratio for IBDverse nominated effector gene\nalready found in TenK10K (log scale)"
  ) +
  theme_classic() +
  scale_fill_viridis_b(
    breaks = log(ybreaks_or),
    labels = ybreaks_or,
    name = "Odds ratio",
    begin = 0,
    limits = range(log(ybreaks_or))) +
  # annotate with slope coeff
  annotate("text", x = -Inf, y = Inf, label = paste0("β = ", slope), hjust = -0.2, vjust = 1.2) +
  scale_y_continuous(
      breaks = log(ybreaks_or),
      labels = function(x) signif(exp(x), 2)) + # Convert log(OR) back to OR for y-axis labels
  theme(legend.position = "right")

out.dir <- "figures"
ggsave(paste0(out.dir,"/Blood_proportion_vs_OR_effector_genes.png"), width = 6.5, height = 5)


# Are IBDverse cat genes RELATIVELY enriched for non-novelty vs other IBDverse categories?

lm_fit <- lm(log(OR) ~ log10(prop_blood), data = tenk_not_novel_rel_res)
slope <- signif(coef(lm_fit)[[2]], 3)
pval <- signif(summary(lm_fit)$coefficients["log10(prop_blood)", "Pr(>|t|)"], 3)
print(paste0("Slope: ", slope, ". p=", pval))

# ybreaks_or = c(1, 1.5, 2, 3, 4, 6)
# ybreaks_log_or = log(ybreaks_or)
tenk_not_novel_rel_res$annotation_type = "Major population"

ggplot(tenk_not_novel_rel_res, aes(x = log10(prop_blood), y = log(OR))) +
  annotate("rect", xmin = -Inf, xmax = -2, ymin = -Inf, ymax = Inf,
           fill = "#4DAF4A", alpha = 0.2) +
  annotate("rect", xmin = -2, xmax = Inf, ymin = -Inf, ymax = Inf,
           fill = "#377EB8", alpha = 0.2) +
  annotate("text", x = -Inf, y = -Inf, label = "Gut-specific (<1% blood-derived cells)", 
           hjust = -0.1, vjust = -0.5, color = "#4DAF4A", size = 3.5, fontface = "plain") +
  annotate("text", x = Inf, y = -Inf, label = "Non-gut-specific", 
           hjust = 1.1, vjust = -0.5, color = "#377EB8", size = 3.5, fontface = "plain") +
  geom_smooth(method = "lm", se = FALSE, color = "black", linetype = "dashed", linewidth = 0.5) +
  # geom_errorbar(aes(ymin = log(CI_low), ymax = log(CI_high), color = category), width = 0.05, linewidth = 0.4) +
  geom_point(aes(fill = category), size = 4, shape = 21, stroke = 0.7) +
  geom_label_repel(aes(label = category), size = 9/.pt, force = 5, min.segment.length = 0,
    label.size = 0, max.overlaps = Inf, label.padding = 0.1, point.padding = unit(5, "lines"), fill = alpha("white", .8)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "lightgrey", linewidth = 0.5) +
  scale_fill_paletteer_d("ggsci::category10_d3", name = NULL) +
  scale_color_paletteer_d("ggsci::category10_d3", name = NULL, guide = "none") +
  labs(
    x = "Fraction of cells derived from blood samples in IBDverse",
    y = "Odds ratio for IBDverse-nominated \neffector gene already found in TenK10K"
  ) +
  theme_classic() +
  scale_y_continuous(
    breaks = log(c(1, 1.5, 2, 3, 4, 5)),
    labels = exp) +
  scale_x_continuous(
    breaks = c(min_x, -4, -3, -2, -1, log10(0.5)), 
    labels = c(0, "0.01%", "0.1%", "1%", "10%", "50%")) +
  annotate("label", x = -Inf, y = Inf, label = paste0("β = ", slope, ", p = ", pval), hjust = -0.2, vjust = 1.2) +
  theme(legend.position = "none")

ggsave(paste0(out.dir, "/Blood_proportion_vs_OR_effector_genes_relative.png"), device = ragg::agg_png, scaling = 1, width = 5.5, height = 5)

# cell-type level results

# min_x = pseudocount for cell type with prop_blood = 0
min_x <- -5

tenk_not_novel_rel_res_celltype <- fread(paste0(input.dir, "/tenk_not_novel_res_relative_celltype.tsv")) |> 
    filter(is.finite(OR), is.finite(prop_blood)) |> 
    mutate(prop_blood_mod = ifelse(prop_blood == 0, 10^min_x, prop_blood))
tenk_not_novel_rel_res_celltype[order(category_new, label_machine),  cell_num := seq_along(label_machine)]

lm_fit <- lm(log(OR) ~ log10(prop_blood_mod), data = tenk_not_novel_rel_res_celltype)
slope <- signif(coef(lm_fit)[[2]], 3)
pval <- signif(summary(lm_fit)$coefficients["log10(prop_blood_mod)", "Pr(>|t|)"], 3)
print(paste0("Slope: ", slope, ". p=", pval))

# ybreaks_or = c(1, 1.5, 2, 3, 4, 6)
# ybreaks_log_or = log(ybreaks_or)
# tenk_not_novel_rel_res$annotation_type = "Major population"


p_main <- tenk_not_novel_rel_res_celltype |> 
  ggplot(aes(x = log10(prop_blood_mod), y = log(OR))) +
  geom_smooth(method = "lm", se = FALSE, color = "black", linetype = "dashed", linewidth = 0.5) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "lightgrey", linewidth = 0.5) +
  geom_text_repel(aes(label = cell_num, color = category_new), size = 6/.pt) +
  geom_point(aes(fill = category_new), size = 4, shape = 21, stroke = 0.7) +
  labs(
    x = "Fraction of cells derived from blood samples in IBDverse",
    y = "Odds ratio for IBDverse-nominated \neffector gene already found in TenK10K"
  ) +
  theme_classic() +
  scale_fill_paletteer_d("ggsci::category10_d3", name = NULL) +
  scale_color_paletteer_d("ggsci::category10_d3", name = NULL) +
  annotate("text", x = -Inf, y = Inf, label = paste0("β = ", slope, ", p = ", pval), hjust = -0.2, vjust = 1.2) +
  scale_y_continuous(
    breaks = log(c(0.1, 0.2, 0.5, 1, 2, 5, 10, 20)),
    labels = c(0.1, 0.2, 0.5, 1, 2, 5, 10, 20)) +
  scale_x_continuous(
    breaks = c(min_x, -4, -3, -2, -1, log10(0.5)), 
    labels = c(0, "0.01%", "0.1%", "1%", "10%", "50%")) +
  theme(legend.position = "none")

# cell_order <- data.table(
#   category =c("T", "B", "Myeloid", "Plasma", "Stem", "Colonocyte", "Enterocyte", "Secretory"),
#   col = c(1, 1, 2, 2, 2, 3, 3, 4)
# )

(p_index <- tenk_not_novel_rel_res_celltype |>
  arrange(category_new, cell_num) |>
  mutate(row = row_number()) |>
  ggplot(aes(x = 0, y = -row * 1.5)) +
  geom_text(aes(label = paste0(cell_num, ": ", label_new), color = category_new),
            hjust = 0, size = 6/.pt, lineheight = 1) +
  scale_color_paletteer_d("ggsci::category10_d3") +
  facet_wrap(~ category_new, ncol = 1, scales = "free_y", space = "free_y", strip.position = "top") +
  theme_void() +
  theme(
    legend.position = "none",
    strip.text = element_text(size = 7, face = "bold", angle = 0, hjust = 0, margin = margin(b = 0)),
    panel.spacing.y = unit(0.25, "lines"),
    strip.clip = "off"
  ) +
  coord_cartesian(clip = "off") +
  scale_x_continuous(expand = expansion(mult = c(0, 0.1))) +
  scale_y_continuous(expand = expansion(add = c(1)))
)

plots_or_celltype <- (p_main | p_index) + plot_layout(widths = c(3, 1))

ggsave(paste0(out.dir, "/Blood_proportion_vs_OR_effector_genes_relative_celltype.png"),
      width = 9, height = 6, device = ragg::agg_png, scaling = 1)
