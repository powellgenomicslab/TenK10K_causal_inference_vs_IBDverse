library(tidyverse)
library(data.table)
library(patchwork)
library(paletteer)
library(marquee)

input.dir <- "pi1_mr_novelty_results"
pi1_combined <- fread(paste0(input.dir, "/pi1_combined.tsv"))
ibdverse_immune = c("B", "Plasma", "Myeloid", "T")
df_annot <- fread(paste0(input.dir, "/IBDverse_annotation_map.tsv")) 

comparison <- list(
    all = expression(TRUE),
    gut = expression(
        prefix != "unannotated" &
        !(suffix %in% c("blood", "ct"))),
    gut_nonimmune = expression(gut & !(prefix %in% ibdverse_immune))
)

pi1_all <- fread(paste0(input.dir, "/pi1_all.txt")) |> 
      mutate(
        label_machine = gsub("dMean__", "", dataset_id),
        label_machine = gsub("_all", "", label_machine),
        # prefix: major cell type
        prefix = str_remove_all(label_machine, "_.*"),
        # suffix: tissue
        suffix = str_remove_all(label_machine, ".*_"),
        # mid: minor cell type
        mid = str_remove_all(label_machine, paste0("_", suffix, "$"))
    )

for (c in names(comparison)) {
  pi1_all[, (c) := eval(comparison[[c]])]
}
median_combined <- fread(paste0(input.dir, "/pi1_median_combined.tsv"))

pi1_summary <- pi1_all |> 
    group_by(discovery) |> 
    summarise(median_pi1 = median(pi1, na.rm = TRUE),
              max_pi1 = max(pi1, na.rm = TRUE),
              max_dataset_id = dataset_id[which.max(pi1)]
              min_pi1 = min(pi1, na.rm = TRUE))

tissue_labels <- c("r" = "Rectum", "ti" = "Terminal ileum", "blood" = "Blood", "ct" = "Cultured cells") |> 
  enframe() |> 
  setDT()
pi1_all[tissue_labels, tissue_label := i.value, on = .(suffix = name)]
pi1_all[df_annot, label_new := i.label_new, on = .(label_machine)]

# heatmap
# heatmap for gut
bold_labels <- df_annot[Level == 1, label_new]
pi1_all[, label_format := ifelse(label_new %in% bold_labels, paste0("{.red3 **", label_new, "**}"), label_new)]
pi1_gut_max <- pi1_all[gut == TRUE][
    , .SD[which.max(pi1)], by = .(discovery)
]

p_gut_heatmap <- pi1_all[gut == TRUE] |> 
    ggplot(aes(x = label_format, y = discovery))  +
    geom_tile(aes(fill = pi1), color = "gray") +
    scale_fill_viridis_c(limits = c(0, 1), direction = -1, option = "magma") +
    facet_grid(
        cols = vars(prefix), rows = vars(tissue_label), scales = "free", space = "free") +
    geom_text(
        data = pi1_gut_max,
        aes(x = label_format, y = discovery, color = "green3"),
        label = "*", size = 5, vjust = 0.75, inherit.aes = FALSE
    ) +
    scale_color_identity(labels = c("green3" = "Max π1"),
        guide = guide_legend(theme = theme(panel.background = element_blank()))) +
    scale_x_discrete(
        expand = expansion(add = 0.05)
    ) +
    scale_y_discrete(expand = expansion(add = 0.05), labels = ~str_replace_all(.x, "_", " ")) +
    labs(x = NULL, y = NULL, fill = "π1", color = NULL,
         caption = "*{.red3 **Red**} cell types include\nall cells from the lineage*") +
    theme(axis.text.x = element_marquee(angle = 90, hjust = 1, vjust = 0.5),
          strip.clip = "off",
          strip.background = element_blank(),
          panel.background = element_rect(fill = NA, color = "black", linewidth = 0.5),
          panel.ontop = TRUE,
          panel.grid = element_blank(),
          axis.ticks = element_line(),
          strip.text.x = element_text(size = 9, face = "bold", angle = 90, hjust = 0),
          strip.text.y = element_text(size = 9, face = "bold"),
          plot.caption = element_marquee(size = 9, margin = margin(t = -10))
    )

# Save figures
ggsave("figures/pi1_gut_heatmap_by_celltype.png", p_gut_heatmap,
       scaling = 1, device = ragg::agg_png, width = 12, height = 10)

# Summarise the max, min and max cell-type per comparison
pi1_combined %>% 
    group_by(comparison) %>% 
    summarise(
        min = min(pi1), 
        max = max(pi1), 
        max_celltype = major_cell_type[which.max(pi1)],
        .groups = "drop"
    ) +
    theme_bw(x )

ggplot(pi1_combined, aes(x = pi1, fill = major_cell_type)) +
    geom_histogram(bins = 50, color = "black") +
    facet_grid(major_cell_type ~ comparison, scales = "free_y") +
    geom_vline(data = median_combined,
               aes(xintercept = median_pi1),
               lty = "dashed", color = "black", size = 1) +
    geom_text(
        data = median_combined,
        aes(x = 0, y = Inf, label = paste0("median: ", signif(median_pi1, 2))),
        inherit.aes = FALSE, hjust = 0, vjust = 1.5, size = 3.5
    ) +
    theme_classic() +
    theme(
        legend.position = "none",
        strip.background = element_blank(),
        strip.text.x = element_text(size = 10, face = "bold"),
        strip.text.y = element_text(size = 10, face = "bold"),
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12),
        axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14)
    ) +
    labs(x = "Maximum π1 per cell-type", y = "Frequency") +
    xlim(c(0, 1))

out.dir <- "replication"
ggsave(paste0(out.dir, "/pi1_TenK10K_vs_IBDverse.png"), width = 12, height = 12)

# Plot this as a boxplot instead
median_by_comparison <- pi1_combined |>
    summarise(median_pi1 = median(pi1, na.rm = TRUE), .by = comparison)

ggplot(pi1_combined, aes(x = major_cell_type, y = pi1, fill = major_cell_type)) +
    geom_boxplot(outlier.shape = NA, alpha = 0.3, whisker.linewidth = 0.2, box.linewidth = 0.2, median.linewidth = 0.2) +
    geom_point(aes(fill = major_cell_type), shape = "circle filled", position = position_dodge(width = 0.5), size = 1, alpha = 0.7) +
    geom_hline(
        data = median_by_comparison,
        aes(yintercept = median_pi1, color = "red3"),
        lty = "dashed", linewidth = 0.5
    ) +
  scale_fill_paletteer_d("ggsci::category10_d3", name = NULL) +
    scale_color_identity(guide = guide_legend(), labels = c("red3" = "Median π1")) +
    theme_classic() +
    theme(
        legend.position = "right",
        strip.background = element_blank(),
        strip.clip = "off",
        strip.text.x = element_text(size = 11, face = "bold"),
        strip.text.y = element_text(size = 12, face = "bold"),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.y = element_text(size = 12),
        axis.title.x = element_blank(),
        panel.spacing.x = unit(1, "lines"),
        axis.title.y = element_text(size = 14)
    ) +
    facet_grid(. ~ comparison, labeller = label_wrap_gen(width = 30)) +
    labs(x = "Major cell-type", y = "Maximum π1 per cell-type", fill = NULL, color = NULL) +
    ylim(c(0, 1))

out.dir <- "figures"
ggsave(paste0(out.dir, "/pi1_TenK10K_vs_IBDverse_boxplot_v2.png"), width = 10, height = 6, scaling = 1.2)

library(ggbeeswarm)

ggplot(pi1_combined, aes(x = major_cell_type, y = pi1, fill = comparison, color = comparison)) +
    geom_beeswarm(dodge.width = 0.8, size = 1.5, alpha = 0.7) +
    theme_classic() +
    theme(
        legend.position = "bottom",
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12),
        axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14)
    ) +
    labs(x = "Major cell-type", y = "Maximum π1 per cell-type", fill = NULL, color = NULL) +
    ylim(c(0, 1))

ggsave(paste0(out.dir, "/pi1_TenK10K_vs_IBDverse_beeswarm.png"), width = 10, height = 8)
  