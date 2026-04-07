library(dplyr)
library(ggplot2)
library(ggpubr)
library(rstatix)
library(tidyr)
library(patchwork)

# ============================================================================
# Load data
# ============================================================================
total_reads <- 24228692

# setwd("/ijc/LABS/STIK/RAW/Hi-C/HiC_marta_Oct25/loops/intersected_loops/chip_coverage")

setwd("R:/R_STIK/Hi-C/HiC_marta_Oct25/loops/intersected_loops/chip_coverage")

col_names <- c("chr", "start", "end", "loop_class",
               "read_count", "bases_covered",
               "anchor_length", "fraction_covered")

cov <- read.table("anchor_coverage_full.bed",
                  sep = "\t", col.names = col_names)

rnd <- read.table("random_coverage_full_30kb.bed",
                  sep = "\t", col.names = col_names) %>%
  mutate(loop_class = "random_30kb")

# ============================================================================
# Combine and compute metrics
# ============================================================================
level_order <- c("REH_specific", "EP1_specific", "ARO_specific",
                 "shared_EP1_REH", "shared_ARO_REH", "shared_EP1_ARO",
                 "shared_all3", "random_30kb")

all_cov <- bind_rows(cov, rnd) %>%
  mutate(
    CPM               = (read_count / total_reads) * 1e6,
    log2_CPM          = log2(CPM + 1),
    mean_depth        = read_count / anchor_length,
    mean_depth_pseudo = mean_depth + 1e-6,
    loop_class        = factor(loop_class, levels = level_order)
  )

# ============================================================================
# Colors
# ============================================================================
class_colors <- c(
  "REH_specific"   = "#e74c3c",
  "EP1_specific"   = "#3498db",
  "ARO_specific"   = "#2ecc71",
  "shared_EP1_REH" = "#9b59b6",
  "shared_ARO_REH" = "#e67e22",
  "shared_EP1_ARO" = "#1abc9c",
  "shared_all3"    = "#34495e",
  "random_30kb"    = "#95a5a6"
)

# ============================================================================
# Wilcoxon: each loop class vs random_30kb only — BH correction
# ============================================================================
loop_classes <- setdiff(level_order, "random_30kb")

wilcox_vs_random <- lapply(loop_classes, function(cls) {
  sub_df <- all_cov %>%
    filter(loop_class %in% c(cls, "random_30kb")) %>%
    mutate(loop_class = droplevels(loop_class))
  
  wt <- wilcox_test(sub_df, log2_CPM ~ loop_class,
                    ref.group = "random_30kb") %>%
    mutate(group1 = cls, group2 = "random_30kb")
  wt
}) %>%
  bind_rows() %>%
  adjust_pvalue(method = "BH") %>%
  add_significance()

sig_vs_random <- wilcox_vs_random %>%
  filter(p.adj < 0.05) %>%
  select(group1, group2, p.adj.signif)

# ============================================================================
# Build brackets: each significant loop class vs random
# ============================================================================
x_pos <- setNames(seq_along(level_order), level_order)

y_max   <- max(all_cov$log2_CPM, na.rm = TRUE)
y_start <- y_max * 0.80

sig_ann <- sig_vs_random %>%
  mutate(
    y.position = seq(y_start,
                     y_start + (y_max * 0.18),
                     length.out = n()),
    x    = x_pos[group1],
    xend = x_pos[group2],
    xmid = (x_pos[group1] + x_pos[group2]) / 2
  )

bracket_segments <- bind_rows(
  sig_ann %>% transmute(x = x,    xend = xend, y = y.position,        yend = y.position),
  sig_ann %>% transmute(x = x,    xend = x,    y = y.position - 0.05, yend = y.position),
  sig_ann %>% transmute(x = xend, xend = xend, y = y.position - 0.05, yend = y.position)
)

bracket_labels <- sig_ann %>%
  transmute(x = xmid, y = y.position + 0.07, label = p.adj.signif)

# ============================================================================
# PLOT A — Boxplot
# ============================================================================
p_box <- ggplot(all_cov, aes(x = loop_class, y = log2_CPM, fill = loop_class)) +
  geom_boxplot(
    outlier.size  = 0.3,
    outlier.alpha = 0.4,
    size          = 0.4,
    width         = 0.65
  ) +
  geom_hline(
    yintercept = median(all_cov$log2_CPM[all_cov$loop_class == "random_30kb"]),
    linetype   = "dashed",
    color      = "gray40",
    size       = 0.5
  ) +
  geom_segment(
    data        = bracket_segments,
    aes(x = x, xend = xend, y = y, yend = yend),
    inherit.aes = FALSE,
    size        = 0.4
  ) +
  geom_text(
    data        = bracket_labels,
    aes(x = x, y = y, label = label),
    inherit.aes = FALSE,
    size        = 4
  ) +
  scale_fill_manual(values = class_colors) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.22))) +
  theme_bw() +
  theme(
    panel.grid      = element_blank(),
    axis.text.x     = element_text(angle = 35, hjust = 1, size = 10, face = "bold"),
    axis.text.y     = element_text(size = 10),
    axis.title      = element_text(size = 12),
    plot.title      = element_text(hjust = 0.5, face = "bold", size = 13),
    plot.subtitle   = element_text(hjust = 0.5, size = 9, color = "gray40"),
    legend.position = "none"
  ) +
  xlab("Loop type") + ylab("Coverage (log2 CPM)") +
  labs(
    title    = "E2A-PBX1 ChIP-seq counts at loop anchors",
    subtitle = "Wilcoxon vs random_30kb, BH-adjusted — significant comparisons only"
  )

# ============================================================================
# PLOT B — Violin + boxplot
# ============================================================================
p_violin <- ggplot(all_cov, aes(x = loop_class, y = mean_depth_pseudo,
                                fill = loop_class)) +
  geom_violin(
    alpha  = 0.65,
    scale  = "width",
    trim   = FALSE,
    size   = 0.3
  ) +
  geom_boxplot(
    width         = 0.12,
    outlier.size  = 0.3,
    outlier.alpha = 0.4,
    size          = 0.4,
    fill          = "white",
    color         = "black"
  ) +
  geom_hline(
    yintercept = median(all_cov$mean_depth_pseudo[all_cov$loop_class == "random_30kb"]),
    linetype   = "dashed",
    color      = "gray40",
    size       = 0.5
  ) +
  scale_y_log10(
    breaks = c(0.001, 0.01, 0.1, 1),
    labels = c("0.001", "0.010", "0.100", "1.000")
  ) +
  scale_fill_manual(values = class_colors) +
  theme_bw() +
  theme(
    panel.grid      = element_blank(),
    axis.text.x     = element_text(angle = 35, hjust = 1, size = 10, face = "bold"),
    axis.text.y     = element_text(size = 10),
    axis.title      = element_text(size = 12),
    plot.title      = element_text(hjust = 0.5, face = "bold", size = 13),
    legend.position = "none"
  ) +
  xlab("Loop type") + ylab("Mean depth (log10 scale)") +
  labs(title = "ChIP-seq Coverage Across Loop Anchor Types")

# ============================================================================
# SAVE
# ============================================================================
pdf("EP1_chip_loop+random_counts_boxplot.pdf", width = 9, height = 7)
print(p_box)
dev.off()

pdf("EP1_chip_loop+random_coverage_violin.pdf", width = 9, height = 6)
print(p_violin)
dev.off()

p_combined <- p_box + p_violin + plot_annotation(tag_levels = "A")

pdf("EP1_chip_loop+random_combined.pdf", width = 18, height = 7)
print(p_combined)
dev.off()