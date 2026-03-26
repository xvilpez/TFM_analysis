library(dplyr)
library(ggplot2)
library(ggpubr)
library(rstatix)
library(tidyr)
library(cowplot)

# ============================================================================
# Load data
# ============================================================================
total_reads <- 24228692

cov <- read.table(
  "anchor_coverage_full.bed",
  sep       = "\t",
  col.names = c("chr", "start", "end", "loop_class",
                "read_count", "bases_covered",
                "anchor_length", "fraction_covered")
) %>%
  mutate(
    CPM              = (read_count / total_reads) * 1e6,
    log2_CPM         = log2(CPM + 1),
    mean_depth       = read_count / anchor_length,
    mean_depth_pseudo = mean_depth + 1e-6,
    loop_class = factor(loop_class, levels = c(
      "REH_specific", "EP1_specific", "ARO_specific",
      "shared_EP1_REH", "shared_ARO_REH", "shared_EP1_ARO", "shared_all3"
    ))
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
  "shared_all3"    = "#34495e"
)

# ============================================================================
# Wilcoxon pairwise — BH correction — keep only significant
# ============================================================================
wilcox_results <- cov %>%
  wilcox_test(log2_CPM ~ loop_class) %>%
  adjust_pvalue(method = "BH") %>%
  add_significance()

sig_comparisons <- wilcox_results %>%
  filter(p.adj < 0.05) %>%
  select(group1, group2, p.adj.signif)

# Stack brackets above the boxes
y_max   <- max(cov$log2_CPM, na.rm = TRUE)
y_start <- y_max * 0.75
sig_comparisons <- sig_comparisons %>%
  mutate(y.position = seq(y_start,
                          y_start + (y_max * 0.22),
                          length.out = n()))

# ============================================================================
# PLOT A — Boxplot: log2 CPM with significant brackets only
# ============================================================================
p_box <- ggplot(cov, aes(x = loop_class, y = log2_CPM, fill = loop_class)) +
  geom_boxplot(
    outlier.size  = 0.3,
    outlier.alpha = 0.4,
    linewidth     = 0.4,
    width         = 0.65
  ) +
  geom_hline(
    yintercept = median(cov$log2_CPM),
    linetype   = "dashed",
    color      = "red",
    linewidth  = 0.5
  ) +
  stat_pvalue_manual(
    sig_comparisons,
    label        = "p.adj.signif",
    tip.length   = 0.01,
    bracket.size = 0.4,
    size         = 4
  ) +
  scale_fill_manual(values = class_colors) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.18))) +
  theme_bw() +
  theme(
    panel.grid      = element_blank(),
    axis.text.x     = element_text(angle = 35, hjust = 1,
                                   size = 10, face = "bold"),
    axis.text.y     = element_text(size = 10),
    axis.title      = element_text(size = 12),
    plot.title      = element_text(hjust = 0.5, face = "bold", size = 13),
    plot.subtitle   = element_text(hjust = 0.5, size = 9, color = "gray40"),
    legend.position = "none"
  ) +
  xlab("Loop type") +
  ylab("Coverage (log2 CPM)") +
  labs(
    title    = "E2A-PBX1 ChIP-seq counts at loop anchors",
    subtitle = "Wilcoxon, BH-adjusted — significant comparisons only"
  )

# ============================================================================
# PLOT B — Violin + boxplot: mean depth log10
# ============================================================================
p_violin <- ggplot(cov, aes(x = loop_class, y = mean_depth_pseudo,
                            fill = loop_class)) +
  geom_violin(
    alpha     = 0.65,
    scale     = "width",
    trim      = FALSE,
    linewidth = 0.3
  ) +
  geom_boxplot(
    width         = 0.12,
    outlier.size  = 0.3,
    outlier.alpha = 0.4,
    linewidth     = 0.4,
    fill          = "white",
    color         = "black"
  ) +
  geom_hline(
    yintercept = median(cov$mean_depth_pseudo),
    linetype   = "dashed",
    color      = "red",
    linewidth  = 0.5
  ) +
  scale_y_log10(
    breaks = c(0.001, 0.01, 0.1, 1),
    labels = c("0.001", "0.010", "0.100", "1.000")
  ) +
  scale_fill_manual(values = class_colors) +
  theme_bw() +
  theme(
    panel.grid      = element_blank(),
    axis.text.x     = element_text(angle = 35, hjust = 1,
                                   size = 10, face = "bold"),
    axis.text.y     = element_text(size = 10),
    axis.title      = element_text(size = 12),
    plot.title      = element_text(hjust = 0.5, face = "bold", size = 13),
    legend.position = "none"
  ) +
  xlab("Loop type") +
  ylab("Coverage (log10 scale)") +
  labs(title = "ChIP-seq Coverage Across Loop Anchor Types")

# ============================================================================
# COMBINED
# ============================================================================
p_combined <- plot_grid(p_box, p_violin,
                        ncol       = 2,
                        labels     = c("A", "B"),
                        rel_widths = c(1, 1))

print(p_combined)

ggsave("EP1_chip_loop_counts_boxplot.pdf",  p_box,      width = 8,  height = 7)
ggsave("EP1_chip_loop_coverage_violin.pdf", p_violin,   width = 8,  height = 6)
ggsave("EP1_chip_loop_combined.pdf",        p_combined, width = 16, height = 7)