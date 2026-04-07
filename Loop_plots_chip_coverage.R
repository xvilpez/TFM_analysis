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
# Shared level order & colors
# ============================================================================
loop_class_levels <- c("REH_specific", "EP1_specific", "ARO_specific",
                       "shared_EP1_REH", "shared_ARO_REH", "shared_EP1_ARO",
                       "shared_all3")

level_order_full <- c(loop_class_levels, "random_30kb")

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
# Compute metrics — loop anchors only (for pairwise plot)
# ============================================================================
cov_loops <- cov %>%
  mutate(
    CPM               = (read_count / total_reads) * 1e6,
    log2_CPM          = log2(CPM + 1),
    mean_depth        = read_count / anchor_length,
    mean_depth_pseudo = mean_depth + 1e-6,
    loop_class        = factor(loop_class, levels = loop_class_levels)
  )

# ============================================================================
# Compute metrics — all classes including random (for vs-random plot)
# ============================================================================
all_cov <- bind_rows(cov, rnd) %>%
  mutate(
    CPM               = (read_count / total_reads) * 1e6,
    log2_CPM          = log2(CPM + 1),
    mean_depth        = read_count / anchor_length,
    mean_depth_pseudo = mean_depth + 1e-6,
    loop_class        = factor(loop_class, levels = level_order_full)
  )

# ============================================================================
# Helper: build bracket segments + labels from sig_ann data frame
# ============================================================================
make_brackets <- function(sig_ann, tick_drop = 0.05, label_lift = 0.07) {
  segments <- bind_rows(
    sig_ann %>% transmute(x = x,    xend = xend, y = y.position,              yend = y.position),
    sig_ann %>% transmute(x = x,    xend = x,    y = y.position - tick_drop,  yend = y.position),
    sig_ann %>% transmute(x = xend, xend = xend, y = y.position - tick_drop,  yend = y.position)
  )
  labels <- sig_ann %>%
    transmute(x = xmid, y = y.position + label_lift, label = p.adj.signif)
  list(segments = segments, labels = labels)
}

# ============================================================================
# PART 1 — Pairwise Wilcoxon between loop classes (no random)
# ============================================================================
wilcox_pairwise <- cov_loops %>%
  wilcox_test(log2_CPM ~ loop_class) %>%
  adjust_pvalue(method = "BH") %>%
  add_significance()

sig_pairwise <- wilcox_pairwise %>%
  filter(p.adj < 0.05) %>%
  select(group1, group2, p.adj.signif)

x_pos_loops <- setNames(seq_along(loop_class_levels), loop_class_levels)

y_max_loops   <- max(cov_loops$log2_CPM, na.rm = TRUE)
y_start_loops <- y_max_loops * 0.75

sig_pairwise <- sig_pairwise %>%
  mutate(
    y.position = seq(y_start_loops,
                     y_start_loops + (y_max_loops * 0.22),
                     length.out = n()),
    x    = x_pos_loops[group1],
    xend = x_pos_loops[group2],
    xmid = (x_pos_loops[group1] + x_pos_loops[group2]) / 2
  )

brk_pairwise <- make_brackets(sig_pairwise)

# --- Plot A1: boxplot, pairwise ---
p_box_pairwise <- ggplot(cov_loops, aes(x = loop_class, y = log2_CPM, fill = loop_class)) +
  geom_boxplot(outlier.size = 0.3, outlier.alpha = 0.4, size = 0.4, width = 0.65) +
  geom_hline(
    yintercept = median(cov_loops$log2_CPM, na.rm = TRUE),
    linetype = "dashed", color = "red", size = 0.5
  ) +
  geom_segment(
    data = brk_pairwise$segments,
    aes(x = x, xend = xend, y = y, yend = yend),
    inherit.aes = FALSE, size = 0.4
  ) +
  geom_text(
    data = brk_pairwise$labels,
    aes(x = x, y = y, label = label),
    inherit.aes = FALSE, size = 4
  ) +
  scale_fill_manual(values = class_colors) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.18))) +
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
    subtitle = "Pairwise Wilcoxon, BH-adjusted — significant comparisons only"
  )

# --- Plot B1: violin + boxplot, pairwise ---
p_violin_pairwise <- ggplot(cov_loops, aes(x = loop_class, y = mean_depth_pseudo, fill = loop_class)) +
  geom_violin(alpha = 0.65, scale = "width", trim = FALSE, size = 0.3) +
  geom_boxplot(width = 0.12, outlier.size = 0.3, outlier.alpha = 0.4,
               size = 0.4, fill = "white", color = "black") +
  geom_hline(
    yintercept = median(cov_loops$mean_depth_pseudo, na.rm = TRUE),
    linetype = "dashed", color = "red", size = 0.5
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
# PART 2 — Wilcoxon: each loop class vs random_30kb, BH correction
# ============================================================================
wilcox_vs_random <- lapply(loop_class_levels, function(cls) {
  sub_df <- all_cov %>%
    filter(loop_class %in% c(cls, "random_30kb")) %>%
    mutate(loop_class = droplevels(loop_class))
  
  wilcox_test(sub_df, log2_CPM ~ loop_class, ref.group = "random_30kb") %>%
    mutate(group1 = cls, group2 = "random_30kb")
}) %>%
  bind_rows() %>%
  adjust_pvalue(method = "BH") %>%
  add_significance()

sig_vs_random <- wilcox_vs_random %>%
  filter(p.adj < 0.05) %>%
  select(group1, group2, p.adj.signif)

x_pos_full <- setNames(seq_along(level_order_full), level_order_full)

y_max_full   <- max(all_cov$log2_CPM, na.rm = TRUE)
y_start_full <- y_max_full * 0.80

sig_vs_random <- sig_vs_random %>%
  mutate(
    y.position = seq(y_start_full,
                     y_start_full + (y_max_full * 0.18),
                     length.out = n()),
    x    = x_pos_full[group1],
    xend = x_pos_full[group2],
    xmid = (x_pos_full[group1] + x_pos_full[group2]) / 2
  )

brk_random <- make_brackets(sig_vs_random)

# --- Plot A2: boxplot, vs random ---
p_box_random <- ggplot(all_cov, aes(x = loop_class, y = log2_CPM, fill = loop_class)) +
  geom_boxplot(outlier.size = 0.3, outlier.alpha = 0.4, size = 0.4, width = 0.65) +
  geom_hline(
    yintercept = median(all_cov$log2_CPM[all_cov$loop_class == "random_30kb"], na.rm = TRUE),
    linetype = "dashed", color = "gray40", size = 0.5
  ) +
  geom_segment(
    data = brk_random$segments,
    aes(x = x, xend = xend, y = y, yend = yend),
    inherit.aes = FALSE, size = 0.4
  ) +
  geom_text(
    data = brk_random$labels,
    aes(x = x, y = y, label = label),
    inherit.aes = FALSE, size = 4
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

# --- Plot B2: violin + boxplot, vs random ---
p_violin_random <- ggplot(all_cov, aes(x = loop_class, y = mean_depth_pseudo, fill = loop_class)) +
  geom_violin(alpha = 0.65, scale = "width", trim = FALSE, size = 0.3) +
  geom_boxplot(width = 0.12, outlier.size = 0.3, outlier.alpha = 0.4,
               size = 0.4, fill = "white", color = "black") +
  geom_hline(
    yintercept = median(all_cov$mean_depth_pseudo[all_cov$loop_class == "random_30kb"], na.rm = TRUE),
    linetype = "dashed", color = "gray40", size = 0.5
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
  labs(title = "ChIP-seq Coverage Across Loop Anchor Types (+ random)")

# ============================================================================
# SAVE — pairwise (loop classes only)
# ============================================================================
pdf("EP1_chip_loop_counts_boxplot.pdf", width = 8, height = 7)
print(p_box_pairwise)
dev.off()

pdf("EP1_chip_loop_coverage_violin.pdf", width = 8, height = 6)
print(p_violin_pairwise)
dev.off()

p_combined_pairwise <- p_box_pairwise + p_violin_pairwise +
  plot_annotation(tag_levels = "A")

pdf("EP1_chip_loop_combined.pdf", width = 16, height = 7)
print(p_combined_pairwise)
dev.off()

# ============================================================================
# SAVE — vs random
# ============================================================================
pdf("EP1_chip_loop+random_counts_boxplot.pdf", width = 9, height = 7)
print(p_box_random)
dev.off()

pdf("EP1_chip_loop+random_coverage_violin.pdf", width = 9, height = 6)
print(p_violin_random)
dev.off()

p_combined_random <- p_box_random + p_violin_random +
  plot_annotation(tag_levels = "A")

pdf("EP1_chip_loop+random_combined.pdf", width = 18, height = 7)
print(p_combined_random)
dev.off()