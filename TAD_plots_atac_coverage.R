library(dplyr)
library(ggplot2)
library(ggpubr)
library(rstatix)
library(tidyr)
library(patchwork)

# ============================================================================
# Load data
# ============================================================================
total_reads <- as.numeric(readLines(
  "R:/R_STIK/Hi-C/HiC_marta_Oct25/TADs/comparison/atac_coverage/total_mapped_reads.txt"
))

setwd("R:/R_STIK/Hi-C/HiC_marta_Oct25/TADs/comparison/atac_coverage")

col_names <- c("chr", "start", "end", "boundary_class",
               "read_count", "bases_covered",
               "boundary_length", "fraction_covered")

cov <- read.table("boundary_coverage_full.bed",
                  sep = "\t", col.names = col_names)

# ============================================================================
# Level order & colors
# ============================================================================
boundary_class_levels <- c("shared", "unique_HiC_REH_EP1", "unique_REH")

class_colors <- c(
  "shared"            = "#3498db",
  "unique_HiC_REH_EP1" = "#9b59b6",
  "unique_REH"        = "#e74c3c"
)

# ============================================================================
# Compute metrics
# ============================================================================
cov_boundaries <- cov %>%
  mutate(
    CPM               = (read_count / total_reads) * 1e6,
    log2_CPM          = log2(CPM + 1),
    mean_depth        = read_count / boundary_length,
    mean_depth_pseudo = mean_depth + 1e-6,
    boundary_class    = factor(boundary_class, levels = boundary_class_levels)
  )

# ============================================================================
# Helper: build bracket segments + labels from sig_ann data frame
# ============================================================================
make_brackets <- function(sig_ann, tick_drop = 0.05, label_lift = 0.07) {
  segments <- bind_rows(
    sig_ann %>% transmute(x = x,    xend = xend, y = y.position,             yend = y.position),
    sig_ann %>% transmute(x = x,    xend = x,    y = y.position - tick_drop, yend = y.position),
    sig_ann %>% transmute(x = xend, xend = xend, y = y.position - tick_drop, yend = y.position)
  )
  labels <- sig_ann %>%
    transmute(x = xmid, y = y.position + label_lift, label = p.adj.signif)
  list(segments = segments, labels = labels)
}

# ============================================================================
# Pairwise Wilcoxon between boundary classes, BH correction
# ============================================================================
wilcox_pairwise <- cov_boundaries %>%
  wilcox_test(log2_CPM ~ boundary_class) %>%
  adjust_pvalue(method = "BH") %>%
  add_significance()

sig_pairwise <- wilcox_pairwise %>%
  filter(p.adj < 0.05) %>%
  select(group1, group2, p.adj.signif)

x_pos <- setNames(seq_along(boundary_class_levels), boundary_class_levels)

y_max   <- max(cov_boundaries$log2_CPM, na.rm = TRUE)
y_start <- y_max * 0.75

if (nrow(sig_pairwise) > 0) {
  sig_pairwise <- sig_pairwise %>%
    mutate(
      y.position = seq(y_start,
                       y_start + (y_max * 0.22),
                       length.out = n()),
      x    = x_pos[group1],
      xend = x_pos[group2],
      xmid = (x_pos[group1] + x_pos[group2]) / 2
    )
  brk_pairwise <- make_brackets(sig_pairwise)
} else {
  brk_pairwise <- list(
    segments = data.frame(x = numeric(), xend = numeric(), y = numeric(), yend = numeric()),
    labels   = data.frame(x = numeric(), y = numeric(), label = character())
  )
}

# ============================================================================
# Plot A — boxplot
# ============================================================================
p_box <- ggplot(cov_boundaries, aes(x = boundary_class, y = log2_CPM, fill = boundary_class)) +
  geom_boxplot(outlier.size = 0.3, outlier.alpha = 0.4, size = 0.4, width = 0.65) +
  geom_hline(
    yintercept = median(cov_boundaries$log2_CPM, na.rm = TRUE),
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
  scale_x_discrete(labels = c(
    "shared"             = "Shared",
    "unique_HiC_REH_EP1" = "Unique HiC-REH-EP1",
    "unique_REH"         = "Unique REH"
  )) +
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
  xlab("Boundary type") + ylab("Coverage (log2 CPM)") +
  labs(
    title    = "REH ATAC-seq signal at TAD boundaries",
    subtitle = "Pairwise Wilcoxon, BH-adjusted — significant comparisons only"
  )

# ============================================================================
# Plot B — violin + boxplot
# ============================================================================
p_violin <- ggplot(cov_boundaries, aes(x = boundary_class, y = mean_depth_pseudo, fill = boundary_class)) +
  geom_violin(alpha = 0.65, scale = "width", trim = FALSE, size = 0.3) +
  geom_boxplot(width = 0.12, outlier.size = 0.3, outlier.alpha = 0.4,
               size = 0.4, fill = "white", color = "black") +
  geom_hline(
    yintercept = median(cov_boundaries$mean_depth_pseudo, na.rm = TRUE),
    linetype = "dashed", color = "red", size = 0.5
  ) +
  scale_y_log10(
    breaks = c(0.001, 0.01, 0.1, 1),
    labels = c("0.001", "0.010", "0.100", "1.000")
  ) +
  scale_fill_manual(values = class_colors) +
  scale_x_discrete(labels = c(
    "shared"             = "Shared",
    "unique_HiC_REH_EP1" = "Unique HiC-REH-EP1",
    "unique_REH"         = "Unique REH"
  )) +
  theme_bw() +
  theme(
    panel.grid      = element_blank(),
    axis.text.x     = element_text(angle = 35, hjust = 1, size = 10, face = "bold"),
    axis.text.y     = element_text(size = 10),
    axis.title      = element_text(size = 12),
    plot.title      = element_text(hjust = 0.5, face = "bold", size = 13),
    legend.position = "none"
  ) +
  xlab("Boundary type") + ylab("Mean depth (log10 scale)") +
  labs(title = "REH ATAC-seq coverage at TAD boundaries")

# ============================================================================
# Plot C — fraction of bases covered (chromatin accessibility breadth)
# ============================================================================
p_fraction <- ggplot(cov_boundaries, aes(x = boundary_class, y = fraction_covered, fill = boundary_class)) +
  geom_boxplot(outlier.size = 0.3, outlier.alpha = 0.4, size = 0.4, width = 0.65) +
  geom_hline(
    yintercept = median(cov_boundaries$fraction_covered, na.rm = TRUE),
    linetype = "dashed", color = "red", size = 0.5
  ) +
  scale_fill_manual(values = class_colors) +
  scale_x_discrete(labels = c(
    "shared"             = "Shared",
    "unique_HiC_REH_EP1" = "Unique HiC-REH-EP1",
    "unique_REH"         = "Unique REH"
  )) +
  scale_y_continuous(limits = c(0, 1), expand = expansion(mult = c(0.02, 0.05))) +
  theme_bw() +
  theme(
    panel.grid      = element_blank(),
    axis.text.x     = element_text(angle = 35, hjust = 1, size = 10, face = "bold"),
    axis.text.y     = element_text(size = 10),
    axis.title      = element_text(size = 12),
    plot.title      = element_text(hjust = 0.5, face = "bold", size = 13),
    legend.position = "none"
  ) +
  xlab("Boundary type") + ylab("Fraction of bases covered") +
  labs(title = "ATAC-seq breadth of coverage at TAD boundaries")

# ============================================================================
# Save individual plots
# ============================================================================
pdf("REH_atac_boundary_boxplot.pdf", width = 6, height = 6)
print(p_box)
dev.off()

pdf("REH_atac_boundary_violin.pdf", width = 6, height = 6)
print(p_violin)
dev.off()

pdf("REH_atac_boundary_fraction.pdf", width = 6, height = 6)
print(p_fraction)
dev.off()

# ============================================================================
# Save combined panel (A | B | C)
# ============================================================================
p_combined <- p_box + p_violin + p_fraction +
  plot_annotation(tag_levels = "A")

pdf("REH_atac_boundary_combined.pdf", width = 18, height = 6)
print(p_combined)
dev.off()

# ============================================================================
# Per-class summary table
# ============================================================================
summary_tbl <- cov_boundaries %>%
  group_by(boundary_class) %>%
  summarise(
    n_boundaries      = n(),
    median_log2_CPM   = median(log2_CPM, na.rm = TRUE),
    mean_log2_CPM     = mean(log2_CPM, na.rm = TRUE),
    median_depth      = median(mean_depth, na.rm = TRUE),
    median_frac_cov   = median(fraction_covered, na.rm = TRUE),
    .groups = "drop"
  )

print(summary_tbl)
write.table(summary_tbl, "per_class_summary_R.tsv",
            sep = "\t", quote = FALSE, row.names = FALSE)

# Print Wilcoxon results
cat("\n--- Wilcoxon pairwise results ---\n")
print(wilcox_pairwise %>% select(group1, group2, p, p.adj, p.adj.signif))