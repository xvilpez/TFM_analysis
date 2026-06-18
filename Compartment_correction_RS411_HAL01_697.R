# ============================================================
#  A/B Chromatin Compartment Analysis – All Samples
#  Cscore bedgraph files: REH replicates + new cell lines
# ============================================================
#
#  Old samples: loaded from pre-corrected bedgraph files
#  New samples: raw bedgraphs with sample-specific corrections applied here
#
#  New sample corrections:
#    s697  – none identified yet
#    HAL01 – none identified yet
#    RS411 – flip chr13 + remove chr4:56,932,708-67,727,403
#
#  NOTE: R variable names cannot start with a digit.
#        The 697 sample is stored as  s697  throughout.
# ============================================================

# ── 0. Libraries ─────────────────────────────────────────────
library(dplyr)
library(ggplot2)
library(reshape2)
library(ggrepel)
library(pheatmap)
library(viridis)
library(tidyr)
library(grid)
library(cowplot)

# ── 1. Paths & output directories ────────────────────────────
data_dir   <- "R:/R_STIK/Hi-C/HiC_marta_Oct25/cscore"
output_dir <- "Comp_Analysis"
files_dir  <- file.path(output_dir, "files")
plots_dir  <- file.path(output_dir, "plots")

for (d in c(files_dir, plots_dir)) {
  if (!dir.exists(d)) dir.create(d, recursive = TRUE)
}

# Helper: significance label
sig_label <- function(p) {
  case_when(
    is.na(p) | p >= 0.05 ~ NA_character_,
    p < 0.001 ~ paste0("p=", formatC(p, format = "e", digits = 1)),
    TRUE ~ paste0("p=", round(p, 3))
  )
}

# ── 2. Data loading ───────────────────────────────────────────
read_cscore <- function(filename, score_colname) {
  read.table(
    file.path(data_dir, filename),
    sep       = "\t",
    col.names = c("chr", "start", "end", score_colname)
  )
}

# ── 2.1 Old samples – already corrected ──────────────────────
REH           <- read_cscore("REH_manual_corrected.bedgraph",           "score_REH")
EP1_1         <- read_cscore("EP1_1_manual_corrected.bedgraph",         "score_EP1_1")
EP1_2         <- read_cscore("EP1_2_manual_corrected.bedgraph",         "score_EP1_2")
EP1_Aro_1     <- read_cscore("EP1_Aro_1_manual_corrected.bedgraph",     "score_EP1_Aro_1")
EP1_Aro_2     <- read_cscore("EP1_Aro_2_manual_corrected.bedgraph",     "score_EP1_Aro_2")
EP1_merge     <- read_cscore("EP1_merge_manual_corrected.bedgraph",     "score_EP1_merge")
EP1_Aro_merge <- read_cscore("EP1_Aro_merge_manual_corrected.bedgraph", "score_EP1_Aro_merge")

# ── 2.2 New samples – raw bedgraphs ──────────────────────────
s697  <- read_cscore("697_tp.bedgraph_correct.bedgraph",   "score_697")
HAL01 <- read_cscore("HALO1_tp.bedgraph_correct.bedgraph", "score_HAL01")
RS411 <- read_cscore("RS411_tp.bedgraph_correct.bedgraph", "score_RS411")

# ── 3. Correct new samples ────────────────────────────────────

# s697: no corrections needed yet
s697_corr <- s697

# HAL01: no corrections needed yet
HAL01_corr <- HAL01

# RS411: flip chr13 + remove chr4:56,932,708-67,727,403
RS411_corr <- RS411

# flip chr13
RS411_corr[RS411_corr$chr == "chr13", "score_RS411"] <-
  -RS411_corr[RS411_corr$chr == "chr13", "score_RS411"]

# remove problematic chr4 region
RS411_corr <- RS411_corr[!(RS411_corr$chr == "chr4" &
                             RS411_corr$end   >= 56932708 &
                             RS411_corr$start <= 67727403), ]

cat("RS411 bins before correction:", nrow(RS411),      "\n")
cat("RS411 bins after  correction:", nrow(RS411_corr), "\n")

# Optionally save corrected new-sample bedgraphs
write.table(s697_corr, "697_manual_corrected.bedgraph",
            sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
write.table(HAL01_corr, "HAL01_manual_corrected.bedgraph",
            sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
write.table(RS411_corr, file.path(files_dir, "RS411_manual_corrected.bedgraph"),
            sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)

# ── 4. Merge all samples ──────────────────────────────────────
# Old samples are already corrected; new samples use their corrected versions
cs_all <- Reduce(
  function(x, y) merge(x, y, by = c("chr", "start", "end")),
  list(REH, EP1_1, EP1_2, EP1_Aro_1, EP1_Aro_2,
       EP1_merge, EP1_Aro_merge,
       s697_corr, HAL01_corr, RS411_corr)
)

cat("Columns:", paste(colnames(cs_all), collapse = ", "), "\n")
cat("Total bins after merge:", nrow(cs_all), "\n")
head(cs_all)

# ── 5. PCA ────────────────────────────────────────────────────
# Sample metadata — order matches columns 4:13 of cs_all
sample_labels <- c("REH", "EP1_1", "EP1_2", "EP1_Aro_1", "EP1_Aro_2",
                   "EP1_merge", "EP1_Aro_merge",
                   "697", "HAL01", "RS411")

conditions <- c("REH", "EP1", "EP1", "EP1_Aro", "EP1_Aro",
                "EP1", "EP1_Aro",
                "697", "HAL01", "RS411")

types <- c("original",  "replicate", "replicate",
           "replicate", "replicate", "merge",     "merge",
           "original",  "original",  "original")

cond_colors <- c(
  "REH"     = "#e74c3c",
  "EP1"     = "#3498db",
  "EP1_Aro" = "#2ecc71",
  "697"     = "#9b59b6",
  "HAL01"   = "#f39c12",
  "RS411"   = "#1abc9c"
)

type_shapes <- c(
  "original"  = 21,
  "replicate" = 22,
  "merge"     = 24
)

cscore_matrix           <- as.matrix(cs_all[, 4:13])
colnames(cscore_matrix) <- sample_labels
cscore_matrix           <- na.omit(cscore_matrix)
cscore_matrix           <- cscore_matrix[apply(cscore_matrix, 1, var) != 0, ]

cat("Bins used for PCA:", nrow(cscore_matrix), "\n")

PCA_all <- prcomp(t(cscore_matrix), scale. = TRUE)
pct_var <- round(100 * PCA_all$sdev^2 / sum(PCA_all$sdev^2), 1)

pca_df           <- as.data.frame(PCA_all$x)
pca_df$sample    <- sample_labels
pca_df$condition <- conditions
pca_df$type      <- types

pca_plot <- ggplot(pca_df,
                   aes(x = PC1, y = PC2,
                       fill  = condition,
                       shape = type,
                       label = sample)) +
  geom_point(color = "black", stroke = 1.3, size = 5) +
  geom_text_repel(size = 4, max.overlaps = 20) +
  scale_fill_manual(values = cond_colors) +
  scale_shape_manual(values = type_shapes) +
  theme_bw() +
  labs(
    title = "Cscore PCA – all samples (corrected)",
    x     = paste0("PC1 (", pct_var[1], "% variance)"),
    y     = paste0("PC2 (", pct_var[2], "% variance)")
  ) +
  theme(
    panel.grid   = element_blank(),
    plot.title   = element_text(size = 15, face = "bold", hjust = 0.5),
    axis.title   = element_text(size = 13),
    axis.text    = element_text(size = 11),
    legend.title = element_text(size = 11),
    legend.text  = element_text(size = 10)
  ) +
  guides(
    fill  = guide_legend(title = "Condition", override.aes = list(shape = 21)),
    shape = guide_legend(title = "Type")
  )

print(pca_plot)
ggsave(file.path(plots_dir, "pca_all_samples_corrected.pdf"),
       pca_plot, width = 10, height = 6)