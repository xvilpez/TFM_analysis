library(VennDiagram)
library(grid)

#############################################
# Utility: count loops
#############################################
count_loops <- function(file) {
  if (!file.exists(file)) {
    warning(paste("File not found:", file))
    return(0)
  }
  nrow(read.table(file, sep = "\t", stringsAsFactors = FALSE))
}

#############################################
# Three-way Venn from classified loop files
#############################################
plot_threeway_venn <- function(
  shared_all3_file,
  only_A_file,
  only_B_file,
  only_C_file,
  shared_AB_file,
  shared_AC_file,
  shared_BC_file,
  label_A,
  label_B,
  label_C,
  output_pdf,
  title_suffix = "±30 kb"
) {
  
  # Count loops in each region
  only_A <- count_loops(only_A_file)
  only_B <- count_loops(only_B_file)
  only_C <- count_loops(only_C_file)
  shared_AB <- count_loops(shared_AB_file)
  shared_AC <- count_loops(shared_AC_file)
  shared_BC <- count_loops(shared_BC_file)
  shared_all3 <- count_loops(shared_all3_file)
  
  # Calculate total areas
  area_A <- only_A + shared_AB + shared_AC + shared_all3
  area_B <- only_B + shared_AB + shared_BC + shared_all3
  area_C <- only_C + shared_AC + shared_BC + shared_all3
  
  # Calculate pairwise intersections (including the triple)
  n12 <- shared_AB + shared_all3  # A ∩ B
  n13 <- shared_AC + shared_all3  # A ∩ C
  n23 <- shared_BC + shared_all3  # B ∩ C
  n123 <- shared_all3             # A ∩ B ∩ C
  
  cat("\n=== Three-way comparison:", label_A, "vs", label_B, "vs", label_C, "===\n")
  cat(label_A, "only:", only_A, "\n")
  cat(label_B, "only:", only_B, "\n")
  cat(label_C, "only:", only_C, "\n")
  cat("Shared", label_A, "&", label_B, "only:", shared_AB, "\n")
  cat("Shared", label_A, "&", label_C, "only:", shared_AC, "\n")
  cat("Shared", label_B, "&", label_C, "only:", shared_BC, "\n")
  cat("Shared all three:", shared_all3, "\n")
  cat("---\n")
  cat("Total", label_A, ":", area_A, "\n")
  cat("Total", label_B, ":", area_B, "\n")
  cat("Total", label_C, ":", area_C, "\n")
  cat("===\n")
  
  venn.plot <- draw.triple.venn(
    area1 = area_A,
    area2 = area_B,
    area3 = area_C,
    n12 = n12,
    n23 = n23,
    n13 = n13,
    n123 = n123,
    category = c(label_A, label_B, label_C),
    fill = c("lightblue", "lightgreen", "lightyellow"),
    cat.col = c("darkblue", "darkgreen", "darkorange"),
    cex = 1.5,
    cat.cex = 1.5,
    cat.dist = c(0.05, 0.05, 0.05),
    main = paste("Loop comparison:", label_A, "vs", label_B, "vs", label_C, title_suffix)
  )
  
  pdf(output_pdf, width = 8, height = 8)
  grid.draw(venn.plot)
  dev.off()
  
  cat("Saved:", output_pdf, "\n\n")
}

#############################################
# Set working directory and plot
#############################################
setwd("/ijc/LABS/STIK/RAW/Hi-C/HiC_marta_Oct25/loops/intersected_loops")

plot_threeway_venn(
  shared_all3_file = "shared_loops_all3_05_10kb_exp30kb.bed",
  only_A_file = "EP1_specific_only_05_10kb_exp30kb.bed",
  only_B_file = "Aro_specific_only_05_10kb_exp30kb.bed",
  only_C_file = "REH_specific_only_05_10kb_exp30kb.bed",
  shared_AB_file = "shared_EP1_Aro_only_05_10kb_exp30kb.bed",
  shared_AC_file = "shared_EP1_REH_only_05_10kb_exp30kb.bed",
  shared_BC_file = "shared_Aro_REH_only_05_10kb_exp30kb.bed",
  label_A = "EP1",
  label_B = "Aro",
  label_C = "REH",
  output_pdf = "venn_EP1_Aro_REH_threeway_exp30kb.pdf",
  title_suffix = "±30 kb"
)