# =============================================================================
# figure1c_sample_overlap_upset.R
#
# Purpose:
#   UpSet plot showing the overlap in sample availability across the three
#   main analysis stages (SD-ASM identification, cis-eQTL mapping, phenotypic
#   association) (Figure 1C).
#
# Input:
#   - intersect.txt : per-sample binary indicator table with columns for
#                      WGBS, WGS, RNA-seq ("gene expression"), and clinical
#                      trait ("Clinical traits") data availability
#
# Output:
#   - Upset_plot.png / Upset_plot.pdf
# =============================================================================

library(UpSetR)
library(ComplexHeatmap)

## Set this to the directory containing intersect.txt
# setwd("/path/to/working_dir")

dat <- read.table("intersect.txt", header = TRUE, check.names = FALSE)[c(-1, -2, -3)]
dat <- dat[rowSums(dat[-2]) >= 1, ][-2]
names(dat)[c(1, 4)] <- c("WGS", "Clinical traits")

## Two additional WGBS samples were sequenced but excluded from WGS-matched
## analyses; flag them as WGBS-available only where WGS data is absent
dat$WGBS[1:4] <- 1
dat[dat$WGBS == 1 & dat$WGS == 0, "WGBS"] <- 0

## Restrict to the three sample-overlap combinations shown in Figure 1C:
## SD-ASM identification (WGBS+WGS), eQTL mapping (+ gene expression),
## phenotypic association (+ clinical traits)
m <- make_comb_mat(dat, mode = "intersect")
m <- m[c("1100", "1010", "1001")]

custom_labels <- c("SD-ASM\nidentification", "eQTL\nmapping", "Phenotypic\n association")

top_anno <- HeatmapAnnotation(
  "Intersection Size" = anno_barplot(
    comb_size(m), border = FALSE,
    gp = gpar(fill = c("#E69F00", "#56B4E9", "#009E73"),
              col  = c("#E69F00", "#56B4E9", "#009E73")),
    add_numbers = TRUE, numbers_rot = 0, height = unit(9, "cm")
  ),
  "Custom Labels" = anno_text(custom_labels, just = "top", rot = 0,
                               gp = gpar(fontsize = 10))
)

right_anno <- upset_right_annotation(m, width = unit(11, "cm"), add_numbers = TRUE)

upset_plot <- function() {
  UpSet(m, comb_order = c(1, 2, 3), set_order = c(1, 2, 3, 4),
        pt_size = unit(5, "mm"), lwd = 3,
        comb_col = c("#E69F00", "#56B4E9", "#009E73"),
        top_annotation = top_anno,
        right_annotation = right_anno,
        left_annotation = NULL,
        show_row_names = TRUE)
}

png("Upset_plot.png", width = 2000, height = 1200, res = 200)
upset_plot()
dev.off()

pdf("Upset_plot.pdf", width = 9, height = 6)
upset_plot()
dev.off()
