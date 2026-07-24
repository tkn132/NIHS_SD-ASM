# =============================================================================
# figure2a_volcano_plot.R
#
# Purpose:
#   Volcano plot of methylation difference vs. -log10(FDR) for all tested
#   heterozygous-SNP CpG windows, highlighting significant SD-ASM loci
#   (Figure 2A).
#
# Input:
#   - all.stats.txt : per-locus SD-ASM test statistics, with columns including
#                      Chr, BP, Methylation.Difference, FDR
#
# Output:
#   - Volcano.pdf / Volcano.png
# =============================================================================

library(ggplot2)
library(tidyverse)
library(ggsci)

# NB: for plotting speed/file size, a random subsample of 500,000 loci is
# drawn from the full test set; significance calls (dashed threshold lines)
# are unaffected since they are based on the fixed FDR/effect-size cutoffs.

dat <- read.table("all.stats.txt", header = TRUE)
names(dat)[1:2] <- c("Chr", "BP")

## Remove the MHC/HLA region (chr6:28,510,120-33,480,577, GRCh38)
## https://www.ncbi.nlm.nih.gov/grc/human/regions/MHC
dat <- dat[!(dat$Chr == "chr6" & dat$BP > 28510120 & dat$BP < 33480577), ]

## Flag SD-ASM significance: |methylation difference| >= 30% and FDR < 0.05
dat$col <- ifelse(abs(dat$Methylation.Difference) > 0.3 & dat$FDR < 0.05,
                   "significant", "non-significant")
table(dat$col)

dat <- dat[dat$FDR != 0, ]
set.seed(1)
sdat <- dat[sample(seq_len(nrow(dat)), 500000), ]

volcano_theme <- theme(
  plot.title = element_text(size = 10),
  axis.title = element_text(size = 7),
  text = element_text(size = 7),
  legend.text = element_text(size = 6),
  legend.title = element_text(size = 7),
  axis.line = element_line(colour = "black"),
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(),
  panel.border = element_blank(),
  panel.background = element_blank()
)

volcano_plot <- ggplot(sdat, aes(y = -log10(FDR), x = Methylation.Difference, color = col)) +
  geom_point() +
  geom_vline(xintercept = c(-0.3, 0.3), col = "gray", linetype = "dashed") +
  geom_hline(yintercept = -log10(0.05), col = "gray", linetype = "dashed") +
  scale_color_manual(values = c("grey", "#bb0c00")) +
  guides(col = "none") +
  volcano_theme

pdf("Volcano.pdf", width = 4, height = 3.7)
print(volcano_plot)
dev.off()

png("Volcano.png", width = 900, height = 700, res = 150)
print(volcano_plot)
dev.off()
