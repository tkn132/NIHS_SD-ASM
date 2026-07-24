# =============================================================================
# figure2bc_gene_feature_enrichment.R
#
# Purpose:
#   Compare the genomic-feature distribution of significant SD-ASM loci
#   against all tested heterozygous SNPs (Figure 2B: all features;
#   Figure 2C: excluding intergenic/intronic regions), and test each feature
#   for enrichment/depletion among SD-ASM loci using Fisher's exact test.
#
# Input:
#   - Gene_Feature_Count.txt : columns "feature", "sig" (SD-ASM loci), "all"
#                              (all tested heterozygous SNPs)
#
# Output:
#   - Gene_Feature_comparison.pdf
#   - Gene_Feature_comparison_nointron_nointergenic.pdf
#   - results_df (printed): per-feature Fisher's exact test results
# =============================================================================

library(ggplot2)
library(tidyr)
library(dplyr)
library(ggsci)

feature_theme <- theme(
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

gf <- read.table("Gene_Feature_Count.txt", header = TRUE)[1:3]

## Collapse minor/composite ncRNA categories and drop overlapping/ambiguous
## feature classes that are not reported individually in the figure
gf <- rbind(gf, c("ncRNA", sum(gf$sig[6:9]), sum(gf$all[6:9])))
gf <- gf[!(gf$feature %in% c(
  "exonic;splicing", "ncRNA_exonic;splicing", "upstream;downstream",
  "UTR5;UTR3", "ncRNA_intronic", "ncRNA_splicing", "ncRNA_exonic"
)), ]

## ---- Figure 2B: all genomic features --------------------------------------
gf_long <- pivot_longer(gf, sig:all, names_to = "Set", values_to = "Frequency") %>%
  as.data.frame() %>%
  group_by(Set) %>%
  mutate(Percentage = Frequency / sum(Frequency) * 100) %>%
  as.data.frame()

fig2b <- ggplot(gf_long, aes(x = Set, y = Percentage, fill = feature)) +
  geom_bar(stat = "identity") +
  scale_x_discrete(labels = c("Heterozygous SNPs", "Significant SD-ASM")) +
  labs(fill = "Gene Feature", x = "", title = "All features") +
  scale_fill_startrek() +
  feature_theme

pdf("Gene_Feature_comparison.pdf", width = 4, height = 3.7)
print(fig2b)
dev.off()

## ---- Figure 2C: excluding intergenic and intronic regions -----------------
gf_long2 <- gf_long[!(gf_long$feature %in% c("intergenic", "intronic")), ] %>%
  group_by(Set) %>%
  mutate(Percentage = Frequency / sum(Frequency) * 100) %>%
  as.data.frame()

fig2c <- ggplot(gf_long2, aes(x = Set, y = Percentage, fill = feature)) +
  geom_bar(stat = "identity") +
  scale_x_discrete(labels = c("Heterozygous SNPs", "Significant SD-ASM")) +
  labs(fill = "Gene Feature", x = "", title = "Without intergenic and intronic regions") +
  scale_fill_startrek() +
  feature_theme

pdf("Gene_Feature_comparison_nointron_nointergenic.pdf", width = 4, height = 3.7)
print(fig2c)
dev.off()

## ---- Feature-level enrichment test (Fisher's exact test) ------------------
total_sig <- sum(gf$sig)
total_non_sig <- sum(gf$all) - total_sig

results <- apply(gf, 1, function(row) {
  feature_sig <- as.numeric(row["sig"])
  feature_all <- as.numeric(row["all"])
  feature_non_sig <- feature_all - feature_sig

  contingency_table <- matrix(
    c(feature_sig, feature_non_sig, total_sig - feature_sig, total_non_sig - feature_non_sig),
    nrow = 2, byrow = TRUE
  )

  fisher_test <- fisher.test(contingency_table)
  c(feature = row["feature"], p_value = fisher_test$p.value,
    odds_ratio = unname(fisher_test$estimate))
})

results_df <- as.data.frame(t(results))
results_df$p_value <- as.numeric(results_df$p_value)
results_df$odds_ratio <- as.numeric(results_df$odds_ratio)

## Bonferroni threshold across features
bonferroni_threshold <- 0.05 / nrow(results_df)

print(results_df[order(results_df$odds_ratio), ])
cat("Bonferroni threshold:", bonferroni_threshold, "\n")
