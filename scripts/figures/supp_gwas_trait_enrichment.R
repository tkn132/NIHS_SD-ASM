# =============================================================================
# supp_gwas_trait_enrichment.R
#
# Purpose:
#   Test whether SD-ASM loci with suggestive GWAS trait associations (p<0.05)
#   are over-represented among specific phenotype categories, using a
#   hypergeometric enrichment test against the total number of trait
#   associations tested across all phenotypes.
#
# Input:
#   - GWASenrich.csv : columns "n" (observed SD-ASM associations per category),
#                      "total" (total associations tested per category)
#
# Output:
#   - dat (printed): per-category enrichment p-values and BH-adjusted q-values
# =============================================================================

dat <- read.csv("GWASenrich.csv", header = TRUE)

background_total <- sum(dat$total)

dat$enrichment_pvalue <- apply(dat, 1, function(row) {
  observed <- as.numeric(row["n"])
  total_in_category <- as.numeric(row["total"])
  phyper(observed - 1, total_in_category, background_total - total_in_category,
         sum(dat$n), lower.tail = FALSE)
})

dat$adjusted_pvalue <- p.adjust(dat$enrichment_pvalue, method = "BH")

print(head(dat[order(dat$adjusted_pvalue), ], 10))
