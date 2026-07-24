# =============================================================================
# supp_maf_comparison.R
#
# Purpose:
#   Compare minor allele frequencies (MAF) of significant SD-ASM variants in
#   the NI cohort against gnomAD reference population frequencies, and inspect
#   allele-frequency patterns within specific eQTL hotspot regions
#   (e.g. the chr12 TAS2R43 cluster, Table 1).
#
# Inputs:
#   - MAF.all.3.csv    : gnomAD (v4.1 genome) allele frequencies for SD-ASM
#                          variants, columns include Chr, Start, Ref, Alt, AF, ...
#   - sigSDASM.afreq    : NI cohort allele frequencies (PLINK2 --freq output)
#   - sigSDASM.bim      : NI cohort variant positions (PLINK .bim)
#
# Output:
#   - MAF.all.2.csv : merged NI + gnomAD frequency table, NI MAF re-oriented
#                      to the gnomAD alt allele where needed
# =============================================================================

library(dplyr)
library(stringr)

## Set this to the directory containing the input files above
# setwd("/path/to/working_dir")

all <- read.csv("MAF.all.3.csv", header = TRUE, na.strings = ".")

ni <- read.table("sigSDASM.afreq")[-6]
bim <- read.table("sigSDASM.bim")
ni <- cbind(ni, bim$V4)
names(ni) <- c("Chr", "SNP", "A1", "A2", "NI", "Start")

## Merge NI and gnomAD frequency tables on chromosome + position
dat <- merge(all, ni, by = c("Chr", "Start"))

## Re-orient the NI allele frequency to match the gnomAD alt allele where the
## reference/alt alleles are flipped between the two datasets
sum(dat$Ref == dat$A2, na.rm = TRUE)
sum(dat$Ref == dat$A1, na.rm = TRUE)

dat$NI <- ifelse(dat$Ref == dat$A2, 1 - dat$NI, as.numeric(dat$NI))
summary(dat$NI)

write.csv(dat, "MAF.all.2.csv", row.names = FALSE)

names(dat) <- str_remove(names(dat), "gnomad41_genome_")

## ---- Inspect specific eQTL hotspot regions (Table 1) -----------------------
## chr22 DDT/MIF cluster
dat[dat$Chr == 22 & dat$Start >= 23943487 & dat$Start < 23944710, ]

## chr12 TAS2R43 cluster (rare-variant enriched, AF<0.05 in gnomAD)
dat[dat$Chr == 12 & dat$Start >= 11070000 & dat$Start < 11150000 & dat$AF < 0.05, ]
dat[dat$Gene.refGene == "TAS2R43", ]

## ---- Genome-wide check for rare-in-gnomAD, common-in-NI loci ---------------
rare_ni_enriched <- dat[dat$AF < 0.01 & dat$NI > 0.05, ]
table(rare_ni_enriched$Chr)
