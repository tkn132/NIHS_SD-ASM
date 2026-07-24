# =============================================================================
# figure3_distance_to_tss.R
#
# Purpose:
#   Compare the distribution of distances to the nearest transcription start
#   site (TSS) between significant SD-ASM loci and all tested heterozygous
#   SNPs, to characterise clustering of SD-ASM loci near TSSs (Figure 3;
#   "77.6% of significant eQTLs were located within 50 Kbp of the TSS").
#
# Inputs:
#   - col3.distance : bedtools-closest output for significant SD-ASM loci
#                      (chr, start, end, distance-to-nearest-TSS)
#   - all.distance   : same, for all tested heterozygous SNPs
#
# Output:
#   - printed summary tables / plots (distance histograms, proportion by bin)
# =============================================================================

library(ggplot2)
library(dplyr)
library(tidyr)
library(stringr)

sig <- read.table("col3.distance")
sig$set <- "sig"
names(sig) <- c("chr", "start", "end", "dis", "set")
sig$chr <- as.integer(str_remove(sig$chr, "chr"))

all <- read.table("all.distance")
all$set <- "all"
names(all) <- c("chr", "start", "end", "dis", "set")
all$chr <- str_remove(all$chr, "chr")

## Remove the MHC/HLA region (chr6:28,510,120-33,480,577, GRCh38)
## https://www.ncbi.nlm.nih.gov/grc/human/regions/MHC
sig <- sig[!(sig$chr == 6 & sig$start > 28510120 & sig$start < 33480577), ]
all <- all[!(all$chr == 6 & all$start > 28510120 & all$start < 33480577), ]

dat <- rbind(sig, all)

## ---- Overlaid histogram, loci within 25 kb of TSS --------------------------
sdat <- dat[abs(dat$dis) < 25000, ]

ggplot() +
  geom_histogram(data = sdat[sdat$set == "sig", ],
                  mapping = aes(x = dis, y = after_stat(count / sum(count))),
                  binwidth = 500, fill = "#E7B800", color = "#E7B800", alpha = 0.4) +
  geom_histogram(data = sdat[sdat$set == "all", ],
                  mapping = aes(x = dis, y = after_stat(count / sum(count))),
                  binwidth = 500, fill = "#00AFBB", color = "#00AFBB", alpha = 0.1) +
  scale_y_continuous(labels = scales::percent, limits = c(0, 0.05)) +
  labs(x = "Distance to TSS (bp)", y = "Proportion of loci")

## ---- Line chart of SD-ASM loci proportion by distance bin ------------------
snp_data <- sig %>%
  mutate(Distance_bin = cut(dis, breaks = seq(-1000000, 1000000, by = 5000),
                             include.lowest = TRUE))

distance_counts <- snp_data %>%
  group_by(Distance_bin) %>%
  summarise(Count = n()) %>%
  mutate(Proportion = Count / sum(Count)) %>%
  filter(!is.na(Distance_bin))

distance_counts$Distance_mid <- as.numeric(gsub("[^0-9.-]", "", levels(distance_counts$Distance_bin)))

ggplot(distance_counts, aes(x = Distance_mid, y = Proportion)) +
  geom_line() +
  labs(x = "Distance to TSS (bp)", y = "Proportion of SD-ASM loci") +
  theme_minimal() +
  scale_x_continuous(breaks = seq(-10000, 10000, by = 5000))

## ---- Bar plot of SD-ASM vs. all loci by 1 kb distance bin ------------------
sig_counts <- sig %>% count(dis, name = "sig")
all_counts <- all %>% count(dis, name = "all")

dat2 <- merge(sig_counts, all_counts, by = "dis", all = TRUE) %>%
  pivot_longer(sig:all, names_to = "set", values_to = "count") %>%
  as.data.frame()

dat2$bin <- cut(dat2$dis, breaks = seq(-200000, 200000, 1000))
dat2 <- dat2[complete.cases(dat2), ]

sdat2 <- dat2 %>%
  group_by(set, bin) %>%
  summarise(count = sum(count, na.rm = TRUE), .groups = "drop") %>%
  group_by(set) %>%
  mutate(Percentage = count / sum(count, na.rm = TRUE) * 100) %>%
  as.data.frame()

ggplot(sdat2[sdat2$set == "sig", ], aes(x = bin, y = Percentage)) +
  geom_bar(stat = "identity") +
  labs(title = "Significant SD-ASM loci")

ggplot(sdat2[sdat2$set == "all", ], aes(x = bin, y = Percentage)) +
  geom_bar(stat = "identity") +
  labs(title = "All tested heterozygous SNPs")
