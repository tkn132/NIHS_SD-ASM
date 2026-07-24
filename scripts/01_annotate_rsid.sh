#!/usr/bin/env bash
# =============================================================================
# 01_annotate_rsid.sh
#
# Purpose:
#   Annotate significant SD-ASM CpG loci (output of the SD-ASM detection step;
#   see Methods "Detection of sequence-dependent allele-specific methylation
#   sites") with dbSNP rsIDs, annotate the full NI WGS/imputed variant set with
#   rsIDs, and intersect the two to identify which SD-ASM loci correspond to
#   well-characterised NI variants used in downstream eQTL and phenotype
#   association steps.
#
#   ANNOVAR's -filter mode requires a specific ref/alt column order; because
#   the SD-ASM table is not phased to a fixed order, annotation is run twice
#   (as-is, then with alleles flipped) and the results are combined to
#   maximise the number of rsID matches recovered.
#
# Requirements:
#   - ANNOVAR (annotate_variation.pl) with the hg38 avsnp150 database
#   - R (dplyr)
#
# Inputs:
#   - $SDASM_CPG_TABLE   CpG-level SD-ASM output table (chr, pos, alleles, ...)
#   - $NI_IMPUTED_BFILE  NI imputed genotype set, PLINK .bed/.bim/.fam
#
# Outputs:
#   - sig_annotation_cpgs_rsid.tsv               SD-ASM loci with dbSNP rsIDs
#   - NI.finalSNPset.hg38.NIKGP.impute.rsid.tsv  NI variant set with dbSNP rsIDs
#   - sig_annotation_cpgs_rsid_NIoverlap.tsv     SD-ASM loci overlapping NI variants
#   - ASM_NIid.txt                               NI variant IDs for overlapping SD-ASM loci
# =============================================================================
set -euo pipefail

# ---- EDIT THESE PATHS FOR YOUR ENVIRONMENT ---------------------------------
WORKDIR=/path/to/working_dir
ANNOVAR_DIR=/path/to/annovar
SDASM_CPG_TABLE=/path/to/sig_annotation_cpgs.tsv
NI_IMPUTED_BFILE=/path/to/NI.finalSNPset.hg38.NIKGP.impute
# -----------------------------------------------------------------------------

cd "$WORKDIR"

## ---------------------------------------------------------------------------
## Step 1: Annotate SD-ASM CpG loci with dbSNP rsIDs
## ---------------------------------------------------------------------------
awk 'NR>1 {print $1,$3,$3,$14,$13}' "$SDASM_CPG_TABLE" | sed -r 's/^chr//' \
  > temp.sig_annotation_cpgs.avinput

perl "$ANNOVAR_DIR/annotate_variation.pl" temp.sig_annotation_cpgs.avinput \
  "$ANNOVAR_DIR/humandb/" -filter -build hg38 -dbtype avsnp150

awk '{print $1, $2, $3, $5, $4}' temp.sig_annotation_cpgs.avinput.hg38_avsnp150_filtered \
  > temp.sig_annotation_cpgs_flip.avinput

perl "$ANNOVAR_DIR/annotate_variation.pl" temp.sig_annotation_cpgs_flip.avinput \
  "$ANNOVAR_DIR/humandb/" -filter -build hg38 -dbtype avsnp150

cat temp*dropped > sig_annotation_cpgs_rsid.tsv
rm temp*

## ---------------------------------------------------------------------------
## Step 2: Annotate the full NI imputed variant set with dbSNP rsIDs
## ---------------------------------------------------------------------------
awk '{print $1,$4,$4,$5,$6,$2}' "${NI_IMPUTED_BFILE}.bim" > temp.avinput

perl "$ANNOVAR_DIR/annotate_variation.pl" temp.avinput \
  "$ANNOVAR_DIR/humandb/" -filter -build hg38 -dbtype avsnp150

awk '{print $1, $2, $3, $5, $4, $6}' temp.avinput.hg38_avsnp150_filtered \
  > temp_flip.avinput

perl "$ANNOVAR_DIR/annotate_variation.pl" temp_flip.avinput \
  "$ANNOVAR_DIR/humandb/" -filter -build hg38 -dbtype avsnp150

cat temp*dropped > NI.finalSNPset.hg38.NIKGP.impute.rsid.tsv
rm temp*

## ---------------------------------------------------------------------------
## Step 3: Intersect SD-ASM rsIDs with the NI variant rsID set
## ---------------------------------------------------------------------------
cut -f2,8 NI.finalSNPset.hg38.NIKGP.impute.rsid.tsv | grep -vw '\.' > NI.rsid.txt

Rscript - <<'EOF'
library(dplyr)

dat <- read.table("NI.rsid.txt")
names(dat) <- c("rsid", "id")

asm <- read.table("sig_annotation_cpgs_rsid.tsv")
names(asm)[2] <- "rsid"

overlap <- merge(dat, asm, by = "rsid")
write.table(overlap[-3], "sig_annotation_cpgs_rsid_NIoverlap.tsv",
            row.names = FALSE, quote = FALSE)
EOF

cut -f2 -d ' ' sig_annotation_cpgs_rsid_NIoverlap.tsv | sort | uniq > ASM_NIid.txt
