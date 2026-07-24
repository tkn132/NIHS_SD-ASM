#!/usr/bin/env bash
# =============================================================================
# 02_eqtl_mapping.sh
#
# Purpose:
#   Prepare genotype and expression data for OSCA and run cis-eQTL mapping
#   restricted to SD-ASM loci (Methods "Association analyses"). Associations
#   are tested between each SD-ASM variant and the expression of transcripts
#   whose TSS lies within +/-2 Mb.
#
# Requirements:
#   - PLINK2
#   - OSCA v0.46.1
#   - R
#
# Inputs:
#   - $NI_IMPUTED_BFILE    NI imputed genotypes, PLINK .bed/.bim/.fam
#   - $ASM_SNP_LIST        SD-ASM variant IDs overlapping NI variants (from step 01)
#   - $EXPR_MATRIX         Normalised transcript expression matrix (samples x probes)
#   - $UUID_MAP            Sample-ID <-> genotype-UUID mapping table
#   - $BIOMART_ANNOTATION  Probe-to-gene/position annotation (Biomart export)
#   - $COVAR / $QCOVAR     Covariate files (sex; age + expression PCs/latent factors)
#
# Outputs:
#   - NI_gene_expression_overlap.{bed,bim,fam}  Genotype subset for eQTL mapping
#   - NI_gene_expression.{bod,oii,opi}          OSCA expression files
#   - NI_ASM_eQTL_age_sex_suggestive.txt        Suggestive (p<5e-5) cis-eQTL results
#
# Note:
#   Final FDR<0.05 cis-eQTL calls reported in the manuscript (3,348 associations
#   involving 2,327 SD-ASM loci and 1,853 genes) are obtained by applying
#   Benjamini-Hochberg correction across all tested SD-ASM variant-gene pairs
#   to this suggestive output.
# =============================================================================
set -euo pipefail

# ---- EDIT THESE PATHS FOR YOUR ENVIRONMENT ---------------------------------
WORKDIR=/path/to/working_dir/gene_expression
NI_IMPUTED_BFILE=/path/to/NI.finalSNPset.hg38.NIKGP.impute
ASM_SNP_LIST=/path/to/ASM_NIid.txt
EXPR_MATRIX=/path/to/NI.exprs.data.fullyProcessed.csv
UUID_MAP=/path/to/UUID_map.txt
BIOMART_ANNOTATION=/path/to/eQTL_map_Biomart_May2020.txt
COVAR=/path/to/covar_sex.txt
QCOVAR=/path/to/qcovar_age_0pc.txt
OSCA=/path/to/osca-0.46.1
PLINK2=/path/to/plink2
# -----------------------------------------------------------------------------

cd "$WORKDIR"

## ---------------------------------------------------------------------------
## Step 1: Identify samples with matched expression + genotype data
## ---------------------------------------------------------------------------
head -1 "$EXPR_MATRIX" | tr ',' '\n' | awk 'NR>1{print $1}' > gene_expression_samples.txt
grep -f gene_expression_samples.txt "$UUID_MAP" | awk '{print $1, $1}' | grep -v 'NA' \
  > gene_expression_samples_UUID.txt

## ---------------------------------------------------------------------------
## Step 2: Extract SD-ASM SNPs for the matched samples
## ---------------------------------------------------------------------------
"$PLINK2" --bfile "$NI_IMPUTED_BFILE" \
  --rm-dup \
  --make-bed \
  --keep gene_expression_samples_UUID.txt \
  --extract "$ASM_SNP_LIST" \
  --out NI_gene_expression_overlap

## ---------------------------------------------------------------------------
## Step 3: Build OSCA expression (.bod), sample (.oii), and probe (.opi) files
## ---------------------------------------------------------------------------
Rscript - <<EOF
map <- read.table("$UUID_MAP", header = TRUE)

trans <- read.csv("$EXPR_MATRIX", header = TRUE, row.names = 1)
trans <- as.data.frame(t(trans))
trans\$XSID <- rownames(trans)
trans <- merge(map, trans, by = "XSID")
trans <- trans[!is.na(trans\$UUID), ]
trans\$XSID <- trans\$UUID
names(trans)[c(1, 2)] <- c("FID", "IID")

write.table(trans, "NI_gene_expression.txt", row.names = FALSE, quote = FALSE)
EOF

"$OSCA" --efile NI_gene_expression.txt --gene-expression --make-bod --out NI_gene_expression

awk '{print $1, $2, 0, 0, "NA"}' NI_gene_expression_overlap.fam > NI_gene_expression.oii

cut -f1 -d, "$EXPR_MATRIX" | tail -n +2 > gene_expression_transcripts.txt

Rscript - <<EOF
trans <- read.table("gene_expression_transcripts.txt")
biomart <- read.table("$BIOMART_ANNOTATION", header = TRUE, sep = "\t")

dat <- biomart[biomart\$ILLUMINA.HumanHT.12.V4.probe %in% trans\$V1, ]
dat\$Strand2 <- ifelse(dat\$Strand == 1, "+", "-")
dat <- dat[dat\$Chromosome.scaffold.name %in% 1:22, ]  # autosomes only

write.table(dat[-5], "annotated.opi", col.names = FALSE, row.names = FALSE,
            quote = FALSE, sep = "\t")
EOF

"$OSCA" --befile NI_gene_expression --update-opi annotated.opi

## ---------------------------------------------------------------------------
## Step 4: Run cis-eQTL mapping (SD-ASM loci vs. transcripts within +/-2Mb of TSS)
## ---------------------------------------------------------------------------
for i in {1..3}; do
  "$OSCA" \
    --eqtl \
    --bfile NI_gene_expression_overlap \
    --befile NI_gene_expression \
    --covar "$COVAR" \
    --qcovar "$QCOVAR" \
    --thread-num 16 \
    --maf 0.05 \
    --task-num 3 \
    --task-id "$i" \
    --out NI_ASM_eQTL_age_sex
done

## ---------------------------------------------------------------------------
## Step 5: Extract suggestive (p < 5e-5) cis-eQTL associations
## ---------------------------------------------------------------------------
for i in {1..3}; do
  "$OSCA" --beqtl-summary "NI_ASM_eQTL_age_sex_3_$i" \
    --query 5e-5 \
    --out "temp${i}.txt"
done

head -1 temp1.txt > NI_ASM_eQTL_age_sex_suggestive.txt
for i in {1..3}; do
  tail -n +2 "temp${i}.txt" >> NI_ASM_eQTL_age_sex_suggestive.txt
done
rm temp*.txt
