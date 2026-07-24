#!/usr/bin/env bash
# =============================================================================
# 03_phenotype_association.sh
#
# Purpose:
#   Test SD-ASM loci (restricted to those acting as eQTLs; see step 02) for
#   association with quantitative phenotypes using a linear mixed model
#   (Methods "Association analyses", phenotype association step).
#
# Requirements:
#   - PLINK2, bcftools, tabix/bgzip
#   - GCTA (v1.94.1)
#   - bedtools
#   - R
#
# Inputs:
#   - $WGS2021_VCF     WGS2021 joint-called VCF (autosomes, core pedigree)
#   - $ASM_RSID_LIST   SD-ASM rsIDs (from step 01)
#   - $PHENO_FILE      Phenotype table (GCTA format: FID IID trait1 trait2 ...)
#   - $COVAR / $QCOVAR Covariate files (sex; age)
#   - $GENE_BED        Gene coordinates (UCSC knownGene, hg38) for nearest-gene annotation
#
# Outputs:
#   - wgs2021.coreped.asm.qc.{bed,bim,fam}     QC'd SD-ASM genotype subset
#   - gcta.out/asmgrm.<trait>.mlma             Per-trait association summary statistics
#   - gcta.out/top.annot/annot.asmgrm.<trait>.txt  Suggestive (p<0.5*) associations
#                                                   annotated with nearest gene
#
#   *Note: the p<0.5 filter below is a loose pre-filter to keep file sizes
#   manageable before nearest-gene annotation; the Bonferroni/suggestive
#   thresholds reported in the manuscript (Table 2) are applied afterwards.
# =============================================================================
set -euo pipefail

# ---- EDIT THESE PATHS FOR YOUR ENVIRONMENT ---------------------------------
WORKDIR=/path/to/working_dir/wgs2021
WGS2021_VCF=/path/to/wgs2021.nochr.autosomes.corepedigree.vcf.gz
ASM_RSID_LIST=/path/to/ASM_rsid.txt
PHENO_FILE=/path/to/pheno.anthro.gctainput.txt
COVAR=/path/to/pheno.anthro.sex.gctainput.txt
QCOVAR=/path/to/pheno.anthro.age.gctainput.txt
GENE_BED=/path/to/knownGene.hg38.Mar20.bed
PLINK2=/path/to/plink2
GCTA=/path/to/gcta
# -----------------------------------------------------------------------------

cd "$WORKDIR"

## ---------------------------------------------------------------------------
## Step 1: Annotate WGS2021 VCF with dbSNP rsIDs
## ---------------------------------------------------------------------------
bcftools query -f '%CHROM\t%INFO/dbsnp_rsid\t%POS\t0\t%ALT\t%REF\n' "$WGS2021_VCF" \
  > rsid.txt
awk -v OFS='\t' '{print $1, $3-1, $3, $2}' rsid.txt > rsid.bed
bgzip -f rsid.bed
tabix -s 1 -b 2 -e 3 rsid.bed.gz

bcftools annotate -c CHROM,FROM,TO,ID -a rsid.bed.gz --threads 32 \
  -Oz -o wgs2021.nochr.autosomes.corepedigree.rsid.vcf.gz "$WGS2021_VCF"

"$PLINK2" --vcf wgs2021.nochr.autosomes.corepedigree.rsid.vcf.gz \
  --make-bed --out wgs2021.nochr.autosomes.corepedigree.rsid

## ---------------------------------------------------------------------------
## Step 2: Extract and QC the SD-ASM variant subset
## ---------------------------------------------------------------------------
"$PLINK2" --bfile wgs2021.nochr.autosomes.corepedigree.rsid \
  --extract "$ASM_RSID_LIST" \
  --make-bed --out wgs2021.coreped.asm

"$PLINK2" --bfile wgs2021.coreped.asm \
  --maf 0.05 --hwe 1e-6 --geno 0.05 --mind 0.05 \
  --make-bed --out wgs2021.coreped.asm.qc

## ---------------------------------------------------------------------------
## Step 3: Build genetic relationship matrices (GRM)
## ---------------------------------------------------------------------------
"$PLINK2" --bfile wgs2021.nochr.autosomes.corepedigree.rsid \
  --maf 0.05 --geno 0.05 --hwe 1e-6 \
  --make-bed --out wgs2021.autosomes.corepedigree.qc

"$GCTA" --bfile wgs2021.autosomes.corepedigree.qc --make-grm \
  --out wgs2021.autosomes.corepedigree.qc

"$GCTA" --bfile wgs2021.coreped.asm --make-grm \
  --out wgs2021.coreped.asm

## ---------------------------------------------------------------------------
## Step 4: Trait-by-trait mixed linear model association (GCTA --mlma)
## ---------------------------------------------------------------------------
mkdir -p gcta.out
N_TRAITS=$(($(head -1 "$PHENO_FILE" | wc -w) - 2))

for i in $(seq 1 "$N_TRAITS"); do
  trait=$(head -1 "$PHENO_FILE" | cut -f$((i + 2)) -d' ')
  "$GCTA" --bfile wgs2021.coreped.asm.qc --mlma \
    --grm wgs2021.coreped.asm \
    --covar "$COVAR" \
    --qcovar "$QCOVAR" \
    --pheno "$PHENO_FILE" \
    --mpheno "$i" \
    --out "gcta.out/asmgrm.$trait"
done

## ---------------------------------------------------------------------------
## Step 5: Pre-filter associations and annotate with nearest gene
## ---------------------------------------------------------------------------
cd gcta.out
ln -sf "$GENE_BED" knownGene.hg38.Mar20.bed
mkdir -p top.annot

for i in $(seq 1 "$N_TRAITS"); do
  trait=$(head -1 "$PHENO_FILE" | cut -f$((i + 2)) -d' ')
  awk '(NR==1) || (($9+0)<0.5)' "asmgrm.$trait.mlma" > "temp.asmgrm.$trait.txt"
  tail -n +2 "temp.asmgrm.$trait.txt" \
    | awk '{print "chr"$1"\t"$3"\t"$3"\t"$2}' \
    | sort -u -k1,1 -k2,2n > "temp.asmgrm.$trait.bed"
  bedtools closest -d -a "temp.asmgrm.$trait.bed" -b knownGene.hg38.Mar20.bed \
    > "temp.asmgrm.$trait.closest"
done

Rscript - <<EOF
pheno <- read.table("$PHENO_FILE", header = TRUE)
traits <- names(pheno)[c(-1, -2)]

for (t in traits) {
  dat  <- read.table(paste0("temp.asmgrm.", t, ".txt"), header = TRUE)
  clst <- read.table(paste0("temp.asmgrm.", t, ".closest"))
  clst <- clst[!duplicated(clst[, c(1, 4, 8)]), c(4, 8, 9)]
  names(clst) <- c("SNP", "gene", "dis2gene")
  dat <- merge(dat, clst, by = "SNP")
  write.table(dat, paste0("top.annot/annot.asmgrm.", t, ".txt"),
              sep = "\t", row.names = FALSE, quote = FALSE)
}
EOF

rm -f temp.asmgrm.*
