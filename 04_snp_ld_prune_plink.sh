#!/bin/bash

# Input VCF and output base name
INPUT_VCF="onlyPruper.miss10.ALL.snp.vcf"
BASE_NAME="snp_dataset"

module load conda
conda activate plink_env

# Step 1: Convert VCF to PLINK format with unique SNP IDs
echo "[1] Converting VCF to PLINK binary format..."
plink \
  --vcf "$INPUT_VCF" \
  --allow-extra-chr \
  --make-bed \
  --out "$BASE_NAME" \
  --set-missing-var-ids @:#:\$1:\$2

# Step 2: LD pruning
echo "[2] Performing LD pruning..."
plink \
  --bfile "$BASE_NAME" \
  --allow-extra-chr \
  --indep-pairwise 50 5 0.96 \
  --out pruned_snps

# Step 3: Extract pruned SNPs
echo "[3] Extracting pruned SNPs..."
plink \
  --bfile "$BASE_NAME" \
  --allow-extra-chr \
  --extract pruned_snps.prune.in \
  --make-bed \
  --out "${BASE_NAME}_pruned"

# Step 4: Convert to VCF
echo "[4] Converting pruned dataset back to VCF..."
plink \
  --bfile "${BASE_NAME}_pruned" \
  --allow-extra-chr \
  --recode vcf \
  --out "${BASE_NAME}_pruned"

echo "[Done] LD pruning complete. Output: ${BASE_NAME}_pruned.*"
