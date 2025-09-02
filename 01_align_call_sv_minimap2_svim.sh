#!/bin/bash -l
#
# Align an assembly to a reference with minimap2 and call SVs with SVIM-asm.
# Usage:
#   sbatch 01_align_call_sv_minimap2_svim.sh <VARIETY> <QUERY_GENOME_FA> <REFERENCE_FA> <OUTPUT_DIR>
# Example:
#   sbatch 01_align_call_sv_minimap2_svim.sh variety /path/variety.genome.fa /path/Prunus_persica_V2-pseudo.fa results/variety

#SBATCH --nodes=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=6G
#SBATCH --time=6:00:00
#SBATCH --job-name=minimap_svim
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err

set -euo pipefail

if [[ $# -lt 4 ]]; then
  echo "Error: missing arguments."
  echo "Usage: $0 <VARIETY> <QUERY_GENOME_FA> <REFERENCE_FA> <OUTPUT_DIR>"
  exit 1
fi

VARIETY="$1"
QUERY_FA="$2"
REF_FA="$3"
OUTPUT_DIR="$4"

THREADS="${SLURM_CPUS_PER_TASK:-4}"
mkdir -p "$OUTPUT_DIR" # optional if output directory is not existing 

# Load tools (adjust module names to your cluster)
module load conda
conda activate svim_asm
module load minimap2
module load samtools

# Outputs
SAM="${OUTPUT_DIR}/mapped_${VARIETY}.sam"
BAM="${OUTPUT_DIR}/mapped_${VARIETY}.bam"
SORT_BAM="${OUTPUT_DIR}/mapped_sorted_${VARIETY}.bam"
SVIM_DIR="${OUTPUT_DIR}/svim_${VARIETY}"
SVLEN_PNG="${OUTPUT_DIR}/sv-lengths_${VARIETY}.png"
VCF="${OUTPUT_DIR}/variants_${VARIETY}.vcf"
LOG="${OUTPUT_DIR}/svim_${VARIETY}.log"

echo "[$(date)] Aligning ${VARIETY} with minimap2..."
minimap2 -t "$THREADS" -ax asm5 "$REF_FA" "$QUERY_FA" > "$SAM"

echo "[$(date)] Converting, sorting, indexing..."
samtools view -@ "$THREADS" -b "$SAM" -o "$BAM"
samtools sort -@ "$THREADS" -o "$SORT_BAM" "$BAM"
samtools index "$SORT_BAM"
rm -f "$SAM"  # save space

echo "[$(date)] Calling SVs with SVIM-asm (haploid)..."
mkdir -p "$SVIM_DIR"
# SVIM-asm writes outputs into a directory; capture stdout/stderr to LOG
svim-asm haploid "$SVIM_DIR" "$SORT_BAM" "$REF_FA" &> "$LOG"

# Normalize filenames for convenience if SVIM writes defaults
if [[ -f "${SVIM_DIR}/sv-lengths.png" ]]; then
  mv "${SVIM_DIR}/sv-lengths.png" "$SVLEN_PNG"
fi
if [[ -f "${SVIM_DIR}/variants.vcf" ]]; then
  mv "${SVIM_DIR}/variants.vcf" "$VCF"
fi

echo "[$(date)] Done. Outputs:"
echo "  Sorted BAM:   $SORT_BAM"
echo "  SVIM log:     $LOG"
echo "  SV lengths:   $SVLEN_PNG"
echo "  Variants VCF: $VCF"
