#!/bin/bash
# Merge single-sample SV (ins/dels) VCFs, collapse near-duplicates with Truvari, and build a VG graph.
# Usage:
#   ./02_merge_collapse_graph.sh <vcf_dir> <output_prefix> <reference.fasta> <conda_env>
# Example:
#   ./02_merge_collapse_graph.sh results/ins_del_vcfs results/pangenome/peach \
#       data/Prunus_persica_V2-pseudo.fa /home/user/miniconda3/envs/svtools
#
# Outputs (given output_prefix=/path/out/prefix):
#   /path/out/prefix.merged.vcf.gz
#   /path/out/prefix.tru-merged.vcf.gz
#   /path/out/prefix.tru-collapsed.vcf
#   /path/out/prefix.graph.vg

set -euo pipefail

if [[ $# -ne 4 ]]; then
  echo "Usage: $0 <vcf_dir> <output_prefix> <reference.fasta> <conda_env>" >&2
  exit 1
fi

VCF_DIR="$1"         # directory containing input VCFs (each 1 sample)
OUT_PREFIX="$2"      # path prefix for outputs (no extension). Its directory will be created.
REFERENCE="$3"       # reference FASTA (indexed if needed by downstream tools)
CONDA_ENV="$4"       # conda env with bcftools, tabix, truvari, vg

# Threads (tweak as needed)
THREADS="${THREADS:-8}"

# Prepare output directory
OUT_DIR="$(dirname "$OUT_PREFIX")"
mkdir -p "$OUT_DIR"

echo "[Setup] Activating conda environment: $CONDA_ENV"
if ! command -v conda >/dev/null 2>&1; then
  echo "conda not found on PATH. Ensure Conda is installed and available." >&2
  exit 1
fi
eval "$(conda shell.bash hook)"
conda activate "$CONDA_ENV"

# Tool checks (fail fast if missing)
for tool in bcftools tabix truvari vg; do
  if ! command -v "$tool" >/dev/null 2>&1; then
    echo "Required tool '$tool' not found. Activate the right environment." >&2
    exit 1
  fi
done

echo "[Step 1] Renaming variant IDs per file to ensure global uniqueness..."
cd "$VCF_DIR"

RENAMED_DIR="${OUT_DIR}/renamed_vcfs"
mkdir -p "$RENAMED_DIR"

# Process .vcf and .vcf.gz
shopt -s nullglob
VCFS=( *.vcf *.vcf.gz )
if [[ ${#VCFS[@]} -eq 0 ]]; then
  echo "No VCFs found in $VCF_DIR" >&2
  exit 1
fi

LIST_FILE="${OUT_DIR}/listofvcffiles.txt"
: > "$LIST_FILE"

for vcf in "${VCFS[@]}"; do
  base="${vcf%.vcf}"
  base="${base%.gz}"            # handle .vcf.gz
  sample="${base}"              # use filename stem as a proxy for sample name
  out_vcf="${RENAMED_DIR}/${sample}_renamed.vcf.gz"

  echo "  * $vcf → $(basename "$out_vcf")"

  # Ensure compressed input for bcftools
  if [[ "$vcf" == *.vcf.gz ]]; then
    in_vcf="$vcf"
  else
    bgzip -f "$vcf"
    tabix -p vcf "${vcf}.gz"
    in_vcf="${vcf}.gz"
  fi

  # Overwrite ID with a unique, reproducible schema including sample name
  # e.g., SAMPLE_CHR_POS_SVTYPE_SVLEN
  bcftools annotate -x ID -I +"${sample}"'_%CHROM_%POS_%INFO/SVTYPE_%SVLEN' \
    -Oz -o "$out_vcf" "$in_vcf"
  tabix -p vcf "$out_vcf"

  echo "$out_vcf" >> "$LIST_FILE"
done

echo "[Step 2] Merging VCFs with bcftools (may take a while)..."
MERGED_VCF="${OUT_PREFIX}.merged.vcf.gz"
bcftools merge -m none --force-samples -l "$LIST_FILE" \
  -Oz -o "$MERGED_VCF" --threads "$THREADS"
tabix -p vcf "$MERGED_VCF"

echo "[Step 3] Collapsing nearby/duplicate SVs with Truvari..."
TRU_MERGED_VCF="${OUT_PREFIX}.tru-merged.vcf.gz"
TRU_COLLAPSED="${OUT_PREFIX}.tru-collapsed.vcf"

truvari collapse \
  -i "$MERGED_VCF" \
  -p 0 -P 0.5 -s 0 \
  -o "$TRU_MERGED_VCF" \
  -c "$TRu_COLLAPSED" \
  --threads "$THREADS" || {
    echo "Warning: Truvari collapse without --threads (older versions). Retrying..."
    truvari collapse \
      -i "$MERGED_VCF" \
      -p 0 -P 0.5 -s 0 \
      -o "$TRU_MERGED_VCF" \
      -c "$TRU_COLLAPSED"
  }

tabix -p vcf "$TRU_MERGED_VCF"

echo "[Step 4] Constructing variation graph with vg..."
VG_OUT="${OUT_PREFIX}.graph.vg"
vg construct -r "$REFERENCE" -v "$TRU_MERGED_VCF" -a -f > "$VG_OUT"

echo "[Done] Outputs:"
echo "  Merged VCF:           $MERGED_VCF"
echo "  Truvari merged VCF:   $TRU_MERGED_VCF"
echo "  Truvari collapsed:    $TRU_COLLAPSED"
echo "  VG graph:             $VG_OUT"
