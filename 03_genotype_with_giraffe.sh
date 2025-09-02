#!/usr/bin/env bash
#SBATCH --job-name=genotyping
#SBATCH --ntasks=1
#SBATCH --mem=12G
#SBATCH --cpus-per-task=8
#SBATCH --partition=all
#SBATCH --time=24:00:00
#SBATCH --output=%x_%j.log
#SBATCH --error=%x_%j.err

set -euo pipefail

# XG for pack/call
XG="/scratch2/074-arabidopsis-MITEs/CACTA_Peach/final_files/graph.xg"

# Samples (add more IDs as needed)
SAMPLES=( SRRXXXXXX )

# Read location template
READ_BASE="/scratch2/074-arabidopsis-MITEs/CACTA_Peach/SRA_Gabino"

# Threads
THREADS="${SLURM_CPUS_PER_TASK:-8}"

# ---------------------------------

# Modules & env
module load samtools
module load conda
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate vg_env

# Sanity checks
for f in "$GBZ" "$DIST" "$MIN" "$ZIPS" "$XG"; do
  [[ -s "$f" ]] || { echo "Missing index: $f" >&2; exit 1; }
done
command -v vg >/dev/null || { echo "vg not found in vg_env" >&2; exit 1; }

for POP in "${SAMPLES[@]}"; do
  echo "==> Processing $POP"
  R1="${READ_BASE}/${POP}/${POP}.sra_1.fastq"
  R2="${READ_BASE}/${POP}/${POP}.sra_2.fastq"
  [[ -s "$R1" && -s "$R2" ]] || { echo "Missing FASTQs for $POP" >&2; exit 1; }

  OUTDIR="${POP}"
  mkdir -p "$OUTDIR"

  GAM="${OUTDIR}/${POP}.gam"
  PACK="${OUTDIR}/${POP}.pack"
  STATS="${OUTDIR}/${POP}.alignments.txt"
  VCF="${OUTDIR}/${POP}.vcf"
  VCFGZ="${VCF}.gz"

  # Map reads to the graph with Giraffe (paired-end)
  vg giraffe -t "$THREADS" \
    -Z "$GBZ" -d "$DIST" -m "$MIN" -z "$ZIPS" \
    -f "$R1" -f "$R2" > "$GAM"

  # Mapping stats
  vg stats -a "$GAM" > "$STATS"

  # Pack coverage and call variants
  vg pack -t "$THREADS" -x "$XG" -g "$GAM" -o "$PACK"
  vg call -t "$THREADS" -a -k "$PACK" "$XG" > "$VCF"

  # Compress & index VCF
  bgzip -f "$VCF"
  tabix -f -p vcf "$VCFGZ"

  if [[ -s "$VCFGZ" ]]; then
    echo "✅ $POP: VCF generated at $VCFGZ"
    # Optional: remove uncompressed VCF to save space
    rm -f "$VCF"
  else
    echo "❌ $POP: VCF generation failed" >&2
    exit 1
  fi

  # Keep only the essentials
  # (GAM + stats + VCF.gz + VCF.gz.tbi are preserved)
  find "$OUTDIR" -type f ! -name "$(basename "$GAM")" \
                         ! -name "$(basename "$STATS")" \
                         ! -name "$(basename "$VCFGZ")" \
                         ! -name "$(basename "$VCFGZ").tbi" \
                         -exec rm -f {} +
done
