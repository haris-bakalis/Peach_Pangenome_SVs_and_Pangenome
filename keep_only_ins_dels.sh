echo "Creating INS/DEL-only VCF..."

INS_DEL_VCF="${OUTPUT_DIR}/${VARIETY}.INS_DEL.vcf.gz"

# If bcftools is available, use it (recommended)
if command -v bcftools >/dev/null 2>&1; then
  RAW_VCF="$variants_vcf"

  # bgzip+index if needed
  if [[ "$RAW_VCF" =~ \.vcf$ ]]; then
    bgzip -c "$RAW_VCF" > "${RAW_VCF}.gz"
    tabix -p vcf "${RAW_VCF}.gz"
    RAW_VCF="${RAW_VCF}.gz"
  fi

  # Keep only DEL and INS
  bcftools view -i 'INFO/SVTYPE="DEL" || INFO/SVTYPE="INS"' \
    -Oz -o "$INS_DEL_VCF" "$RAW_VCF"
  tabix -p vcf "$INS_DEL_VCF"
