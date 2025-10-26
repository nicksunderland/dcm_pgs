#!/bin/bash
# ---- launch_pgs_jobs.sh ----
# Run parallel Swiss Army Knife subset jobs on UK Biobank RAP
# Author: nicholas.sunderland@bristol.ac.uk

# --- Config ---
PROJECT="HERMES3"
COHORT_PREFIX="ukb22828"
RSID_LIST="rsidlist.txt"
BULK_DIR="/Bulk/Imputation/UKB imputation from genotype"
PROJ_DIR="/dcm_pgs"
OUT_PREFIX="dcm_pgs_chr"

# --- Login (if needed) ---
if ! dx whoami &>/dev/null; then
  echo "🔐 Not logged in. Opening RAP login page..."
  dx login
fi

# --- Select project ---
echo "📁 Selecting RAP project: $PROJECT"
dx select "$PROJECT"

# --- Confirm environment ---
echo "✅ Using project:"
dx env

# --- Loop through chromosomes ---
for CHR in {1..22}; do
  INPUT_BGEN="${COHORT_PREFIX}_c${CHR}_b0_v3.bgen"
  INPUT_BGEN_IDX="${COHORT_PREFIX}_c${CHR}_b0_v3.bgen.bgi"
  OUTPUT_BGEN="${OUT_PREFIX}${CHR}_subset.bgen"

  echo "🚀 Launching chromosome $CHR..."

  dx run swiss-army-knife \
    -iin="${PROJ_DIR}/${RSID_LIST}" \
    -iin="${BULK_DIR}/${INPUT_BGEN}" \
    -iin="${BULK_DIR}/${INPUT_BGEN_IDX}" \
    -icmd="bgenix -g ${INPUT_BGEN} -incl-rsids ${RSID_LIST} > ${OUTPUT_BGEN} && \
           bgenix -g ${OUTPUT_BGEN} -index" \
    --tag "chr${CHR}" \
    -imount_inputs=true \
    --instance-type "mem1_ssd2_v2_x2" \
    --priority normal \
    --name "subset_chr${CHR}" \
    --destination "/dcm_pgs/${CHR}/" \
    --yes

  echo "✅ Submitted job for chr${CHR}"
done

echo "🎉 All jobs submitted!"


# PGS_SCORE="PGS004861_clean.txt"
# #  OUTPUT_QC="${OUT_PREFIX}${CHR}_qc"
  ##  OUTPUT_PRS="${OUT_PREFIX}${CHR}_score"
# #  INPUT_SAMP="${COHORT_PREFIX}_c${CHR}_b0_v3.sample"
#    -iin="${BULK_DIR}/${INPUT_SAMP}" \
#    -iin="${PROJ_DIR}/${PGS_SCORE}" \
#plink2 --bgen ${OUTPUT_BGEN} ref-first \
#                  --sample ${INPUT_SAMP} \
#                  --maf 0.01 \
#                  --hwe 1e-6 \
#                  --geno 0.10 \
#                  --mind 0.02 \
#                  --export bgen-1.2 \
#                  --out ${OUTPUT_QC} && \
#           bgenix -g ${OUTPUT_QC}.bgen -index && \
#           plink2 --bgen ${OUTPUT_QC}.bgen ref-first \
#                  --sample ${OUTPUT_QC}.sample \
#                  --score ${PGS_SCORE} 1 4 6 header cols=+scoresums \
#                  --out ${OUTPUT_PRS} && \
#           rm -f ${OUTPUT_BGEN} \
#                 ${OUTPUT_BGEN}.bgi \
#                 ${OUTPUT_QC}.bgen \
#                 ${OUTPUT_QC}.bgen.bgi \
#                 ${OUTPUT_QC}.sample" \