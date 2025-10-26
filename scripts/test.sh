#!/bin/bash
# ---- launch_pgs_jobs.sh ----
# Run parallel Swiss Army Knife subset jobs on UK Biobank RAP
# Author: nicholas.sunderland@bristol.ac.uk

# --- Config ---
PROJECT="HERMES3"
COHORT_PREFIX="ukb22828"
RSID_LIST="rsidlist.txt"
PGS_SCORE="PGS004861_clean.txt"
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
for CHR in {21..21}; do
  echo "🚀 Launching chromosome $CHR..."
  dx run swiss-army-knife \
    -iin="${PROJ_DIR}/dcm_pgs_chr21_qc.bgen" \
    -iin="${PROJ_DIR}/dcm_pgs_chr21_qc.bgen.bgi" \
    -iin="${PROJ_DIR}/dcm_pgs_chr21_qc.sample" \
    -iin="${PROJ_DIR}/${PGS_SCORE}" \
    -icmd="bgenix -g dcm_pgs_chr21_qc.bgen -list | head && \
           plink2 --bgen dcm_pgs_chr21_qc.bgen ref-first \
                  --sample dcm_pgs_chr21_qc.sample \
                  --score ${PGS_SCORE} 1 4 6 header cols=+scoresums \
                  --out dcm_pgs_chr21_score" \
    --tag "chr${CHR}" \
    -imount_inputs=true \
    --instance-type "mem1_ssd2_v2_x2" \
    --priority normal \
    --name "subset_chr${CHR}" \
    --destination "/dcm_pgs/" \
    --yes

  echo "✅ Submitted job for chr${CHR}"
done

echo "🎉 All jobs submitted!"



