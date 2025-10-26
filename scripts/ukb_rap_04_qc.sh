#!/bin/bash
# ---- launch_pgs_jobs.sh ----
# Run Swiss Army Knife index concantenation file on UK Biobank RAP
# Author: nicholas.sunderland@bristol.ac.uk

# --- Config ---
PROJECT="HERMES3"

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

# --- Concat chromosomes ---
echo "🚀 Running plink QC"

dx run swiss-army-knife \
  -iin="/dcm_pgs/combined/combined_ukb_pgs.bgen" \
  -iin="/dcm_pgs/combined/combined_ukb_pgs.bgen.bgi" \
  -iin="/Bulk/Imputation/UKB imputation from genotype/ukb22828_c1_b0_v3.sample" \
  -iin="/dcm_pgs/whitebrit.txt" \
  -icmd="plink2 --bgen combined_ukb_pgs.bgen ref-first \
                --sample ukb22828_c1_b0_v3.sample \
                --keep-fam whitebrit.txt \
                --maf 0.01 \
                --hwe 1e-6 midp\
                --geno 0.10 \
                --mind 0.02 \
                --make-pgen \
                --threads 4 \
                --out qc_ukb_pgs" \
  --tag "qc" \
  -imount_inputs=true \
  --instance-type "mem1_ssd2_v2_x4" \
  --priority normal \
  --name "qc" \
  --destination "/dcm_pgs/qc/" \
  --yes

echo "🎉 Jobs submitted!"