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
echo "🚀 Concatenating chromosomes"

dx run swiss-army-knife \
  -iin="/dcm_pgs/qc/qc_ukb_pgs.pgen" \
  -iin="/dcm_pgs/qc/qc_ukb_pgs.pvar" \
  -iin="/dcm_pgs/qc/qc_ukb_pgs.psam" \
  -iin="/dcm_pgs/PGS004861_clean.txt" \
  -icmd="plink2 --pfile qc_ukb_pgs \
                --score PGS004861_clean.txt 1 4 6 header cols=+scoresums \
                --out score_ukb_pgs" \
  --tag "score" \
  -imount_inputs=true \
  --instance-type "mem1_ssd2_v2_x4" \
  --priority normal \
  --name "score" \
  --destination "/dcm_pgs/score/" \
  --yes

echo "🎉 Jobs submitted!"