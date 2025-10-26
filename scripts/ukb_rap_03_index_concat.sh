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
  -iin="/dcm_pgs/combined/combined_ukb_pgs.bgen" \
  -icmd="bgenix -g combined_ukb_pgs.bgen -index" \
  --tag "idx_cat" \
  -imount_inputs=true \
  --instance-type "mem1_ssd2_v2_x4" \
  --priority normal \
  --name "index_concatenation" \
  --destination "/dcm_pgs/combined/" \
  --yes

echo "🎉 Jobs submitted!"