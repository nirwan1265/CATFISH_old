#!/bin/bash
# Step 2: GWAS, one array element per permutation (whole genome).
# Submit with the vcf2gwas conda environment active:  bsub < 02_gwas_array_LSF.sh
#BSUB -J "bap_gwas[1-100]%50"
#BSUB -q sara
#BSUB -n 1
#BSUB -W 1440
#BSUB -o logs/gwas.%J.%I.out
#BSUB -e logs/gwas.%J.%I.err
set -euo pipefail
source /share/maize/ntanduk/CATFISH/simulation/bap_perm_wy/config.sh
bash "$WORKING_DIR/02_run_gwas_single.sh" "$LSB_JOBINDEX"
