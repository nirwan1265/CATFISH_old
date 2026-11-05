#!/bin/bash
# One permutation, whole genome, with vcf2gwas. The kinship and PCs come from
# 00_prepare_genome_LSF.sh (-k, -cf), so nothing is recomputed per run.
# Runs on a private copy of the genome VCF: vcf2gwas writes and deletes a
# .csi next to the VCF, which breaks other jobs reading the same file.
#   bash 02_run_gwas_single.sh <perm>
set -euo pipefail
source /share/maize/ntanduk/CATFISH/simulation/bap_perm_wy/config.sh
B="$1"; BPAD=$(printf "%03d" "$B")
PHENO="${PHENO_DIR}/pheno_perm_${BPAD}.csv"
OUTDIR="${GWAS_DIR}/perm_${BPAD}"; ASSOC="${OUTDIR}/assoc.txt"
for f in "$PHENO" "$GENOME_VCF" "$KINSHIP" "$PCS_CSV"; do [ -f "$f" ] || { echo "missing $f" >&2; exit 1; }; done
mkdir -p "$OUTDIR"; rm -rf "$OUTDIR/Output"; rm -f "$ASSOC"
JOB_TMP=$(mktemp -d "${TMPDIR:-/tmp}/bap_gwas_${BPAD}_XXXXXX")
trap 'rm -rf "$JOB_TMP"' EXIT
VCF_LOCAL="${JOB_TMP}/genome.vcf.gz"
cp "$GENOME_VCF" "$VCF_LOCAL"
cd "$OUTDIR"
echo "[gwas] perm=$B"
"$VCF2GWAS" -v "$VCF_LOCAL" -pf "$PHENO" -ap -k "$KINSHIP" -cf "$PCS_CSV" -ac -lmm 1 -np -nl -T 1
mapfile -t F < <(find Output -name '*.assoc.txt' -print | sort)
[ "${#F[@]}" -eq 1 ] || { echo "[gwas] expected 1 .assoc.txt, found ${#F[@]}" >&2; printf '  %s\n' "${F[@]}" >&2; exit 1; }
cp "${F[0]}" "$ASSOC"
rm -rf "$OUTDIR/Output"
echo "[gwas] -> $ASSOC ($(($(wc -l < "$ASSOC") - 1)) SNPs)"
