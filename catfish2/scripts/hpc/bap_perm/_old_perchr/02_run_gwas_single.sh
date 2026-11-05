#!/bin/bash
# One (permutation, chromosome) GWAS with vcf2gwas, on a private copy of the
# VCF so parallel jobs never share (and delete) each other's .csi index.
# Same command as the real BAP GWAS: -cf PCA -c 3 -lmm 1.
set -euo pipefail
source /share/maize/ntanduk/CATFISH/simulation/bap_perm_wy/config.sh
B="$1"; CHR="$2"; BPAD=$(printf "%03d" "$B")
PHENO="${PHENO_DIR}/pheno_perm_${BPAD}.csv"
OUTDIR="${GWAS_DIR}/perm_${BPAD}/chr${CHR}"
ASSOC="${OUTDIR}/assoc.txt"
VCF_SRC="${VCF_DIR}/$(echo "$VCF_PREFIX_PATTERN" | sed "s/@CHR@/${CHR}/g")"
[ -f "$PHENO" ] || { echo "missing $PHENO (run step 1)" >&2; exit 1; }
[ -f "$VCF_SRC" ] || { echo "missing $VCF_SRC" >&2; exit 1; }
mkdir -p "$OUTDIR"; rm -rf "$OUTDIR/Output"; rm -f "$ASSOC"
JOB_TMP=$(mktemp -d "${TMPDIR:-/tmp}/bap_gwas_${BPAD}_chr${CHR}_XXXXXX")
trap 'rm -rf "$JOB_TMP"' EXIT
VCF_LOCAL="${JOB_TMP}/$(basename "$VCF_SRC")"
cp "$VCF_SRC" "$VCF_LOCAL"
for ext in csi tbi; do [ -f "${VCF_SRC}.${ext}" ] && cp "${VCF_SRC}.${ext}" "${VCF_LOCAL}.${ext}"; done
cd "$OUTDIR"
echo "[gwas] perm=$B chr=$CHR"
"$VCF2GWAS" -v "$VCF_LOCAL" -pf "$PHENO" -ap -cf 'PCA' -c 3 -lmm 1 -np -nl -T 1
mapfile -t F < <(find Output -name '*.assoc.txt' -print | sort)
[ "${#F[@]}" -eq 1 ] || { echo "[gwas] expected 1 .assoc.txt, found ${#F[@]}" >&2; printf '  %s\n' "${F[@]}" >&2; exit 1; }
cp "${F[0]}" "$ASSOC"
rm -rf "$OUTDIR/Output"            # keep only assoc.txt
echo "[gwas] -> $ASSOC"
