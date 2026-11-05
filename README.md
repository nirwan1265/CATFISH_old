# CATFISH

**CATFISH** (Combining **C**auchy combination test (ACAT), **A**daptive **T**Fisher, **F**isher's method, m**I**n-P, and **S**touffer's method for **H**olistic pathway analysis) is an R package for post-GWAS pathway analysis. It takes MAGMA gene-level results, scores every pathway with five complementary tests combined into one omnibus, and calibrates the omnibus with a Westfall--Young procedure that accounts for gene--gene correlation and controls the family-wise error rate across pathways. It then ranks candidate genes by combining SNP-, gene- and pathway-level evidence.

## What CATFISH does

- Runs five pathway tests on gene-level GWAS statistics (ACAT, Fisher, adaptive soft TFisher, minP, Stouffer) and combines them with ACAT into an omnibus.
- Calibrates the omnibus against a genome-wide correlated null: each null replicate draws one vector of gene z-scores, chromosome by chromosome, from the MAGMA gene--gene correlations, and scores every pathway on that same draw. This gives each pathway a calibrated p-value and a Westfall--Young family-wise p-value.
- Works in log space, so gene p-values far below 1e-300 keep their order.
- Wraps MAGMA (annotation and gene analysis, per chromosome) and reads its output back.
- Ships pathways for maize, sorghum, arabidopsis, PlantCyc, *Drosophila* and human (HumanCyc, Reactome), and reads GMT, BioCyc and PMN files.
- Ranks candidate genes with `gene_prioritize()`.

## Pathway signal archetypes

Different pathway signal shapes favour different tests, and the shape of a real pathway's signal is unknown, so CATFISH runs all five and lets the omnibus pick up signal whatever its shape.

| Archetype | Abbrev | Gene-level p-value pattern | Best methods |
|-----------|--------|---------------------------|-------------|
| Strong Sparse Signal | **SSS** | Few genes with extremely small p-values; most near null | ACAT, minP |
| Moderate Dense Signal | **MDS** | Many genes with modest consistent signal, no dominant gene | Fisher, Stouffer, soft TFisher |
| Weak Diffuse Enrichment | **WDE** | Subtle pathway-wide shift; few or none individually significant | Stouffer, Fisher |
| Mixed Signal Architecture | **MSA** | One/few strong drivers plus broader moderate support | soft TFisher, Fisher, ACAT |
| Single-Gene Proxy | **SGP** | One dominant gene; rest near null | minP, ACAT *(interpret cautiously)* |

## Installation

### 1) MAGMA (external)

Download MAGMA from https://ctg.cncr.nl/software/magma. It is needed for the MAGMA steps only; the pathway and gene-ranking functions run without it.

### 2) CATFISH

CATFISH depends on the ACAT package, which is only on GitHub, so install it first.

```r
# install.packages("remotes")
remotes::install_github("yaowuliu/ACAT")
remotes::install_github("nirwan1265/CATFISH")
```

### 3) MAGMA path

```r
library(CATFISH)
options(magma.path = "/absolute/path/to/magma")
magma_path()
```

## End-to-end workflow

The GWAS itself is run outside CATFISH (GEMMA, GAPIT, PLINK, ...); CATFISH starts from its SNP-level summary statistics. The vignette (`vignette("catfish-workflow", package = "CATFISH")`) has a version of this that runs on simulated data in a few seconds.

### 1. Gene locations

```r
gff3_to_geneloc(gff = "data/reference.gff3", out = "results/species.genes.loc")
```

Built-in gene location files exist for maize, sorghum, arabidopsis and *Drosophila* (`inst/extdata/*.genes.loc`).

### 2. MAGMA annotation

```r
ann <- magma_annotate(
  stats_file     = "data/gwas.tsv",
  rename_columns = c(CHR = "chr", SNP = "rs", POS = "ps", PVALUE = "p_wald"),
  gene_loc       = "results/species.genes.loc",
  out_prefix     = "trait",
  out_dir        = "results/annot",
  window         = c(25, 25),          # kb up/downstream
  nonhuman       = TRUE
)
```

### 3. MAGMA gene analysis, per chromosome

```r
magma_gene(
  bfile          = "data/plink_prefix",
  gene_annot     = ann$gene_annot,
  stats_file     = "data/gwas.tsv",
  rename_columns = c(CHR = "chr", SNP = "rs", PVALUE = "p_wald", NMISS = "n_miss"),
  n_total        = 1000,
  out_prefix     = "trait",
  out_dir        = "results/magma",
  gene_model     = "snp-wise=mean",
  chroms         = 1:10,               # one run per chromosome: trait_chr1, ...
  n_threads      = 10                  # 1 = one after another
)
```

On a cluster, run one chromosome per array job with `chroms = <that chromosome>`.

**Gene model.** Westfall--Young calibration assumes gene p-values are uniform where there is no association. MAGMA's SNP-wise mean model meets this when the GWAS is well calibrated; the SNP-wise top model, and `multi=snp-wise` which includes it, are often inflated in panels with long-range LD, and that inflation is then read as pathway signal. `magma_read_genes()` prints the median of every gene p-value column (about 0.5 expected) so this can be checked.

### 4. Read all chromosomes back

```r
genes <- magma_read_genes(dir = "results/magma", pattern = "^trait_chr.*\\.genes\\.out$")

raw <- list.files("results/magma", "^trait_chr.*\\.genes\\.raw$", full.names = TRUE)
magma_genesraw_to_cor_pairs_banded(raw, "results/magma/trait_gene_cor_pairs.txt",
                                   gene_regex = "^SORBI")   # pattern of your gene IDs
```

MAGMA computes gene--gene correlations only within a chromosome, so all chromosome files are needed; they are written to one file.

### 5. Pathways

```r
pw <- catfish_load_pathways(species = "sorghum")
# built-in: "maize", "sorghum", "arabidopsis", "plant", "fly", "human", "human_reactome"
# your own files: catfish_load_pathways(file = "x_pathways.tsv")  (PMN format),
#                 catfish_load_gmt("x.gmt"), catfish_load_biocyc("x.csv")
```

A data.frame with `pathway_id`, `gene_id` (and optionally `pathway_name`), or a named list of gene vectors, works too.

### 6. Westfall--Young pathway test

```r
res <- catfish_wy_pathways(
  gene_results = genes,
  pathways     = pw,
  gene_col     = "GENE",
  p_col        = "P",                  # P_SNPWISE_MEAN for a multi-model run
  chr_col      = "CHR",
  cor_file     = "results/magma/trait_gene_cor_pairs.txt",
  B            = 10000,
  min_genes    = 5,
  max_genes    = 200,
  seed         = 20260706
)
head(res)
attr(res, "wy_threshold_05")
```

`catfish_omni2_pathways(..., perm_mode = "wy")` runs the same calibration and also reports the five component tests, multiple-testing adjustments and diagnostics.

### 7. Candidate genes

```r
cand <- gene_prioritize(
  magma           = genes,
  magma_p_col     = "P",
  pathway_results = res,
  pathways        = pw,
  gwas            = "data/gwas.tsv",   # optional SNP layer
  gwas_snp_col    = "rs", gwas_p_col = "p_wald",
  gene_annot      = ann$gene_annot,
  targets         = c("GENE_A", "GENE_B")   # optional known genes
)
```

Each gene's score is

```
score = S_gwas + S_magma + 0.5 S_pathway + 0.3 log10(1 + n_pathways)
```

where each S is the rank score -log10(rank / N) of that layer: the strongest SNP in or near the gene, the MAGMA gene p-value, and the strongest tested pathway that contains the gene (by `z_std`).

## Output of `catfish_wy_pathways()`

One row per pathway, sorted by `z_std`:

| Column | Meaning |
|---|---|
| `pathway_id`, `pathway_name` | Pathway |
| `z_std` | Observed omnibus standardised against the pathway's own null (more negative = stronger) |
| `p_fwer_std` | Westfall--Young family-wise p-value; the headline result |
| `p_calibrated` | Per-pathway calibrated p-value |
| `n_genes` | Genes of the pathway present in the gene results |
| `log10_p_omnibus`, `p_omnibus` | Analytic (uncalibrated) omnibus |
| `p_fwer_wy` | Family-wise p-value on the raw omnibus scale |
| `p_bonferroni` | Bonferroni on the analytic omnibus, for comparison |

## Core parameters of `catfish_wy_pathways()`

| Parameter | Meaning | Typical value |
|---|---|---|
| `gene_col`, `p_col`, `chr_col` | Gene ID, gene p-value and chromosome columns | `GENE`, `P`, `CHR` |
| `cor_file` / `cor_pairs` | Gene--gene correlations (gene1, gene2, r) | from `magma_genesraw_to_cor_pairs_banded()` |
| `B` | Null replicates | 10,000 |
| `use_tfisher` | Include soft TFisher in the omnibus | `TRUE` |
| `tau_grid` | Fixed soft TFisher grid; must not be estimated from the data | `c(0.01, 0.05, 0.5)` |
| `min_genes`, `max_genes` | Pathway size limits (genes present) | 5--10, 200 |
| `seed` | Seed for the null draws | any integer |

## Notes

- Gene IDs are matched case-insensitively between gene results, pathways and correlations.
- The built-in sorghum pathways use `SORBI_` gene IDs, the IDs of the built-in `sorghum.genes.loc`.
- Without gene--gene correlations every gene is treated as independent, and the null is no longer LD-aware; CATFISH warns when this happens.
- Gene-length adjustment (`catfish_adjust_gene_p()`) is not needed for MAGMA gene p-values; see the vignette.
