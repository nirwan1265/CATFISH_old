devtools::document()   
devtools::load_all()


library(MAGCAT)
library(metapro)
library(sumFREGAT)
library(metap)

# 1. See where it's installed
.libPaths()

install.packages("TFisher", type = "source")

# 2. Remove TFisher from your main library
remove.packages("TFisher", lib = .libPaths()[1])



# gff3_to_geneloc
gff_path <- "/Users/nirwantandukar/Downloads/Zm-B73-REFERENCE-NAM-5.0_Zm00001eb.1.gff3"

gene_loc_out <- "/Users/nirwantandukar/Downloads/maize.genes.loc"

gff3_to_geneloc(
  gff        = "/Users/nirwantandukar/Downloads/Zm-B73-REFERENCE-NAM-5.0_Zm00001eb.1.gff3",
  out        = "inst/extdata/maize.genes.loc",
  chr_prefix = "chr"   # default; strips "chr"
)

# Gene-loc function
# If it has
magma_annotate(
  stats_file = "/Users/nirwantandukar/Documents/Research/results/GWAS/MLM/nitrogen/GWAS_results/nitrogen_0-5cm_maize_LMM.txt",
  rename_columns = c(
    CHR    = "chr",
    SNP    = "rs",
    POS    = "ps",
    PVALUE = "p_wald"   # not used here but keeps things consistent
  ),
  species    = "maize",        # uses the built-in maize.genes.loc
  out_prefix = "N_maize_MLM",
  out_dir    = "annot",
  window     = c(25, 25)
)

# Run MAGMA on gene-level
# If it has NMISS
magma_gene(
  bfile      = "/Users/nirwantandukar/Documents/Research/data/MAGMA/maize/bed_bim_fam_file/all_maize2",
  gene_annot = "annot/N_maize_MLM.genes.annot",
  stats_file = "/Users/nirwantandukar/Documents/Research/results/GWAS/MLM/nitrogen/GWAS_results/nitrogen_0-5cm_maize_LMM.txt",
  n_total    = 3539,
  rename_columns = c(
    CHR    = "chr",
    SNP    = "rs",
    POS    = "ps",
    PVALUE = "p_wald",
    NMISS  = "n_miss" 
  ),
  out_prefix = "N_maize_MLM",
  out_dir    = "magma_genes",
  gene_model = c("multi=snp-wise"),
)

# Run MAGMA on gene-level per chromosome
magma_gene(
  bfile      = "/Users/nirwantandukar/Documents/Research/data/MAGMA/maize/bed_bim_fam_file/all_maize2",
  gene_annot = "annot/N_maize_MLM.genes.annot",
  stats_file = "/Users/nirwantandukar/Documents/Research/results/GWAS/GAPIT/raw_GWAS_MLM_3PC_N.txt",
  n_total    = 3539,
  rename_columns = c(
    CHR    = "Chr",
    SNP    = "SNP",
    POS    = "Pos",
    PVALUE = "P.value",
    NMISS  = "n_miss"
  ),
  out_prefix = "N_maize_MLM",
  out_dir    = "magma_genes_by_chr",
  gene_model = c("multi=snp-wise"),
  chroms     = 1:10,
  n_threads  = 10       # use up to 10 parallel workers
)


# Run MAGMA on gene-level per chromosome
# If it has total number of observations
magma_gene(
  bfile      = "/Users/nirwantandukar/Documents/Research/data/MAGMA/maize/bed_bim_fam_file/all_maize2",
  gene_annot = "annot/N_maize_MLM.genes.annot",
  stats_file = "/Users/nirwantandukar/Documents/Research/results/GWAS/GAPIT/raw_GWAS_MLM_3PC_N.txt",
  # n_total    = 3539, this is ignored. 
  rename_columns = c(
    CHR    = "Chr",
    SNP    = "SNP",
    POS    = "Pos",
    PVALUE = "P.value",
    NOBS   = "nobs" 
  ),
  out_prefix = "N_maize_MLM",
  out_dir    = "magma_multi_snp_wise_genes_by_chr_N_maize",
  gene_model = c("multi=snp-wise"),
  chroms     = 1:10,
  n_threads  = 10
)



# Combine all the chromosome
## Vector of MAGMA gene files (chr1–chr10)
files <- sprintf("/Users/nirwantandukar/Documents/Research/results/MAGMA/MAGCAT/magma_multi_snp_wise_genes_by_chr_N_maize/N_maize_MLM_chr%d.multi_snp_wise.genes.out", 1:10)

# Optional: if they’re in a subdir, prepend the path, e.g.:
# files <- file.path("magma_genes", sprintf("N_maize_MLM_chr%d.snp_wise_top.genes.out", 1:10))

# Read all files
gene_list <- lapply(files, function(f) {
  if (!file.exists(f)) stop("File not found: ", f)
  utils::read.table(f, header = TRUE, stringsAsFactors = FALSE, check.names = FALSE)
})

# Stack into one big data.frame
genes_all_raw <- do.call(rbind, gene_list)

# ---- (Optional but recommended) deduplicate by GENE: keep smallest P ----
if (!"GENE" %in% names(genes_all_raw) || !"P" %in% names(genes_all_raw)) {
  stop("Combined MAGMA file must have columns 'GENE' and 'P'.")
}

# Order by p-value, then drop duplicates so we keep the best P per gene
o <- order(genes_all_raw$GENE, genes_all_raw$P)
genes_all <- genes_all_raw[o, ]
genes_all <- genes_all[!duplicated(genes_all$GENE), ]

# Nice to have: sort by chromosome/start for sanity
if (all(c("CHR", "START") %in% names(genes_all))) {
  genes_all <- genes_all[order(genes_all$CHR, genes_all$START), ]
}

# Inspect
head(genes_all)
nrow(genes_all)
length(unique(genes_all$GENE))
write.table(
  genes_all,
  file      = "/Users/nirwantandukar/Documents/Research/results/MAGMA/MAGCAT/magma_multi_snp_wise_genes_by_chr_N_maize/magma_N_maize.txt",
  sep       = "\t",      # tab-separated
  quote     = FALSE,     # don't quote strings
  row.names = FALSE      # don't write row numbers
  # col.names = TRUE     # default; can be explicit if you want
)


# ACAT
# 1) MAGMA gene-level results (your merged file, with GENE / P / CHR etc.)
genes_all <- read.table("/Users/nirwantandukar/Documents/Research/results/MAGMA/MAGCAT/magma_multi_snp_wise_genes_by_chr_N_maize/magma_N_maize.txt",
                        header = TRUE, stringsAsFactors = FALSE)
head(genes_all)
# 2) Load maize pathways from CornCyc (inst/extdata/pathway/corncyc_...)
maize_pw <- magcat_load_pathways("maize", gene_col = "Gene-name")
head(maize_pw)

# 3) Run ACAT per pathway (no permutations first):
pw_res <- magcat_acat_pathways(
  gene_results = genes_all,
  species      = "maize",
  gene_col     = "GENE",
  p_col        = "P",
  output       = TRUE,
  out_dir      = "acat_results"
)

head(pw_res)

# 4) If you want permutation-calibrated:
pw_res <- magcat_acat_pathways(
  gene_results = genes_all,
  species      = "maize",
  gene_col     = "GENE",
  #pathways     = maize_pw,
  p_col        = "P",
  B            = 1000,
  seed         = 42,
  output       = TRUE,
  out_dir      = "acat_results"
)

head(pw_res)
# Check the path
attr(pw_res, "file")

# 5.) ordmeta
ord_res <- magcat_ordmeta_pathways(
  gene_results = genes_all,
  species      = "maize",
  gene_col     = "GENE",
  p_col        = "P",
  effect_col   = "ZSTAT",
  is_onetail   = FALSE
)
write.csv(ord_res,"ordmeta.csv")

# wFisher, e.g. weight by NSNPS
wf_res <- magcat_wfisher_pathways(
  gene_results = genes_all,
  species      = "maize",
  gene_col     = "GENE",
  p_col        = "P",
  effect_col   = "ZSTAT",
  weight_col   = "NSNPS",
  is_onetail   = FALSE
)

write.csv(wf_res,"weighted_fisher.csv")


# TFisher
tf_res <- magcat_tfisher_pathways(
  gene_results    = genes_all,
  species         = "maize",
  gene_col        = "GENE",
  p_col           = "P",
  ptrunc          = 0.05,    # truncation threshold
  B_perm          = 5000L,   # permutations for empirical p
  seed            = 123,
  analytic_logical= TRUE,    # chi-square approx (optional)
  output          = TRUE,
  out_dir         = "magcat_tfisher"
)

head(tf_res)
write.csv(tf_res,"tf_res2.csv")

# OMNI
# Full pipeline from MAGMA gene results
omni_res <- magcat_omni_pathways(
  gene_results = magma_genes,   
  species      = "maize",       
  is_onetail   = FALSE,
  weight_col   = "NSNPS",
  B_ordmeta    = 2000,
  B_perm       = 1000,
  seed         = 123,
  output       = TRUE,
  out_dir      = "magcat_omni_full"
)

pw_tab <- read.csv("my_pathway_methods.csv") 
# columns: pathway_id, pathway_name, acat_p, wfisher_p, ordmeta_p, maybe others

omni_from_pvals <- magcat_omni_pathways(
  tab    = pw_tab,
  id_col = "pathway_id",                     # optional, just for sanity
  p_cols = c("acat_p", "wfisher_p", "ordmeta_p"),
  output = TRUE,
  out_dir = "magcat_omni_from_pvals"
)


args(magcat_omni_pathways)

# Directly from a file:
omni_from_pvals <- magcat_omni_pathways(
  tab    = "all_test.csv",
  id_col = "pathway_name",
  p_cols = c("acat_p", "wfisher_p","tpm_p_analytic")
)

write.csv(omni_from_pvals, "omni_from_pvals.csv")
