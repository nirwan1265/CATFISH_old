devtools::document()   
devtools::load_all()


library(MAGCAT)


# gff3_to_geneloc
#gff_path <- "/Users/nirwantandukar/Downloads/Zm-B73-REFERENCE-NAM-5.0_Zm00001eb.1.gff3"

#gene_loc_out <- "/Users/nirwantandukar/Downloads/maize.genes.loc"

#gff3_to_geneloc(
#  gff        = "/Users/nirwantandukar/Downloads/Zm-B73-REFERENCE-NAM-5.0_Zm00001eb.1.gff3",
#  out        = "inst/extdata/maize.genes.loc",
#  chr_prefix = "chr"   # default; strips "chr"
#)

# Gene-loc function
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
  gene_model = c("snp-wise=top")
)


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
  out_dir    = "magma_genes_by_chr",
  gene_model = c("snp-wise=top"),
  chroms     = 1:10,   
  n_threads  = 10       # use up to 10 parallel workers
)


# ACAT
# 1) MAGMA gene-level results (your merged file, with GENE / P / CHR etc.)
genes_all <- read.table("/Users/nirwantandukar/Documents/Research/results/MAGMA/MAGCAT/N_maize_magma_genes_by_chr/N_maize_MLM_chr1.snp_wise_top.genes.out",
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
attr(pw_res, "file")  # path to the CSV we just wrote

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
