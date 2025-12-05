devtools::document()   # if you just added the function
devtools::load_all()   


library(MAGCAT)

gff3_to_geneloc

gff_path <- "/Users/nirwantandukar/Downloads/Zm-B73-REFERENCE-NAM-5.0_Zm00001eb.1.gff3"

gene_loc_out <- "/Users/nirwantandukar/Downloads/maize.genes.loc"

res <- gff3_to_geneloc(
  gff = gff_path,
  out = gene_loc_out
)
