library(genepi.utils)
library(data.table)
dat <- fread("/Users/xx20081/git/dcm_pgs/data/gnomad_pca_variants_chr_pos_ref_alt.tsv")
dat <- chrpos_to_rsid(
  dat,
  chr_col   = "chr",
  pos_col   = "pos",
  ea_col    = "alt",
  nea_col   = "ref",
  flip      = "allow",
  build     = "b37_dbsnp156",
  dbsnp_dir = "/Users/xx20081/Documents/local_data/genome_reference/dbsnp",
  parallel_cores = 10
)
fwrite(dat[,"RSID"], "/Users/xx20081/git/dcm_pgs/data/gnomad_pca_variants_rsid.txt", col.names = FALSE)