# library(tidyverse)
# library(biomaRt)
# library(parallel)

# FUN <- function(i){
#     df <- read_csv(paste("dfre_pd_",as.character(i),".csv",sep=""),show_col_types = FALSE)
#     rs_to_pos = getBM(attributes = c('refsnp_id', 'chr_name', 'chrom_start', 'chrom_end', 'allele'),
#                       filters = 'snp_filter',
#                       values = df$SNP,#rs_alleles %>% pull(rs_id),
#                       mart = snpmart) %>% as.tibble
#     dfb <- inner_join(rs_to_pos, df, by = c('refsnp_id' = 'SNP'))
#     return(dfb)
# }

# codes <- seq(1,11)
# snpmart <- useMart("ENSEMBL_MART_SNP",dataset = "hsapiens_snp",
#                    host = "www.ensembl.org",ensemblRedirect = FALSE)
# cl <- makeCluster(11)
# clusterExport(cl, varlist = 'snpmart') 
# clusterEvalQ(cl, expr = {library(tidyverse)
#                          library(biomaRt)})
# res <- clusterApply(cl = cl, x = codes, fun = FUN)

# dfall <- do.call(bind_rows,res)
# dfall            

library(tidyverse)
list_of_files <- list.files(pattern="dfre_pd_([1-9]|1[0-1]).csv",full.names=FALSE)
#list_of_files

df <- map_dfr(.x=set_names(list_of_files), .f=read_csv, .id="source_file")
#glimpse(df)

library(SNPlocs.Hsapiens.dbSNP144.GRCh38)
snps <- SNPlocs.Hsapiens.dbSNP144.GRCh38
chr_snps <- snpsBySeqname(snps, c(as.character(seq(1,22)),"X"))
chr_snps <- as_tibble(chr_snps)
dfb <- inner_join(chr_snps, df, by = c('RefSNP_id' = 'SNP'))
dfb
write_csv(dfb,"pd_manha.csv") 