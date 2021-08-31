# coloc analysis for COVID GWAS REVISIONS 
library(coloc)
library(tidyverse)
library(cowplot)
theme_set(theme_cowplot())
options(scipen=3)

out.dir <- '~/fusessh/other_projects/covid_damien/coloc/'
setwd(out.dir)

input1  <- 'gwas_locus_200kb_window_sum_stats.txt'
input2  <- 'Slc6a20_pairs.txt'
output  <- '~/fusessh/other_projects/covid_damien/coloc/coloc_out/Slc6a20_coloc_200kb_p12_1e-05'
N1 <- 3795
N2 <- 515
s <- 0.419

a <- read_tsv(input1, col_names = T)
b <- read_tsv(input2, col_names = T)

# b <- b %>% filter(abs(tss_distance) <= 100000)

my.res <- coloc.abf(dataset1=list(beta=a$beta, varbeta=(a$standard_error)^2, 
                                  N=N1,
                                  type="cc",
                                  s=s,
                                  snp = paste("chr",a$hm_variant_id,"_b38",sep="")),
                    dataset2=list(beta=b$slope, 
                                  varbeta=(b$slope_se)^2, 
                                  N=N2,
                                  MAF = b$maf ,
                                  type="quant",
                                  snp = b$variant_id), 
                    p1 = 1e-04, p2 = 1e-04, p12 = 1e-05)


pdf("./coloc_out/coloc_Slc6a20_sensitivity.pdf", height = 9, width = 11)
  sensitivity(my.res,rule="H4 > H3")
dev.off()

write.table(t(my.res$summary), file=paste(output,".summary", sep=""), col.names=T, row.names=F, quote=F, append=F)
write.table(my.res$results, file=paste(output,".results", sep=""), col.names=T, row.names=F, quote=F, append=F)

resdf <- as_tibble(my.res$results)
resdf <- resdf %>% arrange(desc(SNP.PP.H4)) %>% mutate(SNP.PP.H4.rank = row_number())

resdf %>% select(snp, SNP.PP.H4, SNP.PP.H4.rank)

cred.set <- resdf %>% select(snp, SNP.PP.H4, SNP.PP.H4.rank) %>% filter(SNP.PP.H4 > 0.01)
write.table(cred.set, file=paste(output,".credible_set", sep=""), col.names=T, row.names=F, quote=F, append=F)



# For Lztfl1 ----------
# pp12 default 1e-05
# Warning that lowest P.value for eqtl dataset is p = 0.0025575
# PP.H0.abf  PP.H1.abf  PP.H2.abf  PP.H3.abf  PP.H4.abf 
# 0.00004380 0.62900000 0.00000872 0.12500000 0.24600000 
# PP of SNp being causal conditional on H4 being true is > 0.05 whichs is ranked 7 of 2845
# chr3_45818159_G_A_b38    0.0541              7

# pp12 1e-04 > a priori belief that H4 is > H32
# PP.H0.abf  PP.H1.abf  PP.H2.abf  PP.H3.abf  PP.H4.abf 
# 0.00001360 0.19500000 0.00000271 0.03880000 0.76600000 



# For Slc6a20 ----------
# Warning that lowest P.value for eqtl dataset is p = 0.0001
# pp12 default 1e-05
# PP.H0.abf  PP.H1.abf  PP.H2.abf  PP.H3.abf  PP.H4.abf 
# 0.00004010 0.57600000 0.00000953 0.13700000 0.28700000 
# PP of SNP being causal conditional on H4 being true is > 0.05 whichs is ranked 7 of 2845
# chr3_45818159_G_A_b38    0.0539              7

# pp12 default 1e-04
# PP.H0.abf  PP.H1.abf  PP.H2.abf  PP.H3.abf  PP.H4.abf 
# 0.00001120 0.16100000 0.00000266 0.03810000 0.80100000 

# 𝐻0
# : neither trait has a genetic association in the region
# 𝐻1
# : only trait 1 has a genetic association in the region
# 𝐻2
# : only trait 2 has a genetic association in the region
# 𝐻3
# : both traits are associated, but with different causal variants
# 𝐻4
# : both traits are associated and share a single causal variant


# Locus Compare
# install.packages('devtools')
# devtools::install_github("boxiangliu/locuscomparer")
# library(locuscomparer)
# gwas_fn = system.file('extdata','gwas.tsv', package = 'locuscomparer')
# eqtl_fn = system.file('extdata','eqtl.tsv', package = 'locuscomparer')

ca <- a %>% mutate(id = paste("chr",hm_variant_id,"_b38",sep=""), pval.gwas = p_value) %>% select(id, pval.gwas)
cb <- b %>% mutate(id = variant_id, pval.eqtl = pval_nominal) %>% select(id, pval.eqtl)

cc <- left_join(ca, cb)

cc %>%
  ggplot(aes(x=-log10(pval.gwas), y=-log10(pval.eqtl))) +
  geom_point()



