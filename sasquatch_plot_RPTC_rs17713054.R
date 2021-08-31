# 2018/11/30 Plot Selected Kmer comparisons
# Ron Schwessinger

### Works with R/3.5.0-newgcc

library(tidyverse)
library(svglite)

source("/Sasquatch/R_utility/functions_sasq_r_utility.R")

data.dir <- "/database_assembly/idx_correct_assembly/human/DNase/"
pnorm.tag <- "h_ery_1" #identifier for the propensity source used
tissue <- "ENCODE_UW_lung_fetal_merged"
frag.type <- "DNase"
vocab <- PreLoadVocab(data.dir, tissue)
profiles <- PreLoadKmerProfiles(kl = 7, data.dir = data.dir, tissue= tissue, pnorm.tag = pnorm.tag)

# SET YOUR OUTPUT DIRECTORY ---------------------------------------------------------------------
# output directory for tables and plots  
out.dir <- "/covid/3_sasQ/3_plots/"
setwd(out.dir)

# Wrapper for overlap plots from k-mers -----------------------------------------------------------
p1 <- PlotOverlapKmers(
  kmer1="GCAATAC", kmer2="GCAGTAC",
  ymode="merged",
  tissue1=tissue, 
  tissue2=tissue, 
  data.dir=data.dir, pnorm.tag=pnorm.tag, frag.type=frag.type, 
  smooth=TRUE, plot.shoulders = FALSE,
  ylim=c(0,0.01), xlim=c(-75,75),
  preload.profiles = profiles
)

plot(p1)



ggsave(p1, file="rs17713054.svg", width = 8, height = 6)


