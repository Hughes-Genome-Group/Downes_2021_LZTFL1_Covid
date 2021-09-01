---
title: "EMT gene set analysis"
author: "Amy Cross"
date: "02/06/2021"
output: html_document
---

```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE)
knitr::opts_chunk$set(messages = FALSE)
knitr::opts_chunk$set(warnings = FALSE)

library(tidyr)
library(dplyr)
library(circlize)
library(ggplot2)
library(Hmisc)
library(corrplot)

set.seed(42)
setwd("D:/COVID Lung Characterisation Project 2020/DSP analysis/EMT gene set")
```

```{r input files, include = FALSE, echo = FALSE}
#----Dataframe annotations
celltype_names <- c("Activated DC", "AT1", "AT2", "Endothelial cells", "CD8+ T cells", "Dividing macrophages", "Dividing T cells", "Fibroblasts", "MARCO- macrophages", "MARCO+ macrophages", "Mast cells", "Monocytes", "Muscle cells", "NK cells") 

all_celltype_names <- c("Activated DC", "AT1", "AT2", "B cells", "Endothelial cells", "CD4+ T cells", "CD8+ T cells", "Ciliated cells", "DC type 1", "DC type 2", "DC dividing monocyte", "Dividing macrophages", "Dividing NK cells", "Dividing T cells", "Fibroblasts", "Lymph vessel cells", "MARCO- macrophages", "MARCO+ macrophages", "Mast cells", "Monocytes", "Muscle cells", "NK cells", "Plasma cells", "Plasmacytoid DC", "Regulatory T cell", "Neutrophil")

module_names <- c("IFITM2/HSP/ECM", "Apelin/mTOR signalling", "TLR signalling and monocytes", "Interferon responses", "TLR and IL-1 signalling", "Epithelial cells", "Type 2 pneumocytes", "IL-1 response: IL-6/IL-8", "Vasculature", "Cell cycling", "Cytotoxicity and T cells", "Antigen presentation", "Alveolar macrophages", "Macrophage identity", "Fibroblast phenotype", "Hypoxic response", "Chromatin remodeling")


#----Import files
#AOI annotation metadata
annot <- read.csv("D:/COVID Lung Characterisation Project 2020/DSP analysis/EMT gene set/EMT_input/COVID DSP annotation.csv", fileEncoding="UTF-8-BOM") 
annot <-annot[-4,] # remove cun4_004 because no sequencing of this AOI to create a 46x14 df.  AOI in rows and annotations in columns.

#Expression count data for each gene
filtQNcounts<- read.delim("D:/COVID Lung Characterisation Project 2020/DSP analysis/EMT gene set/EMT_input/qn.exprs.corrected.filtered.tsv")
filtQNcounts <- as.data.frame(t(filtQNcounts)) 
rownames(filtQNcounts) <- annot$AOI # 46x1631 df. AOI in rows and genes in columns.

#Module eigengene expression matrix
eigengenes <- read.delim("D:/COVID Lung Characterisation Project 2020/DSP analysis/EMT gene set/EMT_input/wgcna.covid.module.eigengenes.txt", row.names=1)
rownames(eigengenes) <- annot$AOI 
eigengenes$MEgrey <- NULL #Remove the grey module of unassigned genes. 
colnames(eigengenes) <- module_names #46x17 df.  AOI in rows and eigengene modules in cols. 

abundance <- read.csv("D:/COVID Lung Characterisation Project 2020/DSP analysis/EMT gene set/EMT_input/Abundances of lungneut cell types across covid QN data with meanbg.csv", fileEncoding="UTF-8-BOM", row.names = 1)
rownames(abundance) <- all_celltype_names
abundance <- as.data.frame(t(abundance)) # 46x26 df. AOI in rows and the abundance of cell types in columans. 

```

# Input files

annot #AOI annotation metadata

filtQNcounts #Expression count data for each gene

eigengenes #Module eigengene expression matrix

abundance #the relative abundance of cell types estimated using SpatialDecon R package.


## Objective

To explore whether the expression of a group of EMT-associated genes correlates with the presence of certain cell types and/or biological processes. 

## Selected EMT-realted genes derived from the literature

```{r select gene lists, echo=FALSE}
refined_genes <- c("CDH1","DSP","EPCAM","ACTA2","CDH2","AXL","COL1A1","ZEB2","SNAI1","SNAI2","TGFB1","TGFBR1","TGFBR2", "CTNNB1", "FZD6", "FZD7")

"Refined gene list: "
refined_genes
```


## Correlation of gene expression with WGCNA modules (from EMT analysis script 5)

```{r}
coul <- colorRampPalette(c("blue", "white", "red"))

refined_counts_in_covid <- select(filtQNcounts, any_of(refined_genes))

correlation_results <- rcorr(as.matrix(refined_counts_in_covid), as.matrix(eigengenes), type = "spearman") 
corr_p <- correlation_results$P 
corr_p <- corr_p[1:(ncol(corr_p)-17),(ncol(corr_p)-16):(ncol(corr_p))] #matrix of p values
corr_r <- correlation_results$r
corr_r <- corr_r[1:(ncol(corr_r)-17),(ncol(corr_r)-16):(ncol(corr_r))] #matrix of correlation values

refined_corrplot_modvgene <- corrplot(corr_r, tl.col = "black", tl.srt = 45, col = coul(100), p.mat = corr_p, sig.level = 0.05, insig =  "label_sig", tl.cex = 0.8, pch.cex = 1.2, method = "color")

```

```{r echo=FALSE, eval=FALSE, include=FALSE}
setwd("D:/COVID Lung Characterisation Project 2020/DSP analysis/EMT gene set")
write.csv(corr_r, file = "Spearman correlation Rho of EMT genes vs WGCNA modules.csv")
write.csv(corr_p, file = "Spearman correlation pvalue of EMT genes vs WGCNA modules.csv")
```


## Correlation of gene expression with abundances from cell deconvolution
```{r}
refined_abundance <- select(abundance, any_of(celltype_names))

correlation_results <- rcorr(as.matrix(refined_counts_in_covid), as.matrix(refined_abundance), type = "spearman")  
corr_p <- correlation_results$P 
corr_p <- corr_p[1:(ncol(corr_p)-14),(ncol(corr_p)-13):(ncol(corr_p))] #matrix of p values
corr_r <- correlation_results$r
corr_r <- corr_r[1:(ncol(corr_r)-14),(ncol(corr_r)-13):(ncol(corr_r))] #matrix of correlation values

corrplot_refinedgenes_abund_square <- corrplot(corr_r, tl.col = "black", tl.srt = 45, col = coul(100), p.mat = corr_p, sig.level = 0.05, insig =  "label_sig", tl.cex = 0.8, na.label = ".", method = "color", pch.cex = 1.2)

```

```{r echo=FALSE, eval=FALSE, include=FALSE}
setwd("D:/COVID Lung Characterisation Project 2020/DSP analysis/EMT gene set")
write.csv(corr_r, file = "Spearman correlation Rho of EMT genes vs cell abundance.csv")
write.csv(corr_p, file = "Spearman correlation pvalue of EMT genes vs cell abundance.csv")
```
