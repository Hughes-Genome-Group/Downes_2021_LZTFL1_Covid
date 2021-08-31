# visualize deephaem predictions over tissues

library(tidyverse)
library(RColorBrewer)
library(cowplot)
theme_set(theme_cowplot())

# set working directory
out.dir <- "~/fusessh/other_projects/covid_damien/"
setwd(out.dir)

### 4k deepHaem model predicitons ### ==================================
# load prediction and labels | format and save df -----------------
ref.in <- "deepHaem_4k_out/reference_class_scores_deepHaem_4k_predict.bed"
var.in <- "deepHaem_4k_out/variant_class_scores_deepHaem_4k_predict.bed"
logfold.in <- "deepHaem_4k_out/logfold_damage_scores_deepHaem_4k_predict.bed"
total.in <- "deepHaem_4k_out/total_damage_scores_deepHaem_4k_predict.bed"

labels.in <- "labels_full_chrom.txt"

labels <- as_tibble(read.table(labels.in, header = FALSE, colClasses = c("numeric", "character")))
names(labels) <- c("id", "labels")


ref <- as_tibble(read.table(ref.in, header = FALSE))
names(ref) <- c("chr", "pos", "snptag", "ref", "var", labels$labels)
ref <- ref %>% gather(class, value, -c(chr, pos, snptag, ref, var))
ref <- mutate(ref, out = "ref")

var <- as_tibble(read.table(var.in, header = FALSE))
names(var) <- c("chr", "pos", "snptag", "ref", "var", labels$labels)
var <- var %>% gather(class, value, -c(chr, pos, snptag, ref, var))
var <- mutate(var, out = "var")

total <- as_tibble(read.table(total.in, header = FALSE))
names(total) <- c("chr", "pos", "snptag", "ref", "var", labels$labels)
total <- total %>% gather(class, value, -c(chr, pos, snptag, ref, var))
total <- mutate(total, out = "total")

logfold <- as_tibble(read.table(logfold.in, header = FALSE))
names(logfold) <- c("chr", "pos", "snptag", "ref", "var", labels$labels)
logfold <- logfold %>% gather(class, value, -c(chr, pos, snptag, ref, var))
logfold <- mutate(logfold, out = "logfold")

# bind all together into single dataframe
df <- rbind(ref, var, total, logfold)
# remove extra columns 
df <- df %>% select(class, snptag, value, out)

# load detailed labels 
detail.dnase <- as_tibble(read.table('~/fusessh/machine_learning/deepHaem/data/assemble_basenji_compendium/download_data/encode_parsed_dnase.combined.peakcallcheck.txt'))
detail.dnase <- detail.dnase[,c(1:3)]
names(detail.dnase) <- c('desc', 'type', 'id')
detail.dnase$target <- 'open.chrom'
detail.dnase <- detail.dnase %>% select(desc, type, target, id)

detail.chip <- as_tibble(read.table('~/fusessh/machine_learning/deepHaem/data/assemble_basenji_compendium/download_data/encode_parsed_chip_combined.formatted.sorted.available.txt'))
detail.chip <- detail.chip[,c(1:4)]
detail.chip <- detail.chip %>% mutate(V3 = str_replace(V3, '/targets/', '')) %>% mutate(V3 = str_replace(V3, '/', ''))
names(detail.chip) <- c('desc', 'type', 'target', 'id')

details <- rbind(detail.dnase, detail.chip)

# join
df.details <- left_join(df, details, by = c('class' = 'id'))

# check NA example
df.details %>% filter(desc == "fibroblast_of_dermis_female_adult" & snptag == "chr3_45848457_C_T")

df.details <- df.details %>% mutate(desc = as.character(desc), type = as.character(type))

# set for corces and downes
# df.details.sub <- df.details %>% filter(is.na(desc))
df.details.sub <- df.details

df.details.sub <- df.details.sub %>% 
  mutate(desc = if_else(is.na(type), class, desc)) %>%
  mutate(type = case_when(!is.na(type) ~ type,
                         str_detect(class, "CTCF") ~ "TF_ChIP-seq",
                         str_detect(class, "ATAC") ~ "ATAC-seq",
                         str_detect(class, "Don001") ~ "Histone_ChIP-seq",
                         TRUE ~ "ATAC-seq")) %>%
  mutate(target = case_when(!is.na(target) ~ target,
                          str_detect(class, "CTCF") ~ "CTCF",
                          str_detect(class, "ATAC") ~ "open.chrom",
                          str_detect(class, "Don001") ~ "histone.mod",
                          TRUE ~ "open.chrom"))

# df.details[c(1:1536),] <- df.details.sub
df.details <- df.details.sub
df <- df.details 

saveRDS(df , file = "df_dhaem_4k.rds")

# Continue from here =====================================================
df <- readRDS(file = "df_dhaem_4k.rds")

df <- mutate(df, out = factor(out, levels = c("ref", "var", "total", "logfold")))

# RANK DNase/ATAC hits ----------------------------------------------------
odf <- df %>%
  filter(type %in% c("DNase-seq", "ATAC-seq")) %>%
  spread(out, value) %>%
  arrange(desc(abs(total)))

# add rank
odf <- odf %>% 
  arrange(desc(abs(total))) %>%
  mutate(rank = row_number())

odf

# flag tissues of interest 
odf$tissue.of.interest <- "other"

odf <- odf %>% 
  mutate(tissue.of.interest = case_when(str_detect(desc, 'lung') ~ 'lung',
                         str_detect(desc, 'muscle') ~ 'muscle',
                         str_detect(desc, 'fibroblast') ~ 'fibroblast',
                         TRUE ~ 'other'))

# clor by (ORGAN == BLOOD)
blood.desc <- as_tibble(read.table("../covid_snps/tissue_dnase/organ_blood_descriptions.txt", sep= "\t"))
blood.desc$V1 <- as.character(blood.desc$V1)

bodf <- odf %>% mutate(tissue.of.interest = as.character(tissue.of.interest))

for(i in c(1:38)){
 bodf <- bodf %>% mutate(tissue.of.interest = if_else(str_detect(desc, blood.desc$V1[i]), "blood", tissue.of.interest))  
}



# set factor level for plotting
bodf <- bodf %>% mutate(tissue.of.interest =  factor(tissue.of.interest, 
                                            levels = rev(c("lung", "muscle", "fibroblast", "blood", "other"))))

bodf %>% filter(tissue.of.interest == "blood")


# p.scale.free <- odf %>%
#   ggplot(aes(x = rank, y = total)) +
#   geom_point(size=0.5) +
#   facet_wrap(~snptag, scales = "free_y")
# 
# p.scale.free
# 
# ggsave(p.scale.free, filename = "plot_open_chrom_ranking_deepHaem_4k_free_scale.png")

p.scale.fix <- bodf %>%
  ggplot(aes(x = rank, y = total, col = tissue.of.interest)) +
  geom_point(size=1) +
  scale_color_manual(values=rev(brewer.pal(9, "Set1")[c(2,4,3,1,9)])) +
  facet_wrap(~snptag) + 
  theme(strip.background = element_blank())

p.scale.fix

ggsave(p.scale.fix, filename = "plot_open_chrom_ranking_deepHaem_4k_fixed_scale.png")
ggsave(p.scale.fix, filename = "plot_open_chrom_ranking_deepHaem_4k_fixed_scale.pdf")

# save table 
write.table(odf, file = "table_open_chromaint_classifiers_sorted_deephaem_4k.txt", quote = F, row.names = F, col.names = T, sep = "\t")

# TF Ranking -------------------------------------------------------- 
tdf <- df %>%
  filter(type %in% c("TF_ChIP-seq"))

# weird duplicates
tdf <- tdf %>% pivot_wider(names_from = out, values_from = value, values_fn = list(value = mean))

# add rank
tdf <- tdf %>% 
  arrange(desc(abs(total))) %>%
  mutate(rank = row_number())

p.scale.free <- tdf %>%
  ggplot(aes(x = rank, y = total)) +
  geom_point(size=0.5) +
  facet_wrap(~snptag, scales = "free_y")

p.scale.free

ggsave(p.scale.free, filename = "plot_tfb_ranking_deepHaem_4k_free_scale.png")

p.scale.fix <- tdf %>%
  ggplot(aes(x = rank, y = total)) +
  geom_point(size=0.5) +
  facet_wrap(~snptag) + 
  theme(strip.background = element_blank())

p.scale.fix

ggsave(p.scale.fix, filename = "plot_tfb_ranking_deepHaem_4k_fixed_scale.pdf")

# save table 
write.table(tdf, file = "table_tfb_classifiers_sorted_deephaem_4k.txt", quote = F, row.names = F, col.names = T, sep = "\t")


