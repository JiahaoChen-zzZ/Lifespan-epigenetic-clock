rm(list = ls())

library(tidyverse)
library(writexl)
library(readr)
library(IlluminaHumanMethylationEPICv2anno.20a1.hg38)

# data --------------------------------------------------------------------

load('data/GTEx/organ_enriched_genes.Rdata')
load('DMPs.Rdata')
# This tibble named “DMPs” contains 5,636 individuals, a “Sample_Name” column, and 3,237 CpG columns.

DMP = colnames(dt_beta_DMP)[-1]

v2 = read_rds('EPICv2_reannotated_manifest_v1.0.rds')
# This data is the publicly available, re-annotated EPIC v2.0 annotation manifest mentioned in this article.

dt_anno = distinct(v2, Name, GENCODEv47_Gene_ID)
dtc_anno = tibble(CpG = DMP) %>% 
  left_join(dt_anno, by = c('CpG' = 'Name'))

anno_raw = getAnnotation(IlluminaHumanMethylationEPICv2anno.20a1.hg38)
dt_anno2 = anno_raw@listData %>% 
  as.data.frame() %>% 
  mutate(Name = str_extract(Name, '^.*(?=_)')) %>% 
  distinct(Name, UCSC_RefGene_Name)
dtc_anno2 = tibble(CpG = DMP) %>% 
  left_join(dt_anno2, by = c('CpG' = 'Name'))

anno = dtc_anno %>% 
  left_join(dtc_anno2, by = 'CpG') %>% 
  mutate(id1 = str_extract(GENCODEv47_Gene_ID, '^.*?(?=;)'),
         id2 = str_extract(UCSC_RefGene_Name, '^.*?(?=;)'))

library(org.Hs.eg.db)
library(AnnotationDbi)

ID1 = mapIds(org.Hs.eg.db,
             keys = anno$id1,
             column = "SYMBOL",
             keytype = "ENSEMBL",
             multiVals = "first")

ID2 = mapIds(org.Hs.eg.db,
             keys = anno$id2,
             column = "ENSEMBL",
             keytype = "SYMBOL",
             multiVals = "first")

anno = anno %>% 
  mutate(ID1 = as.character(ID1),
         ID2 = as.character(ID2),
         ID1 = ifelse(ID1 == 'NULL', NA, ID1),
         ID2 = ifelse(ID2 == 'NULL', NA, ID2)) %>% 
  mutate(`Ensembl ID` = ifelse(is.na(id1), ID2, id1))

ID3 = mapIds(org.Hs.eg.db,
             keys = anno$`Ensembl ID`,
             column = "SYMBOL",
             keytype = "ENSEMBL",
             multiVals = "first")
anno = anno %>% 
  mutate(`Gene Symbol` = as.character(ID3),
         `Gene Symbol` = ifelse(`Gene Symbol` == 'NULL', NA, `Gene Symbol`)) %>% 
  select(CpG, `Ensembl ID`, `Gene Symbol`)

CpG_organ = anno %>% 
  left_join(organ_enriched_genes, by = c('Ensembl ID' = 'gene'))

write_xlsx(CpG_organ, path = 'CpG_organ.xlsx')
# The data recorded in “Supplementary Table 23. Details of DMPs with Annotations of Organ Enrichment”
