rm(list = ls())

library(tidyverse)
library(fs)
library(writexl)
library(readr)
library(readxl)

# GTEx organ --------------------------------------------------------------

dir = dir_ls('data/GTEx/tissue/') %>% 
  as.character()
# The “tissue” folder contains 54 publicly available GTEx v10 tissue datasets.
names = dir %>% 
  str_replace('data/GTEx/tissue/gene_reads_v10_', '') %>% 
  str_replace('.gct.gz', '')
GTEx = tibble(dir = dir, tissue = names)
write_xlsx(GTEx, path = 'data/GTEx/GTEx.xlsx')

organ = read_excel('data/GTEx/GTEx_organ.xlsx') %>% 
  pull(organ)
# The file “GTEx_organ.xlsx” is a table of tissue-to-organ mappings provided by Oh et al., which corresponds to “Supplementary Table 22. Tissue to organ mapping in GTEx (v10)” in this paper.

for(i in 31:length(names)){
  assign(names[i], read.table(dir[i], skip = 2, header = TRUE, sep = "\t") %>% 
           mutate(organ = organ[i],
                  genes = str_replace(Name, '\\..*', '')) %>% 
           relocate(Name, genes, organ))
}

GTEx_v10 = tibble(sample_id = character(0), 
                  donor_id = character(0),
                  tissue_id = character(0),
                  organ = character(0))
for(i in 1:length(names)){
  dti = get(names[i])
  i_sample_id = colnames(dti)[5:ncol(dti)]
  i_donor_id = str_extract(i_sample_id, '^.*(?=\\..*\\..*\\..*)')
  i_tissue_id = names[i]
  i_organ = organ[i]
  out = tibble(sample_id = i_sample_id, 
               donor_id = i_donor_id,
               tissue_id = i_tissue_id,
               organ = i_organ)
  GTEx_v10 = bind_rows(GTEx_v10, out)
}
write_xlsx(GTEx_v10, path = 'data/GTEx/GTEx_Analysis_v10_RNAseq_samples.xlsx')


# DESeq2 ------------------------------------------------------------------

rm(list = ls())
gc()

GTEx_data = read.table('data/GTEx/GTEx_Analysis_2022-06-06_v10_RNASeQCv2.4.2_gene_reads.gct', 
                       skip = 2, header = TRUE, sep = "\t")

GTEx_data = GTEx_data %>% 
  mutate(gene = str_replace(Name, '\\..*', '')) %>% 
  relocate(Name, gene)

GTEx_data = GTEx_data %>% 
  filter(!str_detect(Name, '_PAR_Y'))

#
library(DESeq2)

gtex_data = GTEx_data %>% 
  select(-c(Name, gene, Description))
rownames(gtex_data) = GTEx_data$gene

GTEx_v10 = read_excel('data/GTEx/GTEx_Analysis_v10_RNAseq_samples.xlsx')
sum(!GTEx_v10$sample_id %in% colnames(gtex_data))

gtex_data = GTEx_data %>% 
  select(all_of(GTEx_v10$sample_id))
rownames(gtex_data) = GTEx_data$gene

sum(colnames(gtex_data) != GTEx_v10$sample_id)

sample_annotation = GTEx_v10 %>% 
  mutate(tissue_type = factor(tissue_id)) %>% 
  select(tissue_type) %>% 
  as.data.frame() 
rownames(sample_annotation) = GTEx_v10$sample_id

# 
dds = DESeqDataSetFromMatrix(
  countData = gtex_data,
  colData = sample_annotation,
  design = ~ tissue_type
)
dds = estimateSizeFactors(dds)
normalized_counts = counts(dds, normalized = TRUE)


# Identifying enriched genes ----------------------------------------------

organ_groups = select(GTEx_v10, sample_id, organ)

organ_expression = normalized_counts %>%
  as.data.frame() %>%
  mutate(gene = rownames(.))
rownames(organ_expression) = NULL

organ_expression = organ_expression %>%
  pivot_longer(cols = -gene, names_to = "sample_id", values_to = "expression") %>%
  left_join(organ_groups, by = "sample_id")

count(organ_expression, is.na(organ))
organ_expression_rm = drop_na(organ_expression, organ)

organ_expression_out = organ_expression_rm %>%
  group_by(gene, organ) %>%
  summarise(max_expression = max(expression)) %>% 
  ungroup()

#
organ_enriched_genes = organ_expression_out %>%
  group_by(gene) %>%
  mutate(
    top_organ = organ[which.max(max_expression)],
    top_expression = max(max_expression),
    second_expression = sort(max_expression, decreasing = TRUE)[2],
    enrichment_ratio = ifelse(second_expression == 0, NA, top_expression / second_expression)
  ) %>%
  filter(top_expression != 0 & max_expression == top_expression)

organ_enriched_genes = organ_enriched_genes %>%
  mutate(GTExOrgan = organ, 
         GTExFoldChange = enrichment_ratio,
         GTExEnriched = enrichment_ratio >= 4) %>% 
  select(gene, GTExOrgan, GTExFoldChange, GTExEnriched)

nrow(distinct(organ_enriched_genes, gene)) == nrow(organ_enriched_genes)

save(organ_enriched_genes, file = 'data/GTEx/organ_enriched_genes.Rdata')


