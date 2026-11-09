rm(list = ls())

library(tidyverse)
library(ENmix)
library(minfi)
library(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
library(IlluminaHumanMethylationEPICv2manifest)
library(IlluminaHumanMethylationEPICv2anno.20a1.hg38)
library(readr)
library(readxl)
library(writexl)

load('data.Rdata')
# This tibble named “data” contains all variables for 5,636 individuals, except for the β values of CpG sites.
dt = data

beta = read_rds('DNAm.rds')
# This matrix named “beta” contains 889,074 rows of CpG sites, one “Sample_Name” column, and 5,636 columns of individuals.

# Horvath, Hannum, PhenoAge ----
pheno = data.frame(SampleID = dt$Sample_Name,
                   Age = dt$calendar_age,
                   Female = as.numeric(dt$sex) - 1)

mscore = methscore(datMeth = beta, datPheno = pheno, normalize = F)
write_xlsx(mscore, path = 'Other_DNAm_clock.xlsx')


# GrimAge, DunedinPACE, Systems Age ----
library(methylCIPHERv2)

beta_t = t(beta)

pheno2 = tibble(Sample_Name = rownames(beta_t)) %>% 
  left_join(select(dt, Sample_Name, calendar_age, sex), by = 'Sample_Name') %>% 
  rename(ID = Sample_Name, Age = calendar_age) %>% 
  mutate(Female = ifelse(sex == 2, 1, 0)) %>% 
  select(-sex) %>% 
  as.data.frame()

res = calc_clocks(
  DNAm = beta_t,
  clocks = c("GrimAgeV1", "GrimAgeV2", "DunedinPACE", "SystemsAge"),
  pheno = pheno2,
  min_clocks_coverage = 0.7,
  min_samples_coverage = 0.7
)

cover = clocks_coverage(res, all_columns = TRUE) %>% 
  filter(clock_id %in% c("GrimAgeV2", 'DunedinPACE', 'Brain', 'Heart', 'Immune', 'Kidney', 'Liver', 'MusculoSkeletal', 'SystemsAge'))

sample = samples_coverage(res) %>% 
  filter(clock_id %in% c("GrimAgeV2", 'DunedinPACE', 'Brain', 'Heart', 'Immune', 'Kidney', 'Liver', 'MusculoSkeletal', 'SystemsAge'))

clock_SystemsAge = res$scores %>% 
  as_tibble(rownames = 'Sample_Name') %>% 
  select(all_of(c("Sample_Name", "GrimAgeV2", 'DunedinPACE', 'Brain', 'Heart', 'Immune', 'Kidney', 'Liver', 'MusculoSkeletal', 'SystemsAge')))

save(clock_SystemsAge, file = 'clock_SystemsAge.Rdata')
write_xlsx(cover, path = 'coverage.xlsx')
