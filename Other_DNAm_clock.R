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

# 
pheno = data.frame(SampleID = dt$Sample_Name,
                   Age = dt$calendar_age,
                   Female = as.numeric(dt$sex) - 1)

mscore = methscore(datMeth = beta, datPheno = pheno, normalize = F)
write_xlsx(mscore, path = 'Other_DNAm_clock.xlsx')