rm(list = ls())

library(tidyverse)
library(readr)
library(broom)
library(writexl)
library(tictoc)
library(future)
library(furrr)
library(purrr)

# load --------------------------------------------------------------------

load('data.Rdata')
# This tibble named “data” contains all variables for 5,636 individuals, except for the β values of CpG sites.
dt = data

beta = read_rds('DNAm.rds')
# This matrix named “beta” contains 889,074 rows of CpG sites, one “Sample_Name” column, and 5,636 columns of individuals.

#
lm_fun_optimized = function(cpg_row, cpg_name, model_data){
  # cpg_row: beta value
  # cpg_name: CpG name
  # model_data: covariates
  
  df = bind_cols(cpg = cpg_row, model_data)
  
  fit = lm(
    cpg ~ calendar_age + sex + CD8T + CD4T + NK + Bcell + Mono + Neu,
    data = df
  )
  
  out = broom::tidy(fit)[2, ] %>%
    mutate(CpG = cpg_name) %>%
    relocate(CpG)
  
  return(out)
}

model_cols = c("calendar_age", "sex", "CD8T", "CD4T", "NK", "Bcell", "Mono", "Neu")
model_data = dt[, model_cols]

cpg_names = rownames(beta)

chunk_size = 100
n_total = nrow(beta)
chunks = split(1:n_total, ceiling(seq_along(1:n_total)/chunk_size))

tic()

results = map_dfr(chunks, function(chunk){
  
  beta_chunk = beta[chunk, , drop = FALSE]
  names_chunk = cpg_names[chunk]
  
  future_map2_dfr(
    .x = split(beta_chunk, seq_len(nrow(beta_chunk))),
    .y = names_chunk,
    ~ lm_fun_optimized(
      cpg_row = .x,
      cpg_name = .y,
      model_data = model_data
    )
  )
})

toc()

# FDR
out = results %>%
  mutate(FDR = p.adjust(p.value, method = "BH"))

DMP = out %>% 
  filter(abs(estimate) > 0.002 & FDR < 0.01)

beta_DMP = beta[rownames(beta) %in% DMP$CpG, ]
name = colnames(beta_DMP)

DMPs = beta_DMP %>% 
  t() %>% 
  as_tibble() %>% 
  mutate(Sample_Name = name) %>% 
  relocate(Sample_Name)

# save --------------------------------------------------------------------
write_rds(DMPs, file = 'DMPs.Rdata')
# This tibble named “DMPs” contains 5,636 individuals, a “Sample_Name” column, and 3,237 CpG columns.
