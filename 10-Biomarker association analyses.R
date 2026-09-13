rm(list = ls())

library(tidyverse)
library(writexl)
library(rms)
library(survival)
library(survminer)
library(broom)
library(lubridate)
library(ggsurvfit)
library(ggplotify)
library(janitor)
library(data.table)
library(aplot)
library(patchwork)
library(rstatix)
library(janitor)
library(readxl)

# data --------------------------------------------------------------------

load('data.Rdata')
# This tibble named “data” contains all variables for 5,636 individuals, except for the β values of CpG sites.

dt = data

# function ----------------------------------------------------------------

mul_lm_fun = function(data, data_label, x, x_label, outcome, 
                      label_outcome, cov, cov_plus, cov_label){
  dt = data
  dt$x = data[[x]]
  dt$outcome = data[[outcome]]
  
  if(is.null(cov) & is.null(cov_plus)){
    fom = as.formula(str_c('outcome ~ x'))
  } else if(!is.null(cov) & is.null(cov_plus)){
    fom = as.formula(str_c('outcome ~ x', '+', str_c(cov, collapse = '+')))
  } else if(is.null(cov) & !is.null(cov_plus)) {
    fom = as.formula(str_c('outcome ~ x', '+', str_c(cov_plus, collapse = '+')))
  } else {
    fom = as.formula(str_c('outcome ~ x', '+', str_c(cov, collapse = '+'), '+', str_c(cov_plus, collapse = '+')))
  }
  
  fit = lm(fom, data = dt)
  con = confint(fit)
  
  tidy = tidy(fit) %>% 
    mutate(beta = estimate, up = con[, '97.5 %'], 
           low = con[, '2.5 %'], 
           `β (95% CI)` = str_c(str_trim(format(round(beta, 4), nsmall = 4)), ' (', 
                                str_trim(format(round(low, 4), nsmall = 4)), ', ', 
                                str_trim(format(round(up, 4), nsmall = 4)), ')'),
           p_value = ifelse(p.value < 0.001, '<0.001', str_trim(format(round(p.value, 3), nsmall = 3)))) %>% 
    rename(var = term) %>% 
    select(var, `β (95% CI)`, p_value, p.value, beta, low, up)
  if(class(dt$x) == 'factor'){
    tidy = tidy %>% 
      slice(2:(length(levels(dt$x)))) %>% 
      mutate(var = str_replace(var, '^x', ''))
    lm_dt = tibble(var = c(x_label, str_c(levels(dt$x)[1])), 
                   `β (95% CI)` = c(NA, '0 (Reference)'), p_value = NA, p.value = NA) %>% 
      bind_rows(tidy) %>% 
      mutate(data = data_label, N = nrow(dt), outcome = label_outcome, model = cov_label) %>% 
      rename(group = var) %>% 
      select(data, outcome, model, everything())
  }
  else {
    lm_dt = tidy %>% 
      slice(2) %>% 
      mutate(var = x_label) %>% 
      mutate(data = data_label, N = nrow(dt), outcome = label_outcome, model = cov_label) %>% 
      rename(group = var) %>% 
      select(data, N, outcome, model, everything())
  }
  
  return(lm_dt)
}


# 3-18 --------------------------------------------------------------------

dt1 = dt %>% 
  filter(calendar_age < 18)

dt1 = dt1 %>% 
  mutate(age_dev = case_when(age_dev > quantile(dt1$age_dev, 0.97) ~ quantile(dt1$age_dev, 0.97),
                             TRUE ~ age_dev))

outcomes1 = c("zHeight", "zWeight", "zBMI",
              "T3", "T4",
              "E2", "Te")
labels1 = c("zHeight", "zWeight", "zBMI",
            "T3", "T4", "Estradiol",
            "Testosterone")


x = c("age_dev")
x_label = c("Lifespan DNAmAgeDev")
cov0 = NULL
cov1 = c("sex")
cov2 = c("calendar_age", "sex")

cov_label = str_c('Model ', 0:2)


lm_map_dt1 = expand_grid(d1 = tibble(data = map(outcomes1, ~ drop_na(dt1, all_of(.x))),
                                     data_label = rep('Total', length(c(outcomes1))),
                                     cov_plus = map(1:length(outcomes1),
                                                    ~ c("CD8T", "CD4T", "NK", "Bcell", "Mono", "Neu")),
                                     outcome = c(outcomes1),
                                     label_outcome = c(labels1)),
                         tibble(x, x_label),
                         tibble(cov = list(cov0, cov1, cov2), cov_label)) %>%
  unnest(cols = c(d1)) %>%
  mutate(N = map_int(data, ~ nrow(.)))


mul_lm_dt1 = pmap_dfr(lm_map_dt1, ~ mul_lm_fun(data = ..1, data_label = ..2, cov_plus = ..3,
                                               x = ..6, x_label = ..7, outcome = ..4,
                                               label_outcome = ..5,
                                               cov = ..8, cov_label = ..9))

output1 = mul_lm_dt1 %>%
  group_by(model) %>% 
  mutate(FDR = p.adjust(p.value, method = 'BH'),
         FDR_g = ifelse(FDR < 0.001, '<0.001', str_trim(format(round(FDR, 3), nsmall = 3)))) %>% 
  ungroup() %>% 
  select(data, N, outcome, model, group, `β (95% CI)`, 
         p_value, FDR_g, beta, low, up, p.value, FDR) %>% 
  pivot_wider(id_cols = c(data, N, outcome, group),
              names_from = model, 
              values_from = c(`β (95% CI)`, p_value, FDR_g, beta, low, up, p.value, FDR),
              names_glue = '{model}_{.value}') %>% 
  relocate(data, N, outcome, group,
           contains('Model 0'),
           contains('Model 1'), contains('Model 2'),
           contains('Model 3'), contains('Model 4'))

output10 = mul_lm_dt1 %>%
  group_by(model) %>% 
  mutate(FDR = p.adjust(p.value, method = 'BH'),
         FDR_g = ifelse(FDR < 0.001, '<0.001', str_trim(format(round(FDR, 3), nsmall = 3)))) %>% 
  ungroup() %>% 
  select(data, N, outcome, model, group, `β (95% CI)`, 
         p_value, FDR_g, beta, low, up, p.value, FDR) %>% 
  pivot_wider(id_cols = c(data, N, outcome, group),
              names_from = model, 
              values_from = c(`β (95% CI)`, p_value, FDR_g),
              names_glue = '{model}_{.value}') %>% 
  relocate(data, N, outcome, group,
           contains('Model 0'),
           contains('Model 1'), contains('Model 2'),
           contains('Model 3'), contains('Model 4'))

write_xlsx(list(output1, output10), path = 'lm_3-18.xlsx')


# 18-65 -------------------------------------------------------------------

dt21 = filter(dt, calendar_age < 65 & calendar_age >= 18)

dt21 = dt21 %>% 
  mutate(age_dev = case_when(age_dev > quantile(dt21$age_dev, 0.99) ~ quantile(dt21$age_dev, 0.99),
                             age_dev < quantile(dt21$age_dev, 0.01) ~ quantile(dt21$age_dev, 0.01),
                             TRUE ~ age_dev))

outcomes2 = c("SBP", "DBP", 
              "TCHOL", "TRIG", "LDL_CH", "HDL_CH",
              "GLU")
labels2 = c("SBP", "DBP",
            "TC", "TG", "LDL‒C", "HDL‒C",
            "FBG")

x = c("age_dev")
x_label = c("Lifespan DNAmAgeDev")
cov0 = NULL
cov1 = c("calendar_age", "sex")
cov2 = c("calendar_age", "sex", "education", "marriage")
cov3 = c("calendar_age", "sex", "education", "marriage", 
         "BMI")
cov4 = c("calendar_age", "sex", "education", "marriage", 
         "BMI", "smoke", "hypertension")

cov_label = str_c('Model ', 0:4)


lm_map_dt21 = expand_grid(d1 = tibble(data = map(outcomes2, ~ drop_na(dt21, all_of(.x))), 
                                      data_label = rep('Total', length(c(outcomes2))),
                                      cov_plus = map(1:length(outcomes2), 
                                                     ~ c("CD8T", "CD4T", "NK", "Bcell", "Mono", "Neu")),
                                      outcome = c(outcomes2), 
                                      label_outcome = c(labels2)),
                          tibble(x, x_label),
                          tibble(cov = list(cov0, cov1, cov2, cov3, cov4), cov_label)) %>% 
  unnest(cols = c(d1)) %>% 
  mutate(N = map_int(data, ~ nrow(.)))

mul_lm_dt21 = pmap_dfr(lm_map_dt21, ~ mul_lm_fun(data = ..1, data_label = ..2, cov_plus = ..3,
                                                 x = ..6, x_label = ..7, outcome = ..4,
                                                 label_outcome = ..5,
                                                 cov = ..8, cov_label = ..9))

output21 = mul_lm_dt21 %>% 
  select(data, outcome, model, group, `β (95% CI)`, 
         p_value, beta, low, up, p.value) %>% 
  pivot_wider(id_cols = c(data, outcome, group),
              names_from = model, 
              values_from = c(`β (95% CI)`, p_value, beta, low, up, p.value),
              names_glue = '{model}_{.value}') %>% 
  relocate(data, outcome, group,
           contains('Model 0'),
           contains('Model 1'), contains('Model 2'),
           contains('Model 3'), contains('Model 4'))

output210 = mul_lm_dt21 %>% 
  select(data, outcome, model, group, `β (95% CI)`, 
         p_value) %>% 
  pivot_wider(id_cols = c(data, outcome, group),
              names_from = model, 
              values_from = c(`β (95% CI)`, p_value),
              names_glue = '{model}_{.value}') %>% 
  relocate(data, outcome, group,
           contains('Model 0'),
           contains('Model 1'), contains('Model 2'),
           contains('Model 3'), contains('Model 4'))

write_xlsx(list(output21, output210), path = 'lm_18-65.xlsx')


# 65-90 -------------------------------------------------------------------

dt22 = filter(dt, calendar_age >= 65 & calendar_age < 90)

dt22 = dt22 %>% 
  mutate(age_dev = case_when(age_dev > quantile(dt22$age_dev, 0.99) ~ quantile(dt22$age_dev, 0.99),
                             age_dev < quantile(dt22$age_dev, 0.01) ~ quantile(dt22$age_dev, 0.01),
                             TRUE ~ age_dev))

outcomes2 = c("BADL", "IADL",
              "MMSE",
              "CREA", "hsCRP")
labels2 = c("BADL", "IADL",
            "MMSE",
            "Plasma creatinine", "hs‒CRP")

x = c("age_dev")
x_label = c("Lifespan DNAmAgeDev")
cov0 = NULL
cov1 = c("calendar_age", "sex")
cov2 = c("calendar_age", "sex", "education", "marriage")
cov3 = c("calendar_age", "sex", "education", "marriage", 
         "BMI")
cov4 = c("calendar_age", "sex", "education", "marriage", 
         "BMI", "smoke", "hypertension")

cov_label = str_c('Model ', 0:4)


lm_map_dt22 = expand_grid(d1 = tibble(data = map(outcomes2, ~ drop_na(dt22, all_of(.x))), 
                                      data_label = rep('Total', length(c(outcomes2))),
                                      cov_plus = map(1:length(outcomes2), 
                                                     ~ c("CD8T", "CD4T", "NK", "Bcell", "Mono", "Neu")),
                                      outcome = c(outcomes2), 
                                      label_outcome = c(labels2)),
                          tibble(x, x_label),
                          tibble(cov = list(cov0, cov1, cov2, cov3, cov4), cov_label)) %>% 
  unnest(cols = c(d1)) %>% 
  mutate(N = map_int(data, ~ nrow(.)))


mul_lm_dt22 = pmap_dfr(lm_map_dt22, ~ mul_lm_fun(data = ..1, data_label = ..2, cov_plus = ..3,
                                                 x = ..6, x_label = ..7, outcome = ..4,
                                                 label_outcome = ..5,
                                                 cov = ..8, cov_label = ..9))

output22 = mul_lm_dt22 %>% 
  select(data, outcome, model, group, `β (95% CI)`, 
         p_value, beta, low, up, p.value) %>% 
  pivot_wider(id_cols = c(data, outcome, group),
              names_from = model, 
              values_from = c(`β (95% CI)`, p_value, beta, low, up, p.value),
              names_glue = '{model}_{.value}') %>% 
  relocate(data, outcome, group,
           contains('Model 0'),
           contains('Model 1'), contains('Model 2'),
           contains('Model 3'), contains('Model 4'))

output220 = mul_lm_dt22 %>% 
  select(data, outcome, model, group, `β (95% CI)`, 
         p_value) %>% 
  pivot_wider(id_cols = c(data, outcome, group),
              names_from = model, 
              values_from = c(`β (95% CI)`, p_value),
              names_glue = '{model}_{.value}') %>% 
  relocate(data, outcome, group,
           contains('Model 0'),
           contains('Model 1'), contains('Model 2'),
           contains('Model 3'), contains('Model 4'))

write_xlsx(list(output22, output220), path = 'lm_65-90.xlsx')


# 90-118 ------------------------------------------------------------------

dt3 = dt %>% 
  filter(calendar_age >= 90)

dt3 = dt3 %>% 
  mutate(age_dev = case_when(age_dev > quantile(dt3$age_dev, 0.99) ~ quantile(dt3$age_dev, 0.99),
                             age_dev < quantile(dt3$age_dev, 0.01) ~ quantile(dt3$age_dev, 0.01),
                             TRUE ~ age_dev))

outcomes3 = c("MMSE", "MMSE_f", "MMSE_dev", 
              "Lym", "NLR", "PLR", 
              "ALB", "ALT", "AST", "GGT", "TBIL", 
              "CREA", "Ucr", "UA",
              "HGS_L", "HGS_R", "s4", "s6", "FI")
labels3 = c("MMSE (baseline)", "MMSE (follow-up)", "MMSE (decline rate, per year)",
            "Lymphocyte", "NLR", "PLR",
            "ALB", "ALT", "AST", 
            "GGT", "TBIL", "Plasma creatinine", 
            "Urine creatinine", "UA",
            "Hand grip strength (left)", "Hand grip strength (right)",
            "4‒meter gait speed", "6‒meter gait speed", "FI")

x = c("age_dev")
x_label = c("Lifespan DNAmAgeDev")
cov0 = NULL
cov1 = c("calendar_age", "sex")
cov2 = c("calendar_age", "sex", "education", "marriage")
cov3 = c("calendar_age", "sex", "education", "marriage", 
         "BMI")
cov4 = c("calendar_age", "sex", "education", "marriage", 
         "BMI", "smoke", "hypertension")

cov_label = str_c('Model ', 0:4)

cov_plus = map(1:length(outcomes3), 
               ~ c("CD8T", "CD4T", "NK", "Bcell", "Mono", "Neu"))
cov_plus[[2]] = c("CD8T", "CD4T", "NK", "Bcell", "Mono", "Neu", 'MMSE')
cov_plus[[3]] = c("CD8T", "CD4T", "NK", "Bcell", "Mono", "Neu", 'MMSE')

lm_map_dt3 = expand_grid(d1 = tibble(data = map(outcomes3, ~ drop_na(dt3, all_of(.x))), 
                                     data_label = rep('Total', length(c(outcomes3))),
                                     cov_plus = cov_plus,
                                     outcome = c(outcomes3), 
                                     label_outcome = c(labels3)),
                         tibble(x, x_label),
                         tibble(cov = list(cov0, cov1, cov2, cov3, cov4), cov_label)) %>% 
  unnest(cols = c(d1)) %>% 
  mutate(N = map_int(data, ~ nrow(.)))


mul_lm_dt3 = pmap_dfr(lm_map_dt3, ~ mul_lm_fun(data = ..1, data_label = ..2, cov_plus = ..3,
                                               x = ..6, x_label = ..7, outcome = ..4,
                                               label_outcome = ..5,
                                               cov = ..8, cov_label = ..9))

output3 = mul_lm_dt3 %>% 
  select(data, outcome, N, model, group, `β (95% CI)`, 
         p_value, beta, low, up, p.value) %>% 
  pivot_wider(id_cols = c(data, outcome, N, group),
              names_from = model, 
              values_from = c(`β (95% CI)`, p_value, beta, low, up, p.value),
              names_glue = '{model}_{.value}') %>% 
  relocate(data, outcome, N, group,
           contains('Model 0'),
           contains('Model 1'), contains('Model 2'),
           contains('Model 3'), contains('Model 4'))

output30 = mul_lm_dt3 %>% 
  select(data, outcome, N, model, group, `β (95% CI)`, 
         p_value) %>% 
  pivot_wider(id_cols = c(data, outcome, N, group),
              names_from = model, 
              values_from = c(`β (95% CI)`, p_value),
              names_glue = '{model}_{.value}') %>% 
  relocate(data, outcome, N, group,
           contains('Model 0'),
           contains('Model 1'), contains('Model 2'),
           contains('Model 3'), contains('Model 4'))

write_xlsx(list(output3, output30), path = 'lm_90-118.xlsx')
