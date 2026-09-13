rm(list = ls())

library(tidyverse)
library(tidyplots)
library(writexl)
library(patchwork)
library(survival)
library(survminer)
library(broom)
library(lubridate)
library(ggsurvfit)
library(ggplotify)
library(lemon)
library(janitor)
library(aplot)
library(readxl)
library(cowplot)
library(ggupset)
library(viridis)
library(ggcor)
library(haven)
library(foreign)
library(stringi)

# data --------------------------------------------------------------------

load('data.Rdata')
# This tibble named “data” contains all variables for 5,636 individuals, except for the β values of CpG sites.
dt = data

# Cox function ------------------------------------------------------------

rate_fun = function(data, data_label, outcome, outtime, group, label, group_label){
  dt = data
  dt$y = data[[outcome]]
  dt$time = data[[outtime]]
  
  if(str_detect(group, 'Total')){
    event_dt = count(dt, y) %>%
      mutate(N = sum(n), event_rate = round(n/N*100, 2), outcome = label,
             `Event/Total (%)` = str_c(n, '/', N, ' (', str_trim(format(event_rate, nsmall = 2)), ')')) %>% 
      filter(y == 1) %>%
      mutate(group = group_label) %>% 
      select(outcome, group, `Event/Total (%)`)
  } else {
    dt$group = data[[group]]
    
    event_dt = count(dt, group, y) %>%
      group_by(group) %>% 
      mutate(N = sum(n), event_rate = round(n/N*100, 2), outcome = label,
             `Event/Total (%)` = str_c(n, '/', N, ' (', str_trim(format(event_rate, nsmall = 2)), ')')) %>% 
      filter(y == 1) %>%
      ungroup %>% 
      select(outcome, group, `Event/Total (%)`)
    if(nrow(event_dt) != length(levels(dt$group))){
      event_dt = event_dt %>% 
        bind_rows(map_dfr(rep(1, length(levels(dt$group)) - nrow(event_dt)), ~ tibble(NA))) %>% 
        select(-`NA`)
    }
  }
  event_dt = event_dt %>% 
    mutate(data = data_label) %>% 
    select(data, everything())
  
  return(event_dt)
}
mul_cox_fun = function(data, data_label, x, x_label, outcome, outtime, label_outcome, cov, cov_plus, cov_label, weight){
  dt = data
  dt$x = data[[x]]
  dt$outcome = data[[outcome]]
  dt$outtime = data[[outtime]]
  
  Surv_obj = Surv(dt$outtime, dt$outcome)
  if(is.null(cov) & is.null(cov_plus)){
    fom = as.formula(str_c('Surv_obj ~ x'))
  } else if(!is.null(cov) & is.null(cov_plus)){
    fom = as.formula(str_c('Surv_obj ~ x', '+', str_c(cov, collapse = '+')))
  } else if(is.null(cov) & !is.null(cov_plus)) {
    fom = as.formula(str_c('Surv_obj ~ x', '+', str_c(cov_plus, collapse = '+')))
  } else {
    fom = as.formula(str_c('Surv_obj ~ x', '+', str_c(cov, collapse = '+'), '+', str_c(cov_plus, collapse = '+')))
  }
  
  if(is.null(weight)){
    fit = coxph(fom, data = dt)
  } else {
    dt$weight = data[[weight]]
    fit = coxph(fom, data = dt, weights = weight)
  }
  con = confint(fit)
  
  fit_ph = cox.zph(fit)
  tidy = tidy(fit) %>% 
    mutate(HR = exp(estimate), up = exp(con[, '97.5 %']), 
           low = exp(con[, '2.5 %']), 
           `HR (95% CI)` = str_c(str_trim(format(round(HR, 4), nsmall = 4)), ' (', 
                                 str_trim(format(round(low, 4), nsmall = 4)), ', ', 
                                 str_trim(format(round(up, 4), nsmall = 4)), ')'),
           p_value = ifelse(p.value < 0.001, '<0.001', str_trim(format(round(p.value, 3), nsmall = 3)))) %>% 
    rename(var = term) %>% 
    select(var, `HR (95% CI)`, p_value, p.value, HR, low, up)
  if(class(dt$x) == 'factor'){
    tidy = tidy %>% 
      slice(1:(length(levels(dt$x))-1)) %>% 
      mutate(var = str_replace(var, '^x', ''))
    cox_dt = tibble(var = c(x_label, str_c(levels(dt$x)[1])), 
                    `HR (95% CI)` = c(NA, '1.00 (Reference)'), p_value = NA, p.value = NA) %>% 
      bind_rows(tidy) %>% 
      mutate(data = data_label, outcome = label_outcome, model = cov_label) %>% 
      rename(group = var) %>% 
      select(data, outcome, model, everything()) %>% 
      mutate(ph_p_value = c(fit_ph[[1]][1, 3], rep(NA, length(levels(dt$x)))))
  }
  else {
    cox_dt = tidy %>% 
      slice(1) %>% 
      mutate(var = x_label) %>% 
      mutate(data = data_label, outcome = label_outcome, model = cov_label) %>% 
      rename(group = var) %>% 
      select(data, outcome, model, everything()) %>% 
      mutate(ph_p_value = fit_ph[[1]][1, 3])
  }
  return(cox_dt)
}

# Disease associations ----------------------------------------------------

dta = dt 

outcomes = c("Hypertension2", "Diabetes2", "comCVD2",
             "MetS2", "comCEE2", "COG2", "comCAE2",
             "INF2", "CKD2", "LC2", "Sarcopenia2")
outtimes = rep("time", 11)
labels_outcomes = c("Hypertension", "Diabetes", "Cardiometabolic events",
                    "Metabolic syndrome", "Arterial events", "Cognitive impairment", "Cardiac events",
                    "Inflammaging", "Kidney function decline", "Liver cirrhosis", "Sarcopenia")

rate_map_dt = expand_grid(d1 = tibble(data = map(outcomes, ~ dta %>% 
                                                   mutate(X1 = .[[str_replace(.x, '2', '1')]],
                                                          X2 = .[[.x]]) %>% 
                                                   filter(X1 == 0 & !is.na(X2))
                                                 ), data_label = rep('Total', 11),
                                      outcomes, outtimes, labels_outcomes),
                          group = c('Total'), 
                          tibble('continuous')
                          ) %>% 
  unnest(d1)
# X1 represents the disease status at baseline, and X2 represents the disease status at follow-up.

rate_dt = pmap_dfr(rate_map_dt, ~ rate_fun(data = ..1, data_label = ..2, 
                                           outcome = ..3, outtime = ..4, 
                                           group = ..6, label = ..5,
                                           group_label = ..7))

### association 
organ_names = c("Organismal", 
                "Adipose", "Artery", "Brain", "Heart", "Immune", "Kidney", 
                "Liver", "Muscle")

x = c('age_dev', str_c('age_dev_', organ_names[-1]))
x_label = organ_names
cov0 = NULL
cov1 = c("calendar_age", "sex")
cov2 = c("calendar_age", "sex", "education", "marriage")
cov3 = c("calendar_age", "sex", "education", "marriage", "BMI", "smoke")
cov4 = c("calendar_age", "sex", "education", "marriage", "BMI", "smoke", 
         "SBP", "DBP", "GLU")

cov_label = str_c('Model ', 0:4)

cox_map_dt = expand_grid(d1 = tibble(data = map(outcomes, ~ dta %>% 
                                                  mutate(X1 = .[[str_replace(.x, '2', '1')]],
                                                         X2 = .[[.x]]) %>% 
                                                  filter(X1 == 0 & !is.na(X2))), 
                                     data_label = rep('Total', 11),
                                     cov_plus = map(1, ~ c("CD8T", "CD4T", "NK", "Bcell", "Mono", "Neu")),
                                     outcomes, outtimes, labels_outcomes),
                         tibble(x, x_label),
                         tibble(cov = list(cov0, cov1, cov2, cov3, cov4), cov_label),
                         tibble(weight = list(NULL))) %>% 
  unnest(cols = c(d1))

cox_map_dt$cov_plus[226:270] = map(cox_map_dt$cov_plus[226:270], ~ 
                                     c("CD8T", "CD4T", "NK", "Bcell", "Mono", "Neu", "MMSE"))


mul_cox_dt = pmap_dfr(cox_map_dt, ~ mul_cox_fun(data = ..1, data_label = ..2, cov_plus = ..3,
                                                x = ..7, x_label = ..8, outcome = ..4,
                                                outtime = ..5, label_outcome = ..6,
                                                cov = ..9, cov_label = ..10, weigh = ..11))

rate_dtc = rate_dt %>% 
  select(-group)

output_age_dev = full_join(mul_cox_dt, rate_dtc, by = c('data', 'outcome')) %>% 
  select(data, outcome, model, group, `Event/Total (%)`, `HR (95% CI)`, 
         p_value, ph_p_value, HR, low, up, p.value) %>% 
  pivot_wider(id_cols = c(data, outcome, group, `Event/Total (%)`),
              names_from = model, 
              values_from = c(`HR (95% CI)`, p_value, ph_p_value, HR, low, up, p.value),
              names_glue = '{model}_{.value}') %>% 
  relocate(data, outcome, group, `Event/Total (%)`, 
           contains('Model 0'),
           contains('Model 1'), contains('Model 2'),
           contains('Model 3'), contains('Model 4'))

output_age_dev0 = full_join(mul_cox_dt, rate_dtc, by = c('data', 'outcome')) %>% 
  select(data, outcome, model, group, `Event/Total (%)`, `HR (95% CI)`, 
         p_value, ph_p_value) %>% 
  pivot_wider(id_cols = c(data, outcome, group, `Event/Total (%)`),
              names_from = model, 
              values_from = c(`HR (95% CI)`, p_value, ph_p_value),
              names_glue = '{model}_{.value}') %>% 
  relocate(data, outcome, group, `Event/Total (%)`, 
           contains('Model 0'),
           contains('Model 1'), contains('Model 2'),
           contains('Model 3'), contains('Model 4'))

write_xlsx(list(output_age_dev, output_age_dev0), path = 'disease_organ_cox.xlsx')




