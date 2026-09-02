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


# data --------------------------------------------------------------------

load('data.Rdata')
# This tibble named “data” contains all variables for 5,636 individuals, except for the β values of CpG sites.
dt = data

# RCS ---------------------------------------------------------------------

dtd = dt %>% 
  mutate(death = case_when(survival_time <= 0 ~ NA_real_,
                           survival_time > 0 ~ death)) %>% 
  filter(calendar_age >= 45 & !is.na(death)) %>% 
  as.data.frame()

## Set datadist ----
dd = datadist(dtd)
options(datadist = "dd")

## all-cause death ----

## Determine optimal number of knots for rcs by minimizing AIC
nk = c(3, 4, 5)
ls_aic = vector("list", 3)
Surv_obj = Surv(dtd$survival_time, dtd$death)

for (i in nk) {
  mod_tmp = cph(
    formula = Surv_obj ~ rcs(age_dev, i) + calendar_age + sex + education + 
      marriage + BMI + smoke + hypertension + CD8T + CD4T + NK + Bcell + Mono + Neu,
    data = dtd,
    x = TRUE, y = TRUE, surv = TRUE
  )
  ls_aic[[i]] = str_c("AIC of model with ", i, " knots: ", round(AIC(mod_tmp), 3))
  ls_aic = ls_aic[!sapply(ls_aic, is.null)]
  rm(mod_tmp)
}

ls_aic


## Model using rcs with 3 knots
ref_value = 0 
dd$limits["Adjust to", "age_dev"] = ref_value

mod1 = cph(
  formula = Surv_obj ~ rcs(age_dev, 3) + calendar_age + sex + education + 
    marriage + BMI + smoke + hypertension + CD8T + CD4T + NK + Bcell + Mono + Neu,
  data = dtd,
  x = TRUE, y = TRUE, surv = TRUE
)

pred_mod1 = rms::Predict(mod1, age_dev, 
                         fun = exp, ref.zero = TRUE)

# Reset datadist
options(datadist = NULL)


# Cox ---------------------------------------------------------------------
rm(list = ls())

load('data.Rdata')
# This tibble named “data” contains all variables for 5,636 individuals, except for the β values of CpG sites.

dt = data %>% 
  mutate(death = case_when(survival_time <= 0 ~ NA_real_,
                           survival_time > 0 ~ death) %>% as_factor()) %>% 
  filter(calendar_age >= 45 & !is.na(death))

MAE_DNAmAge = read_excel('MAE_DNAmAge.xlsx')
# This is data recording the MAE of the clock model.


## Cox function ----
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
           `HR (95% CI)` = str_c(str_trim(format(round(HR, 2), nsmall = 2)), ' (', 
                                 str_trim(format(round(low, 2), nsmall = 2)), ', ', 
                                 str_trim(format(round(up, 2), nsmall = 2)), ')'),
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

## Age ----
dtt = dt %>% 
  mutate(death = as.numeric(death) - 1)

dtt2 = dtt %>% 
  filter(!is.na(CV_death))

outcomes = c("death", "CV_death", "nonCV_death")
outtimes = c("survival_time", "survival_time", "survival_time")
labels_outcomes = c("All-cause mortality", "Cardiovascular mortality", 
                    "Non-cardiovascular mortality")

rate_map_dt = expand_grid(d1 = tibble(data = list(dtt, dtt2, dtt2), 
                                      data_label = rep('Total', 3),
                                      outcomes, outtimes, labels_outcomes),
                          group = c('Total'), 
                          tibble('continuous')) %>% 
  unnest(d1)

rate_dt = pmap_dfr(rate_map_dt, ~ rate_fun(data = ..1, data_label = ..2, 
                                           outcome = ..3, outtime = ..4, 
                                           group = ..6, label = ..5,
                                           group_label = ..7)) %>% 
  select(-group)


### calendar_age, DNAmAge ----
x = c("calendar_age", 'pred_calendar_age')
x_label = c("Chronological age, per year", "DNAmAge, per year")
cov0 = NULL

cov_label = str_c('Model ', 0)

cox_map_dt = expand_grid(d1 = tibble(data = list(dtt, dtt2, dtt2), 
                                     data_label = rep('Total', 3),
                                     cov_plus = map(1:3, ~ c("CD8T", "CD4T", "NK", "Bcell", "Mono", "Neu")),
                                     outcomes, outtimes, labels_outcomes),
                         tibble(x, x_label),
                         tibble(cov = list(cov0), cov_label),
                         tibble(weight = list(NULL))) %>% 
  unnest(cols = c(d1))

mul_cox_dt = pmap_dfr(cox_map_dt, ~ mul_cox_fun(data = ..1, data_label = ..2, cov_plus = ..3,
                                                x = ..7, x_label = ..8, outcome = ..4,
                                                outtime = ..5, label_outcome = ..6,
                                                cov = ..9, cov_label = ..10, weigh = ..11))

output_age_con = full_join(mul_cox_dt, rate_dt, by = c('data', 'outcome')) %>% 
  select(data, outcome, model, group, `Event/Total (%)`, `HR (95% CI)`, 
         p_value, ph_p_value, HR, low, up, p.value)

output_age_con0 = full_join(mul_cox_dt, rate_dt, by = c('data', 'outcome')) %>% 
  select(data, outcome, model, group, `Event/Total (%)`, `HR (95% CI)`, 
         p_value, ph_p_value)

write_xlsx(list(output_age_con, output_age_con0), path = 'cox_age.xlsx')


## DNAmAgeDev ---- 
### rate ----
dtt = dt %>% 
  mutate(death = as.numeric(death) - 1) %>% 
  mutate(dev_c = case_when(abs(age_dev) <= MAE_DNAmAge$MAE ~ 'Normal agers',
                           age_dev < -MAE_DNAmAge$MAE ~ 'Slow agers',
                           age_dev > MAE_DNAmAge$MAE ~ 'Fast agers'),
         dev_c = factor(dev_c, levels = c('Normal agers', 'Slow agers', 'Fast agers')))

dtt2 = dtt %>% 
  filter(!is.na(CV_death))

outcomes = c("death", "CV_death", "nonCV_death")
outtimes = c("survival_time", "survival_time", "survival_time")
labels_outcomes = c("All-cause mortality", "Cardiovascular mortality", 
                    "Non-cardiovascular mortality")

rate_map_dt = expand_grid(d1 = tibble(data = list(dtt, dtt2, dtt2), 
                                      data_label = rep('Total', 3),
                                      outcomes, outtimes, labels_outcomes),
                          group = c('Total', 'dev_c'), 
                          tibble('continuous')) %>% 
  unnest(d1)

rate_dt = pmap_dfr(rate_map_dt, ~ rate_fun(data = ..1, data_label = ..2, 
                                           outcome = ..3, outtime = ..4, 
                                           group = ..6, label = ..5,
                                           group_label = ..7))

### association ----

x = c("age_dev", "dev_c")
x_label = c("DNAmAgeDev, per year", "DNAmAgeDev category")
cov0 = NULL
cov1 = c("calendar_age", "sex")
cov2 = c("calendar_age", "sex", "education", "marriage")
cov3 = c("calendar_age", "sex", "education", "marriage", 
         "BMI")
cov4 = c("calendar_age", "sex", "education", "marriage", 
         "BMI", "smoke", "hypertension")

cov_label = str_c('Model ', 0:4)

cox_map_dt = expand_grid(d1 = tibble(data = list(dtt, dtt2, dtt2), 
                                     data_label = rep('Total', 3),
                                     cov_plus = map(1:3, ~ c("CD8T", "CD4T", "NK", "Bcell", "Mono", "Neu")),
                                     outcomes, outtimes, labels_outcomes),
                         tibble(x, x_label),
                         tibble(cov = list(cov0, cov1, cov2, cov3, cov4), cov_label),
                         tibble(weight = list(NULL))) %>% 
  unnest(cols = c(d1))

mul_cox_dt = pmap_dfr(cox_map_dt, ~ mul_cox_fun(data = ..1, data_label = ..2, cov_plus = ..3,
                                                x = ..7, x_label = ..8, outcome = ..4,
                                                outtime = ..5, label_outcome = ..6,
                                                cov = ..9, cov_label = ..10, weigh = ..11))

rate_dtc = rate_dt %>% 
  mutate(group = case_when(group == 'continuous' ~ x_label[1],
                           TRUE ~ group)) %>% 
  bind_rows(rate_dt %>% 
              mutate(group = case_when(group == 'continuous' ~ x_label[2],
                                       TRUE ~ group))) %>% 
  distinct()

output_age_dev = full_join(mul_cox_dt, rate_dtc, by = c('data', 'outcome', 'group')) %>% 
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

output_age_dev0 = full_join(mul_cox_dt, rate_dtc, by = c('data', 'outcome', 'group')) %>% 
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

write_xlsx(list(output_age_dev, output_age_dev0), path = 'cox_age_dev.xlsx')
