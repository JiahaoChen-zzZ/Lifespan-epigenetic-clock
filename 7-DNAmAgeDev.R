rm(list = ls())

library(tidyverse)
library(broom)

load('data.Rdata')
# This tibble named “data” contains all variables for 5,636 individuals, except for the β values of CpG sites.
dt = data

fit_Organismal = lm(pred_calendar_age ~ calendar_age, data = dt)
fit_Brain = lm(pred_calendar_age_Brain ~ calendar_age, data = dt)
fit_Immune = lm(pred_calendar_age_Immune ~ calendar_age, data = dt)
fit_Muscle = lm(pred_calendar_age_Muscle ~ calendar_age, data = dt)
fit_Artery = lm(pred_calendar_age_Artery ~ calendar_age, data = dt)
fit_Liver = lm(pred_calendar_age_Liver ~ calendar_age, data = dt)
fit_Heart = lm(pred_calendar_age_Heart ~ calendar_age, data = dt)
fit_Adipose = lm(pred_calendar_age_Adipose ~ calendar_age, data = dt)
fit_Kidney = lm(pred_calendar_age_Kidney ~ calendar_age, data = dt)

dt = dt %>% 
  mutate(age_dev = augment(fit_Organismal)$.resid,
         age_dev_Brain = augment(fit_Brain)$.resid,
         age_dev_Immune = augment(fit_Immune)$.resid,
         age_dev_Muscle = augment(fit_Muscle)$.resid,
         age_dev_Artery = augment(fit_Artery)$.resid,
         age_dev_Liver = augment(fit_Liver)$.resid,
         age_dev_Heart = augment(fit_Heart)$.resid,
         age_dev_Adipose = augment(fit_Adipose)$.resid,
         age_dev_Kidney = augment(fit_Kidney)$.resid) %>% 
  relocate(id, cohort, sex, calendar_age, 
           age_dev,
           age_dev_Brain,
           age_dev_Immune,
           age_dev_Muscle,
           age_dev_Artery,
           age_dev_Liver,
           age_dev_Heart,
           age_dev_Adipose,
           age_dev_Kidney)

data = dt
save(data, file = 'data.Rdata')
