rm(list = ls())

library(tidyverse)
library(readr)
library(writexl)
library(broom)
library(DEswan)
library(locfit)
library(patchwork)
library(ggcor)
library(readxl)

# load --------------------------------------------------------------------

load('DMPs.Rdata')
# This tibble named “DMPs” contains 5,636 individuals, a “Sample_Name” column, and 3,237 CpG columns.
dt_DMPs = DMPs
rm(DMPs)

load('data.Rdata')
# This tibble named “data” contains all variables for 5,636 individuals, except for the β values of CpG sites.
dt = data

# weight ------------------------------------------------------------------

dens = density(dt$calendar_age, from = min(dt$calendar_age), to = max(dt$calendar_age), n = 512)
dt$age_density = approx(dens$x, dens$y, xout = dt$calendar_age)$y
dt$weight = 1 / dt$age_density
dt$weight = dt$weight / mean(dt$weight)
dt = relocate(dt, weight, calendar_age)


# Aging peaks --------------------------------------------------------

dtm_DE = dt %>% 
  select(Sample_Name, age_dev, calendar_age, sex, CD8T, CD4T, NK, Bcell, Mono, Neu, weight) %>% 
  left_join(dt_DMPs, by = 'Sample_Name') %>% 
  as.data.frame()

DMPs = colnames(dtm_DE)[12:3248]

library(furrr)
library(progressr)
plan(multisession, workers = 40)

# Create a vector of window center ages
window_centers = 12:100

## overall ----
### window size = 20 ----

handlers(global = TRUE)
handlers("txtprogressbar")

# 
Mod_DEswan_20 = future_map_dfr(window_centers, function(j) {

  dtm = dtm_DE %>% 
    filter(calendar_age >= (j-10) & calendar_age <= (j+10)) %>% 
    mutate(window = case_when(calendar_age >= (j-10) & calendar_age < j ~ 0,
                              calendar_age <= (j+10) & calendar_age > j ~ 1) %>% as_factor()
    ) %>% 
    drop_na(window)
  
  if (nrow(dtm) < 50) return(NULL)
  
  Mod_DEswan_j = map_dfr(DMPs, function(i) {
    fit = glm(dtm[[i]] ~ age_dev*window + sex + CD8T + CD4T + NK + Bcell + Mono + Neu,
              data = dtm, family = gaussian, weights = weight)
    fit.res = car::Anova(fit, type = "2")
    
    tidy = broom::tidy(fit) %>% 
      mutate(beta = estimate, 
             p.value = c(NA, fit.res$`Pr(>Chisq)`),
             `β` = str_trim(format(round(beta, 4), nsmall = 4)),
             p_value = ifelse(p.value < 0.001, '<0.001', str_trim(format(round(p.value, 3), nsmall = 3)))) %>% 
      rename(var = term) %>% 
      select(var, `β`, p_value, p.value, beta)
    
    lm_dt = tidy %>% 
      filter(var %in% c('age_dev', 'window1', 'age_dev:window1')) %>% 
      mutate(var = 1:3) %>% 
      pivot_wider(names_from = var, values_from = c(`β`, p_value, p.value, beta),
                  names_sep = '') %>% 
      mutate(CpG = i) %>% 
      mutate(data = j, N = nrow(dtm)) %>% 
      select(data, CpG, N, everything())
    
    return(lm_dt)
  }, .progress = FALSE)
  
  Mod_DEswan_j = Mod_DEswan_j %>%
    mutate(FDR1 = p.adjust(p.value1, method = 'BH'),
           FDR2 = p.adjust(p.value2, method = 'BH'),
           FDR3 = p.adjust(p.value3, method = 'BH'))
  
  return(Mod_DEswan_j)
}, .progress = TRUE, .options = furrr_options(packages = c("broom", "car", "dplyr", "tidyr", "purrr"), seed = TRUE))

write_rds(Mod_DEswan_20, file = 'Mod_weight_DEswan_20.rds')


### window size = 15 ----

handlers(global = TRUE)
handlers("txtprogressbar")

# 
Mod_DEswan_15 = future_map_dfr(window_centers, function(j) {

  dtm = dtm_DE %>% 
    filter(calendar_age >= (j-7.5) & calendar_age <= (j+7.5)) %>% 
    mutate(window = case_when(calendar_age >= (j-7.5) & calendar_age < j ~ 0,
                              calendar_age <= (j+7.5) & calendar_age > j ~ 1) %>% as_factor()
    ) %>% 
    drop_na(window)
  
  if (nrow(dtm) < 50) return(NULL)
  
  Mod_DEswan_j = map_dfr(DMPs, function(i) {
    fit = glm(dtm[[i]] ~ age_dev*window + sex + CD8T + CD4T + NK + Bcell + Mono + Neu,
              data = dtm, family = gaussian, weights = weight)
    fit.res = car::Anova(fit, type = "2")
    
    tidy = broom::tidy(fit) %>% 
      mutate(beta = estimate, 
             p.value = c(NA, fit.res$`Pr(>Chisq)`),
             `β` = str_trim(format(round(beta, 4), nsmall = 4)),
             p_value = ifelse(p.value < 0.001, '<0.001', str_trim(format(round(p.value, 3), nsmall = 3)))) %>% 
      rename(var = term) %>% 
      select(var, `β`, p_value, p.value, beta)
    
    lm_dt = tidy %>% 
      filter(var %in% c('age_dev', 'window1', 'age_dev:window1')) %>% 
      mutate(var = 1:3) %>% 
      pivot_wider(names_from = var, values_from = c(`β`, p_value, p.value, beta),
                  names_sep = '') %>% 
      mutate(CpG = i) %>% 
      mutate(data = j, N = nrow(dtm)) %>% 
      select(data, CpG, N, everything())
    
    return(lm_dt)
  }, .progress = FALSE)
  
  Mod_DEswan_j = Mod_DEswan_j %>%
    mutate(FDR1 = p.adjust(p.value1, method = 'BH'),
           FDR2 = p.adjust(p.value2, method = 'BH'),
           FDR3 = p.adjust(p.value3, method = 'BH'))
  
  return(Mod_DEswan_j)
}, .progress = TRUE, .options = furrr_options(packages = c("broom", "car", "dplyr", "tidyr", "purrr"), seed = TRUE))

write_rds(Mod_DEswan_15, file = 'Mod_weight_DEswan_15.rds')


### window size = 25 ----

handlers(global = TRUE)
handlers("txtprogressbar")

# 
Mod_DEswan_25 = future_map_dfr(window_centers, function(j) {

  dtm = dtm_DE %>% 
    filter(calendar_age >= (j-12.5) & calendar_age <= (j+12.5)) %>% 
    mutate(window = case_when(calendar_age >= (j-12.5) & calendar_age < j ~ 0,
                              calendar_age <= (j+12.5) & calendar_age > j ~ 1) %>% as_factor()
    ) %>% 
    drop_na(window)
  
  if (nrow(dtm) < 50) return(NULL)
  
  Mod_DEswan_j = map_dfr(DMPs, function(i) {
    fit = glm(dtm[[i]] ~ age_dev*window + sex + CD8T + CD4T + NK + Bcell + Mono + Neu,
              data = dtm, family = gaussian, weights = weight)
    fit.res = car::Anova(fit, type = "2")
    
    tidy = broom::tidy(fit) %>% 
      mutate(beta = estimate, 
             p.value = c(NA, fit.res$`Pr(>Chisq)`),
             `β` = str_trim(format(round(beta, 4), nsmall = 4)),
             p_value = ifelse(p.value < 0.001, '<0.001', str_trim(format(round(p.value, 3), nsmall = 3)))) %>% 
      rename(var = term) %>% 
      select(var, `β`, p_value, p.value, beta)
    
    lm_dt = tidy %>% 
      filter(var %in% c('age_dev', 'window1', 'age_dev:window1')) %>% 
      mutate(var = 1:3) %>% 
      pivot_wider(names_from = var, values_from = c(`β`, p_value, p.value, beta),
                  names_sep = '') %>% 
      mutate(CpG = i) %>% 
      mutate(data = j, N = nrow(dtm)) %>% 
      select(data, CpG, N, everything())
    
    return(lm_dt)
  }, .progress = FALSE)
  
  Mod_DEswan_j = Mod_DEswan_j %>%
    mutate(FDR1 = p.adjust(p.value1, method = 'BH'),
           FDR2 = p.adjust(p.value2, method = 'BH'),
           FDR3 = p.adjust(p.value3, method = 'BH'))
  
  return(Mod_DEswan_j)
}, .progress = TRUE, .options = furrr_options(packages = c("broom", "car", "dplyr", "tidyr", "purrr"), seed = TRUE))

write_rds(Mod_DEswan_25, file = 'Mod_weight_DEswan_25.rds')


## male ----
dtm_male_DE = dtm_DE %>% 
  filter(sex == 1)

### window size = 20 ----

handlers(global = TRUE)
handlers("txtprogressbar")

# 
Mod_DEswan_male_20 = future_map_dfr(window_centers, function(j) {

  dtm = dtm_male_DE %>% 
    filter(calendar_age >= (j-10) & calendar_age <= (j+10)) %>% 
    mutate(window = case_when(calendar_age >= (j-10) & calendar_age < j ~ 0,
                              calendar_age <= (j+10) & calendar_age > j ~ 1) %>% as_factor()
    ) %>% 
    drop_na(window)
  
  if (nrow(dtm) < 50) return(NULL)
  
  Mod_DEswan_male_j = map_dfr(DMPs, function(i) {
    fit = glm(dtm[[i]] ~ age_dev*window + CD8T + CD4T + NK + Bcell + Mono + Neu,
              data = dtm, family = gaussian, weights = weight)
    fit.res = car::Anova(fit, type = "2")
    
    tidy = broom::tidy(fit) %>% 
      mutate(beta = estimate, 
             p.value = c(NA, fit.res$`Pr(>Chisq)`),
             `β` = str_trim(format(round(beta, 4), nsmall = 4)),
             p_value = ifelse(p.value < 0.001, '<0.001', str_trim(format(round(p.value, 3), nsmall = 3)))) %>% 
      rename(var = term) %>% 
      select(var, `β`, p_value, p.value, beta)
    
    lm_dt = tidy %>% 
      filter(var %in% c('age_dev', 'window1', 'age_dev:window1')) %>% 
      mutate(var = 1:3) %>% 
      pivot_wider(names_from = var, values_from = c(`β`, p_value, p.value, beta),
                  names_sep = '') %>% 
      mutate(CpG = i) %>% 
      mutate(data = j, N = nrow(dtm)) %>% 
      select(data, CpG, N, everything())
    
    return(lm_dt)
  }, .progress = FALSE)
  
  Mod_DEswan_male_j = Mod_DEswan_male_j %>%
    mutate(FDR1 = p.adjust(p.value1, method = 'BH'),
           FDR2 = p.adjust(p.value2, method = 'BH'),
           FDR3 = p.adjust(p.value3, method = 'BH'))
  
  return(Mod_DEswan_male_j)
}, .progress = TRUE, .options = furrr_options(packages = c("broom", "car", "dplyr", "tidyr", "purrr"), seed = TRUE))

write_rds(Mod_DEswan_male_20, file = 'Mod_weight_DEswan_male_20.rds')


### window size = 15 ----

handlers(global = TRUE)
handlers("txtprogressbar")

# 
Mod_DEswan_male_15 = future_map_dfr(window_centers, function(j) {

  dtm = dtm_male_DE %>% 
    filter(calendar_age >= (j-7.5) & calendar_age <= (j+7.5)) %>% 
    mutate(window = case_when(calendar_age >= (j-7.5) & calendar_age < j ~ 0,
                              calendar_age <= (j+7.5) & calendar_age > j ~ 1) %>% as_factor()
    ) %>% 
    drop_na(window)
  
  if (nrow(dtm) < 50) return(NULL)
  
  Mod_DEswan_male_j = map_dfr(DMPs, function(i) {
    fit = glm(dtm[[i]] ~ age_dev*window + CD8T + CD4T + NK + Bcell + Mono + Neu,
              data = dtm, family = gaussian, weights = weight)
    fit.res = car::Anova(fit, type = "2")
    
    tidy = broom::tidy(fit) %>% 
      mutate(beta = estimate, 
             p.value = c(NA, fit.res$`Pr(>Chisq)`),
             `β` = str_trim(format(round(beta, 4), nsmall = 4)),
             p_value = ifelse(p.value < 0.001, '<0.001', str_trim(format(round(p.value, 3), nsmall = 3)))) %>% 
      rename(var = term) %>% 
      select(var, `β`, p_value, p.value, beta)
    
    lm_dt = tidy %>% 
      filter(var %in% c('age_dev', 'window1', 'age_dev:window1')) %>% 
      mutate(var = 1:3) %>% 
      pivot_wider(names_from = var, values_from = c(`β`, p_value, p.value, beta),
                  names_sep = '') %>% 
      mutate(CpG = i) %>% 
      mutate(data = j, N = nrow(dtm)) %>% 
      select(data, CpG, N, everything())
    
    return(lm_dt)
  }, .progress = FALSE)
  
  Mod_DEswan_male_j = Mod_DEswan_male_j %>%
    mutate(FDR1 = p.adjust(p.value1, method = 'BH'),
           FDR2 = p.adjust(p.value2, method = 'BH'),
           FDR3 = p.adjust(p.value3, method = 'BH'))
  
  return(Mod_DEswan_male_j)
}, .progress = TRUE, .options = furrr_options(packages = c("broom", "car", "dplyr", "tidyr", "purrr"), seed = TRUE))

write_rds(Mod_DEswan_male_15, file = 'Mod_weight_DEswan_male_15.rds')


### window size = 25 ----

handlers(global = TRUE)
handlers("txtprogressbar")

# 
Mod_DEswan_male_25 = future_map_dfr(window_centers, function(j) {

  dtm = dtm_male_DE %>% 
    filter(calendar_age >= (j-12.5) & calendar_age <= (j+12.5)) %>% 
    mutate(window = case_when(calendar_age >= (j-12.5) & calendar_age < j ~ 0,
                              calendar_age <= (j+12.5) & calendar_age > j ~ 1) %>% as_factor()
    ) %>% 
    drop_na(window)
  
  if (nrow(dtm) < 50) return(NULL)
  
  Mod_DEswan_male_j = map_dfr(DMPs, function(i) {
    fit = glm(dtm[[i]] ~ age_dev*window + CD8T + CD4T + NK + Bcell + Mono + Neu,
              data = dtm, family = gaussian, weights = weight)
    fit.res = car::Anova(fit, type = "2")
    
    tidy = broom::tidy(fit) %>% 
      mutate(beta = estimate, 
             p.value = c(NA, fit.res$`Pr(>Chisq)`),
             `β` = str_trim(format(round(beta, 4), nsmall = 4)),
             p_value = ifelse(p.value < 0.001, '<0.001', str_trim(format(round(p.value, 3), nsmall = 3)))) %>% 
      rename(var = term) %>% 
      select(var, `β`, p_value, p.value, beta)
    
    lm_dt = tidy %>% 
      filter(var %in% c('age_dev', 'window1', 'age_dev:window1')) %>% 
      mutate(var = 1:3) %>% 
      pivot_wider(names_from = var, values_from = c(`β`, p_value, p.value, beta),
                  names_sep = '') %>% 
      mutate(CpG = i) %>% 
      mutate(data = j, N = nrow(dtm)) %>% 
      select(data, CpG, N, everything())
    
    return(lm_dt)
  }, .progress = FALSE)
  
  Mod_DEswan_male_j = Mod_DEswan_male_j %>%
    mutate(FDR1 = p.adjust(p.value1, method = 'BH'),
           FDR2 = p.adjust(p.value2, method = 'BH'),
           FDR3 = p.adjust(p.value3, method = 'BH'))
  
  return(Mod_DEswan_male_j)
}, .progress = TRUE, .options = furrr_options(packages = c("broom", "car", "dplyr", "tidyr", "purrr"), seed = TRUE))

write_rds(Mod_DEswan_male_25, file = 'Mod_weight_DEswan_male_25.rds')


## female ----
dtm_female_DE = dtm_DE %>% 
  filter(sex == 2)

### window size = 20 ----

handlers(global = TRUE)
handlers("txtprogressbar")

# 
Mod_DEswan_female_20 = future_map_dfr(window_centers, function(j) {

  dtm = dtm_female_DE %>% 
    filter(calendar_age >= (j-10) & calendar_age <= (j+10)) %>% 
    mutate(window = case_when(calendar_age >= (j-10) & calendar_age < j ~ 0,
                              calendar_age <= (j+10) & calendar_age > j ~ 1) %>% as_factor()
    ) %>% 
    drop_na(window)
  
  if (nrow(dtm) < 50) return(NULL)
  
  Mod_DEswan_female_j = map_dfr(DMPs, function(i) {
    fit = glm(dtm[[i]] ~ age_dev*window + CD8T + CD4T + NK + Bcell + Mono + Neu,
              data = dtm, family = gaussian, weights = weight)
    fit.res = car::Anova(fit, type = "2")
    
    tidy = broom::tidy(fit) %>% 
      mutate(beta = estimate, 
             p.value = c(NA, fit.res$`Pr(>Chisq)`),
             `β` = str_trim(format(round(beta, 4), nsmall = 4)),
             p_value = ifelse(p.value < 0.001, '<0.001', str_trim(format(round(p.value, 3), nsmall = 3)))) %>% 
      rename(var = term) %>% 
      select(var, `β`, p_value, p.value, beta)
    
    lm_dt = tidy %>% 
      filter(var %in% c('age_dev', 'window1', 'age_dev:window1')) %>% 
      mutate(var = 1:3) %>% 
      pivot_wider(names_from = var, values_from = c(`β`, p_value, p.value, beta),
                  names_sep = '') %>% 
      mutate(CpG = i) %>% 
      mutate(data = j, N = nrow(dtm)) %>% 
      select(data, CpG, N, everything())
    
    return(lm_dt)
  }, .progress = FALSE)
  
  Mod_DEswan_female_j = Mod_DEswan_female_j %>%
    mutate(FDR1 = p.adjust(p.value1, method = 'BH'),
           FDR2 = p.adjust(p.value2, method = 'BH'),
           FDR3 = p.adjust(p.value3, method = 'BH'))
  
  return(Mod_DEswan_female_j)
}, .progress = TRUE, .options = furrr_options(packages = c("broom", "car", "dplyr", "tidyr", "purrr"), seed = TRUE))

write_rds(Mod_DEswan_female_20, file = 'Mod_weight_DEswan_female_20.rds')


### window size = 15 ----

handlers(global = TRUE)
handlers("txtprogressbar")

# 
Mod_DEswan_female_15 = future_map_dfr(window_centers, function(j) {

  dtm = dtm_female_DE %>% 
    filter(calendar_age >= (j-7.5) & calendar_age <= (j+7.5)) %>% 
    mutate(window = case_when(calendar_age >= (j-7.5) & calendar_age < j ~ 0,
                              calendar_age <= (j+7.5) & calendar_age > j ~ 1) %>% as_factor()
    ) %>% 
    drop_na(window)
  
  if (nrow(dtm) < 50) return(NULL)
  
  Mod_DEswan_female_j = map_dfr(DMPs, function(i) {
    fit = glm(dtm[[i]] ~ age_dev*window + CD8T + CD4T + NK + Bcell + Mono + Neu,
              data = dtm, family = gaussian, weights = weight)
    fit.res = car::Anova(fit, type = "2")
    
    tidy = broom::tidy(fit) %>% 
      mutate(beta = estimate, 
             p.value = c(NA, fit.res$`Pr(>Chisq)`),
             `β` = str_trim(format(round(beta, 4), nsmall = 4)),
             p_value = ifelse(p.value < 0.001, '<0.001', str_trim(format(round(p.value, 3), nsmall = 3)))) %>% 
      rename(var = term) %>% 
      select(var, `β`, p_value, p.value, beta)
    
    lm_dt = tidy %>% 
      filter(var %in% c('age_dev', 'window1', 'age_dev:window1')) %>% 
      mutate(var = 1:3) %>% 
      pivot_wider(names_from = var, values_from = c(`β`, p_value, p.value, beta),
                  names_sep = '') %>% 
      mutate(CpG = i) %>% 
      mutate(data = j, N = nrow(dtm)) %>% 
      select(data, CpG, N, everything())
    
    return(lm_dt)
  }, .progress = FALSE)
  
  Mod_DEswan_female_j = Mod_DEswan_female_j %>%
    mutate(FDR1 = p.adjust(p.value1, method = 'BH'),
           FDR2 = p.adjust(p.value2, method = 'BH'),
           FDR3 = p.adjust(p.value3, method = 'BH'))
  
  return(Mod_DEswan_female_j)
}, .progress = TRUE, .options = furrr_options(packages = c("broom", "car", "dplyr", "tidyr", "purrr"), seed = TRUE))

write_rds(Mod_DEswan_female_15, file = 'Mod_weight_DEswan_female_15.rds')


### window size = 25 ----

handlers(global = TRUE)
handlers("txtprogressbar")

# 
Mod_DEswan_female_25 = future_map_dfr(window_centers, function(j) {

  dtm = dtm_female_DE %>% 
    filter(calendar_age >= (j-12.5) & calendar_age <= (j+12.5)) %>% 
    mutate(window = case_when(calendar_age >= (j-12.5) & calendar_age < j ~ 0,
                              calendar_age <= (j+12.5) & calendar_age > j ~ 1) %>% as_factor()
    ) %>% 
    drop_na(window)
  
  if (nrow(dtm) < 50) return(NULL)
  
  Mod_DEswan_female_j = map_dfr(DMPs, function(i) {
    fit = glm(dtm[[i]] ~ age_dev*window + CD8T + CD4T + NK + Bcell + Mono + Neu,
              data = dtm, family = gaussian, weights = weight)
    fit.res = car::Anova(fit, type = "2")
    
    tidy = broom::tidy(fit) %>% 
      mutate(beta = estimate, 
             p.value = c(NA, fit.res$`Pr(>Chisq)`),
             `β` = str_trim(format(round(beta, 4), nsmall = 4)),
             p_value = ifelse(p.value < 0.001, '<0.001', str_trim(format(round(p.value, 3), nsmall = 3)))) %>% 
      rename(var = term) %>% 
      select(var, `β`, p_value, p.value, beta)
    
    lm_dt = tidy %>% 
      filter(var %in% c('age_dev', 'window1', 'age_dev:window1')) %>% 
      mutate(var = 1:3) %>% 
      pivot_wider(names_from = var, values_from = c(`β`, p_value, p.value, beta),
                  names_sep = '') %>% 
      mutate(CpG = i) %>% 
      mutate(data = j, N = nrow(dtm)) %>% 
      select(data, CpG, N, everything())
    
    return(lm_dt)
  }, .progress = FALSE)
  
  Mod_DEswan_female_j = Mod_DEswan_female_j %>%
    mutate(FDR1 = p.adjust(p.value1, method = 'BH'),
           FDR2 = p.adjust(p.value2, method = 'BH'),
           FDR3 = p.adjust(p.value3, method = 'BH'))
  
  return(Mod_DEswan_female_j)
}, .progress = TRUE, .options = furrr_options(packages = c("broom", "car", "dplyr", "tidyr", "purrr"), seed = TRUE))

write_rds(Mod_DEswan_female_25, file = 'Mod_weight_DEswan_female_25.rds')

# FDR ---------------------------------------------------------------------

### overall ----

#
Mod_weight_DEswan_p_20_1 = Mod_weight_DEswan_20 %>% 
  rename(calendar_age = data) %>% 
  group_by(calendar_age) %>% 
  filter(FDR1 < 0.05 & FDR2 < 0.05 & FDR3 < 0.05) %>% 
  summarise(`0.05` = n())

Mod_weight_DEswan_p_20_2 = Mod_weight_DEswan_20 %>% 
  rename(calendar_age = data) %>% 
  group_by(calendar_age) %>% 
  filter(FDR1 < 0.01 & FDR2 < 0.01 & FDR3 < 0.01) %>% 
  summarise(`0.01` = n())

Mod_weight_DEswan_p_20_3 = Mod_weight_DEswan_20 %>% 
  rename(calendar_age = data) %>% 
  group_by(calendar_age) %>% 
  filter(FDR1 < 0.001 & FDR2 < 0.001 & FDR3 < 0.001) %>% 
  summarise(`0.001` = n())

Mod_weight_DEswan_p_20 = tibble(calendar_age = 12:100) %>% 
  left_join(Mod_weight_DEswan_p_20_1, by = 'calendar_age') %>% 
  left_join(Mod_weight_DEswan_p_20_2, by = 'calendar_age') %>% 
  left_join(Mod_weight_DEswan_p_20_3, by = 'calendar_age') %>% 
  map_dfc(~ replace_na(.x, 0))

#
Mod_weight_DEswan_p_15_1 = Mod_weight_DEswan_15 %>% 
  rename(calendar_age = data) %>% 
  group_by(calendar_age) %>% 
  filter(FDR1 < 0.05 & FDR2 < 0.05 & FDR3 < 0.05) %>% 
  summarise(`0.05` = n())

Mod_weight_DEswan_p_15_2 = Mod_weight_DEswan_15 %>% 
  rename(calendar_age = data) %>% 
  group_by(calendar_age) %>% 
  filter(FDR1 < 0.01 & FDR2 < 0.01 & FDR3 < 0.01) %>% 
  summarise(`0.01` = n())

Mod_weight_DEswan_p_15_3 = Mod_weight_DEswan_15 %>% 
  rename(calendar_age = data) %>% 
  group_by(calendar_age) %>% 
  filter(FDR1 < 0.001 & FDR2 < 0.001 & FDR3 < 0.001) %>% 
  summarise(`0.001` = n())

Mod_weight_DEswan_p_15 = tibble(calendar_age = 12:100) %>% 
  left_join(Mod_weight_DEswan_p_15_1, by = 'calendar_age') %>% 
  left_join(Mod_weight_DEswan_p_15_2, by = 'calendar_age') %>% 
  left_join(Mod_weight_DEswan_p_15_3, by = 'calendar_age') %>% 
  map_dfc(~ replace_na(.x, 0))

#
Mod_weight_DEswan_p_25_1 = Mod_weight_DEswan_25 %>% 
  rename(calendar_age = data) %>% 
  group_by(calendar_age) %>% 
  filter(FDR1 < 0.05 & FDR2 < 0.05 & FDR3 < 0.05) %>% 
  summarise(`0.05` = n())

Mod_weight_DEswan_p_25_2 = Mod_weight_DEswan_25 %>% 
  rename(calendar_age = data) %>% 
  group_by(calendar_age) %>% 
  filter(FDR1 < 0.01 & FDR2 < 0.01 & FDR3 < 0.01) %>% 
  summarise(`0.01` = n())

Mod_weight_DEswan_p_25_3 = Mod_weight_DEswan_25 %>% 
  rename(calendar_age = data) %>% 
  group_by(calendar_age) %>% 
  filter(FDR1 < 0.001 & FDR2 < 0.001 & FDR3 < 0.001) %>% 
  summarise(`0.001` = n())

Mod_weight_DEswan_p_25 = tibble(calendar_age = 12:100) %>% 
  left_join(Mod_weight_DEswan_p_25_1, by = 'calendar_age') %>% 
  left_join(Mod_weight_DEswan_p_25_2, by = 'calendar_age') %>% 
  left_join(Mod_weight_DEswan_p_25_3, by = 'calendar_age') %>% 
  map_dfc(~ replace_na(.x, 0))

### male ----

#
Mod_weight_DEswan_male_p_20_1 = Mod_weight_DEswan_male_20 %>% 
  rename(calendar_age = data) %>% 
  group_by(calendar_age) %>% 
  filter(FDR1 < 0.05 & FDR2 < 0.05 & FDR3 < 0.05) %>% 
  summarise(`0.05` = n())

Mod_weight_DEswan_male_p_20_2 = Mod_weight_DEswan_male_20 %>% 
  rename(calendar_age = data) %>% 
  group_by(calendar_age) %>% 
  filter(FDR1 < 0.01 & FDR2 < 0.01 & FDR3 < 0.01) %>% 
  summarise(`0.01` = n())

Mod_weight_DEswan_male_p_20_3 = Mod_weight_DEswan_male_20 %>% 
  rename(calendar_age = data) %>% 
  group_by(calendar_age) %>% 
  filter(FDR1 < 0.001 & FDR2 < 0.001 & FDR3 < 0.001) %>% 
  summarise(`0.001` = n())

Mod_weight_DEswan_male_p_20 = tibble(calendar_age = 12:100) %>% 
  left_join(Mod_weight_DEswan_male_p_20_1, by = 'calendar_age') %>% 
  left_join(Mod_weight_DEswan_male_p_20_2, by = 'calendar_age') %>% 
  left_join(Mod_weight_DEswan_male_p_20_3, by = 'calendar_age') %>% 
  map_dfc(~ replace_na(.x, 0))

#
Mod_weight_DEswan_male_p_15_1 = Mod_weight_DEswan_male_15 %>% 
  rename(calendar_age = data) %>% 
  group_by(calendar_age) %>% 
  filter(FDR1 < 0.05 & FDR2 < 0.05 & FDR3 < 0.05) %>% 
  summarise(`0.05` = n())

Mod_weight_DEswan_male_p_15_2 = Mod_weight_DEswan_male_15 %>% 
  rename(calendar_age = data) %>% 
  group_by(calendar_age) %>% 
  filter(FDR1 < 0.01 & FDR2 < 0.01 & FDR3 < 0.01) %>% 
  summarise(`0.01` = n())

Mod_weight_DEswan_male_p_15_3 = Mod_weight_DEswan_male_15 %>% 
  rename(calendar_age = data) %>% 
  group_by(calendar_age) %>% 
  filter(FDR1 < 0.001 & FDR2 < 0.001 & FDR3 < 0.001) %>% 
  summarise(`0.001` = n())

Mod_weight_DEswan_male_p_15 = tibble(calendar_age = 12:100) %>% 
  left_join(Mod_weight_DEswan_male_p_15_1, by = 'calendar_age') %>% 
  left_join(Mod_weight_DEswan_male_p_15_2, by = 'calendar_age') %>% 
  left_join(Mod_weight_DEswan_male_p_15_3, by = 'calendar_age') %>% 
  map_dfc(~ replace_na(.x, 0))

#
Mod_weight_DEswan_male_p_25_1 = Mod_weight_DEswan_male_25 %>% 
  rename(calendar_age = data) %>% 
  group_by(calendar_age) %>% 
  filter(FDR1 < 0.05 & FDR2 < 0.05 & FDR3 < 0.05) %>% 
  summarise(`0.05` = n())

Mod_weight_DEswan_male_p_25_2 = Mod_weight_DEswan_male_25 %>% 
  rename(calendar_age = data) %>% 
  group_by(calendar_age) %>% 
  filter(FDR1 < 0.01 & FDR2 < 0.01 & FDR3 < 0.01) %>% 
  summarise(`0.01` = n())

Mod_weight_DEswan_male_p_25_3 = Mod_weight_DEswan_male_25 %>% 
  rename(calendar_age = data) %>% 
  group_by(calendar_age) %>% 
  filter(FDR1 < 0.001 & FDR2 < 0.001 & FDR3 < 0.001) %>% 
  summarise(`0.001` = n())

Mod_weight_DEswan_male_p_25 = tibble(calendar_age = 12:100) %>% 
  left_join(Mod_weight_DEswan_male_p_25_1, by = 'calendar_age') %>% 
  left_join(Mod_weight_DEswan_male_p_25_2, by = 'calendar_age') %>% 
  left_join(Mod_weight_DEswan_male_p_25_3, by = 'calendar_age') %>% 
  map_dfc(~ replace_na(.x, 0))

### female ----

#
Mod_weight_DEswan_female_p_20_1 = Mod_weight_DEswan_female_20 %>% 
  rename(calendar_age = data) %>% 
  group_by(calendar_age) %>% 
  filter(FDR1 < 0.05 & FDR2 < 0.05 & FDR3 < 0.05) %>% 
  summarise(`0.05` = n())

Mod_weight_DEswan_female_p_20_2 = Mod_weight_DEswan_female_20 %>% 
  rename(calendar_age = data) %>% 
  group_by(calendar_age) %>% 
  filter(FDR1 < 0.01 & FDR2 < 0.01 & FDR3 < 0.01) %>% 
  summarise(`0.01` = n())

Mod_weight_DEswan_female_p_20_3 = Mod_weight_DEswan_female_20 %>% 
  rename(calendar_age = data) %>% 
  group_by(calendar_age) %>% 
  filter(FDR1 < 0.001 & FDR2 < 0.001 & FDR3 < 0.001) %>% 
  summarise(`0.001` = n())

Mod_weight_DEswan_female_p_20 = tibble(calendar_age = 12:100) %>% 
  left_join(Mod_weight_DEswan_female_p_20_1, by = 'calendar_age') %>% 
  left_join(Mod_weight_DEswan_female_p_20_2, by = 'calendar_age') %>% 
  left_join(Mod_weight_DEswan_female_p_20_3, by = 'calendar_age') %>% 
  map_dfc(~ replace_na(.x, 0))

#
Mod_weight_DEswan_female_p_15_1 = Mod_weight_DEswan_female_15 %>% 
  rename(calendar_age = data) %>% 
  group_by(calendar_age) %>% 
  filter(FDR1 < 0.05 & FDR2 < 0.05 & FDR3 < 0.05) %>% 
  summarise(`0.05` = n())

Mod_weight_DEswan_female_p_15_2 = Mod_weight_DEswan_female_15 %>% 
  rename(calendar_age = data) %>% 
  group_by(calendar_age) %>% 
  filter(FDR1 < 0.01 & FDR2 < 0.01 & FDR3 < 0.01) %>% 
  summarise(`0.01` = n())

Mod_weight_DEswan_female_p_15_3 = Mod_weight_DEswan_female_15 %>% 
  rename(calendar_age = data) %>% 
  group_by(calendar_age) %>% 
  filter(FDR1 < 0.001 & FDR2 < 0.001 & FDR3 < 0.001) %>% 
  summarise(`0.001` = n())

Mod_weight_DEswan_female_p_15 = tibble(calendar_age = 12:100) %>% 
  left_join(Mod_weight_DEswan_female_p_15_1, by = 'calendar_age') %>% 
  left_join(Mod_weight_DEswan_female_p_15_2, by = 'calendar_age') %>% 
  left_join(Mod_weight_DEswan_female_p_15_3, by = 'calendar_age') %>% 
  map_dfc(~ replace_na(.x, 0))

#
Mod_weight_DEswan_female_p_25_1 = Mod_weight_DEswan_female_25 %>% 
  rename(calendar_age = data) %>% 
  group_by(calendar_age) %>% 
  filter(FDR1 < 0.05 & FDR2 < 0.05 & FDR3 < 0.05) %>% 
  summarise(`0.05` = n())

Mod_weight_DEswan_female_p_25_2 = Mod_weight_DEswan_female_25 %>% 
  rename(calendar_age = data) %>% 
  group_by(calendar_age) %>% 
  filter(FDR1 < 0.01 & FDR2 < 0.01 & FDR3 < 0.01) %>% 
  summarise(`0.01` = n())

Mod_weight_DEswan_female_p_25_3 = Mod_weight_DEswan_female_25 %>% 
  rename(calendar_age = data) %>% 
  group_by(calendar_age) %>% 
  filter(FDR1 < 0.001 & FDR2 < 0.001 & FDR3 < 0.001) %>% 
  summarise(`0.001` = n())

Mod_weight_DEswan_female_p_25 = tibble(calendar_age = 12:100) %>% 
  left_join(Mod_weight_DEswan_female_p_25_1, by = 'calendar_age') %>% 
  left_join(Mod_weight_DEswan_female_p_25_2, by = 'calendar_age') %>% 
  left_join(Mod_weight_DEswan_female_p_25_3, by = 'calendar_age') %>% 
  map_dfc(~ replace_na(.x, 0))