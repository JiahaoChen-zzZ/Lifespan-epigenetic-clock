rm(list = ls())

library(tidyverse)
library(haven)
library(lubridate)
library(readxl)
library(janitor)
library(writexl)
library(aplot)
library(ggpointdensity)
library(tidyplots)
library(aplot)
library(patchwork)
library(broom)
library(cowplot)
library(car)
library(glmnet)


# data --------------------------------------------------------------------

load('data.Rdata')
# This tibble named “data” contains all variables for 5,636 individuals, except for the β values of CpG sites.
dt = data

outcomes = c("smoke_num", 
             "drink_freq", 
             'FV',
             'Redmeat', 'Poultry', 
             'EM', 
             'VA', 'MA', 'LA', 'SA')
labels = c("Smoking", 
           "Drinking", 
           'Fresh fruits and vegetables',
           'Red meat', 'Poultry', 
           'Eggs and milk', 
           'VPA', 'MPA', 'LPA', 'SED')

# Overall importance ------------------------------------------------------

organ_names = c("Adipose", "Artery", "Brain", "Heart", "Immune", "Kidney", 
                "Liver", "Muscle")

clocks = c('age_dev', str_c('age_dev_', organ_names))
clocks_labels = c('Organismal', organ_names)

df = dt %>% 
  mutate(across(all_of(clocks), ~ case_when(.x > quantile(.x, 0.99) ~ quantile(.x, 0.99),
                                            .x < quantile(.x, 0.01) ~ quantile(.x, 0.01),
                                            TRUE ~ .x))
  )

dto = drop_na(df,
              smoke_num,
              drink_freq,
              FV,
              Redmeat, Poultry, 
              EM, 
              VA, MA, LA, SA) %>%
  arrange(calendar_age)

## Complete sample combined model ----

covariates = c('calendar_age', 'sex', 'BMI', 
               'CD8T', 'CD4T', 'NK', 'Bcell', 'Mono', 'Neu')

formula_full = as.formula(
  str_c("age_dev ~", 
        str_c(outcomes, collapse = " + "), '+',
        str_c(covariates, collapse = " + ")
  ))

full_model = lm(formula_full, data = dto)
full_summary = tidy(full_model)
vif_values = vif(full_model)
vif_df = tibble(variable = names(vif_values), vif = vif_values)

## Calculate the independent R² contribution for each variable ----

calculate_incremental_r2 = function(data, full_model, var_to_remove) {
  dt = data
  
  full_r2 = glance(full_model)$r.squared
  
  full_formula = formula(full_model)
  all_vars = all.vars(full_formula)
  
  remaining_vars = setdiff(all_vars[-1], var_to_remove)
  reduced_formula = as.formula(
    str_c(all_vars[1], " ~ ", str_c(remaining_vars, collapse = " + "))
  )
  
  reduced_model = lm(reduced_formula, data = dt)
  reduced_r2 = glance(reduced_model)$r.squared
  
  delta_r2 = full_r2 - reduced_r2
  
  anova_test = anova(reduced_model, full_model)
  incremental_p_value = anova_test$`Pr(>F)`[2]
  
  out = tibble(r2 = delta_r2, 
               p_value = incremental_p_value)
  
  return(out)
}

# Calculate the incremental R² for each lifestyle variable (excluding covariates)
r2_contributions = map_dfr(outcomes, ~ calculate_incremental_r2(dto, full_model, .x))

# Calculate the total R² for lifestyle variables
base_formula = as.formula(str_c("age_dev ~ ", str_c(covariates, collapse = " + ")))
base_model = lm(base_formula, data = dto)
base_r2 = glance(base_model)$r.squared
lifestyle_total_r2 = glance(full_model)$r.squared - base_r2

r2_percentages = (r2_contributions$r2 / lifestyle_total_r2) * 100

importance_df = tibble(
  variable = outcomes,
  r2_percentage = r2_percentages,
  beta = full_summary$estimate[(1:length(outcomes)) + 1],
  p_value = r2_contributions$p_value
)

importance_df = importance_df[order(-importance_df$r2_percentage), ]

out_imp = tibble(variable = outcomes, label = labels) %>% 
  left_join(importance_df, by = 'variable')
write_xlsx(out_imp, path = 'lifestyle_imp_total.xlsx')


# Age group analyses ------------------------------------------------------


dto = df %>%
  mutate(
    age_group = case_when(
      calendar_age >= 3 & calendar_age < 18 ~ "3-18",
      calendar_age >= 18 & calendar_age < 65 ~ "18-65", 
      calendar_age >= 65 & calendar_age < 90 ~ "65-90",
      calendar_age >= 90 ~ "90-118",
    ),
    age_group = factor(age_group, levels = c("3-18", "18-65", "65-90", "90-118"))
  )

## Calculate the incremental R² for a variable within a specific age group ---- 

age_group_vars = list(
  "3-18" = setdiff(outcomes, c('smoke_num', 'drink_freq')),
  "18-65" = outcomes,
  "65-90" = outcomes,
  "90-118" = outcomes
)

cov_group_vars = list(
  "3-18" = setdiff(covariates, c('hypertension')),
  "18-65" = covariates,
  "65-90" = covariates,
  "90-118" = covariates
)

# One step function
simple_one_step_screening = function(data, covariates, lifestyle_vars, 
                                     outcome, top_n = 3) {
  dt = data %>% 
    select(all_of(c(outcome, lifestyle_vars, covariates))) %>% 
    drop_na()
  
  # base model
  base_formula = as.formula(str_c(outcome, "~", str_c(covariates, collapse = " + ")))
  base_model = lm(base_formula, data = dt)
  base_r2 = glance(base_model)$r.squared
  
  # Calculate the independent contribution of each lifestyle variable
  results = list()
  
  for(var in lifestyle_vars) {

    test_formula = as.formula(str_c(outcome, "~", var, '+', str_c(covariates, collapse = " + ")))
    test_model = lm(test_formula, data = dt)
    
    test_r2 = glance(test_model)$r.squared
    delta_r2 = test_r2 - base_r2
    
    anova_test = anova(base_model, test_model)
    p_value = anova_test$`Pr(>F)`[2]
    
    coef_summary = tidy(test_model)
    beta = coef_summary$estimat[which(coef_summary$term == var)]
    
    results[[var]] = tibble(
      variable = var,
      incremental_r2 = delta_r2,
      p_value = p_value,
      beta = beta,
      stringsAsFactors = FALSE
    )
  }
  
  results_df = bind_rows(results)
  
  # Filtering and sorting
  filtered_df = results_df %>%
    arrange(desc(incremental_r2))
  
  # Select the top_n most important variables
  if(nrow(filtered_df) > top_n) {
    selected_df = head(filtered_df, top_n)
  } else {
    selected_df = filtered_df
  }
  
  return(selected_df$variable)
}

# A list of stored results
importance_results = list()

# Run through each age group
for(current_age_group in c("3-18", "18-65", "65-90", "90-118")) {
  
  cat("Age group being processed:", current_age_group, "\n")
  
  current_vars = age_group_vars[[current_age_group]]
  current_covs = cov_group_vars[[current_age_group]]
  
  df_sub = dto %>% 
    filter(age_group == current_age_group) %>%
    select(all_of(c('age_dev', current_vars, current_covs))) %>%
    drop_na()
  
  # Variable filtering
  current_vars = simple_one_step_screening(df_sub, current_covs, current_vars, 'age_dev')
  
  formula_str = str_c(
    "age_dev ~",
    str_c(current_vars, collapse = " + "),
    "+",
    str_c(current_covs, collapse = " + ")
  )
  
  formula_full = as.formula(formula_str)
  full_model = lm(formula_full, data = df_sub)
  full_r2 = glance(full_model)$r.squared
  
  # base model
  base_formula_str = str_c("age_dev ~", str_c(current_covs, collapse = " + "))
  base_formula = as.formula(base_formula_str)
  base_model = lm(base_formula, data = df_sub)
  base_r2 = glance(base_model)$r.squared
  
  lifestyle_total_r2 = max(0, full_r2 - base_r2)
  
  # Calculate the incremental R² for each variable
  var_contributions = map_dfr(current_vars, function(var) {

    reduced_formula_str = str_c(
      "age_dev ~",
      str_c(setdiff(current_vars, var), collapse = " + "),
      "+",
      str_c(current_covs, collapse = " + ")
    )
    reduced_formula = as.formula(reduced_formula_str)
    reduced_model = lm(reduced_formula, data = df_sub)
    reduced_r2 = glance(reduced_model)$r.squared
    
    delta_r2 = full_r2 - reduced_r2
    
    anova_test = anova(reduced_model, full_model)
    incremental_p_value = anova_test$`Pr(>F)`[2]
    
    out = tibble(r2 = delta_r2, 
                 p_value = incremental_p_value)
    
    return(out)
  })
  
  var_percentages = (var_contributions$r2 / lifestyle_total_r2) * 100
  
  coef_summary = tidy(full_model)
  
  result_df = tibble(
    age_group = current_age_group,
    variable = current_vars,
    r2_contribution = var_contributions$r2,
    r2_percentage = var_percentages,
    beta = coef_summary$estimate[(1:length(current_vars)) + 1],
    p_value = var_contributions$p_value,
    n = nrow(df_sub)
  )
  
  result_df = result_df %>% arrange(desc(r2_percentage))
  
  importance_results[[current_age_group]] = result_df
}

all_importance = bind_rows(importance_results)

out_importance = all_importance %>% 
  left_join(tibble(variable = outcomes, label = labels), by = 'variable')
write_xlsx(out_importance, path = 'lifestyle_imp_group.xlsx')


# Organ -------------------------------------------------------------------


organ_names = c("Adipose", "Artery", "Brain", "Heart", "Immune", "Kidney", 
                "Liver", "Muscle")
organ_clock = str_c('age_dev_', organ_names)

organ_importance = tibble()

for(i in organ_clock){
  
  importance_results = list()
  
  for(current_age_group in c("3-18", "18-65", "65-90", "90-118")) {
    
    cat("Age group being processed:", current_age_group, "\n")
    
    current_vars = age_group_vars[[current_age_group]]
    current_covs = cov_group_vars[[current_age_group]]
    
    df_sub = dto %>% 
      filter(age_group == current_age_group) %>%
      select(all_of(c(i, current_vars, current_covs))) %>%
      drop_na()
    
    # Variable filtering
    current_vars = simple_one_step_screening(df_sub, current_covs, current_vars, i)
    
    formula_str = str_c(
      i, ' ~ ',
      str_c(current_vars, collapse = " + "),
      "+",
      str_c(current_covs, collapse = " + ")
    )
    
    formula_full = as.formula(formula_str)
    full_model = lm(formula_full, data = df_sub)
    full_r2 = glance(full_model)$r.squared
    
    # base model
    base_formula_str = str_c(i, ' ~ ', str_c(current_covs, collapse = " + "))
    base_formula = as.formula(base_formula_str)
    base_model = lm(base_formula, data = df_sub)
    base_r2 = glance(base_model)$r.squared
    
    lifestyle_total_r2 = max(0, full_r2 - base_r2)
    
    # Calculate the incremental R² for each variable
    var_contributions = map_dfr(current_vars, function(var) {

      reduced_formula_str = str_c(
        i, ' ~ ',
        str_c(setdiff(current_vars, var), collapse = " + "),
        "+",
        str_c(current_covs, collapse = " + ")
      )
      reduced_formula = as.formula(reduced_formula_str)
      reduced_model = lm(reduced_formula, data = df_sub)
      reduced_r2 = glance(reduced_model)$r.squared
      
      delta_r2 = full_r2 - reduced_r2
      
      anova_test = anova(reduced_model, full_model)
      incremental_p_value = anova_test$`Pr(>F)`[2]
      
      out = tibble(r2 = delta_r2, 
                   p_value = incremental_p_value)
      
      return(out)
    })
    
    var_percentages = (var_contributions$r2 / lifestyle_total_r2) * 100
    
    coef_summary = tidy(full_model)
    
    result_df = tibble(
      organ = organ_names[organ_clock == i],
      age_group = current_age_group,
      variable = current_vars,
      r2_contribution = var_contributions$r2,
      r2_percentage = var_percentages,
      beta = coef_summary$estimate[(1:length(current_vars)) + 1],
      p_value = var_contributions$p_value,
      n = nrow(df_sub)
    )
    
    result_df = result_df %>% arrange(desc(r2_percentage))
    
    importance_results[[current_age_group]] = result_df
  }
  
  i_importance = bind_rows(importance_results)
  
  organ_importance = bind_rows(organ_importance, i_importance)
}

organ_importance_s = organ_importance %>% 
  group_by(organ, age_group) %>% 
  filter(r2_percentage == max(r2_percentage)) %>% 
  ungroup() %>% 
  left_join(tibble(variable = outcomes, label = labels), by = 'variable') %>% 
  mutate(sign = ifelse(beta > 0, 1, -1),
         LF = ifelse(p_value < 0.05, '*', '')) %>% 
  select(1:3, 9, everything())
write_xlsx(organ_importance_s, path = 'lifestyle_imp_organ.xlsx')

