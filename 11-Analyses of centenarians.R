rm(list = ls())

library(tidyverse)
library(writexl)
library(broom)
library(ggplotify)
library(janitor)
library(aplot)
library(patchwork)
library(rstatix)
library(janitor)
library(readxl)
library(tidyplots)
library(survival)
library(survminer)


# data --------------------------------------------------------------------

load('data.Rdata')
# This tibble named “data” contains all variables for 5,636 individuals, except for the β values of CpG sites.

organ_names = c("Adipose", "Artery", "Brain", "Heart", "Immune", "Kidney", 
                "Liver", "Muscle")

dt = data %>% 
  filter(calendar_age >= 100)

# MAE ---------------------------------------------------------------------

dt1 = dt %>% 
  select(id, age_dev, 37:44)
colnames(dt1) = c("id", "Lifespan DNAmAge", "Horvath1", "Horvath2", "Hannum", "PCHorvath1", 
                  "PCHorvath2", "PCHannum", "PhenoAge", "PCPhenoAge")
dt1 = dt1 %>% 
  pivot_longer(cols = -id, names_to = 'clock', values_to = 'age_dev') %>% 
  group_nest(clock) %>% 
  {
    left_join(tibble(clock = rev(c("Lifespan DNAmAge", "Horvath1", "PCHorvath1", 
                                   "Horvath2", "PCHorvath2", "Hannum", 
                                   "PCHannum", "PhenoAge", "PCPhenoAge"))), ., by = 'clock')
  } %>% 
  mutate(MAE = map_dbl(data, ~ mean(abs(.x$age_dev)))) %>% 
  select(-data) %>% 
  mutate(clock = factor(clock, levels = rev(c("Lifespan DNAmAge", "Horvath1", "PCHorvath1", 
                                              "Horvath2", "PCHorvath2", "Hannum", 
                                              "PCHannum", "PhenoAge", "PCPhenoAge"))))


# AgeDev distribution of organ clocks -------------------------------------

## 100+ ----
dt31 = dt %>% 
  select(id, calendar_age, 6:13) %>% 
  arrange(calendar_age)
colnames(dt31) = c("id", "calendar_age", 
                   "Brain", "Immune", "Muscle", 
                   "Artery", "Liver", "Heart", "Adipose", "Kidney")

dt31 = dt31 %>% 
  pivot_longer(cols = -c(id, calendar_age), 
               names_to = 'clock', values_to = 'age_dev') %>% 
  group_by(clock) %>% 
  summarise(mean = mean(age_dev), sd = sd(age_dev), 
            median = median(age_dev), IQR = IQR(age_dev)) %>% 
  ungroup() %>% 
  mutate(cate = 'Centenarians')

## 90+ ----
dt32 = data %>% 
  filter(calendar_age >= 90 & calendar_age < 100) %>% 
  select(id, calendar_age, 6:13) %>% 
  arrange(calendar_age)
colnames(dt32) = c("id", "calendar_age", 
                   "Brain", "Immune", "Muscle", 
                   "Artery", "Liver", "Heart", "Adipose", "Kidney")

dt32 = dt32 %>% 
  pivot_longer(cols = -c(id, calendar_age), 
               names_to = 'clock', values_to = 'age_dev') %>% 
  group_by(clock) %>% 
  summarise(mean = mean(age_dev), sd = sd(age_dev), 
            median = median(age_dev), IQR = IQR(age_dev)) %>% 
  ungroup() %>% 
  mutate(cate = 'Nonagenarians')

## 80+ ----
dt33 = data %>% 
  filter(calendar_age >= 80 & calendar_age < 90) %>% 
  select(id, calendar_age, 6:13) %>% 
  arrange(calendar_age)
colnames(dt33) = c("id", "calendar_age", 
                   "Brain", "Immune", "Muscle", 
                   "Artery", "Liver", "Heart", "Adipose", "Kidney")

dt33 = dt33 %>% 
  pivot_longer(cols = -c(id, calendar_age), 
               names_to = 'clock', values_to = 'age_dev') %>% 
  group_by(clock) %>% 
  summarise(mean = mean(age_dev), sd = sd(age_dev), 
            median = median(age_dev), IQR = IQR(age_dev)) %>% 
  ungroup() %>% 
  mutate(cate = 'Octogenarians')

## 65+ ----
dt34 = data %>% 
  filter(calendar_age >= 65 & calendar_age < 80) %>% 
  select(id, calendar_age, 6:13) %>% 
  arrange(calendar_age)
colnames(dt34) = c("id", "calendar_age", 
                   "Brain", "Immune", "Muscle", 
                   "Artery", "Liver", "Heart", "Adipose", "Kidney")

dt34 = dt34 %>% 
  pivot_longer(cols = -c(id, calendar_age), 
               names_to = 'clock', values_to = 'age_dev') %>% 
  group_by(clock) %>% 
  summarise(mean = mean(age_dev), sd = sd(age_dev), 
            median = median(age_dev), IQR = IQR(age_dev)) %>% 
  ungroup() %>% 
  mutate(cate = 'Young-old adults \n (65-80)')


## 25-65 ----
dt35 = data %>% 
  filter(calendar_age >= 25 & calendar_age < 65) %>% 
  select(id, calendar_age, 6:13) %>% 
  arrange(calendar_age)
colnames(dt35) = c("id", "calendar_age", 
                   "Brain", "Immune", "Muscle", 
                   "Artery", "Liver", "Heart", "Adipose", "Kidney")

dt35 = dt35 %>% 
  pivot_longer(cols = -c(id, calendar_age), 
               names_to = 'clock', values_to = 'age_dev') %>% 
  group_by(clock) %>% 
  summarise(mean = mean(age_dev), sd = sd(age_dev), 
            median = median(age_dev), IQR = IQR(age_dev)) %>% 
  ungroup() %>% 
  mutate(cate = 'Midlife adults \n (25-65)')


## 18-25 ----
dt36 = data %>% 
  filter(calendar_age >= 18 & calendar_age < 25) %>% 
  select(id, calendar_age, 6:13) %>% 
  arrange(calendar_age)
colnames(dt36) = c("id", "calendar_age", 
                   "Brain", "Immune", "Muscle", 
                   "Artery", "Liver", "Heart", "Adipose", "Kidney")

dt36 = dt36 %>% 
  pivot_longer(cols = -c(id, calendar_age), 
               names_to = 'clock', values_to = 'age_dev') %>% 
  group_by(clock) %>% 
  summarise(mean = mean(age_dev), sd = sd(age_dev), 
            median = median(age_dev), IQR = IQR(age_dev)) %>% 
  ungroup() %>% 
  mutate(cate = 'Emerging adults \n (18-25)')

## 3-18 ----
dt37 = data %>% 
  filter(calendar_age >= 3 & calendar_age < 18) %>% 
  select(id, calendar_age, 6:13) %>% 
  arrange(calendar_age)
colnames(dt37) = c("id", "calendar_age", 
                   "Brain", "Immune", "Muscle", 
                   "Artery", "Liver", "Heart", "Adipose", "Kidney")

dt37 = dt37 %>% 
  pivot_longer(cols = -c(id, calendar_age), 
               names_to = 'clock', values_to = 'age_dev') %>% 
  group_by(clock) %>% 
  summarise(mean = mean(age_dev), sd = sd(age_dev), 
            median = median(age_dev), IQR = IQR(age_dev)) %>% 
  ungroup() %>% 
  mutate(cate = 'Children and \n adolescents (3-18)') 

## Merge ----
dt3 = bind_rows(dt31, dt32, dt33, dt34, dt35, dt36, dt37) %>% 
  mutate(clock = factor(clock, levels = rev(c("Adipose", "Artery", 
                                              "Brain", "Heart", "Immune", "Kidney", 
                                              "Liver", "Muscle"))),
         cate = factor(cate, levels = c("Children and \n adolescents (3-18)",
                                        "Emerging adults \n (18-25)",
                                        "Midlife adults \n (25-65)",
                                        "Young-old adults \n (65-80)",
                                        "Octogenarians",
                                        "Nonagenarians", 
                                        "Centenarians")))

# Hierarchical clustering ----

## Compute the Euclidean distance matrix
pred_matrix = dt31 %>% 
  select(all_of(organ_names)) %>% 
  as.matrix()
hclust_result = hclust(dist(t(pred_matrix), method = "manhattan"))


library(ggdendro)

ggdendrogram = function (data, segments = TRUE, labels = TRUE, leaf_labels = TRUE, 
                         rotate = FALSE, theme_dendro = TRUE, ...) 
{
  dataClass <- if (inherits(data, "dendro")) 
    data$class
  else class(data)
  angle <- if (dataClass %in% c("dendrogram", "hclust")) {
    ifelse(rotate, 0, 90)
  }
  else {
    ifelse(rotate, 90, 0)
  }
  hjust <- if (dataClass %in% c("dendrogram", "hclust")) {
    ifelse(rotate, 1, 1)
  }
  else {
    0.5
  }
  if (!is.dendro(data)) 
    data <- dendro_data(data)
  p <- ggplot() + geom_blank()
  if (segments && !is.null(data$segments)) {
    p <- p + geom_segment(data = segment(data), aes(x = .data[["x"]], 
                                                    y = .data[["y"]], xend = .data[["xend"]], yend = .data[["yend"]]),
                          linewidth = 0.1)
  }
  if (leaf_labels && !is.null(data$leaf_labels)) {
    p <- p + geom_text(data = leaf_label(data), aes(x = .data[["x"]], 
                                                    y = .data[["y"]], label = .data[["label"]]), hjust = hjust, 
                       angle = angle, ...)
  }
  if (labels) {
    p <- p + scale_x_continuous(breaks = seq_along(data$labels$label), 
                                labels = data$labels$label)
  }
  if (rotate) {
    p <- p + coord_flip()
    p <- p + scale_y_continuous()
  }
  else {
    p <- p + scale_y_continuous()
  }
  if (theme_dendro) 
    p <- p + theme_dendro()
  p <- p + theme(axis.text.x = element_text(angle = angle, 
                                            hjust = 1, vjust = 0.5)) + theme(axis.text.y = element_text(angle = angle, 
                                                                                                        hjust = 1))
  p
}


dendr = dendro_data(hclust_result, type = "rectangle")

hclust_plot = ggdendrogram(dendr, rotate = TRUE) + 
  labs(title = 'Hierarchical clustering') + 
  theme(
    axis.text.x = element_blank(), 
    axis.text.y = element_text(color = 'black', family = 'sans', size = 5),
    plot.title = element_text(color = 'black', family = 'sans', size = 5, hjust = 0.5)
  )

# The youngest organ ------------------------------------------------------

## 100+ ----
dt41 = dt %>% 
  select(id, calendar_age, 6:13) %>% 
  arrange(calendar_age)
colnames(dt41) = c("id", "calendar_age", 
                   "Brain", "Immune", "Muscle", 
                   "Artery", "Liver", "Heart", "Adipose", "Kidney")

dt41 = dt41 %>% 
  pivot_longer(cols = -c(id, calendar_age), 
               names_to = 'clock', values_to = 'age_dev') %>% 
  group_by(id) %>% 
  mutate(most = clock[which.min(age_dev)]) %>% 
  ungroup() %>% 
  mutate(cate = 'Centenarians') %>% 
  distinct(id, calendar_age, most, cate)


## 90+ ----
dt42 = data %>% 
  filter(calendar_age >= 90 & calendar_age < 100) %>% 
  select(id, calendar_age, 6:13) %>% 
  arrange(calendar_age)
colnames(dt42) = c("id", "calendar_age", 
                   "Brain", "Immune", "Muscle", 
                   "Artery", "Liver", "Heart", "Adipose", "Kidney")

dt42 = dt42 %>% 
  pivot_longer(cols = -c(id, calendar_age), 
               names_to = 'clock', values_to = 'age_dev') %>% 
  group_by(id) %>% 
  mutate(most = clock[which.min(age_dev)]) %>% 
  ungroup() %>% 
  mutate(cate = 'Nonagenarians') %>% 
  distinct(id, calendar_age, most, cate)

## 80+ ----
dt43 = data %>% 
  filter(calendar_age >= 80 & calendar_age < 90) %>% 
  select(id, calendar_age, 6:13) %>% 
  arrange(calendar_age)
colnames(dt43) = c("id", "calendar_age", 
                   "Brain", "Immune", "Muscle", 
                   "Artery", "Liver", "Heart", "Adipose", "Kidney")

dt43 = dt43 %>% 
  pivot_longer(cols = -c(id, calendar_age), 
               names_to = 'clock', values_to = 'age_dev') %>% 
  group_by(id) %>% 
  mutate(most = clock[which.min(age_dev)]) %>% 
  ungroup() %>% 
  mutate(cate = 'Octogenarians') %>% 
  distinct(id, calendar_age, most, cate)

## 65+ ----
dt44 = data %>% 
  filter(calendar_age >= 65 & calendar_age < 80) %>% 
  select(id, calendar_age, 6:13) %>% 
  arrange(calendar_age)
colnames(dt44) = c("id", "calendar_age", 
                   "Brain", "Immune", "Muscle", 
                   "Artery", "Liver", "Heart", "Adipose", "Kidney")

dt44 = dt44 %>% 
  pivot_longer(cols = -c(id, calendar_age), 
               names_to = 'clock', values_to = 'age_dev') %>% 
  group_by(id) %>% 
  mutate(most = clock[which.min(age_dev)]) %>% 
  ungroup() %>% 
  mutate(cate = 'Young-old adults \n (65-80)') %>% 
  distinct(id, calendar_age, most, cate)

## 25-65 ----
dt45 = data %>% 
  filter(calendar_age >= 25 & calendar_age < 65) %>% 
  select(id, calendar_age, 6:13) %>% 
  arrange(calendar_age)
colnames(dt45) = c("id", "calendar_age", 
                   "Brain", "Immune", "Muscle", 
                   "Artery", "Liver", "Heart", "Adipose", "Kidney")

dt45 = dt45 %>% 
  pivot_longer(cols = -c(id, calendar_age), 
               names_to = 'clock', values_to = 'age_dev') %>% 
  group_by(id) %>% 
  mutate(most = clock[which.min(age_dev)]) %>% 
  ungroup() %>% 
  mutate(cate = 'Midlife adults \n (25-65)') %>% 
  distinct(id, calendar_age, most, cate)


## 18-25 ----
dt46 = data %>% 
  filter(calendar_age >= 18 & calendar_age < 25) %>% 
  select(id, calendar_age, 6:13) %>% 
  arrange(calendar_age)
colnames(dt46) = c("id", "calendar_age", 
                   "Brain", "Immune", "Muscle", 
                   "Artery", "Liver", "Heart", "Adipose", "Kidney")

dt46 = dt46 %>% 
  pivot_longer(cols = -c(id, calendar_age), 
               names_to = 'clock', values_to = 'age_dev') %>% 
  group_by(id) %>% 
  mutate(most = clock[which.min(age_dev)]) %>% 
  ungroup() %>% 
  mutate(cate = 'Emerging adults \n (18-25)') %>% 
  distinct(id, calendar_age, most, cate)


## 3-18 ----
dt47 = data %>% 
  filter(calendar_age < 18) %>% 
  select(id, calendar_age, 6:13) %>% 
  arrange(calendar_age)
colnames(dt47) = c("id", "calendar_age", 
                   "Brain", "Immune", "Muscle", 
                   "Artery", "Liver", "Heart", "Adipose", "Kidney")

dt47 = dt47 %>% 
  pivot_longer(cols = -c(id, calendar_age), 
               names_to = 'clock', values_to = 'age_dev') %>% 
  group_by(id) %>% 
  mutate(most = clock[which.min(age_dev)]) %>% 
  ungroup() %>% 
  mutate(cate = 'Children and \n adolescents (3-18)') %>% 
  distinct(id, calendar_age, most, cate)


## Merge ----
dt4 = bind_rows(dt41, dt42, dt43, dt44, dt45, dt46, dt47) %>% 
  mutate(most = factor(most, levels = c("Adipose", "Artery", 
                                        "Brain", "Heart", "Immune", "Kidney", 
                                        "Liver", "Muscle")),
         cate = factor(cate, levels = c("Children and \n adolescents (3-18)",
                                        "Emerging adults \n (18-25)",
                                        "Midlife adults \n (25-65)",
                                        "Young-old adults \n (65-80)",
                                        "Octogenarians",
                                        "Nonagenarians", 
                                        "Centenarians"))) %>% 
  count(cate, most) %>% 
  group_by(cate) %>% 
  mutate(N = sum(n),
         per = n/N*100) %>% 
  ungroup()


# Adipose-Kidney-Muscle ---------------------------------------------------

# The “XXX_ageotype” variables represent slow agers, normal agers, and fast agers in different organs, as classified based on the MAE of their respective models.
dt51 = dt %>% 
  mutate(triad = case_when(Adipose_ageotype == 'Slow agers' & 
                             Kidney_ageotype == 'Slow agers' & 
                             Muscle_ageotype == 'Slow agers' ~ 1,
                           TRUE ~ 0)) %>% 
  mutate(cate = 'Centenarians')
count(dt51, triad)

dt52 = data %>% 
  filter(calendar_age >= 90 & calendar_age < 100) %>% 
  mutate(triad = case_when(Adipose_ageotype == 'Slow agers' & 
                             Kidney_ageotype == 'Slow agers' & 
                             Muscle_ageotype == 'Slow agers' ~ 1,
                           TRUE ~ 0)) %>% 
  mutate(cate = 'Nonagenarians')
count(dt52, triad)

dt53 = data %>% 
  filter(calendar_age >= 80 & calendar_age < 90) %>% 
  mutate(triad = case_when(Adipose_ageotype == 'Slow agers' & 
                             Kidney_ageotype == 'Slow agers' & 
                             Muscle_ageotype == 'Slow agers' ~ 1,
                           TRUE ~ 0)) %>% 
  mutate(cate = 'Octogenarians')
count(dt53, triad)

dt54 = data %>% 
  filter(calendar_age >= 65 & calendar_age < 80) %>% 
  mutate(triad = case_when(Adipose_ageotype == 'Slow agers' & 
                             Kidney_ageotype == 'Slow agers' & 
                             Muscle_ageotype == 'Slow agers' ~ 1,
                           TRUE ~ 0)) %>% 
  mutate(cate = 'Young-old adults \n (65-80)')
count(dt54, triad)

dt55 = data %>% 
  filter(calendar_age >= 25 & calendar_age < 65) %>% 
  mutate(triad = case_when(Adipose_ageotype == 'Slow agers' & 
                             Kidney_ageotype == 'Slow agers' & 
                             Muscle_ageotype == 'Slow agers' ~ 1,
                           TRUE ~ 0)) %>% 
  mutate(cate = 'Midlife adults \n (25-65)')
count(dt55, triad)

dt56 = data %>% 
  filter(calendar_age >= 18 & calendar_age < 25) %>% 
  mutate(triad = case_when(Adipose_ageotype == 'Slow agers' & 
                             Kidney_ageotype == 'Slow agers' & 
                             Muscle_ageotype == 'Slow agers' ~ 1,
                           TRUE ~ 0)) %>% 
  mutate(cate = 'Emerging adults \n (18-25)')
count(dt56, triad)

dt57 = data %>% 
  filter(calendar_age < 18) %>% 
  mutate(triad = case_when(Adipose_ageotype == 'Slow agers' & 
                             Kidney_ageotype == 'Slow agers' & 
                             Muscle_ageotype == 'Slow agers' ~ 1,
                           TRUE ~ 0)) %>% 
  mutate(cate = 'Children and \n adolescents (3-18)')
count(dt57, triad)

dt5 = bind_rows(dt51, dt52, dt53, dt54, dt55, dt56, dt57) %>% 
  mutate(cate = factor(cate, levels = c("Children and \n adolescents (3-18)",
                                        "Emerging adults \n (18-25)",
                                        "Midlife adults \n (25-65)",
                                        "Young-old adults \n (65-80)",
                                        "Octogenarians",
                                        "Nonagenarians", 
                                        "Centenarians")))

dtt5 = dt5 %>% 
  select(id, triad, cate) %>% 
  group_by(cate) %>% 
  summarise(n = sum(triad),
            N = n(),
            per = n/N*100) %>% 
  ungroup()


# Association analyses ----------------------------------------------------

dt6 = dt5 %>% 
  mutate(death = case_when(survival_time <= 0 ~ NA_real_,
                           survival_time > 0 ~ death)) %>% 
  filter(calendar_age >= 45 & !is.na(death)) %>% 
  mutate(triad = as_factor(triad)) %>% 
  mutate(death = as.numeric(death) - 1)

## Cox ----
dt60 = dt6 %>% 
  mutate(cate = 'Overall \n (45-118)')

dt61 = dt6 %>%
  filter(cate == 'Centenarians')

dt62 = dt6 %>%
  filter(cate == 'Nonagenarians')

dt63 = dt6 %>%
  filter(cate == 'Octogenarians')

dt64 = dt6 %>%
  filter(cate == 'Young-old adults \n (65-80)')

dt65 = dt6 %>%
  filter(cate == 'Midlife adults \n (25-65)') %>% 
  mutate(cate = 'Midlife adults \n (45-65)')


### function ----
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
      complete(group, y, fill = list(n = 0)) %>% 
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
    select(data, everything()) %>% 
    mutate(x = group_label)
  
  return(event_dt)
}
mul_cox_fun = function(data, data_label, x, x_label, outcome, outtime, 
                       label_outcome, cov, cov_plus, cov_label, weight){
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
  cox_dt = cox_dt %>% 
    mutate(x = x_label)
  return(cox_dt)
}


### rate ---- 

outcomes = c("death")
outtimes = c("survival_time")
labels_outcomes = c("All-cause mortality")

rate_map_dt = expand_grid(d1 = tibble(data = list(dt60, dt61, dt62, dt63, dt64, dt65), 
                                      data_label = rev(c("Midlife adults \n (45-65)",
                                                         "Young-old adults \n (65-80)",
                                                         "Octogenarians",
                                                         "Nonagenarians", 
                                                         "Centenarians",
                                                         "Overall \n (45-118)"))),
                          tibble(outcomes, outtimes, labels_outcomes,
                                 group = c("triad"), 
                                 label = c('Adipose-Kidney-Muscle slow agers'))
) %>% 
  unnest(d1)

rate_dt = pmap_dfr(rate_map_dt, ~ rate_fun(data = ..1, data_label = ..2, 
                                           outcome = ..3, outtime = ..4, 
                                           group = ..6, label = ..5,
                                           group_label = ..7))

### association ----
x = c("triad")
x_label = c('Adipose-Kidney-Muscle slow agers')
cov0 = NULL
cov1 = c("calendar_age", "sex")
cov2 = c("calendar_age", "sex", "education", "marriage")
cov3 = c("calendar_age", "sex", "education", "marriage", 
         "BMI")
cov4 = c("calendar_age", "sex", "education", "marriage", 
         "BMI", "smoke", "hypertension")

cov_label = str_c('Model ', 0:4)

cox_map_dt = expand_grid(d1 = tibble(data = list(dt60, dt61, dt62, dt63, dt64, dt65), 
                                     data_label = rev(c("Midlife adults \n (45-65)",
                                                        "Young-old adults \n (65-80)",
                                                        "Octogenarians",
                                                        "Nonagenarians", 
                                                        "Centenarians",
                                                        "Overall \n (45-118)")),
                                     cov_plus = map(1:6, ~ c("CD8T", "CD4T", "NK", "Bcell", "Mono", "Neu")),
),
tibble(outcomes, outtimes, labels_outcomes, x, x_label),
tibble(cov = list(cov0, cov1, cov2, cov3, cov4), cov_label),
tibble(weight = list(NULL))) %>% 
  unnest(cols = c(d1))

mul_cox_dt = pmap_dfr(cox_map_dt[1:20,], ~ mul_cox_fun(data = ..1, data_label = ..2, cov_plus = ..3,
                                                       x = ..7, x_label = ..8, outcome = ..4,
                                                       outtime = ..5, label_outcome = ..6,
                                                       cov = ..9, cov_label = ..10, weigh = ..11))

output_triad = full_join(mul_cox_dt, rate_dt, by = c('data', 'outcome', 'group', 'x')) %>% 
  select(data, outcome, x, model, group, `Event/Total (%)`, `HR (95% CI)`, 
         p_value, ph_p_value, HR, low, up, p.value) %>% 
  pivot_wider(id_cols = c(data, outcome, x, group, `Event/Total (%)`),
              names_from = model, 
              values_from = c(`HR (95% CI)`, p_value, ph_p_value, HR, low, up, p.value),
              names_glue = '{model}_{.value}') %>% 
  relocate(data, outcome, x, group, `Event/Total (%)`, 
           contains('Model 0'),
           contains('Model 1'), contains('Model 2'),
           contains('Model 3'), contains('Model 4'))

output_triad0 = full_join(mul_cox_dt, rate_dt, by = c('data', 'outcome', 'group', 'x')) %>% 
  select(data, outcome, x, model, group, `Event/Total (%)`, `HR (95% CI)`, 
         p_value, ph_p_value) %>% 
  pivot_wider(id_cols = c(data, outcome, x, group, `Event/Total (%)`),
              names_from = model, 
              values_from = c(`HR (95% CI)`, p_value, ph_p_value),
              names_glue = '{model}_{.value}') %>% 
  relocate(data, outcome, x, group, `Event/Total (%)`, 
           contains('Model 0'),
           contains('Model 1'), contains('Model 2'),
           contains('Model 3'), contains('Model 4'))


write_xlsx(list(output_triad, output_triad0), path = 'Centenarians_cox.xlsx')


## lm ----
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

dt66 = dt5 %>% 
  filter(calendar_age >= 65)


outcomes = c("BADL", "IADL", "MMSE",
             "HGS_L", "HGS_R", "s4", "s6", 
             "FI")
labels = c("BADL", "IADL", "MMSE",
           "Hand grip strength (left)", "Hand grip strength (right)",
           "4‒meter gait speed", "6‒meter gait speed", 
           "FI")

x = c("triad")
x_label = c("Adipose-Kidney-Muscle slow agers")
cov0 = NULL
cov1 = c("calendar_age", "sex")
cov2 = c("calendar_age", "sex", "education", "marriage")
cov3 = c("calendar_age", "sex", "education", "marriage", 
         "BMI")
cov4 = c("calendar_age", "sex", "education", "marriage", 
         "BMI", "smoke", "hypertension")

cov_label = str_c('Model ', 0:4)

cov_plus = map(1:length(outcomes), 
               ~ c("CD8T", "CD4T", "NK", "Bcell", "Mono", "Neu"))

lm_map_dt = expand_grid(d1 = tibble(data = map(outcomes, ~ drop_na(dt66, all_of(.x))), 
                                    data_label = rep('Total', length(c(outcomes))),
                                    cov_plus = cov_plus,
                                    outcome = c(outcomes), 
                                    label_outcome = c(labels)),
                        tibble(x, x_label),
                        tibble(cov = list(cov0, cov1, cov2, cov3, cov4), cov_label)) %>% 
  unnest(cols = c(d1)) %>% 
  mutate(N = map_int(data, ~ nrow(.)))


mul_lm_dt = pmap_dfr(lm_map_dt, ~ mul_lm_fun(data = ..1, data_label = ..2, cov_plus = ..3,
                                             x = ..6, x_label = ..7, outcome = ..4,
                                             label_outcome = ..5,
                                             cov = ..8, cov_label = ..9))

output = mul_lm_dt %>% 
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

output0 = mul_lm_dt %>% 
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

write_xlsx(list(output, output0), path = 'Centenarians_lm.xlsx')
