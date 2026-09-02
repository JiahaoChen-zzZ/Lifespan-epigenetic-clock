rm(list = ls())

library(tidyverse)
library(readr)
library(readxl)
library(writexl)
library(mlr3verse)
library(future)
library(broom)
library(rstatix)
library(modelr)
library(mlr3mbo)
library(bbotk)
library(mlr3extralearners)
library(DiceKriging)
library(randomForest)
library(lightgbm)
library(fastshap)
library(ggrepel)
library(tidyplots)
library(ggstar)
library(patchwork)
library(Hmisc)
library(compareGroups)

# load --------------------------------------------------------------------

load('DMPs.Rdata')
# This tibble named “DMPs” contains 5,636 individuals, a “Sample_Name” column, and 3,237 CpG columns.
load('data.Rdata')
# This tibble named “data” contains all variables for 5,636 individuals, except for the β values of CpG sites.


# function ----------------------------------------------------------------

boruta_shap_lgbm = function(X, y, n_trials = 20, percentile = 100, p_value = 0.05, seed = 123,
                            Model) {
  set.seed(seed)
  
  tentative = colnames(X)
  out = tibble(Trial = integer(0), `Number of CpG` = integer(0), MAE = numeric(0))
  
  # MAE
  mae_metric = function(preds, dtrain) {
    labels = get_field(dtrain, "label")
    mae = mean(abs(labels - preds))
    return(list(
      name = "mae",
      value = mae,
      higher_better = FALSE
    ))
  }
  
  # 80% train, 20% test
  task = X %>% 
    mutate(Y = y) %>% 
    as_task_regr(target = "Y", id = 'clock')
  split = partition(task, ratio = 0.8)
  train_idx = split$train
  
  for (trial in 1:n_trials) {
    X_train = X[train_idx, tentative]
    y_train = y[train_idx]
    X_valid = X[-train_idx, tentative]
    y_valid = y[-train_idx]
    
    # 
    shadow_data = X_train %>% 
      mutate(across(everything(), ~ sample(.x, size = length(.x)))) %>% 
      setNames(paste0("shadow_", names(X_train)))
    
    # 
    dtrain = lgb.Dataset(
      data = as.matrix(bind_cols(X_train, shadow_data)),
      label = y_train
    )
    
    # 
    dvalid = lgb.Dataset.create.valid(
      dataset = dtrain,
      data = as.matrix(bind_cols(X_valid, shadow_data[1:nrow(X_valid), ])),
      label = y_valid
    )
    
    # 
    params = list(
      num_leaves = Model$param_set$values$num_leaves,
      max_depth = Model$param_set$values$max_depth,
      min_data_in_leaf = Model$param_set$values$min_data_in_leaf,
      
      learning_rate = Model$param_set$values$learning_rate,
      
      
      lambda_l1 = Model$param_set$values$lambda_l1,
      lambda_l2 = Model$param_set$values$lambda_l2,
      
      feature_fraction = Model$param_set$values$feature_fraction,
      bagging_fraction = Model$param_set$values$bagging_fraction,
      bagging_freq = Model$param_set$values$bagging_freq,
      
      objective = "regression",
      metric = 'custom',
      boosting = "gbdt"
      
    )
    model = lgb.train(
      params = params,
      data = dtrain,
      valids = list(valid = dvalid),
      eval = mae_metric,  
      nrounds = Model$param_set$values$num_iterations,
      verbose = -1,
      early_stopping_rounds = 20
    )
    rm(dtrain, dvalid)
    
    # SHAP
    shap_values = predict(
      model,
      as.matrix(bind_cols(X_train, shadow_data)),
      type = 'contrib',
      num_iteration = model$best_iter
    )
    rm(shadow_data)
    
    # 
    shap_importance = colMeans(abs(shap_values))
    rm(shap_values)
    
    # 
    real_importance = shap_importance[1:ncol(X_train)]
    shadow_importance = shap_importance[(ncol(X_train) + 1):(length(shap_importance)-1)]
    rm(shap_importance)
    
    # 
    threshold = quantile(shadow_importance, probs = percentile / 100)
    hits = real_importance > threshold
    rm(real_importance)
    
    # 
    rejected_new = names(X)[!hits]
    
    tentative = setdiff(tentative, rejected_new)
    
    #
    cat(sprintf("Trial %d: Best MAE = %.4f (Iteration %d)\n",
                trial, model$best_score, model$best_iter))
    cat("tentative: ", tentative, "\n", "features number = ", length(tentative), '\n')
    out = bind_rows(out, tibble(Trial = trial, `Number of CpG` = length(tentative),
                                MAE = model$best_score))
    rm(model)
    gc()
    # 
    if (length(tentative) == 0) break
  }
  
  # 
  list(tentative = tentative, out = out)
}

# data --------------------------------------------------------------------
dt = DMPs %>% 
  left_join(select(data, Sample_Name, calendar_age), by = 'Sample_Name') %>% 
  relocate(Sample_Name, calendar_age)
rm(DMPs)

colnames(dt) = str_replace_all(colnames(dt), '-', '_')

dtm = select(dt, Sample_Name) %>% 
  mutate(row_ids = 1:nrow(.))
class(dtm$row_ids)

dt = dt %>% 
  select(-Sample_Name)

# task --------------------------------------------------------------------

task_all = as_task_regr(dt, target = "calendar_age", id = 'clock')
task = as_task_regr(dt, target = "calendar_age", id = 'clock')
task

min_age = min(dt$calendar_age)
max_age = max(dt$calendar_age)

# Tuner ---------------------------------------------------------------------

## Bayesian optimization ----
surrogate = srlrn(lrn("regr.km", covtype = "matern5_2", 
                      optim.method = "BFGS", 
                      control = list(trace = FALSE)))
acq_optimizer = acqo(opt("nloptr", algorithm = "NLOPT_GN_ORIG_DIRECT"),  
                     terminator = trm("stagnation", iters = 100, threshold = 1e-5))
tuner = tnr("mbo", 
            loop_function = bayesopt_ego, 
            surrogate = surrogate, 
            acq_function = acqf("ei"), 
            acq_optimizer = acq_optimizer)


# sample ------------------------------------------------------------------

set.seed(123)
split = partition(task, ratio = 0.7)
task$row_roles$use = split$train
task

dt_set = select(data, id, Sample_Name, calendar_age) %>% 
  left_join(dtm, by = 'Sample_Name') %>% 
  mutate(set = case_when(row_ids %in% split$train ~ 'train set',
                         row_ids %in% split$test ~ 'test set')) %>% 
  arrange(row_ids)


# organ -------------------------------------------------------------------

CpG_organ = read_xlsx('CpG_organ.xlsx')
# The data recorded in “Supplementary Table 23. Details of DMPs with Annotations of Organ Enrichment”
count(CpG_organ, GTExOrgan) %>% 
  arrange(desc(n))

# Brain -------------------------------------------------------------------
task = as_task_regr(dt, target = "calendar_age", id = 'clock')
task$row_roles$use = split$train

Brain_CpG = CpG_organ %>% 
  filter(GTExOrgan == 'Brain') %>% 
  pull(CpG)

task_Brain_all = task_all
task_Brain_all$select(cols = Brain_CpG)

task_Brain = task
task_Brain$select(cols = Brain_CpG)
task_Brain
task_Brai$row_roles$use

## benchmark ----
learner = lrn("regr.lightgbm",

              num_leaves = to_tune(20, 300),          
              max_depth = to_tune(3, 12),             
              min_data_in_leaf = to_tune(10, 100),    
              
              learning_rate = to_tune(0.01, 0.3, logscale = TRUE),  
              num_iterations = to_tune(100, 1000),    
              
              lambda_l1 = to_tune(0, 10),      
              lambda_l2 = to_tune(0, 10),      
              
              feature_fraction = to_tune(0.5, 1.0),  
              bagging_fraction = to_tune(0.5, 1.0),   
              bagging_freq = to_tune(1, 10)           
)

learner$param_set
para_lightgbm = c('num_leaves', 'max_depth', 'min_data_in_leaf',
                  'learning_rate', 'num_iterations', 
                  'lambda_l1', 'lambda_l2',
                  'feature_fraction', 'bagging_fraction', 'bagging_freq')
set.seed(123)
instance = tune(tuner = tuner, task = task_Brain, learner = learner, 
                resampling = rsmp("cv", folds = 5),  
                measure = msr("regr.mae"), term_evals = 25)
instance$result
instance$archive
learner$param_set$values = instance$result_learner_param_vals
benchmarking_best_lightgbm_Brain = learner
benchmarking_best_lightgbm_Brain$train(task_Brain, row_ids = split$train)

# predict
lightgbm_Brain_pred_train = benchmarking_best_lightgbm_Brain$predict(task_Brain, row_ids = split$train)
lightgbm_Brain_pred_test = benchmarking_best_lightgbm_Brain$predict(task_Brain_all, row_ids = split$test)
train_lightgbm_Brain = tibble(row_ids = lightgbm_Brain_pred_train$row_ids, 
                              calendar_age = lightgbm_Brain_pred_train$truth,
                              pred_calendar_age = lightgbm_Brain_pred_train$response)
test_lightgbm_Brain = tibble(row_ids = lightgbm_Brain_pred_test$row_ids, 
                             calendar_age = lightgbm_Brain_pred_test$truth,
                             pred_calendar_age = lightgbm_Brain_pred_test$response)
fit_lightgbm_Brain = lm(pred_calendar_age ~ calendar_age, data = test_lightgbm_Brain)

optimal_lightgbm_Brain = instance$result %>% 
  select(all_of(para_lightgbm)) %>% 
  mutate(selection = 'optimal')

test_prediction_Brain = test_lightgbm_Brain %>% 
  left_join(dtm, by = 'row_ids')

out_lightgbm_Brain = list(train_tuning = instance$archive$data %>% 
                            select(all_of(para_lightgbm), regr.mae) %>% 
                            left_join(optimal_lightgbm_Brain, by = para_lightgbm) %>%
                            mutate(learning_rate = exp(learning_rate)) %>% 
                            set_names(c('Max number of leaves in one tree', 
                                        'Max depth for tree model',
                                        'Minimal number of data in one leaf',
                                        'Shrinkage rate', 
                                        'Number of boosting iterations', 
                                        'L1 regularization', 
                                        'L2 regularization',
                                        'Fraction of random feature selection in each iteration', 
                                        'Fraction of random sample selection in each iteration', 
                                        'Frequency for bagging'), 
                                      'MAE', 'selection'),
                          train_prediction = train_lightgbm_Brain %>% left_join(dtm, by = 'row_ids'),
                          test_prediction = test_prediction_Brain,
                          Extrapolation = tibble(Min = min(test_prediction_Brain$pred_calendar_age),
                                                 Max = max(test_prediction_Brain$pred_calendar_age),
                                                 `n (<Min)` = sum(test_prediction_Brain$pred_calendar_age < min_age),
                                                 `n (>Max)` = sum(test_prediction_Brain$pred_calendar_age > max_age)),
                          test_evaluate = tibble(r = cor(test_lightgbm_Brain$calendar_age, test_lightgbm_Brain$pred_calendar_age, method = 'pearson'),
                                                 R2 = rsquare(fit_lightgbm_Brain, test_lightgbm_Brain),
                                                 RMSE = rmse(fit_lightgbm_Brain, test_lightgbm_Brain),
                                                 MAE = mae(fit_lightgbm_Brain, test_lightgbm_Brain)))

## feature selection ----
benchmarking_best_lightgbm_Brain$param_set$values

# data 
X = dt %>% 
  slice(split$train) %>% 
  select(all_of(Brain_CpG)) %>% 
  as.data.frame()
y = dt %>% 
  slice(split$train) %>% 
  pull(calendar_age)

result_Brain = boruta_shap_lgbm(X, y, n_trials = 200, percentile = 100, p_value = 0.05, 
                                Model = benchmarking_best_lightgbm_Brain)
fs_Brain = list(feature = tibble(CpG = result_Brain$tentative), 'feature selection process' = result_Brain$out)

## 特征选择后调参 ----
fs_Brain = fs_Brain$feature

task_Brain$select(cols = fs_Brain$CpG)
task_Brain$row_roles$use

learner = lrn("regr.lightgbm",
              
              num_leaves = to_tune(20, 300),          
              max_depth = to_tune(3, 12),             
              min_data_in_leaf = to_tune(10, 100),    
              
              learning_rate = to_tune(0.01, 0.3, logscale = TRUE),  
              num_iterations = to_tune(100, 1000),    
              
              lambda_l1 = to_tune(0, 10),      
              lambda_l2 = to_tune(0, 10),      
              
              feature_fraction = to_tune(0.5, 1.0),  
              bagging_fraction = to_tune(0.5, 1.0),   
              bagging_freq = to_tune(1, 10)           
)

learner$param_set
para_lightgbm = c('num_leaves', 'max_depth', 'min_data_in_leaf',
                  'learning_rate', 'num_iterations', 
                  'lambda_l1', 'lambda_l2',
                  'feature_fraction', 'bagging_fraction', 'bagging_freq')
set.seed(1234)
instance = tune(tuner = tuner, task = task_Brain, learner = learner, 
                resampling = rsmp("cv", folds = 5),  
                measure = msr("regr.mae"), term_evals = 25)
instance$result
instance$archive
learner$param_set$values = instance$result_learner_param_vals
fs_best_lightgbm_Brain = learner
fs_best_lightgbm_Brain$train(task_Brain, row_ids = split$train)

# predict
lightgbm_Brain_pred_train = fs_best_lightgbm_Brain$predict(task_Brain, row_ids = split$train)
lightgbm_Brain_pred_test = fs_best_lightgbm_Brain$predict(task_Brain_all, row_ids = split$test)
train_lightgbm_Brain = tibble(row_ids = lightgbm_Brain_pred_train$row_ids, 
                              calendar_age = lightgbm_Brain_pred_train$truth,
                              pred_calendar_age = lightgbm_Brain_pred_train$response)
test_lightgbm_Brain = tibble(row_ids = lightgbm_Brain_pred_test$row_ids, 
                             calendar_age = lightgbm_Brain_pred_test$truth,
                             pred_calendar_age = lightgbm_Brain_pred_test$response)
fit_lightgbm_Brain = lm(pred_calendar_age ~ calendar_age, data = test_lightgbm_Brain)

optimal_lightgbm_Brain = instance$result %>% 
  select(all_of(para_lightgbm)) %>% 
  mutate(selection = 'optimal')

test_prediction_Brain = test_lightgbm_Brain %>% 
  left_join(dtm, by = 'row_ids')

out_lightgbm_Brain = list(train_tuning = instance$archive$data %>% 
                            select(all_of(para_lightgbm), regr.mae) %>% 
                            left_join(optimal_lightgbm_Brain, by = para_lightgbm) %>%
                            mutate(learning_rate = exp(learning_rate)) %>% 
                            set_names(c('Max number of leaves in one tree', 
                                        'Max depth for tree model',
                                        'Minimal number of data in one leaf',
                                        'Shrinkage rate', 
                                        'Number of boosting iterations', 
                                        'L1 regularization', 
                                        'L2 regularization',
                                        'Fraction of random feature selection in each iteration', 
                                        'Fraction of random sample selection in each iteration', 
                                        'Frequency for bagging'), 
                                      'MAE', 'selection'),
                          train_prediction = train_lightgbm_Brain %>% left_join(dtm, by = 'row_ids'),
                          test_prediction = test_prediction_Brain,
                          Extrapolation = tibble(Min = min(test_prediction_Brain$pred_calendar_age),
                                                 Max = max(test_prediction_Brain$pred_calendar_age),
                                                 `n (<Min)` = sum(test_prediction_Brain$pred_calendar_age < min_age),
                                                 `n (>Max)` = sum(test_prediction_Brain$pred_calendar_age > max_age)),
                          test_evaluate = tibble(r = cor(test_lightgbm_Brain$calendar_age, test_lightgbm_Brain$pred_calendar_age, method = 'pearson'),
                                                 R2 = rsquare(fit_lightgbm_Brain, test_lightgbm_Brain),
                                                 RMSE = rmse(fit_lightgbm_Brain, test_lightgbm_Brain),
                                                 MAE = mae(fit_lightgbm_Brain, test_lightgbm_Brain)))

## Overall prediction ----

dt_Brain = dt %>% 
  select(calendar_age, all_of(fs_Brain$CpG))

# sample --
task_Brain = as_task_regr(dt_Brain, target = "calendar_age", id = 'clock')

# learner --
learner = lrn("regr.lightgbm",
              
              num_leaves = fs_best_lightgbm_Brain$param_set$values$num_leaves,
              max_depth = fs_best_lightgbm_Brain$param_set$values$max_depth,
              min_data_in_leaf = fs_best_lightgbm_Brain$param_set$values$min_data_in_leaf,
              
              learning_rate = fs_best_lightgbm_Brain$param_set$values$learning_rate,
              num_iterations = fs_best_lightgbm_Brain$param_set$values$num_iterations,
              
              lambda_l1 = fs_best_lightgbm_Brain$param_set$values$lambda_l1,
              lambda_l2 = fs_best_lightgbm_Brain$param_set$values$lambda_l2,
              
              feature_fraction = fs_best_lightgbm_Brain$param_set$values$feature_fraction,
              bagging_fraction = fs_best_lightgbm_Brain$param_set$values$bagging_fraction,
              bagging_freq = fs_best_lightgbm_Brain$param_set$values$bagging_freq
)

# prediction --
set.seed(123)
cv10 = rsmp("cv", folds = 5)
cv10$instantiate(task_Brain)

rr = mlr3::resample(task = task_Brain, learner = learner, resampling = cv10, store_models = TRUE)

pred = tibble(row_ids = rr$prediction()$row_ids, 
              calendar_age = rr$prediction()$truth,
              pred_calendar_age = rr$prediction()$response)
pred = pred %>% 
  left_join(dtm, by = 'row_ids') %>% 
  left_join(select(data, Sample_Name,
                   CD8T, CD4T, NK, Bcell, Mono, Neu), by = 'Sample_Name')

fit_Brain = lm(pred_calendar_age ~ calendar_age, data = pred)

out_evaluate_Brain = tibble(r = cor(pred$calendar_age, pred$pred_calendar_age, method = 'pearson'),
                            R2 = rsquare(fit_Brain, pred),
                            RMSE = rmse(fit_Brain, pred),
                            MAE = mae(fit_Brain, pred))

## save ----
pred_Brain = select(pred, Sample_Name, pred_calendar_age) %>% 
  rename(pred_calendar_age_Brain = pred_calendar_age)

rm(list = c("cv10", "learner", "pre", "pred", 
            "benchmarking_best_lightgbm_Brain", 
            "fit_Brain", "fit_lightgbm_Brain", 
            "fs_best_lightgbm_Brain", "lightgbm_Brain_pred_test", 
            "lightgbm_Brain_pred_train", "lightgbm_plot_Brain", 
            "out_evaluate_Brain", "out_lightgbm_Brain", 
            "test_lightgbm_Brain", "train_lightgbm_Brain"))


# Other organ-specific clocks are constructed in the same way.
# ......

# save --------------------------------------------------------------------
data = data %>% 
  left_join(pred_Brain, by = 'Sample_Name') %>% 
  left_join(pred_Immune, by = 'Sample_Name') %>% 
  left_join(pred_Muscle, by = 'Sample_Name') %>% 
  left_join(pred_Artery, by = 'Sample_Name') %>% 
  left_join(pred_Liver, by = 'Sample_Name') %>% 
  left_join(pred_Heart, by = 'Sample_Name') %>% 
  left_join(pred_Adipose, by = 'Sample_Name') %>% 
  left_join(pred_Kidney, by = 'Sample_Name') %>% 
  relocate(contains("pred_calendar_age"))

save(data, file = 'data.Rdata')
