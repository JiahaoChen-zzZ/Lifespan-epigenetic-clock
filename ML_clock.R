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


# data --------------------------------------------------------------------
dt = DMPs %>% 
  left_join(select(data, Sample_Name, calendar_age), by = 'Sample_Name') %>% 
  relocate(Sample_Name, calendar_age)
rm(DMPs)

colnames(dt) = str_replace_all(colnames(dt), '-', '_')

dtm = select(dt, Sample_Name) %>% 
  mutate(row_ids = 1:nrow(.))

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


# benchmark ---------------------------------------------------------------

## sample ----
set.seed(123)
split = partition(task, ratio = 0.7)
task$row_roles$use = split$train
task

dt_set = select(dt_combined_imp, id, Sample_Name, calendar_age) %>% 
  left_join(dtm, by = 'Sample_Name') %>% 
  mutate(set = case_when(row_ids %in% split$train ~ 'train set',
                         row_ids %in% split$test ~ 'test set')) %>% 
  arrange(row_ids)

## lasso ----
learner = lrn("regr.glmnet", alpha = 1,
              s = to_tune(0.001, 10))
learner$param_set
para_lasso = c('s')
set.seed(123)
instance = tune(tuner = tuner, task = task, learner = learner, 
                resampling = rsmp("cv", folds = 5),  
                measure = msr("regr.mae"), term_evals = 25)
instance$result
instance$archive
learner$param_set$values = instance$result_learner_param_vals
benchmarking_best_lasso = learner
benchmarking_best_lasso$train(task, row_ids = split$train)

# predict
lasso_pred_train = benchmarking_best_lasso$predict(task, row_ids = split$train)
lasso_pred_test = benchmarking_best_lasso$predict(task_all, row_ids = split$test)
train_lasso = tibble(row_ids = lasso_pred_train$row_ids, 
                     calendar_age = lasso_pred_train$truth,
                     pred_calendar_age = lasso_pred_train$response)
test_lasso = tibble(row_ids = lasso_pred_test$row_ids, 
                    calendar_age = lasso_pred_test$truth,
                    pred_calendar_age = lasso_pred_test$response)
fit_lasso = lm(pred_calendar_age ~ calendar_age, data = test_lasso)

optimal_lasso = instance$result %>% 
  select(all_of(para_lasso)) %>% 
  mutate(selection = 'optimal')

test_prediction = test_lasso %>% 
  left_join(dtm, by = 'row_ids')

out_lasso = list(train_tuning = instance$archive$data %>% 
                   select(all_of(para_lasso), regr.mae) %>% 
                   left_join(optimal_lasso, by = para_lasso) %>% 
                   set_names(c('𝜆'), 'MAE', 'selection'),
                 train_prediction = train_lasso %>% left_join(dtm, by = 'row_ids'),
                 test_prediction = test_prediction,
                 Extrapolation = tibble(Min = min(test_prediction$pred_calendar_age),
                                        Max = max(test_prediction$pred_calendar_age),
                                        `n (<Min)` = sum(test_prediction$pred_calendar_age < min_age),
                                        `n (>Max)` = sum(test_prediction$pred_calendar_age > max_age)),
                 test_evaluate = tibble(r = cor(test_lasso$calendar_age, test_lasso$pred_calendar_age, method = 'pearson'),
                                        R2 = rsquare(fit_lasso, test_lasso),
                                        RMSE = rmse(fit_lasso, test_lasso),
                                        MAE = mae(fit_lasso, test_lasso)))

## Elastic Net ----
learner = lrn("regr.glmnet", 
              alpha = to_tune(0, 1),
              s = to_tune(0.001, 10))
learner$param_set
para_elastic = c('alpha', 's')
set.seed(123)
instance = tune(tuner = tuner, task = task, learner = learner, 
                resampling = rsmp("cv", folds = 5),  
                measure = msr("regr.mae"), term_evals = 25)
instance$result
instance$archive
learner$param_set$values = instance$result_learner_param_vals
benchmarking_best_elastic = learner
benchmarking_best_elastic$train(task, row_ids = split$train)

# predict
elastic_pred_train = benchmarking_best_elastic$predict(task, row_ids = split$train)
elastic_pred_test = benchmarking_best_elastic$predict(task_all, row_ids = split$test)
train_elastic = tibble(row_ids = elastic_pred_train$row_ids, 
                       calendar_age = elastic_pred_train$truth,
                       pred_calendar_age = elastic_pred_train$response)
test_elastic = tibble(row_ids = elastic_pred_test$row_ids, 
                      calendar_age = elastic_pred_test$truth,
                      pred_calendar_age = elastic_pred_test$response)
fit_elastic = lm(pred_calendar_age ~ calendar_age, data = test_elastic)

optimal_elastic = instance$result %>% 
  select(all_of(para_elastic)) %>% 
  mutate(selection = 'optimal')

test_prediction = test_elastic %>% 
  left_join(dtm, by = 'row_ids')

out_elastic = list(train_tuning = instance$archive$data %>% 
                     select(all_of(para_elastic), regr.mae) %>% 
                     left_join(optimal_elastic, by = para_elastic) %>% 
                     set_names(c('α', '𝜆'), 'MAE', 'selection'),
                   train_prediction = train_elastic %>% left_join(dtm, by = 'row_ids'),
                   test_prediction = test_prediction,
                   Extrapolation = tibble(Min = min(test_prediction$pred_calendar_age),
                                          Max = max(test_prediction$pred_calendar_age),
                                          `n (<Min)` = sum(test_prediction$pred_calendar_age < min_age),
                                          `n (>Max)` = sum(test_prediction$pred_calendar_age > max_age)),
                   test_evaluate = tibble(r = cor(test_elastic$calendar_age, test_elastic$pred_calendar_age, method = 'pearson'),
                                          R2 = rsquare(fit_elastic, test_elastic),
                                          RMSE = rmse(fit_elastic, test_elastic),
                                          MAE = mae(fit_elastic, test_elastic)))

## KNN ----
learner = lrn("regr.kknn",
              k = to_tune(p_int(3, 30)), 
              distance = to_tune(p_dbl(1, 5)))
learner$param_set
para_knn = c('k', 'distance')
set.seed(123)
instance = tune(tuner = tuner, task = task, learner = learner, 
                resampling = rsmp("cv", folds = 5),  
                measure = msr("regr.mae"), term_evals = 25)
instance$result
instance$archive
learner$param_set$values = instance$result_learner_param_vals
benchmarking_best_knn = learner
benchmarking_best_knn$train(task, row_ids = split$train)

# predict
knn_pred_train = benchmarking_best_knn$predict(task, row_ids = split$train)
knn_pred_test = benchmarking_best_knn$predict(task_all, row_ids = split$test)
train_knn = tibble(row_ids = knn_pred_train$row_ids, 
                   calendar_age = knn_pred_train$truth,
                   pred_calendar_age = knn_pred_train$response)
test_knn = tibble(row_ids = knn_pred_test$row_ids, 
                  calendar_age = knn_pred_test$truth,
                  pred_calendar_age = knn_pred_test$response)
fit_knn = lm(pred_calendar_age ~ calendar_age, data = test_knn)

optimal_knn = instance$result %>% 
  select(all_of(para_knn)) %>% 
  mutate(selection = 'optimal')

test_prediction = test_knn %>% 
  left_join(dtm, by = 'row_ids')

out_knn = list(train_tuning = instance$archive$data %>% 
                 select(all_of(para_knn), regr.mae) %>% 
                 left_join(optimal_knn, by = para_knn) %>% 
                 set_names(c('Number of neighbors', 'Minkowski distance'), 
                           'MAE', 'selection'),
               train_prediction = train_knn %>% left_join(dtm, by = 'row_ids'),
               test_prediction = test_prediction,
               Extrapolation = tibble(Min = min(test_prediction$pred_calendar_age),
                                      Max = max(test_prediction$pred_calendar_age),
                                      `n (<Min)` = sum(test_prediction$pred_calendar_age < min_age),
                                      `n (>Max)` = sum(test_prediction$pred_calendar_age > max_age)),
               test_evaluate = tibble(r = cor(test_knn$calendar_age, test_knn$pred_calendar_age, method = 'pearson'),
                                      R2 = rsquare(fit_knn, test_knn),
                                      RMSE = rmse(fit_knn, test_knn),
                                      MAE = mae(fit_knn, test_knn)))

## SVM ----
learner = lrn("regr.svm",
              type = "eps-regression",
              kernel = "radial",
              cost = to_tune(0.1, 10),
              gamma = to_tune(0, 5))
learner$param_set
para_svm = c('cost', 'gamma')
set.seed(123)
instance = tune(tuner = tuner, task = task, learner = learner, 
                resampling = rsmp("cv", folds = 5),  
                measure = msr("regr.mae"), term_evals = 25)
instance$result
instance$archive
learner$param_set$values = instance$result_learner_param_vals
benchmarking_best_svm = learner
benchmarking_best_svm$train(task, row_ids = split$train)

# predict
svm_pred_train = benchmarking_best_svm$predict(task, row_ids = split$train)
svm_pred_test = benchmarking_best_svm$predict(task_all, row_ids = split$test)
train_svm = tibble(row_ids = svm_pred_train$row_ids, 
                   calendar_age = svm_pred_train$truth,
                   pred_calendar_age = svm_pred_train$response)
test_svm = tibble(row_ids = svm_pred_test$row_ids, 
                  calendar_age = svm_pred_test$truth,
                  pred_calendar_age = svm_pred_test$response)
fit_svm = lm(pred_calendar_age ~ calendar_age, data = test_svm)

optimal_svm = instance$result %>% 
  select(all_of(para_svm)) %>% 
  mutate(selection = 'optimal')

test_prediction = test_svm %>% 
  left_join(dtm, by = 'row_ids')

out_svm = list(train_tuning = instance$archive$data %>% 
                 select(all_of(para_svm), regr.mae) %>% 
                 left_join(optimal_svm, by = para_svm) %>% 
                 set_names(c('Cost of constraints violation', 
                             'γ'), 
                           'MAE', 'selection'),
               train_prediction = train_svm %>% left_join(dtm, by = 'row_ids'),
               test_prediction = test_prediction,
               Extrapolation = tibble(Min = min(test_prediction$pred_calendar_age),
                                      Max = max(test_prediction$pred_calendar_age),
                                      `n (<Min)` = sum(test_prediction$pred_calendar_age < min_age),
                                      `n (>Max)` = sum(test_prediction$pred_calendar_age > max_age)),
               test_evaluate = tibble(r = cor(test_svm$calendar_age, test_svm$pred_calendar_age, method = 'pearson'),
                                      R2 = rsquare(fit_svm, test_svm),
                                      RMSE = rmse(fit_svm, test_svm),
                                      MAE = mae(fit_svm, test_svm)))

## Random Forest ----
learner = lrn("regr.randomForest", 
              ntree = to_tune(p_int(100, 300)),
              mtry = to_tune(p_int(6, 20)),
              nodesize = to_tune(p_int(2, 5)),
              maxnodes = to_tune(p_int(3, 7)))
learner$param_set
para_random = c('ntree', 'mtry', 'nodesize', 'maxnodes')
set.seed(123)
instance = tune(tuner = tuner, task = task, learner = learner, 
                resampling = rsmp("cv", folds = 5),  
                measure = msr("regr.mae"), term_evals = 25)
instance$result
instance$archive
learner$param_set$values = instance$result_learner_param_vals
benchmarking_best_random = learner
benchmarking_best_random$train(task, row_ids = split$train)

# predict
random_pred_train = benchmarking_best_random$predict(task, row_ids = split$train)
random_pred_test = benchmarking_best_random$predict(task_all, row_ids = split$test)
train_random = tibble(row_ids = random_pred_train$row_ids, 
                      calendar_age = random_pred_train$truth,
                      pred_calendar_age = random_pred_train$response)
test_random = tibble(row_ids = random_pred_test$row_ids, 
                     calendar_age = random_pred_test$truth,
                     pred_calendar_age = random_pred_test$response)
fit_random = lm(pred_calendar_age ~ calendar_age, data = test_random)

optimal_random = instance$result %>% 
  select(all_of(para_random)) %>% 
  mutate(selection = 'optimal')

test_prediction = test_random %>% 
  left_join(dtm, by = 'row_ids')

out_random = list(train_tuning = instance$archive$data %>% 
                    select(all_of(para_random), regr.mae) %>% 
                    left_join(optimal_random, by = para_random) %>% 
                    set_names(c('Number of trees', 
                                'Number of variables randomly sampled as candidates',
                                'Minimum size of terminal nodes',
                                'Maximum number of terminal nodes trees'), 
                              'MAE', 'selection'),
                  train_prediction = train_random %>% left_join(dtm, by = 'row_ids'),
                  test_prediction = test_prediction,
                  Extrapolation = tibble(Min = min(test_prediction$pred_calendar_age),
                                         Max = max(test_prediction$pred_calendar_age),
                                         `n (<Min)` = sum(test_prediction$pred_calendar_age < min_age),
                                         `n (>Max)` = sum(test_prediction$pred_calendar_age > max_age)),
                  test_evaluate = tibble(r = cor(test_random$calendar_age, test_random$pred_calendar_age, method = 'pearson'),
                                         R2 = rsquare(fit_random, test_random),
                                         RMSE = rmse(fit_random, test_random),
                                         MAE = mae(fit_random, test_random)))

## LightGBM ----
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
              bagging_freq = to_tune(1, 10))

learner$param_set
para_lightgbm = c('num_leaves', 'max_depth', 'min_data_in_leaf',
                  'learning_rate', 'num_iterations', 
                  'lambda_l1', 'lambda_l2',
                  'feature_fraction', 'bagging_fraction', 'bagging_freq')
set.seed(123)
instance = tune(tuner = tuner, task = task, learner = learner, 
                resampling = rsmp("cv", folds = 5),  
                measure = msr("regr.mae"), term_evals = 25)
instance$result
instance$archive
learner$param_set$values = instance$result_learner_param_vals
benchmarking_best_lightgbm = learner
benchmarking_best_lightgbm$train(task, row_ids = split$train)

# predict
lightgbm_pred_train = benchmarking_best_lightgbm$predict(task, row_ids = split$train)
lightgbm_pred_test = benchmarking_best_lightgbm$predict(task_all, row_ids = split$test)
train_lightgbm = tibble(row_ids = lightgbm_pred_train$row_ids, 
                        calendar_age = lightgbm_pred_train$truth,
                        pred_calendar_age = lightgbm_pred_train$response)
test_lightgbm = tibble(row_ids = lightgbm_pred_test$row_ids, 
                       calendar_age = lightgbm_pred_test$truth,
                       pred_calendar_age = lightgbm_pred_test$response)
fit_lightgbm = lm(pred_calendar_age ~ calendar_age, data = test_lightgbm)

optimal_lightgbm = instance$result %>% 
  select(all_of(para_lightgbm)) %>% 
  mutate(selection = 'optimal')

test_prediction = test_lightgbm %>% 
  left_join(dtm, by = 'row_ids')

out_lightgbm = list(train_tuning = instance$archive$data %>% 
                      select(all_of(para_lightgbm), regr.mae) %>% 
                      left_join(optimal_lightgbm, by = para_lightgbm) %>%
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
                    train_prediction = train_lightgbm %>% left_join(dtm, by = 'row_ids'),
                    test_prediction = test_prediction,
                    Extrapolation = tibble(Min = min(test_prediction$pred_calendar_age),
                                           Max = max(test_prediction$pred_calendar_age),
                                           `n (<Min)` = sum(test_prediction$pred_calendar_age < min_age),
                                           `n (>Max)` = sum(test_prediction$pred_calendar_age > max_age)),
                    test_evaluate = tibble(r = cor(test_lightgbm$calendar_age, test_lightgbm$pred_calendar_age, method = 'pearson'),
                                           R2 = rsquare(fit_lightgbm, test_lightgbm),
                                           RMSE = rmse(fit_lightgbm, test_lightgbm),
                                           MAE = mae(fit_lightgbm, test_lightgbm)))

# feature selection -------------------------------------------------------
benchmarking_best_lightgbm$param_set$values

# 
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

# data 
X = dt %>% 
  slice(split$train) %>% 
  select(-calendar_age) %>% 
  as.data.frame()
y = dt %>% 
  slice(split$train) %>% 
  pull(calendar_age)

# BorutaShap
result = boruta_shap_lgbm(X, y, n_trials = 200, percentile = 100, p_value = 0.05, 
                          Model = benchmarking_best_lightgbm)
fs = list(feature = tibble(CpG = result$tentative), 'feature selection process' = result$out)

# Feature selection followed by hyperparameter tuning -----------------------------------------------------------------
fs = fs$feature

task$select(cols = fs$CpG)
task$row_roles$use

## LightGBM ----
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
              bagging_freq = to_tune(1, 10))

learner$param_set
para_lightgbm = c('num_leaves', 'max_depth', 'min_data_in_leaf',
                  'learning_rate', 'num_iterations', 
                  'lambda_l1', 'lambda_l2',
                  'feature_fraction', 'bagging_fraction', 'bagging_freq')
set.seed(1234)
instance = tune(tuner = tuner, task = task, learner = learner, 
                resampling = rsmp("cv", folds = 5),  
                measure = msr("regr.mae"), term_evals = 25)
instance$result
instance$archive
learner$param_set$values = instance$result_learner_param_vals
fs_best_lightgbm = learner
fs_best_lightgbm$train(task, row_ids = split$train)

# predict
lightgbm_pred_train = fs_best_lightgbm$predict(task, row_ids = split$train)
lightgbm_pred_test = fs_best_lightgbm$predict(task_all, row_ids = split$test)
train_lightgbm = tibble(row_ids = lightgbm_pred_train$row_ids, 
                        calendar_age = lightgbm_pred_train$truth,
                        pred_calendar_age = lightgbm_pred_train$response)
test_lightgbm = tibble(row_ids = lightgbm_pred_test$row_ids, 
                       calendar_age = lightgbm_pred_test$truth,
                       pred_calendar_age = lightgbm_pred_test$response)
fit_lightgbm = lm(pred_calendar_age ~ calendar_age, data = test_lightgbm)

optimal_lightgbm = instance$result %>% 
  select(all_of(para_lightgbm)) %>% 
  mutate(selection = 'optimal')

test_prediction = test_lightgbm %>% 
  left_join(dtm, by = 'row_ids')

out_lightgbm = list(train_tuning = instance$archive$data %>% 
                      select(all_of(para_lightgbm), regr.mae) %>% 
                      left_join(optimal_lightgbm, by = para_lightgbm) %>%
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
                    train_prediction = train_lightgbm %>% left_join(dtm, by = 'row_ids'),
                    test_prediction = test_prediction,
                    Extrapolation = tibble(Min = min(test_prediction$pred_calendar_age),
                                           Max = max(test_prediction$pred_calendar_age),
                                           `n (<Min)` = sum(test_prediction$pred_calendar_age < min_age),
                                           `n (>Max)` = sum(test_prediction$pred_calendar_age > max_age)),
                    test_evaluate = tibble(r = cor(test_lightgbm$calendar_age, test_lightgbm$pred_calendar_age, method = 'pearson'),
                                           R2 = rsquare(fit_lightgbm, test_lightgbm),
                                           RMSE = rmse(fit_lightgbm, test_lightgbm),
                                           MAE = mae(fit_lightgbm, test_lightgbm)))

## SHAP ----
library(kernelshap)
library(shapviz)

task = as_task_regr(dt, target = "calendar_age", id = 'clock')
task

task$select(cols = fs$CpG)
task$row_roles$use

X = dt %>% 
  select(all_of(fs$CpG)) %>% 
  as.data.frame()
y = dt %>% 
  pull(calendar_age)

# 
dtrain = lgb.Dataset(
  data = as.matrix(X),
  label = y
)

# 
params = list(
  num_leaves = fs_best_lightgbm$param_set$values$num_leaves,
  max_depth = fs_best_lightgbm$param_set$values$max_depth,
  min_data_in_leaf = fs_best_lightgbm$param_set$values$min_data_in_leaf,
  
  learning_rate = fs_best_lightgbm$param_set$values$learning_rate,
  
  
  lambda_l1 = fs_best_lightgbm$param_set$values$lambda_l1,
  lambda_l2 = fs_best_lightgbm$param_set$values$lambda_l2,
  
  feature_fraction = fs_best_lightgbm$param_set$values$feature_fraction,
  bagging_fraction = fs_best_lightgbm$param_set$values$bagging_fraction,
  bagging_freq = fs_best_lightgbm$param_set$values$bagging_freq,
  
  objective = "regression",
  boosting = "gbdt"
)

model = lgb.train(
  params = params,
  data = dtrain,
  nrounds = fs_best_lightgbm$param_set$values$num_iterations,
  verbose = -1
)
rm(dtrain)

sv_lgb = shapviz(model, X_pred = as.matrix(X))
p_bar = sv_importance(sv_lgb, max_display = 26, kind = "bar") + 
  theme_classic()
p_beeswarm = sv_importance(sv_lgb, max_display = 26, kind = "beeswarm") + 
  theme_classic() + 
  theme(axis.text.y = element_blank())

# Overall prediction --------------------------------------------------------------------

dts = dt %>% 
  select(calendar_age, all_of(fs$CpG))

# sample --
task = as_task_regr(dts, target = "calendar_age", id = 'clock')

# learner --
learner = lrn("regr.lightgbm",
              
              num_leaves = fs_best_lightgbm$param_set$values$num_leaves,
              max_depth = fs_best_lightgbm$param_set$values$max_depth,
              min_data_in_leaf = fs_best_lightgbm$param_set$values$min_data_in_leaf,
              
              learning_rate = fs_best_lightgbm$param_set$values$learning_rate,
              num_iterations = fs_best_lightgbm$param_set$values$num_iterations,
              
              lambda_l1 = fs_best_lightgbm$param_set$values$lambda_l1,
              lambda_l2 = fs_best_lightgbm$param_set$values$lambda_l2,
              
              feature_fraction = fs_best_lightgbm$param_set$values$feature_fraction,
              bagging_fraction = fs_best_lightgbm$param_set$values$bagging_fraction,
              bagging_freq = fs_best_lightgbm$param_set$values$bagging_freq
)


# prediction --
set.seed(123)
cv10 = rsmp("cv", folds = 5)
cv10$instantiate(task)

rr = mlr3::resample(task = task, learner = learner, resampling = cv10, store_models = TRUE)

pred = tibble(row_ids = rr$prediction()$row_ids, 
              calendar_age = rr$prediction()$truth,
              pred_calendar_age = rr$prediction()$response)
pred = pred %>% 
  left_join(dtm, by = 'row_ids') %>% 
  left_join(select(data, Sample_Name,
                   CD8T, CD4T, NK, Bcell, Mono, Neu), by = 'Sample_Name')

fit = lm(pred_calendar_age ~ calendar_age, data = pred)

out_evaluate = tibble(r = cor(pred$calendar_age, pred$pred_calendar_age, method = 'pearson'),
                      R2 = rsquare(fit, pred),
                      RMSE = rmse(fit, pred),
                      MAE = mae(fit, pred))

