test_that("only AUROC and AUPRC are accepted as classification metrics", {
  skip_on_cran()
  local_temp_wd()
  d <- sim_classification()
  expect_error(compute_features.training.ML(features_train = d$X, target_var = d$y, task_type = "classification",
                                            trait.positive = "R", metric = "Accuracy", k_folds = 2, n_rep = 1),
               "Choose either \"AUROC\" or \"AUPRC\"")
})

test_that("invalid inputs of compute_features.training.ML() stop with a clear message", {
  d <- sim_classification()
  s <- sim_survival()
  train <- function(...) compute_features.training.ML(features_train = d$X, k_folds = 2, n_rep = 1, ...)

  expect_error(train(target_var = d$y, trait.positive = "R"), "`task_type` must be provided")
  expect_error(train(task_type = "classification", target_var = d$y, trait.positive = "r"),
               "with and without `trait.positive`")
  expect_error(train(task_type = "classification", target_var = replace(d$y, 1, NA), trait.positive = "R"),
               "missing values")
  expect_error(train(task_type = "classification", target_var = d$y, trait.positive = "R", LODO = TRUE),
               "`batch_var`")
  expect_error(compute_features.training.ML(features_train = s$X, task_type = "survival", time_var = s$time,
                                            event_var = factor(ifelse(s$event == 1, "dead", "alive"))),
               "1 = event, 0 = censored")

  # an event indicator given as a factor keeps its 0/1 values
  testthat::local_mocked_bindings(compute_k_fold_CV_survival = function(df_features, df_outcome, ...) df_outcome)
  outcome <- compute_features.training.ML(features_train = s$X, task_type = "survival", time_var = s$time,
                                          event_var = factor(s$event))
  expect_equal(outcome$event, s$event)
})

test_that("invalid inputs of compute_features.ML() stop with a clear message", {
  d <- sim_classification()
  s <- sim_survival()
  coldata <- data.frame(resp = d$y, cohort = "A", row.names = rownames(d$X))
  coldata_surv <- data.frame(time = s$time, event = s$event, row.names = rownames(s$X))
  ml <- function(...) compute_features.ML(features_train = d$X[1:60, ], features_test = d$X[61:80, ], coldata = coldata,
                                          k_folds = 2, n_rep = 1, ...)
  ml_surv <- function(coldata, ...) compute_features.ML(features_train = s$X[1:40, ], features_test = s$X[41:60, ],
                                                        coldata = coldata, task_type = "survival", ...)

  expect_error(ml(trait = "resp", trait.positive = "R"), "`task_type` must be provided")
  expect_error(ml(task_type = "regression", trait = "resp", trait.positive = "R"), "`task_type` must be provided")
  expect_error(ml(task_type = "classification", trait.positive = "R"), "both `trait` and `trait.positive`")
  expect_error(ml(task_type = "classification", trait = "resp", trait.positive = "r"), "with and without `trait.positive`")
  expect_error(ml(task_type = "classification", trait = "resp", trait.positive = "R", LODO = TRUE), "`batch_id`")
  expect_error(ml_surv(coldata_surv), "both `time_var` and `event_var`")
  expect_error(ml_surv(coldata_surv[-(1:3), ], time_var = "time", event_var = "event"), "missing from rownames\\(coldata\\)")
  expect_error(ml_surv(transform(coldata_surv, event = ifelse(event == 1, "dead", "alive")), time_var = "time",
                       event_var = "event"), "1 = event, 0 = censored")

  # an event indicator given as a factor keeps its 0/1 values
  testthat::local_mocked_bindings(compute_k_fold_CV_survival = function(df_features, df_outcome, ...) stop("event: ", paste(df_outcome$event, collapse = "")))
  expect_error(ml_surv(transform(coldata_surv, event = factor(event)), time_var = "time", event_var = "event"),
               paste0("event: ", paste(s$event[1:40], collapse = "")))
})

test_that("custom fold function without tunable arguments and the default metric (AUROC) trains", {
  skip_on_cran()
  local_temp_wd()
  d <- sim_classification()
  build <- function(x) data.frame(f_sum = x$g1 + x$g2, f_g3 = x$g3, row.names = rownames(x))
  fold_fun <- function(data, folds = NULL, bestune = NULL, ...) {
    if (!is.null(bestune)) return(list(cbind(build(data), target = data$target), list(), bestune))
    for (nm in names(folds)) {
      tr <- folds[[nm]]; te <- setdiff(seq_len(nrow(data)), tr)
      saveRDS(list(train_data = cbind(build(data[tr, ]), target = data$target[tr]), test_data = build(data[te, ]),
                   rowIndex = te, fold_name = nm, obs_test = data$target[te]),
              file.path("Results", paste0("fold_", nm, ".rds")))
    }
  }

  # metric = NULL: AUROC is used
  res <- suppressWarnings(compute_features.training.ML(features_train = d$X, target_var = d$y,
                                                       task_type = "classification", trait.positive = "R",
                                                       k_folds = 2, n_rep = 1, fold_construction_fun = fold_fun))
  expect_s3_class(res$Model, "train")
  expect_setequal(setdiff(colnames(res$Model$trainingData), ".outcome"), c("f_sum", "f_g3"))
})

test_that("compute_cv_AUC() returns the median per model, sorted, and saves the plots", {
  local_temp_wd()
  fake <- function(auroc) list(resample = data.frame(AUROC = auroc, AUPRC = auroc - 0.1),
                               trainingData = data.frame(.outcome = factor(c("no", "yes", "yes"))))
  models <- list(A = fake(c(0.6, 0.7, 0.9)), B = fake(c(0.8, 0.9, 0.95)))

  res <- compute_cv_AUC(models, file_name = "test", AUC_type = "AUROC", return = TRUE)
  expect_named(res$AUROC, c("model", "Median_AUROC", "MAD_AUROC"))
  expect_named(res$AUPRC, c("model", "Median_AUPRC", "MAD_AUPRC"))
  expect_equal(res$AUROC$model, c("B", "A"))
  expect_equal(res$AUROC$Median_AUROC, c(0.9, 0.7))
  expect_equal(res$Top_model, "B")
  expect_true(all(file.exists(file.path("Results", c("AUROC_CV_methods_test.pdf", "AUPRC_CV_methods_test.pdf")))))
})

test_that("AUROC and AUPRC of a resample with a single class are NA, with a warning", {
  expect_warning(expect_true(is.na(calculate_auc_roc_resample(rep("yes", 5), runif(5)))), "single class")
  expect_warning(expect_true(is.na(calculate_auc_prc_resample(rep("no", 5), runif(5)))), "single class")
  expect_equal(calculate_auc_roc_resample(c("yes", "yes", "no", "no"), c(0.9, 0.8, 0.2, 0.1)), 1)
})

test_that("AUROC and AUPRC do not depend on the order of samples with the same predicted probability", {
  obs <- rep(c("yes", "no"), each = 10)
  constant <- rep(0.5, 20)
  # a model predicting the same probability for all samples: AUROC 0.5, AUPRC = proportion of positive samples
  expect_equal(calculate_auc_roc_resample(obs, constant), 0.5)
  expect_equal(calculate_auc_roc_resample(rev(obs), constant), 0.5)
  expect_equal(calculate_auc_prc_resample(obs, constant), 0.5)
  expect_equal(calculate_auc_prc_resample(rev(obs), constant), 0.5)

  set.seed(1)
  obs <- sample(c("yes", "no"), 60, replace = TRUE)
  prob <- round(runif(60) * 4) / 4 # 5 distinct probabilities
  yes_first <- order(obs == "no"); no_first <- order(obs == "yes")
  expect_equal(calculate_auc_roc_resample(obs[yes_first], prob[yes_first]), calculate_auc_roc_resample(obs[no_first], prob[no_first]))
  expect_equal(calculate_auc_prc_resample(obs[yes_first], prob[yes_first]), calculate_auc_prc_resample(obs[no_first], prob[no_first]))
  sens_spec <- function(i) get_sensitivity_specificity(data.frame(yes = prob[i]), obs[i], "model")
  expect_equal(calculate_auroc(sens_spec(yes_first)$fpr, sens_spec(yes_first)$Sensitivity),
               calculate_auroc(sens_spec(no_first)$fpr, sens_spec(no_first)$Sensitivity))

  # a model ranking every positive sample above every negative one: AUROC and AUPRC are 1
  obs <- rep(c("yes", "no"), c(5, 15)); prob <- seq(1, 0, length.out = 20)
  expect_equal(calculate_auc_roc_resample(obs, prob), 1)
  expect_equal(calculate_auc_prc_resample(obs, prob), 1)
})

test_that("LODO folds stop with a clear message when a cohort has fewer samples than k_folds", {
  d <- data.frame(x = rnorm(23), target = rep(c("no", "yes"), length.out = 23), dataset = rep(c("A", "B", "C"), c(10, 10, 3)))
  expect_error(construct_stratified_cohort_folds(d, "dataset", "target", k_folds = 5, n_rep = 1),
               "fewer samples than `k_folds` \\(5\\): C \\(3\\)")
  folds <- construct_stratified_cohort_folds(d, "dataset", "target", k_folds = 3, n_rep = 2)
  expect_length(folds, 6)
})

test_that("invalid inputs of compute_prediction() stop with a clear message", {
  d <- sim_classification()
  train_data <- cbind(d$X[1:60, ], target = factor(ifelse(d$y[1:60] == "R", "yes", "no"), levels = c("no", "yes")))
  fit <- function(data) suppressWarnings(caret::train(target ~ ., data = data, method = "glm",
                                                      trControl = caret::trainControl(method = "none", classProbs = TRUE)))
  model <- fit(train_data)
  X_test <- d$X[61:80, ]; y_test <- d$y[61:80]
  predict_test <- function(...) utils::capture.output(res <- compute_prediction(...)) # the function prints its progress

  expect_error(compute_prediction(model, X_test, y_test, trait.positive = "R", task_type = "regression"), "`task_type` must be")
  expect_error(compute_prediction(list(Model = model), X_test, y_test, trait.positive = "R"), "caret train object")
  expect_error(compute_prediction(model, X_test, task_type = "survival", time_var = 1:20, event_var = rep(1, 20)),
               "`\\$Model_object`")
  expect_error(compute_prediction(model, X_test, y_test, trait.positive = "r"), "with and without `trait.positive`")
  expect_error(compute_prediction(model, X_test[y_test == "R", ], y_test[y_test == "R"], trait.positive = "R"),
               "with and without `trait.positive`")
  expect_error(compute_prediction(model, X_test, replace(y_test, 1, NA), trait.positive = "R"), "missing values")

  # a model with a single feature can be used for prediction
  one_feature <- fit(train_data[, c("g1", "target")])
  utils::capture.output(pred <- compute_prediction(one_feature, X_test, y_test, trait.positive = "R"))
  expect_true(pred$AUC$AUROC$estimate >= 0 && pred$AUC$AUROC$estimate <= 1)
  expect_equal(nrow(pred$Predictions), 20)
})

test_that("get_curves() saves the ROC and precision-recall plots, for one curve and for several cohorts", {
  local_temp_wd()
  set.seed(1)
  metrics <- function(n) {
    obs <- factor(sample(c("no", "yes"), n, replace = TRUE), levels = c("no", "yes"))
    get_sensitivity_specificity(data.frame(yes = stats::runif(n)), obs, "model")
  }
  auc <- list(estimate = 0.7, lower = 0.6, upper = 0.8)
  m <- metrics(40)

  # one curve, also when the table is a tibble and without a file name
  suppressWarnings(get_curves(data = tibble::as_tibble(m), color = "model", auc_roc = auc, auc_prc = auc))
  expect_true(all(file.exists(file.path("Results", c("ROC_curve.pdf", "PRC_curve.pdf")))))

  # several cohorts need LODO = TRUE and vectors named with the cohorts
  stacked <- rbind(cbind(metrics(30), cohort = "A"), cbind(metrics(30), cohort = "B"))
  per_cohort <- lapply(auc, function(v) c(A = v, B = v))
  suppressWarnings(get_curves(data = stacked, color = "cohort", auc_roc = per_cohort, auc_prc = per_cohort, LODO = TRUE,
                              file.name = "cohorts"))
  expect_true(file.exists(file.path("Results", "ROC_curve_cohorts.pdf")))
  expect_error(get_curves(data = stacked, color = "cohort", auc_roc = per_cohort, auc_prc = per_cohort), "several curves")
  expect_error(get_curves(data = stacked, color = "cohort", auc_roc = lapply(per_cohort, unname), auc_prc = per_cohort,
                          LODO = TRUE), "named with the values")
})

test_that("bootstrap_auc() does not reset the random numbers of the session and is reproducible", {
  set.seed(1)
  pred <- data.frame(yes = stats::runif(40))
  target <- factor(sample(c("no", "yes"), 40, replace = TRUE), levels = c("no", "yes"))
  set.seed(5); expected <- stats::runif(1)
  set.seed(5); b1 <- pipeML:::bootstrap_auc(pred, target, "m", B = 50)
  expect_equal(stats::runif(1), expected)
  expect_equal(pipeML:::bootstrap_auc(pred, target, "m", B = 50), b1) # same seed, same intervals
})

test_that("roc_at_grid() and prc_at_grid() read the curves with straight lines between points, as they are drawn", {
  # ROC: a tied group of both classes gives a diagonal from (0, 0) to (0.5, 0.5); then a vertical step at 0.5
  roc <- pipeML:::roc_at_grid(fpr = c(0.5, 0.5, 1), sensitivity = c(0.5, 1, 1), grid = c(0, 0.25, 0.5, 0.75, 1))
  expect_equal(roc, c(0, 0.25, 1, 1, 1)) # a step function gave 0 at 0.25

  # PR: start (0, 1) -> (0.5, 1) -> vertical drop to (0.5, 0.5) -> (1, 0.75)
  prc <- pipeML:::prc_at_grid(recall = c(0.5, 0.5, 1), precision = c(1, 0.5, 0.75), grid = c(0, 0.25, 0.5, 0.75, 1))
  expect_equal(prc, c(1, 1, 1, 0.625, 0.75)) # 0.625: halfway between (0.5, 0.5) and (1, 0.75)

  # without ties the ROC curve is a staircase: same values as reading it as a step function
  set.seed(1)
  y <- factor(sample(c("no", "yes"), 40, replace = TRUE), levels = c("no", "yes"))
  ss <- pipeML:::get_sensitivity_specificity(data.frame(yes = stats::runif(40)), y, "m")
  grid <- seq(0, 1, length.out = 101)
  step <- c(0, cummax(ss$Sensitivity))[findInterval(grid, c(0, ss$fpr))]
  expect_equal(pipeML:::roc_at_grid(ss$fpr, ss$Sensitivity, grid), step)
})

test_that("bootstrap_auc() resamples each class separately, so a small test set with few positives works", {
  set.seed(2)
  target <- factor(c(rep("yes", 2), rep("no", 18)), levels = c("no", "yes")) # ~12% of plain resamples had no positive
  pred <- data.frame(yes = c(stats::runif(2, 0.5, 1), stats::runif(18, 0, 0.7)))
  b <- pipeML:::bootstrap_auc(pred, target, "m", B = 200)
  expect_false(anyNA(b$AUROC$values))
  expect_true(b$AUROC$lower <= b$AUROC$upper)
})
