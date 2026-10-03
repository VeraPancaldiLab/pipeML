# End-to-end tests of training, prediction and SHAP. They train several models, so they are skipped on CRAN;
# survival training (6 models with hyperparameter grids) takes several minutes and only runs locally.

test_that("classification: standard training, prediction and SHAP of the final model", {
  skip_on_cran()
  local_temp_wd()
  d <- sim_classification()
  train <- 1:60; test <- 61:80

  res <- suppressWarnings(compute_features.training.ML(features_train = d$X[train, ], target_var = d$y[train],
                                                       task_type = "classification", trait.positive = "R",
                                                       metric = "AUROC", k_folds = 2, n_rep = 1))
  expect_s3_class(res$Model, "train")
  expect_equal(nrow(res$Model$trainingData), length(train))
  expect_length(list.files("Results", pattern = "^fold_model_", recursive = TRUE), 0)

  # treebag (no hyperparameters) keeps the AUROC of each resample
  expect_gt(length(unique(res$ML_Models$BAG$resample$AUROC)), 1)

  # every final model is fitted with the hyperparameters selected by `metric` (reported in $bestTune)
  for (m in res$ML_Models) {
    if (isS4(m$finalModel)) next # kernlab models do not store the hyperparameters they were fitted with
    fitted_tune <- vapply(as.data.frame(m$finalModel$tuneValue), as.character, character(1))
    expect_equal(vapply(as.data.frame(m$bestTune)[names(fitted_tune)], as.character, character(1)), fitted_tune)
  }

  # $pred: each row keeps the probability of its own sample (the AUROC of a resample is reproduced from its rows)
  p <- res$ML_Models$KNN$pred
  p$obs <- res$ML_Models$KNN$trainingData$.outcome[p$rowIndex]
  for (g in split(p, list(p$Resample, p$k), drop = TRUE)) {
    expect_equal(pipeML:::calculate_auc_roc_resample(g$obs, g$yes), g$AUROC[1])
  }

  pred <- compute_prediction(model = res$Model, test_data = d$X[test, ], target_var = d$y[test],
                             task_type = "classification", trait.positive = "R")
  expect_true(pred$AUC$AUROC$estimate >= 0 && pred$AUC$AUROC$estimate <= 1)

  shap <- compute_shap_values(res$Model, task_type = "classification")
  expect_equal(dim(shap), c(length(train), 3))
  expect_equal(rownames(shap), rownames(d$X)[train])
  expect_false(anyNA(shap))

  # baseline + SHAP values = predicted probability of each training sample
  X_model <- res$Model$trainingData[, colnames(shap)]
  expect_equal(unname(attr(shap, "baseline") + rowSums(shap)),
               unname(predict(res$Model, X_model, type = "prob")[, "yes"]), tolerance = 1e-6)

  # reproducible with the same seed
  expect_equal(compute_shap_values(res$Model, task_type = "classification"), shap)

  # wrong object
  expect_error(compute_shap_values(res, task_type = "classification"), "caret train object")
})

test_that("classification standard path: redundant features are removed before the CV unless preprocess = FALSE", {
  skip_on_cran()
  local_temp_wd()
  d <- sim_classification()
  d$X$g1_copy <- d$X$g1 * 2 # correlation 1 with g1
  train <- 1:60; test <- 61:80

  res <- suppressWarnings(compute_features.training.ML(features_train = d$X[train, ], target_var = d$y[train],
                                                       task_type = "classification", trait.positive = "R",
                                                       metric = "AUROC", k_folds = 2, n_rep = 1))
  # every final model is trained on the preprocessed features
  for (m in res$ML_Models) {
    kept <- setdiff(colnames(m$trainingData), ".outcome")
    expect_length(kept, 3)
    expect_false(all(c("g1", "g1_copy") %in% kept))
  }
  expect_true(all(res$ML_Models$RF$results$mtry <= 3))

  # the test set keeps all its features: compute_prediction() uses the ones of the model
  pred <- compute_prediction(model = res$Model, test_data = d$X[test, ], target_var = d$y[test],
                             task_type = "classification", trait.positive = "R")
  expect_true(pred$AUC$AUROC$estimate >= 0 && pred$AUC$AUROC$estimate <= 1)

  res_raw <- suppressWarnings(compute_features.training.ML(features_train = d$X[train, ], target_var = d$y[train],
                                                           task_type = "classification", trait.positive = "R",
                                                           metric = "AUROC", k_folds = 2, n_rep = 1,
                                                           preprocess = FALSE))
  expect_setequal(setdiff(colnames(res_raw$Model$trainingData), ".outcome"), colnames(d$X))
})

test_that("classification: custom fold construction function without and with tunable parameters", {
  skip_on_cran()
  local_temp_wd()
  d <- sim_classification()
  # f_dup is redundant (correlation 1 with f_sum): removed inside each fold and in the final training set
  build <- function(x, scale) data.frame(f_sum = (x$g1 + x$g2) * scale, f_g3 = x$g3 * scale,
                                         f_dup = (x$g1 + x$g2) * scale * 2, row.names = rownames(x))

  fold_fun <- function(data, folds = NULL, bestune = NULL, scale = 1, ...) {
    tunable <- length(scale) > 1 || !is.null(bestune$scale)
    if (!is.null(bestune)) {
      s <- if (!is.null(bestune$scale)) bestune$scale[1] else scale
      return(list(cbind(build(data, s), target = data$target), list(scale_used = s),
                  if (tunable) data.frame(scale = s) else bestune))
    }
    for (nm in names(folds)) {
      tr <- folds[[nm]]; te <- setdiff(seq_len(nrow(data)), tr)
      one <- function(s) {
        out <- list(train_data = cbind(build(data[tr, ], s), target = data$target[tr]),
                    test_data = build(data[te, ], s), rowIndex = te, fold_name = nm, obs_test = data$target[te])
        if (tunable) out$params <- data.frame(scale = s)
        out
      }
      saveRDS(if (tunable) lapply(scale, one) else one(scale), file.path("Results", paste0("fold_", nm, ".rds")))
    }
  }

  for (args in list(list(fold_construction_args_fixed = list(scale = 2)),
                    list(fold_construction_args_tunable = list(scale = c(1, 2))),
                    list(fold_construction_args_fixed = list(scale = 2), preprocess = FALSE),
                    list(fold_construction_args_tunable = list(scale = c(1, 2)), preprocess = FALSE))) {
    n_features <- if (isFALSE(args$preprocess)) 3 else 2
    # a fold file left by an interrupted run must not be read as an extra resample
    saveRDS(list(), file.path("Results", "fold_Fold9.Rep9.rds"))

    res <- suppressWarnings(do.call(compute_features.training.ML,
                                    c(list(features_train = d$X, target_var = d$y, task_type = "classification",
                                           trait.positive = "R", metric = "AUROC", k_folds = 2, n_rep = 1,
                                           fold_construction_fun = fold_fun), args)))
    # final model trained on the custom features (preprocessed unless preprocess = FALSE); fold files removed after use
    kept <- setdiff(colnames(res$Model$trainingData), ".outcome")
    expect_length(kept, n_features)
    expect_true(all(kept %in% c("f_sum", "f_g3", "f_dup")) && "f_g3" %in% kept)
    expect_length(list.files("Results", pattern = "^fold_.*\\.rds$"), 0)
    expect_length(unique(res$Model$pred$Resample), 2)

    # the rf grid is sized from the features the models are trained on (custom features after preprocessing)
    expect_true(all(res$ML_Models$RF$results$mtry <= n_features))
    # the final model of every algorithm is trained on those features
    expect_true(all(vapply(res$ML_Models, function(m) length(setdiff(colnames(m$trainingData), ".outcome")), 1) == n_features))

    # every final model is fitted with the hyperparameters selected by `metric` (the AUROC-best row of $results)
    for (m in res$ML_Models) {
      if (isS4(m$finalModel)) next # kernlab models do not store the hyperparameters they were fitted with
      fitted_tune <- vapply(as.data.frame(m$finalModel$tuneValue), as.character, character(1))
      cols <- intersect(names(fitted_tune), names(m$results))
      if (length(cols) == 0) next
      best <- as.data.frame(m$results)[which.max(m$results$AUROC), cols, drop = FALSE]
      expect_equal(vapply(best, as.character, character(1)), fitted_tune[cols])
    }

    shap <- compute_shap_values(res$Model, task_type = "classification")
    expect_setequal(colnames(shap), kept)
    expect_equal(nrow(shap), nrow(d$X))
  }
})

test_that("survival: standard training (parallel), prediction and SHAP of the final model", {
  skip_on_cran()
  skip_on_ci()
  skip_if_not_installed("censored")
  local_temp_wd()
  d <- sim_survival(n = 60)
  d$X$g3_copy <- d$X$g3 * 2 # correlation 1 with g3: one of the two is removed before the CV

  res <- suppressWarnings(compute_features.training.ML(features_train = d$X, task_type = "survival",
                                                       time_var = d$time, event_var = d$event,
                                                       k_folds = 2, n_rep = 1, ncores = 2))
  expect_false(is.null(res$Model$Model_object))
  expect_length(colnames(res$Model$trainingData), 5)
  expect_true(all(c("g1", "g2", "time", "event") %in% colnames(res$Model$trainingData)))
  expect_false(all(c("g3", "g3_copy") %in% colnames(res$Model$trainingData)))

  pred <- compute_prediction(model = res$Model, test_data = d$X, task_type = "survival",
                             time_var = d$time, event_var = d$event)
  expect_true(pred$c_index >= 0 && pred$c_index <= 1)

  # the final model does not depend on ncores (the per-fold seeds run in the session only when sequential)
  res_seq <- suppressWarnings(compute_features.training.ML(features_train = d$X, task_type = "survival",
                                                           time_var = d$time, event_var = d$event,
                                                           k_folds = 2, n_rep = 1, ncores = NULL))
  risk_of <- function(r) as.numeric(unlist(pipeML:::predict_and_evaluate_survival(r$Model$Model_object, d$X)$preds))
  expect_equal(risk_of(res_seq), risk_of(res))

  shap <- compute_shap_values(res$Model, task_type = "survival")
  expect_equal(dim(shap), c(nrow(d$X), 3))
  risk <- as.numeric(unlist(pipeML:::predict_and_evaluate_survival(res$Model$Model_object, d$X, NULL, NULL)$preds))
  expect_equal(unname(attr(shap, "baseline") + rowSums(shap)), risk, tolerance = 1e-6)
})

test_that("survival LODO: folds use the cohort, and the cohort label is not a predictor", {
  skip_on_cran()
  skip_on_ci()
  skip_if_not_installed("censored")
  local_temp_wd()
  d <- sim_survival(n = 90)
  cohort <- rep(c("CohortA", "CohortB", "CohortC"), each = 30)

  res <- suppressWarnings(compute_features.training.ML(features_train = d$X, task_type = "survival",
                                                       time_var = d$time, event_var = d$event,
                                                       LODO = TRUE, batch_var = cohort,
                                                       k_folds = 3, n_rep = 1, ncores = 2))

  predictors <- colnames(workflows::extract_mold(res$Model$Model_object)$predictors)
  expect_false(any(c("dataset", "strata") %in% predictors))
  expect_setequal(colnames(res$Model$trainingData), c("g1", "g2", "g3", "time", "event"))

  # every held-out fold contains samples of every cohort
  rm_ <- res$Model$Resample_matrix
  keep <- !duplicated(paste(rm_$Resample, rm_$rowIndex))
  counts <- table(rm_$Resample[keep], cohort[rm_$rowIndex[keep]])
  expect_true(all(counts > 0))
})

test_that("survival one-step (compute_features.ML) uses the custom fold construction function", {
  skip_on_cran()
  skip_on_ci()
  skip_if_not_installed("censored")
  local_temp_wd()
  d <- sim_survival(n = 80)
  train <- 1:60; test <- 61:80
  coldata <- data.frame(time = d$time, event = d$event, row.names = rownames(d$X))
  build <- function(x) data.frame(f_sum = x$g1 + x$g2, f_g3 = x$g3, row.names = rownames(x))

  # the fold function receives the outcome in data ("time" and "event" columns), as "target" for classification
  seen_outcome <- new.env()
  fold_fun <- function(data, folds = NULL, bestune = NULL, ...) {
    seen_outcome$ok <- all(c("time", "event") %in% colnames(data))
    time <- data$time; event <- data$event
    data <- data[, setdiff(colnames(data), c("time", "event")), drop = FALSE]
    add_outcome <- function(f, idx) { f$time <- time[idx]; f$event <- event[idx]; f }
    if (!is.null(bestune)) return(list(add_outcome(build(data), seq_len(nrow(data))), list(), bestune))
    for (nm in names(folds)) {
      tr <- folds[[nm]]; te <- setdiff(seq_len(nrow(data)), tr)
      saveRDS(list(train_data = add_outcome(build(data[tr, ]), tr), test_data = add_outcome(build(data[te, ]), te),
                   rowIndex = te), file.path("Results", paste0("fold_", nm, ".rds")))
    }
  }

  res <- suppressWarnings(compute_features.ML(features_train = d$X[train, ], features_test = build(d$X[test, ]),
                                              coldata = coldata, task_type = "survival",
                                              time_var = "time", event_var = "event", k_folds = 2, n_rep = 1,
                                              ncores = 2, fold_construction_fun = fold_fun))

  expect_true(seen_outcome$ok)
  predictors <- colnames(workflows::extract_mold(res$Model$Model$Model_object)$predictors)
  expect_setequal(setdiff(predictors, "(Intercept)"), c("f_sum", "f_g3"))
  expect_true(is.numeric(res$C_index) && res$C_index >= 0 && res$C_index <= 1)

  # the forest's mtry grid is sized from the custom features (2), not from the input columns (3)
  forest <- res$Model$ML_Models$rand_forest_aorsf
  expect_false(is.null(forest))
  expect_true(all(forest$Results_folds$mtry <= 2))
})

test_that("survival custom folds with tunable arguments run with the default ncores (sequential)", {
  skip_on_cran()
  skip_on_ci()
  skip_if_not_installed("censored")
  local_temp_wd()
  d <- sim_survival(n = 60)
  build <- function(x, s) data.frame(f_sum = (x$g1 + x$g2) * s, f_g3 = x$g3, time = x$time, event = x$event,
                                     row.names = rownames(x))
  fold_fun <- function(data, folds = NULL, bestune = NULL, scale = 1, ...) {
    if (!is.null(bestune)) return(list(build(data, bestune$scale[1]), list(), data.frame(scale = bestune$scale[1])))
    for (nm in names(folds)) {
      tr <- folds[[nm]]; te <- setdiff(seq_len(nrow(data)), tr)
      out <- lapply(scale, function(s) list(train_data = build(data[tr, ], s), test_data = build(data[te, ], s),
                                             rowIndex = te, fold_name = nm, params = data.frame(scale = s)))
      saveRDS(out, file.path("Results", paste0("fold_", nm, ".rds")))
    }
  }

  # ncores = NULL used to create a cluster with 0 workers ("subscript out of bounds")
  res <- suppressWarnings(compute_features.training.ML(features_train = d$X, task_type = "survival",
                                                       time_var = d$time, event_var = d$event, k_folds = 2, n_rep = 1,
                                                       fold_construction_fun = fold_fun,
                                                       fold_construction_args_tunable = list(scale = c(1, 2))))
  expect_true(res$Custom_output$Parameters$scale %in% c(1, 2))
  expect_setequal(setdiff(colnames(res$Model$trainingData), c("time", "event")), c("f_sum", "f_g3"))
  expect_length(list.files("Results", pattern = "^fold_.*\\.rds$"), 0)
})
