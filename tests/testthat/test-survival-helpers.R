test_that("compute_ml_survival() returns the error of a failed fit instead of a silent NULL", {
  skip_if_not_installed("censored")
  d <- sim_survival()
  df <- cbind(d$X, event = d$event)   # no "time" column: the fit must fail

  res <- pipeML:::compute_ml_survival(df_train = df, outcome_col = "time", event_col = "event",
                                      model = "cox_ph_survival", models_hyperparameters = NULL)

  expect_s3_class(res, "pipeML_fit_error")
  expect_equal(res$model, "cox_ph_survival")
  expect_true(nzchar(res$error))
})

test_that("compute_ml_survival() fits a model without hyperparameters when given an empty bestTune", {
  skip_if_not_installed("censored")
  d <- sim_survival()
  df <- cbind(d$X, time = d$time, event = d$event)

  res <- pipeML:::compute_ml_survival(df_train = df, outcome_col = "time", event_col = "event",
                                      model = "cox_ph_survival",
                                      models_hyperparameters = list(tibble::tibble(.rows = 1)))

  expect_false(inherits(res, "pipeML_fit_error"))
  expect_s3_class(res, "workflow")
})

test_that("predict_and_evaluate_survival() returns risk scores without outcome columns", {
  skip_if_not_installed("censored")
  d <- sim_survival()
  df <- cbind(d$X, time = d$time, event = d$event)
  fit <- pipeML:::compute_ml_survival(df_train = df, outcome_col = "time", event_col = "event",
                                      model = "cox_ph_survival", models_hyperparameters = NULL)

  pred <- pipeML:::predict_and_evaluate_survival(fit, d$X, NULL, NULL)

  expect_equal(length(unlist(pred$preds)), nrow(d$X))
  expect_true(is.na(pred$c_index))
})

test_that("aggregate_results() excludes configurations with a failed fit and keeps failed models in place", {
  ok_rows <- function(resample, penalty) {
    tibble::tibble(predictions = c(0.1, 0.2), model = "m1", Resample = resample, rowIndex = 1:2,
                   penalty = penalty, c_index = c(0.7, 0.7))
  }
  failed_row <- function(model, resample, penalty) {
    tibble::tibble(model = model, Resample = resample, fit_error = "boom", penalty = penalty)
  }
  # 2 folds x 2 models. Model 1: penalty 0.5 fails in fold 2. Model 2: fails everywhere
  all_loaded <- list(
    list(dplyr::bind_rows(ok_rows("Fold1", 0.1), ok_rows("Fold1", 0.5)), failed_row("m2", "Fold1", 0.1)),
    list(dplyr::bind_rows(ok_rows("Fold2", 0.1), failed_row("m1", "Fold2", 0.5)), failed_row("m2", "Fold2", 0.1))
  )

  msgs <- testthat::capture_messages(res <- pipeML:::aggregate_results(all_loaded, task = "survival"))
  expect_length(msgs, 2)
  expect_match(msgs[1], "Model m1: 1 fit\\(s\\) failed; 1 hyperparameter configuration\\(s\\) excluded")
  expect_match(msgs[2], "Model m2: .*model excluded")

  expect_length(res, 2)
  expect_equal(unique(res[[1]]$Prediction_folds$penalty), 0.1)
  expect_equal(res[[1]]$bestTune$penalty, 0.1)
  expect_null(res[[2]])
})

test_that("aggregate_results() excludes configurations whose C-index is NA in a resample, with a message", {
  rows <- function(resample, penalty, c_index) {
    tibble::tibble(predictions = c(0.1, 0.2), model = "m1", Resample = resample, rowIndex = 1:2,
                   penalty = penalty, c_index = c_index)
  }
  # penalty 0.1 has the best C-index in Fold1 but an NA C-index in Fold2 (e.g. no event in that held-out fold)
  all_loaded <- list(
    list(dplyr::bind_rows(rows("Fold1", 0.1, 0.9), rows("Fold1", 0.5, 0.7))),
    list(dplyr::bind_rows(rows("Fold2", 0.1, NA_real_), rows("Fold2", 0.5, 0.6)))
  )
  expect_message(res <- pipeML:::aggregate_results(all_loaded, task = "survival"),
                 "Model m1: 1 fit\\(s\\) failed; 1 hyperparameter configuration\\(s\\) excluded. Error: C-index is NA")
  expect_equal(res[[1]]$bestTune$penalty, 0.5)
  expect_equal(unique(res[[1]]$Prediction_folds$penalty), 0.5)
  expect_false("fit_error" %in% names(res[[1]]$Prediction_folds))
})

test_that("survival predictions are risk scores: higher for samples with a risk factor, C-index above 0.5", {
  skip_if_not_installed("censored")
  set.seed(1); n <- 150
  X <- data.frame(g1 = rnorm(n), g2 = rnorm(n))
  df <- cbind(X, time = stats::rexp(n, exp(1.2 * X$g1)), event = 1) # g1 raises the hazard
  # a model predicting a linear predictor (Cox) and one predicting a survival time (tree)
  for (model in c("cox_ph_survival", "decision_tree_partykit")) {
    fit <- compute_ml_survival(df_train = df, outcome_col = "time", event_col = "event", model = model,
                               models_hyperparameters = list(NULL))
    pred <- predict_and_evaluate_survival(fit, df, "time", "event")
    expect_gt(cor(pred$preds$.pred, df$g1, method = "spearman"), 0.5)
    expect_gt(pred$c_index, 0.6)
  }
})

test_that("plot_survival_performance() forms the risk groups, also with few distinct predictions", {
  local_temp_wd()
  set.seed(1); n <- 120
  obs <- data.frame(time = stats::rexp(n), event = stats::rbinom(n, 1, 0.8))
  pred <- function(v) list(preds = tibble::tibble(.pred = v), c_index = 0.6, c_index_lower = 0.5, c_index_upper = 0.7)

  p <- suppressWarnings(suppressMessages(plot_survival_performance(obs, pred(rnorm(n)), n_groups = 3, file_name = "t")))
  expect_s3_class(p, "ggsurvplot")
  expect_true(file.exists(file.path("Results", "Survival_KM_t.pdf")))

  # two distinct predictions: the samples with the same prediction stay together, so only 2 groups
  expect_message(suppressWarnings(plot_survival_performance(obs, pred(rep(c(-2, -1), c(80, 40))), n_groups = 3)),
                 "Only 2 risk groups")
  expect_true(file.exists(file.path("Results", "Survival_KM.pdf")))

  expect_error(plot_survival_performance(obs, pred(rep(1, n)), n_groups = 2), "same risk for all samples")
  expect_error(plot_survival_performance(obs, pred(rnorm(n)), n_groups = 1), "at least 2")

  # default: 2 groups, as in compute_prediction()
  p2 <- suppressWarnings(plot_survival_performance(obs, pred(rnorm(n))))
  expect_equal(levels(p2$data.survplot$risk_group), c("Low risk", "High risk"))

  # samples without a prediction are left out with a message
  v <- rnorm(n); v[c(3, 9)] <- NA
  expect_message(suppressWarnings(plot_survival_performance(obs, pred(v))), "2 sample\\(s\\) without a prediction")

  # clear errors: more groups than samples, rows that do not match the predictions
  expect_error(plot_survival_performance(obs[1:5, ], pred(rnorm(5)), n_groups = 8), "must not be larger than the number of samples")
  expect_error(plot_survival_performance(obs[1:3, ], pred(rnorm(n))), "one row per prediction")
})

test_that("wrapper_train_best_hyperparams_survival() returns the model hyperparameters, and NULL when the fit fails", {
  skip_if_not_installed("censored")
  skip_if_not_installed("aorsf")
  set.seed(1); n <- 60
  train <- data.frame(g1 = rnorm(n), g2 = rnorm(n), time = stats::rexp(n), event = stats::rbinom(n, 1, 0.8))
  # final mode of a fold function with one tunable feature parameter (`scale`)
  fold_fun <- function(data, bestune = NULL, ...) {
    feats <- data.frame(f1 = data$g1 * bestune$scale, f2 = data$g2, time = data$time, event = data$event)
    list(feats, list(), data.frame(scale = bestune$scale))
  }
  optimized <- list(bestTune = data.frame(trees = 10, min_n = 5, mtry = 2, scale = 2))

  res <- pipeML:::wrapper_train_best_hyperparams_survival(train, optimized, "rand_forest_aorsf", fold_fun, NULL)
  expect_setequal(colnames(res$Model$bestTune), c("trees", "min_n", "mtry")) # feature parameter only in Parameters
  expect_equal(res$custom_output$Parameters$scale, 2)

  # a feature with only missing values: the fit fails, the model is excluded with a warning
  na_fun <- function(data, bestune = NULL, ...) {
    list(data.frame(f1 = NA_real_, f2 = data$g2, time = data$time, event = data$event), list(), data.frame(scale = 1))
  }
  expect_warning(out <- pipeML:::wrapper_train_best_hyperparams_survival(train, list(bestTune = data.frame(scale = 1)),
                                                                         "cox_ph_survival", na_fun, NULL,
                                                                         preprocess = FALSE),
                 "could not be fitted on all training samples")
  expect_null(out)
})

test_that("the number of bags of bag_tree_rpart (trees) reaches the engine as nbagg", {
  skip_if_not_installed("censored")
  set.seed(1)
  d <- data.frame(a = rnorm(80), b = rnorm(80), time = stats::rexp(80), event = stats::rbinom(80, 1, 0.8))
  fit <- pipeML:::compute_ml_survival(d, outcome_col = "time", event_col = "event", model = "bag_tree_rpart",
                                      models_hyperparameters = list(tibble::tibble(trees = 7)))
  expect_length(parsnip::extract_fit_engine(fit)$mtrees, 7)

  fold_fun <- function(data, bestune = NULL, ...) list(data[c("a", "b", "time", "event")], list(), data.frame(scale = 1))
  res <- pipeML:::wrapper_train_best_hyperparams_survival(d, list(bestTune = data.frame(trees = 12, scale = 1)),
                                                          "bag_tree_rpart", fold_fun, NULL)
  expect_length(parsnip::extract_fit_engine(res$Model$fitted)$mtrees, 12)
})

test_that("predict_and_evaluate_survival(ci = FALSE) gives the same C-index without bootstrap, and the bootstrap keeps the session RNG", {
  skip_if_not_installed("censored")
  set.seed(1)
  d <- data.frame(a = rnorm(40), b = rnorm(40), time = stats::rexp(40), event = stats::rbinom(40, 1, 0.8))
  fit <- pipeML:::compute_ml_survival(d, outcome_col = "time", event_col = "event", model = "cox_ph_survival",
                                      models_hyperparameters = NULL)
  with_ci <- pipeML:::predict_and_evaluate_survival(fit, d, "time", "event")
  no_ci <- pipeML:::predict_and_evaluate_survival(fit, d, "time", "event", ci = FALSE)
  expect_equal(no_ci$c_index, with_ci$c_index)
  expect_true(is.na(no_ci$c_index_lower) && is.na(no_ci$c_index_upper))
  expect_equal(no_ci$preds, with_ci$preds)

  # the bootstrap (seed 123) does not reset the random numbers of the session
  set.seed(999); expected <- stats::runif(1)
  set.seed(999); invisible(pipeML:::predict_and_evaluate_survival(fit, d, "time", "event"))
  expect_equal(stats::runif(1), expected)
})

test_that("compute_cindex_ci() resamples events and censored samples separately, so every resample has events", {
  set.seed(1)
  d <- data.frame(time = stats::rexp(20), event = c(1, 1, rep(0, 18)), .pred = stats::rnorm(20))
  ci <- pipeML:::compute_cindex_ci(d, n_boot = 200)
  expect_false(anyNA(c(ci$c_index, ci$CI_lower, ci$CI_upper)))
  expect_true(ci$CI_lower <= ci$CI_upper)
})
