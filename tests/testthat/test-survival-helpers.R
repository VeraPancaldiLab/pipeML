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
