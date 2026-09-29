test_that("only AUROC and AUPRC are accepted as classification metrics", {
  skip_on_cran()
  local_temp_wd()
  d <- sim_classification()
  expect_error(compute_features.training.ML(features_train = d$X, target_var = d$y, task_type = "classification",
                                            trait.positive = "R", metric = "Accuracy", k_folds = 2, n_rep = 1),
               "Choose either \"AUROC\" or \"AUPRC\"")
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
