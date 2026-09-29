test_that("preprocess_features removes near-constant and correlated features and keeps the outcome", {
  set.seed(1)
  d <- data.frame(a = rnorm(50), const = 1, target = factor(rep(c("no", "yes"), 25)))
  d$a_copy <- d$a + rnorm(50, sd = 1e-3)   # |r| > 0.9 with a
  d$b <- rnorm(50)

  out <- pipeML:::preprocess_features(d, target_col = "target", cor_thresh = 0.9)

  expect_false("const" %in% colnames(out))
  expect_true(xor("a" %in% colnames(out), "a_copy" %in% colnames(out)))
  expect_true(all(c("b", "target") %in% colnames(out)))
})

test_that("default survival hyperparameter grids have no duplicated configurations", {
  set.seed(1)
  X <- data.frame(g1 = rnorm(60), g2 = rnorm(60), g3 = rnorm(60))  # small data: values are capped
  for (m in c("proportional_hazards_glmnet", "decision_tree_partykit", "bag_tree_rpart", "rand_forest_aorsf")) {
    grid <- tidyr::expand_grid(!!!pipeML:::get_default_hyperparams(m, train_x = X, v = 2))
    expect_equal(nrow(grid), nrow(dplyr::distinct(grid)), info = m)
  }
})

test_that("data_example_survival codes the event as 0/1", {
  expect_true(all(pipeML::data_example_survival$status %in% c(0, 1)))
  expect_true(all(c("time", "status") %in% colnames(pipeML::data_example_survival)))
})

test_that("ensure_caret() attaches caret", {
  pipeML:::ensure_caret()
  expect_true("package:caret" %in% search())
})
