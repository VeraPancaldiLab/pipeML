test_that("preprocess_features removes near-constant and correlated features", {
  set.seed(1)
  target <- rep(c("no", "yes"), 25)
  d <- data.frame(a = rnorm(50), const = 1, row.names = paste0("S", 1:50))
  d$a_copy <- d$a + rnorm(50, sd = 1e-3)   # |r| > 0.9 with a
  d$b <- rnorm(50)

  out <- pipeML:::preprocess_features(d, cor_thresh = 0.9)

  expect_false("const" %in% colnames(out))
  expect_true(xor("a" %in% colnames(out), "a_copy" %in% colnames(out)))
  expect_true("b" %in% colnames(out))
  expect_equal(rownames(out), rownames(d))

  # a feature that is constant in one class and varies in the other separates the classes: it is kept
  d$marker <- ifelse(target == "yes", rnorm(50, mean = 5), 0)
  expect_true("marker" %in% colnames(pipeML:::preprocess_features(d)))
})

test_that("default survival hyperparameter grids have no duplicated configurations", {
  set.seed(1)
  X <- data.frame(g1 = rnorm(60), g2 = rnorm(60), g3 = rnorm(60))  # small data: values are capped
  for (m in c("proportional_hazards_glmnet", "decision_tree_partykit", "bag_tree_rpart", "rand_forest_aorsf")) {
    grid <- tidyr::expand_grid(!!!pipeML:::get_default_hyperparams(m, train_x = X, v = 2))
    expect_equal(nrow(grid), nrow(dplyr::distinct(grid)), info = m)
  }
})

test_that("default survival grids only tune arguments that reach the engine, and work without train_x", {
  # parsnip does not pass cost_complexity to partykit, nor learn_rate / sample_size / stop_iter to mboost
  expect_setequal(names(pipeML:::get_default_hyperparams("decision_tree_partykit")), c("tree_depth", "min_n"))
  expect_setequal(names(pipeML:::get_default_hyperparams("boost_tree_mboost")), c("trees", "min_n", "tree_depth"))
  # without train_x the forests have no mtry (engine default) instead of an error
  expect_setequal(names(pipeML:::get_default_hyperparams("rand_forest_aorsf")), c("trees", "min_n"))
  X <- data.frame(g1 = rnorm(60), g2 = rnorm(60), g3 = rnorm(60))
  expect_true(all(pipeML:::get_default_hyperparams("rand_forest_aorsf", train_x = X, v = 2)$mtry <= 3))
})

test_that("data_example_survival codes the event as 0/1", {
  expect_true(all(pipeML::data_example_survival$status %in% c(0, 1)))
  expect_true(all(c("time", "status") %in% colnames(pipeML::data_example_survival)))
})

test_that("ensure_caret() attaches caret", {
  pipeML:::ensure_caret()
  expect_true("package:caret" %in% search())
})

test_that("preprocess_features stops with a clear message on invalid inputs", {
  set.seed(1)
  d <- data.frame(a = rnorm(40), b = rnorm(40))
  expect_error(preprocess_features(cbind(d, grp = rep(c("x", "z"), each = 20))),
               "All features must be numeric. Non-numeric features: grp")
  expect_error(preprocess_features(data.frame(const = rep(1, 40))), "No feature is left after preprocessing")
})

test_that("get_tune_grid() uses log-spaced lambda and estimates the svmRadial sigma from the features", {
  lambda <- pipeML:::get_tune_grid("lasso", NULL)$lambda
  expect_length(lambda, 20)
  expect_equal(range(lambda), c(0.001, 1))
  expect_equal(sum(lambda < 0.05), 11) # log-spaced (2 with the former evenly spaced grid): values spread over the small lambdas too

  set.seed(1)
  few  <- data.frame(matrix(rnorm(60 * 3), 60), target = "no")
  many <- data.frame(matrix(rnorm(60 * 300), 60), target = "no")
  g_few <- pipeML:::get_tune_grid("svmRadial", few)
  expect_equal(nrow(g_few), 9)
  expect_identical(pipeML:::get_tune_grid("svmRadial", few), g_few) # deterministic: same grid in every fold
  # the estimated sigma shrinks as features are added
  expect_lt(max(pipeML:::get_tune_grid("svmRadial", many)$sigma), min(g_few$sigma))
})
