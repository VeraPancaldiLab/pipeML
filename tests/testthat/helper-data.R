# Small simulated datasets shared by the tests

# Classification: the response depends on g1 + g2; g3 is noise
sim_classification <- function(n = 80, seed = 1) {
  set.seed(seed)
  X <- data.frame(g1 = rnorm(n), g2 = rnorm(n), g3 = rnorm(n))
  rownames(X) <- paste0("P", seq_len(n))
  y <- ifelse(X$g1 + X$g2 + rnorm(n, sd = 0.7) > 0, "R", "NR")
  list(X = X, y = y)
}

# Survival: the risk depends on g1 + g2; g3 is noise
sim_survival <- function(n = 60, seed = 1) {
  set.seed(seed)
  X <- data.frame(g1 = rnorm(n), g2 = rnorm(n), g3 = rnorm(n))
  rownames(X) <- paste0("P", seq_len(n))
  list(X = X, time = stats::rexp(n, exp(0.8 * (X$g1 + X$g2))), event = stats::rbinom(n, 1, 0.8))
}

# Run the code of a test in a temporary working directory. pipeML writes to Results/, which it creates in the
# working directory when the package is loaded, so it is created here too
local_temp_wd <- function(env = parent.frame()) {
  dir <- tempfile("pipeML_test_")
  dir.create(file.path(dir, "Results"), recursive = TRUE)
  old <- setwd(dir)
  do.call(on.exit, list(bquote(setwd(.(old))), add = TRUE), envir = env)
  invisible(dir)
}
