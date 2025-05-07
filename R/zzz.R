.onLoad <- function(libname, pkgname) {
  results_dir <- file.path("Results")
  dir.create(results_dir, showWarnings = FALSE, recursive = TRUE)
}
