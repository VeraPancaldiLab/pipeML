
utils::globalVariables(c(
  "obs", "yes", "AUROC", "AUPRC", "Precision", "F1", "MCC", "Kappa", "AccuracySD",
  "KappaSD", "features", "Overall", "importance", "weighted_importance", "Category",
  "Trait", "dataset", "Models", "binomial", "Actual", "Prediction", ".outcome",
  "Mean_AUROC", "Mean_AUPRC", "SD_AUROC", "SD_AUPRC", "Mean_Accuracy", "SD_Accuracy",
  "resample", "rowIndex", "Samples", "color.roc", "color.prc", "Folder", "AUC_roc",
  "AUC_prc", "Cohort", "medianAUROC", "medianAUPRC", "fpr", ".", "meanImp", "Variable",
  "Value", "Decision", "feature", "shap_value", "mean_shap", "direction", "seed",
  "tp", "fp", "fn", "is_yes", "calibrated_yes", "model", "pred", "no", "Resample", "Sensitivity",
  "Specificity", "Accuracy", "value", "job",
  "metrics",
  "Accuracy_resample",
  "Kappa_resample",
  "MAD_AUROC",
  "MAD_AUPRC",
  "MAD_Accuracy",
  "c_index",
  "c_index_median",
  "Median_CINDEX",
  "MAD_CINDEX",
  "time",
  "event",
  ".config_id",
  ".pred",
  "n_resamples",
  "parameter_i"
))

#' Compute Boruta algorithm
#'
#' @param data A data frame with the column of the variable to predict named "target" and the predictor features as additional columns.
#' @param seed A numeric value used to set the random seed for reproducibility.
#' @param fix Logical. If TRUE, applies TentativeRoughFix() from the Boruta package to resolve tentative features.
#'
#' @return A list containing:
#' \itemize{
#'   \item A data frame with feature importance statistics.
#'   \item A character vector indicating the Boruta decision for each feature (Confirmed, Tentative, or Rejected).
#' }
compute_boruta <- function(data, seed, fix = TRUE) {

  set.seed(seed)
  boruta_output <- Boruta::Boruta(target ~ ., data = data, doTrace = 0)

  if (fix) {
    roughFixMod <- Boruta::TentativeRoughFix(boruta_output)
    boruta_output <- roughFixMod
  }

  imps <- Boruta::attStats(boruta_output)
  decision <- as.character(imps$decision)

  res <- imps %>%
    data.frame() %>%
    tibble::rownames_to_column("Variable") %>%
    dplyr::select(-decision)


  return(list(res, decision))
}

#' Merge Boruta Results
#'
#' Merge results from multiple Boruta runs to identify robust feature selections.
#'
#' @param importance_values A list of data frames with feature importance values from each iteration.
#' @param decisions A list of character vectors with the decision labels from each iteration.
#' @param file_name A string used for naming the output plot file.
#' @param iterations Integer. The number of Boruta iterations performed.
#' @param threshold A numeric value between 0 and 1. Features labeled as 'Confirmed' or 'Tentative' in more than \code{threshold * iterations} will be retained.
#' @param return Logical. Whether to save plots in the "Results/" directory.
#'
#' @return A list containing:
#' \itemize{
#'   \item A vector of confirmed features.
#'   \item A vector of tentative features.
#'   \item A data frame with median importance values and final decisions.
#' }
merge_boruta_results = function(importance_values, decisions, file_name, iterations, threshold, return = TRUE){

  ### Construct matrix of importance
  combined_importance <- do.call(rbind, importance_values)
  combined_results_long <- combined_importance %>% #Matrix for plotting
    tidyr::pivot_longer(cols = meanImp, names_to = "Measure", values_to = "Value")

  median_df <- combined_importance %>% #Calculate the median for each column, grouped by the variable name
    dplyr::group_by(Variable) %>%
    dplyr::summarize(dplyr::across(dplyr::everything(), \(x) stats::median(x, na.rm = TRUE)))

  ### Retrieve important and tentatives variables

  combined_results <- do.call(cbind, decisions)
  rownames(combined_results) = median_df$Variable
  decisions_summary <- apply(combined_results, 1, function(x) {
    table(factor(x, levels = c("Confirmed", "Tentative", "Rejected")))
  })
  confirmed_vars <- names(which(decisions_summary["Confirmed",] >= round(threshold*iterations)))
  tentative_vars <- names(which(decisions_summary["Tentative",] >= round(threshold*iterations)))

  # For plotting
  combined_results_long$Decision = "Rejected"
  combined_results_long$Decision[which(combined_results_long$Variable %in% confirmed_vars)] = "Confirmed"
  combined_results_long$Decision[which(combined_results_long$Variable %in% tentative_vars)] = "Tentative"

  mean_order <- median_df %>% #Extract the order of variables for plotting
    dplyr::arrange(meanImp) %>%
    dplyr::pull(Variable)

  # For result
  median_df$Decision = "Rejected"
  median_df$Decision[which(median_df$Variable %in% confirmed_vars)] = "Confirmed"
  median_df$Decision[which(median_df$Variable %in% tentative_vars)] = "Tentative"

  # Plot variable importance boxplots
  if(return){
    grDevices::pdf(paste0("Results/Boruta_variable_importance_", file_name, ".pdf"), width = 8, height = 12)
    print(ggplot2::ggplot(combined_results_long, ggplot2::aes(x = factor(Variable, levels = mean_order), y = Value, fill = Decision)) +
            ggplot2::geom_bar(stat = "identity", position = "dodge") +
            ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, hjust = 1)) +
            ggplot2::coord_flip() +
            ggplot2::labs(x = "Features", y = "Importance", title = paste0("Variable Importance by Boruta after ", iterations, " bootstraps\n", file_name)) +
            ggplot2::scale_fill_manual(values = c("Confirmed" = "green", "Tentative" = "yellow", "Rejected" = "red")) +
            ggplot2::facet_wrap(~ Measure, scales = "free_y"))
    grDevices::dev.off()
  }

  return(list(Confirmed = confirmed_vars, Tentative = tentative_vars, Matrix_Importance = median_df))
}

#' Compute Feature Selection Using Repeated Boruta Algorithm
#'
#' Repeatedly applies the Boruta feature selection algorithm and aggregates results to determine consistently selected features.
#'
#' @param data A data frame with the column "target" (factor) as the response and other columns as features.
#' @param iterations Integer. The number of Boruta iterations to perform.
#' @param fix Logical. If TRUE, applies TentativeRoughFix() to resolve tentative features after each iteration.
#' @param tentative Logical. Whether to include tentative features as confirmed in the training dataset.
#' @param doParallel Logical. Whether to use parallel processing.
#' @param workers Integer. Number of CPU cores to use for parallel execution. If NULL, uses all available cores minus one.
#' @param file_name A string for naming output plots and CSV files saved in the "Results/" directory.
#' @param threshold A numeric value between 0 and 1. A feature must be confirmed in more than \code{threshold * iterations} to be finally labeled as confirmed.
#' @param return Logical. Whether to save the resulting plots in the "Results/" directory.
#' @param verbose Boolen value to whether print or no the function messages
#'
#' @return A list containing:
#' \itemize{
#'   \item A vector of confirmed features.
#'   \item A vector of tentative features.
#'   \item A data frame with median importance values and final decisions.
#' }
#'
#' @export
#'
feature.selection.boruta <- function(data, iterations = NULL, fix = FALSE, tentative = FALSE, doParallel = F, workers=NULL, file_name = NULL, threshold = NULL, return, verbose = T) {
  if(doParallel){
    if(is.null(iterations) == T){
      stop("No iterations specified for running in parallel, please set a number. If you want to run feature selection once consider setting doParallel = F")
    }else{
      if(is.null(workers)==T){
        num_cores <- parallel::detectCores() - 1
      }else{
        num_cores <- workers
      }

      cl = parallel::makeCluster(num_cores) #Forking just copy the R session in its current state. - makeCluster() all must be exported (copied) to the clusters, which can add some overhead
      doParallel::registerDoParallel(cl)

      if(verbose){
        message("Running ", iterations, " iterations of the Boruta algorithm using ", num_cores, " cores")
      }

      res <- foreach::foreach(seed = sample.int(100000, iterations)) %dopar% {

        tryCatch({
          # If successful, return the result and the seed
          list(result = compute_boruta(data, seed, fix),
               error = NULL,
               seed = seed)
        }, error = function(e) {
          # If an error occurs, return the error message and the seed for debugging
          list(result = NULL, error = e$message, seed = seed)
        })
      }

      parallel::stopCluster(cl)
      unregister_dopar() #Stop Dopar from running in the background

    }

    # Extract the first sublist of each element
    matrix_of_importance <- lapply(res, function(x) x[[1]])

    # Extract the second sublist of each element
    features_labels <- lapply(res, function(x) x[[2]])

    res = merge_boruta_results(matrix_of_importance, features_labels, file_name = file_name, iterations = iterations, threshold = threshold, return = return)
  }else{
    if(is.null(iterations) == T){
      stop("No iterations specified for running Boruta algorithm for feature selection")
    }else{
      if(verbose){
        message("Running ", iterations, " iterations of the Boruta algorithm")
      }
      res = list()
      for (i in 1:iterations) {
        res[[i]] = compute_boruta(data, seed = sample.int(100000, 1), fix)
      }

      # Extract the first sublist of each element
      matrix_of_importance <- lapply(res, function(x) x[[1]])

      # Extract the second sublist of each element
      features_labels <- lapply(res, function(x) x[[2]])

      res = merge_boruta_results(matrix_of_importance, features_labels, file_name = file_name, iterations = iterations, threshold = threshold, return = return)
    }

  }

  if(tentative == F){
    if(length(res$Confirmed) <= 1){ #No enough features selected for training model
      if(verbose){
        message("No enough features selected for training a model")
      }
      return(NULL)
    }else{
      if(verbose){
        message("\nKeeping only features confirmed in more than ", threshold," of the times for training...............................................................\n\n")
        cat("If you want to consider also tentative features, please specify tentative = TRUE in the parameters.\n\n")
      }
      train_data = data[,colnames(data)%in%res$Confirmed, drop = F] %>%
        dplyr::mutate(target = data$target)
    }
  }else{
    sum_features = length(res$Confirmed) + length(res$Tentative)
    if(sum_features <= 1){
      if(verbose){
        message("No enough features selected for training a model")
      }
      return(NULL)
    }else{
      if(verbose){
        cat("Keeping features confirmed and tentative in more than 80% of the times for training...............................................................\n\n")
      }
      train_data = data[,colnames(data)%in%c(res$Confirmed, res$Tentative), drop = F] %>%
        dplyr::mutate(target = data$target)
    }
  }

  return(train_data)

}

#' Perform repeated stratified k-fold cross-validation for model training and tuning
#'
#' This function performs repeated stratified k-fold cross-validation on a dataset to train and tune hyperparameters for 13 machine learning methods. Optionally, it can also perform model stacking and Boruta-based feature selection. Performance is evaluated using user-specified metrics such as Accuracy, AUROC, or AUPRC.
#'
#' @param train_data A data frame containing features and a target column named 'target' corresponding to the response variable to predict.
#' @param k_folds Integer. Number of folds for k-fold cross-validation. Default is 5.
#' @param n_rep Integer. Number of repetitions of the k-fold cross-validation. Default is 100.
#' @param stacking Logical. Whether to perform model stacking. Default is FALSE.
#' @param metric Character. Metric used for hyperparameter tuning and model evaluation. Supported values are "Accuracy", "AUROC", and "AUPRC".
#' @param file_name Character. File name used for saving output plots in the `Results/` directory.
#' @param LODO Logical. If TRUE, performs Leave-One-Dataset-Out (LODO) cross-validation by stratifying folds based on cohort membership.
#' @param ncores Integer. Number of cores to use for parallelization. If not given, detectCores() - 1 will be used.
#' @param return Logical. Whether to return the results and generated plots.
#' @param fold_construction_fun Function. A custom function used to construct the cross-validation folds.
#' This function must accept a \code{bestune} argument, which is used internally to inject optimized parameters
#' after hyperparameter tuning. If \code{bestune = NULL}, the function will explore a parameter grid across folds
#' (parallelized with \code{foreach}); if \code{bestune} is provided, the optimized parameters will be applied
#' to rebuild the features on the full training data.
#' @param fold_construction_args_fixed List. A list of arguments passed to \code{fold_construction_fun}
#' that remain fixed during both cross-validation and final training.
#' @param fold_construction_args_tunable List. A list of arguments passed to \code{fold_construction_fun}
#' that define the hyperparameters to be tuned during cross-validation. Each element should contain
#' candidate values for tuning.
#'
#' @return A list containing:
#' \itemize{
#'   \item Features used during training
#'   \item The selected machine learning model
#'   \item All trained machine learning models
#' }
#'
#' If \code{stacking = TRUE}, the list will also include:
#' \itemize{
#'   \item Base models
#'   \item Meta-learner
#'   \item Matrix of weighted feature importance (see \code{calculate_feature_importance_stacking()})
#' }
#'
#'
compute_k_fold_CV = function(train_data, k_folds, n_rep, stacking = FALSE, metric = "Accuracy", file_name = NULL, LODO = FALSE,
                             ncores = NULL, return = FALSE, fold_construction_fun = NULL,
                             fold_construction_args_fixed = NULL,
                             fold_construction_args_tunable = NULL){

  ml_methods_names <- c("treebag", "rf", "C5.0",
                        "glmnet", "knn", "rpart", "lasso", "ridge",
                        "svmRadial", "svmLinear", "xgbTree")

  if(!(metric %in% c("AUROC", "AUPRC","Accuracy"))){
    stop("The metric assigned is not supported. Choose either accuracy or AUC.")
  }

  ######### Machine Learning models

  ######### Stratify K fold cross-validation
  if(LODO == T){
    multifolds = construct_stratified_cohort_folds(train_data, 'dataset', 'target', k_folds = k_folds, n_rep = n_rep)
    train_data = train_data %>% dplyr::select(-dataset) #After creating multifolds we remove this variable to be able to train
  }else{
    multifolds = caret::createMultiFolds(train_data[,'target'], k = k_folds, times = n_rep) #repeated folds
  }

  ### Implement parallelization
  if(is.null(ncores) == TRUE){
    ncores = parallel::detectCores() - 2
  }

  if(is.null(fold_construction_fun)){ #No custom function provided, using normal CV

    custom_outputs = NULL

    cl <- parallel::makeCluster(ncores)
    doParallel::registerDoParallel(cl)

    trainControl <- caret::trainControl(index = multifolds, method="repeatedcv", number=k_folds, repeats=n_rep, verboseIter = F, allowParallel = T, classProbs = TRUE, savePredictions=T)

    ##################################################### ML models
    #To do: Re-calculate accuracy values based on tuning parameters optimized by the cv AUC - now the values are based on accuracy! be careful

    ################## Bagged CART
    fit.treebag <- caret::train(target~., data = train_data, method = "treebag", metric = "Accuracy",trControl = trainControl)

    ################## RF
    fit.rf <- caret::train(target~., data = train_data, method = "rf", metric = "Accuracy",trControl = trainControl)

    ################## C5.0
    fit.c50 <- caret::train(target~., data = train_data, method = "C5.0", metric = "Accuracy",trControl = trainControl)

    ################## GLMNET - Regularized Logistic Regression (Elastic net)
    fit.glmnet <- caret::train(target~., data = train_data, method="glmnet", metric="Accuracy",trControl=trainControl)

    ################## KNN - k-Nearest Neighbors
    fit.knn <- caret::train(target~., data = train_data, method="knn", metric="Accuracy",trControl=trainControl)

    ################## CART - Classification and Regression Trees (CART),
    fit.cart <- caret::train(target~., data = train_data, method="rpart", metric="Accuracy",trControl=trainControl)

    ################## Regularized Lasso
    fit.lasso <- caret::train(target~., data = train_data, method="glmnet", metric="Accuracy",trControl=trainControl, tuneGrid = expand.grid(alpha = 1, lambda = seq(0.001, 1, length = 20)))

    ################## Ridge regression
    fit.ridge <- caret::train(target~., data = train_data, method="glmnet", metric="Accuracy",trControl=trainControl, tuneGrid = expand.grid(alpha = 0, lambda = seq(0.001, 1, length = 20)))

    ################## Support Vector Machine with Radial Kernel
    fit.svm_radial <- caret::train(target ~ ., data = train_data, method = "svmRadial", metric = "Accuracy", trControl = trainControl)

    ################## Support Vector Machine with Linear Kernel
    fit.svm_linear <- caret::train(target ~ ., data = train_data, method = "svmLinear", metric = "Accuracy", trControl = trainControl)

    ################## XGboost
    # param_grid <- expand.grid(
    #   nrounds = seq(from = 200, to = nrounds, by = 50),
    #   eta = c(0.025, 0.05, 0.1, 0.3),
    #   max_depth = c(2, 3, 4, 5, 6),
    #   gamma = 0,
    #   colsample_bytree = 1,
    #   min_child_weight = 1,
    #   subsample = 1
    # )

    #### Set allowParallel = F
    #Disable parallelization in xgbTree cause the xgboost algorithm has its own internal parallelization, controlled by the nthread parameter — it uses this to speed up tree construction within a single model.
    #The caret package, on the other hand, can parallelize between models — for example, it can train different cross-validation folds or hyperparameter combinations at the same time if a parallel backend is registered
    #https://stackoverflow.com/questions/39528392/parallel-processing-with-xgboost-and-caret
    #If both are ON it can slower performance (lead to over-parallelization and CPU contention)
    trainControl <- caret::trainControl(index = multifolds, method="repeatedcv", number=k_folds, repeats=n_rep, verboseIter = F, allowParallel = F, classProbs = TRUE, savePredictions=T)

    invisible(utils::capture.output({fit.xgbTree <- caret::train(target~., data=train_data, method="xgbTree", metric = "Accuracy", trControl=trainControl)}, type = "output"))

    parallel::stopCluster(cl)  # stop the cluster after parallel execution
    unregister_dopar() #Stop Dopar from running in the background

    # Store models in a named list
    models <- list(
      BAG = fit.treebag,
      RF = fit.rf,
      C50 = fit.c50,
      GLMNET = fit.glmnet,
      KNN = fit.knn,
      CART = fit.cart,
      LASSO = fit.lasso,
      RIDGE = fit.ridge,
      SVM_radial = fit.svm_radial,
      SVM_linear = fit.svm_linear,
      XGboost = fit.xgbTree
    )

  }else{

    # Custom fold construction (is running in parallel)
    do.call(fold_construction_fun, c(list(data = train_data, folds = multifolds), fold_construction_args_fixed, fold_construction_args_tunable))

    ### Extract the file names of the folds
    result_files <- list.files("Results", pattern = "^fold_.*\\.rds$", full.names = TRUE)
    fold_data = vector("list", length(result_files))

    # Initialize master list to store everything in memory
    models_all_folds <- vector("list", length(result_files))

    # Iterate across folds and inside each subfold corresponding to each param combination (if exist)
    for (fold_i in seq_along(result_files)) { ### number of folds (k_fold x n_rep)
      cat("\nRunning ML models with fold", fold_i, "\n")

      result = readRDS(result_files[[fold_i]]) ## per resample

      # Each fold contains multiple parameter sets (list of lists) --> fold_construction_args_tunable != NULL
      if (!is.null(fold_construction_args_tunable)) {

        models_all_params <- vector("list", length(result))

        for (parameter_i in seq_along(result)) {

          train_data_i <- result[[parameter_i]][["train_data"]]
          test_data_i  <- result[[parameter_i]][["test_data"]]

          # Preprocessing features (remove collinear variables and no-variance)
          train_data_i <- preprocess_features(train_data_i, cor_thresh = 0.9, target_col = "target")

          # Replace in original train/test datasets
          result[[parameter_i]][["train_data"]] <- train_data_i
          result[[parameter_i]][["test_data"]]  <- test_data_i[, setdiff(colnames(train_data_i), "target")]

          # Run all ML methods for this parameter configuration
          models <- lapply(
             ml_methods_names,
             function(method) {
               tune_grid <- get_tune_grid(method, train_data)

               model_name <- if (method %in% c("lasso", "ridge")) "glmnet" else method

               do.call(compute_custom_k_fold_CV,
                         list(processed_folds = result[[parameter_i]],
                              ml_method = model_name,
                              tuneGrid = tune_grid))
            }
        )

        models_all_params[[parameter_i]] <- models
      }

      # Store all parameter results for this fold
      models_all_folds[[fold_i]] <- models_all_params

      }else { # Custom function does not have hyperparams to tune --> fold_construction_args_tunable != NULL
        train_data_i <- result[["train_data"]]
        test_data_i <- result[["test_data"]]

        # Preprocessing
        train_data_i <- preprocess_features(train_data_i, cor_thresh = 0.9, target_col = "target")

        # Replace
        result[["train_data"]] = train_data_i
        result[["test_data"]] = test_data_i[, setdiff(colnames(train_data_i), "target")]

        models = lapply(
          ml_methods_names,
          function(method){

            tune_grid <- get_tune_grid(method, train_data)

            model_name <- if (method %in% c("lasso", "ridge")) "glmnet" else method

            # Custom CV validation and hyperparameter tuning
            do.call(compute_custom_k_fold_CV,
                    list(processed_folds = result,
                         ml_method = model_name,
                         tuneGrid = tune_grid))


          }
        )

        models_all_folds[[fold_i]] <- models
      }
    }

    file.remove(result_files) ## Delete the files after using them
    agg <- aggregate_results(models_all_folds, task = 'classification')

    ## Sanity check (each param conf has to be evaluated in all resamples)
    for(i in 1:length(agg)){
      hp_cols_all = names(agg[[i]][["bestTune"]]) ### Hyperparameter names
      x = agg[[i]][["Prediction_folds"]] %>%
        dplyr::distinct(Resample, dplyr::across(all_of(hp_cols_all))) %>%
        dplyr::count(dplyr::across(all_of(hp_cols_all)), name = "n_resamples") %>%
        dplyr::arrange(desc(n_resamples))

      # Expected number of resamples (folds × repeats)
      expected_resamples <- length(unique(agg[[i]][["Prediction_folds"]]$Resample))

      # Sanity check
      if (any(x$n_resamples != expected_resamples)) {
        stop("Inconsistent number of resamples detected for parameter configuration\n",
                    paste0(hp_cols_all, collapse = " "))
      }
    }

    model_names <- c(
      BAG = "treebag",
      RF = "rf",
      C50 = "C5.0",
      GLMNET = "glmnet",
      KNN = "knn",
      CART = "rpart",
      LASSO = "lasso",
      RIDGE = "ridge",
      SVM_radial = "svmRadial",
      SVM_linear = "svmLinear",
      XGboost = "xgbTree"
    )

    ################################ Train model with optimized hyperparameters

    if (!is.null(fold_construction_args_tunable)) {
      optimized_models <- lapply(seq_along(ml_methods_names), function(i) {
        wrapper_train_best_hyperparams_classification(
          train_data,
          agg[[i]],                     # model-specific aggregated results
          ml_methods_names[i],
          fold_construction_fun,
          fold_construction_args_fixed
        )
      })

      # Split components across lists
      training_sets  <- lapply(optimized_models, `[[`, "training_set")
      custom_outputs <- lapply(optimized_models, `[[`, "custom_output")
      models         <- lapply(optimized_models, `[[`, "Model")

      # Assign pretty names
      names(training_sets) <- names(model_names)
      names(custom_outputs) <- names(model_names)
      names(models) <- names(model_names)
    }else{
      models = agg
      names(models) <- names(model_names)
    }

  }

  ####### Optimized based on metric (only AUC or Accuracy available)
  if(metric == "AUROC" || metric == "AUPRC"){

    if(is.null(fold_construction_fun)){

      hyperparams <- list(
        BAG = NULL,
        RF = "mtry",
        C50 = c("trials", "model", "winnow"),
        GLMNET = c("alpha", "lambda"),
        KNN = "k",
        CART = "cp",
        LASSO = c("alpha", "lambda"),
        RIDGE = c("alpha", "lambda"),
        SVM_radial = c("sigma", "C"),
        SVM_linear = "C",
        XGboost = c("nrounds", "max_depth", "eta", "gamma", "colsample_bytree", "min_child_weight", "subsample")
      )

      # Iterate over models
      models <- lapply(names(models), function(name) {
        model <- models[[name]]
        names(model)[names(model) == "pred"] <- "Prediction_folds"
        names(model)[names(model) == "results"] <- "Results_folds"
        names(model)[names(model) == "resample"] <- "Resample_matrix"

        res <- calculate_cv_metrics(model, metric, hyperparams[[name]])

        # Remove transition objects (might change after with better code)
        model$Prediction_folds = NULL
        model$Results_folds = NULL
        model$Resample_matrix = NULL

        # Replace model fields (we need the original names cause at one point when we run default functions from caret it uses these names)
        model$pred <- res$Prediction_folds
        model$resample <- res$Resample_matrix
        model$results <- res$Results_folds
        model$bestTune = res$bestTune

        return(model)
      })

      custom_outputs = NULL

    }else{ #### There is custom fold function

      if(!is.null(fold_construction_args_tunable)){
        hyperparams <- list(
          BAG = names(fold_construction_args_tunable),
          RF = c("mtry",names(fold_construction_args_tunable)),
          C50 = c("trials", "model", "winnow", names(fold_construction_args_tunable)),
          GLMNET = c("alpha", "lambda", names(fold_construction_args_tunable)),
          KNN = c("k", names(fold_construction_args_tunable)),
          CART = c("cp",names(fold_construction_args_tunable)),
          LASSO = c("alpha", "lambda", names(fold_construction_args_tunable)),
          RIDGE = c("alpha", "lambda", names(fold_construction_args_tunable)),
          SVM_radial = c("sigma", "C", names(fold_construction_args_tunable)),
          SVM_linear = c("C", names(fold_construction_args_tunable)),
          XGboost = c("nrounds", "max_depth", "eta", "gamma", "colsample_bytree", "min_child_weight", "subsample", names(fold_construction_args_tunable))
        )

        # Iterate over models and update everything
        updated_results <- lapply(names(models), function(name){
          res <- calculate_cv_metrics(models[[name]], metric, hyperparams[[name]]) #Re-tuned hyperparams based on a different metric than 'Accuracy'

          # Re-calculate training set based on new optimized hyperparams
          training_all <- do.call(
            fold_construction_fun,
            c(list(data = train_data, bestune = res$bestTune), fold_construction_args_fixed)
          )

          # Preprocess features
          training_set <- preprocess_features(training_all[[1]], cor_thresh = 0.9, target_col = "target")

          # Retrieve custom_output
          custom_output <- training_all[[2]]

          model <- models[[name]]
          model$Prediction_folds <- res$Prediction_folds
          model$Resample_matrix <- res$Resample_matrix
          model$Results_folds <- res$Results_folds

          if("none" %in% model$bestTune){
            model$bestTune <- tibble::tibble(parameter = "none")
          }else{
            model$bestTune <- res$bestTune %>% select(-colnames(training_all[[3]]))
          }

          list(
            Model = model,
            training_set = training_set,
            custom_output = c(training_all[[2]], list("Parameters" = training_all[[3]]))
          )
        })

        names(updated_results) <- names(models)

        # Split into three lists
        models <- lapply(updated_results, `[[`, "Model")
        training_sets <- lapply(updated_results, `[[`, "training_set")
        custom_outputs <- lapply(updated_results, `[[`, "custom_output")

      }else{
        hyperparams <- list(
          BAG = NULL,
          RF = "mtry",
          C50 = c("trials", "model", "winnow"),
          GLMNET = c("alpha", "lambda"),
          KNN = "k",
          CART = "cp",
          LASSO = c("alpha", "lambda"),
          RIDGE = c("alpha", "lambda"),
          SVM_radial = c("sigma", "C"),
          SVM_linear = "C",
          XGboost = c("nrounds", "max_depth", "eta", "gamma", "colsample_bytree", "min_child_weight", "subsample")
        )

        # Iterate over models and update everything
        updated_results <- lapply(names(models), function(name){
          res <- calculate_cv_metrics(models[[name]], metric, hyperparams[[name]]) #Re-tuned hyperparams based on a different metric than 'Accuracy' (e.g. AUROC)

          # Re-calculate training set
          training_all <- do.call(
            fold_construction_fun,
            c(list(data = train_data, bestune = res$bestTune), fold_construction_args_fixed)
          )

          # Preprocess features
          training_set <- preprocess_features(training_all[[1]], cor_thresh = 0.9, target_col = "target")

          model <- models[[name]]
          model$Prediction_folds <- res$Prediction_folds
          model$Resample_matrix <- res$Resample_matrix
          model$Results_folds <- res$Results_folds


          list(
            Model = model,
            training_set = training_set,
            custom_output = training_all[[2]]
          )
        })

        names(updated_results) <- names(models)

        # Split into three lists
        models <- lapply(updated_results, `[[`, "Model")
        training_sets <- lapply(updated_results, `[[`, "training_set")
        custom_outputs <- lapply(updated_results, `[[`, "custom_output")
      }
    }
  }

  ############### Collect ML models
  names(models) = names(hyperparams)

  fit.treebag <- models$BAG
  fit.rf <- models$RF
  fit.c50 <- models$C50
  fit.glmnet <- models$GLMNET
  fit.knn <- models$KNN
  fit.cart <- models$CART
  fit.lasso <- models$LASSO
  fit.ridge <- models$RIDGE
  fit.svm_radial <- models$SVM_radial
  fit.svm_linear <- models$SVM_linear
  fit.xgbTree <- models$XGboost

  ############################################# These predictions are use for the meta-learner because it needs the predictions from the models in the complete dataset (might change in the future)

  if(is.null(fold_construction_fun)){

    ###Prediction with best tuned hyper-parameters (Missing to add platt scaling to calibrated probabilities (when tested it didnt converge, need to be checked)) See https://www.cs.cornell.edu/~alexn/papers/calibration.icml05.crc.rev3.pdf

    # -------------------------------------> Missing: Train models with bestTune from CV (only for meta-learner: stacking)

    ###Bagged CART

    predictions.bag <- data.frame(stats::predict(fit.treebag, newdata = train_data, type = "prob")) %>% #Predictions using tuned model
      dplyr::select(yes) %>%
      dplyr::rename(BAG = yes)

    ###Random Forest

    predictions.rf = data.frame(stats::predict(fit.rf, newdata = train_data, type = "prob")) %>%
      dplyr::select(yes) %>%
      dplyr::rename(RF = yes) #Predictions of model (already ordered)

    ###C5.0

    predictions.c50 = data.frame(stats::predict(fit.c50$finalModel, newdata = train_data, type = "prob")) %>%
      dplyr::select(yes) %>%
      dplyr::rename(C50 = yes)  #Predictions of model (already ordered)

    ### LG

    # predictions.glm = stats::predict(fit.glm, newdata = train_data, type = "prob") %>%
    #   data.frame() %>%
    #   dplyr::select(yes) %>%
    #   dplyr::rename(GLM = yes)  #Predictions of model (already ordered)

    ### LDA

    # predictions.lda = stats::predict(fit.lda, newdata = train_data, type = "prob") %>%
    #   data.frame() %>%
    #   dplyr::select(yes) %>%
    #   dplyr::rename(LDA = yes)  #Predictions of model (already ordered)

    ### GLMNET

    predictions.glmnet = stats::predict(fit.glmnet, newdata = train_data, type = "prob") %>%
      data.frame() %>%
      dplyr::select(yes) %>%
      dplyr::rename(GLMNET = yes)  #Predictions of model (already ordered)

    ### KNN

    predictions.knn = stats::predict(fit.knn, newdata = train_data, type = "prob") %>%
      data.frame() %>%
      dplyr::select(yes) %>%
      dplyr::rename(KNN = yes) #Predictions of model (already ordered)

    ## CART

    predictions.cart = stats::predict(fit.cart, newdata = train_data, type = "prob") %>%
      data.frame() %>%
      dplyr::select(yes) %>%
      dplyr::rename(CART = yes)  #Predictions of model (already ordered)

    ## Regularized Lasso

    predictions.lasso = stats::predict(fit.lasso, newdata = train_data, type = "prob") %>%
      data.frame() %>%
      dplyr::select(yes) %>%
      dplyr::rename(LASSO = yes)  #Predictions of model (already ordered)

    ## Ridge regression

    predictions.ridge = stats::predict(fit.ridge, newdata = train_data, type = "prob") %>%
      data.frame() %>%
      dplyr::select(yes) %>%
      dplyr::rename(RIDGE = yes)  #Predictions of model (already ordered)

    ## SVM radial

    predictions.svm_radial = stats::predict(fit.svm_radial, newdata = train_data, type = "prob") %>%
      data.frame() %>%
      dplyr::select(yes) %>%
      dplyr::rename(SVM_radial = yes)  #Predictions of model (already ordered)

    ## SVM linear

    predictions.svm_linear = stats::predict(fit.svm_linear, newdata = train_data, type = "prob") %>%
      data.frame() %>%
      dplyr::select(yes) %>%
      dplyr::rename(SVM_linear = yes)  #Predictions of model (already ordered)

    ## XGboost

    predictions.xgboost = stats::predict(fit.xgbTree, newdata = train_data, type = "prob") %>%
      data.frame() %>%
      dplyr::select(yes) %>%
      dplyr::rename(XGboost = yes)  #Predictions of model (already ordered)

  }else{ # Training data is calculated from custom functions

    cl <- parallel::makeCluster(ncores)
    doParallel::registerDoParallel(cl)

    ######## Bagged CART
    cat("\nRunning BAG....................\n")

    # Train model with bestTune from CV
    temp = fit.treebag
    fit.treebag <- caret::train(
      target ~ .,
      data = training_sets$BAG,
      method = "treebag",
      trControl = trainControl(method = "none", classProbs = TRUE, allowParallel = TRUE),
      tuneGrid = temp$bestTune
    )

    # Return caret-like object
    fit.treebag$results = temp$Results_folds
    fit.treebag$pred = temp$Prediction_folds
    fit.treebag$resample = temp$Resample_matrix
    fit.treebag$bestTune = temp$bestTune

    # Predictions in trained model
    predictions.bag = predict(fit.treebag, newdata = training_sets$BAG, type = "prob") %>%
      data.frame() %>%
      dplyr::select(yes) %>%
      dplyr::rename(BAG = yes) #Predictions of model (already ordered)

    ######## RF
    cat("\nRunning RF....................\n")

    # Train model with bestTune from CV
    temp = fit.rf
    fit.rf <- caret::train(
      target ~ .,
      data = training_sets$RF,
      method = "rf",
      trControl = trainControl(method = "none", classProbs = TRUE, allowParallel = TRUE),
      tuneGrid = temp$bestTune
    )

    # Return caret-like object
    fit.rf$results = temp$Results_folds
    fit.rf$pred = temp$Prediction_folds
    fit.rf$resample = temp$Resample_matrix
    fit.rf$bestTune = temp$bestTune

    # Predictions in trained model
    predictions.rf = predict(fit.rf, newdata = training_sets$RF, type = "prob") %>%
      data.frame() %>%
      dplyr::select(yes) %>%
      dplyr::rename(RF = yes) #Predictions of model (already ordered)

    ######## C5.0
    cat("\nRunning C50....................\n")

    # Train model with bestTune from CV
    temp = fit.c50
    fit.c50 <- caret::train(
      target ~ .,
      data = training_sets$C50,
      method = "C5.0",
      trControl = trainControl(method = "none", classProbs = TRUE, allowParallel = TRUE),
      tuneGrid = temp$bestTune
    )

    # Return caret-like object
    fit.c50$results = temp$Results_folds
    fit.c50$pred = temp$Prediction_folds
    fit.c50$resample = temp$Resample_matrix
    fit.c50$bestTune = temp$bestTune

    # Predictions in trained model
    predictions.c50 = predict(fit.c50, newdata = training_sets$C50, type = "prob") %>%
      data.frame() %>%
      dplyr::select(yes) %>%
      dplyr::rename(C50 = yes) #Predictions of model (already ordered)

    ### GLMNET

    cat("\nRunning GLMNET....................\n")
    # Train model with bestTune from CV
    temp = fit.glmnet
    fit.glmnet <- caret::train(
      target ~ .,
      data = training_sets$GLMNET,
      method = "glmnet",
      trControl = trainControl(method = "none", classProbs = TRUE, allowParallel = TRUE),
      tuneGrid = temp$bestTune
    )

    # Return caret-like object
    fit.glmnet$results = temp$Results_folds
    fit.glmnet$pred = temp$Prediction_folds
    fit.glmnet$resample = temp$Resample_matrix
    fit.glmnet$bestTune = temp$bestTune

    # Predictions in trained model
    predictions.glmnet = predict(fit.glmnet, newdata = training_sets$GLMNET, type = "prob") %>%
      data.frame() %>%
      dplyr::select(yes) %>%
      dplyr::rename(GLMNET = yes) #Predictions of model (already ordered)

    ### KNN

    cat("\nRunning KNN....................\n")
    # Train model with bestTune from CV
    temp = fit.knn
    fit.knn <- caret::train(
      target ~ .,
      data = training_sets$KNN,
      method = "knn",
      trControl = trainControl(method = "none", classProbs = TRUE, allowParallel = TRUE),
      tuneGrid = temp$bestTune
    )

    # Return caret-like object
    fit.knn$results = temp$Results_folds
    fit.knn$pred = temp$Prediction_folds
    fit.knn$resample = temp$Resample_matrix
    fit.knn$bestTune = temp$bestTune

    # Predictions in trained model
    predictions.knn = predict(fit.knn, newdata = training_sets$KNN, type = "prob") %>%
      data.frame() %>%
      dplyr::select(yes) %>%
      dplyr::rename(KNN = yes) #Predictions of model (already ordered)

    ### CART

    cat("\nRunning CART....................\n")
    # Train model with bestTune from CV
    temp = fit.cart
    fit.cart <- caret::train(
      target ~ .,
      data = training_sets$CART,
      method = "rpart",
      trControl = trainControl(method = "none", classProbs = TRUE, allowParallel = TRUE),
      tuneGrid = temp$bestTune
    )

    # Return caret-like object
    fit.cart$results = temp$Results_folds
    fit.cart$pred = temp$Prediction_folds
    fit.cart$resample = temp$Resample_matrix
    fit.cart$bestTune = temp$bestTune

    # Predictions in trained model
    predictions.cart = predict(fit.cart, newdata = training_sets$CART, type = "prob") %>%
      data.frame() %>%
      dplyr::select(yes) %>%
      dplyr::rename(CART = yes) #Predictions of model (already ordered)

    ## Regularized Lasso

    cat("\nRunning Lasso....................\n")
    # Train model with bestTune from CV
    temp = fit.lasso
    fit.lasso <- caret::train(
      target ~ .,
      data = training_sets$LASSO,
      method = "glmnet",
      trControl = trainControl(method = "none", classProbs = TRUE, allowParallel = TRUE),
      tuneGrid = temp$bestTune
    )

    # Return caret-like object
    fit.lasso$results = temp$Results_folds
    fit.lasso$pred = temp$Prediction_folds
    fit.lasso$resample = temp$Resample_matrix
    fit.lasso$bestTune = temp$bestTune

    # Predictions in trained model
    predictions.lasso = predict(fit.lasso, newdata = training_sets$LASSO, type = "prob") %>%
      data.frame() %>%
      dplyr::select(yes) %>%
      dplyr::rename(LASSO = yes) #Predictions of model (already ordered)

    ## Ridge regression

    cat("\nRunning Ridge....................\n")
    # Train model with bestTune from CV
    temp = fit.ridge
    fit.ridge <- caret::train(
      target ~ .,
      data = training_sets$RIDGE,
      method = "glmnet",
      trControl = trainControl(method = "none", classProbs = TRUE, allowParallel = TRUE),
      tuneGrid = temp$bestTune
    )

    # Return caret-like object
    fit.ridge$results = temp$Results_folds
    fit.ridge$pred = temp$Prediction_folds
    fit.ridge$resample = temp$Resample_matrix
    fit.ridge$bestTune = temp$bestTune

    # Predictions in trained model
    predictions.ridge = predict(fit.ridge, newdata = training_sets$RIDGE, type = "prob") %>%
      data.frame() %>%
      dplyr::select(yes) %>%
      dplyr::rename(RIDGE = yes) #Predictions of model (already ordered)

    ## SVM radial

    cat("\nRunning SVM radial....................\n")
    # Train model with bestTune from CV
    temp = fit.svm_radial
    fit.svm_radial <- caret::train(
      target ~ .,
      data = training_sets$SVM_radial,
      method = "svmRadial",
      trControl = trainControl(method = "none", classProbs = TRUE, allowParallel = TRUE),
      tuneGrid = temp$bestTune
    )

    # Return caret-like object
    fit.svm_radial$results = temp$Results_folds
    fit.svm_radial$pred = temp$Prediction_folds
    fit.svm_radial$resample = temp$Resample_matrix
    fit.svm_radial$bestTune = temp$bestTune

    # Predictions in trained model
    predictions.svm_radial = predict(fit.svm_radial, newdata = training_sets$SVM_radial, type = "prob") %>%
      data.frame() %>%
      dplyr::select(yes) %>%
      dplyr::rename(SVM_radial = yes) #Predictions of model (already ordered)

    ### SVM linear

    cat("\nRunning SVM linear....................\n")
    # Train model with bestTune from CV
    temp = fit.svm_linear
    fit.svm_linear <- caret::train(
      target ~ .,
      data = training_sets$SVM_linear,
      method = "svmLinear",
      trControl = trainControl(method = "none", classProbs = TRUE, allowParallel = TRUE),
      tuneGrid = temp$bestTune
    )

    # Return caret-like object
    fit.svm_linear$results = temp$Results_folds
    fit.svm_linear$pred = temp$Prediction_folds
    fit.svm_linear$resample = temp$Resample_matrix
    fit.svm_linear$bestTune = temp$bestTune

    # Predictions in trained model
    predictions.svm_linear = predict(fit.svm_linear, newdata = training_sets$SVM_linear, type = "prob") %>%
      data.frame() %>%
      dplyr::select(yes) %>%
      dplyr::rename(SVM_linear = yes) #Predictions of model (already ordered)

    ## XGboost

    cat("\nRunning XGboost....................\n")
    # Train model with bestTune from CV
    temp = fit.xgbTree
    fit.xgbTree <- caret::train(
      target ~ .,
      data = training_sets$XGboost,
      method = "xgbTree",
      trControl = trainControl(method = "none", classProbs = TRUE, allowParallel = FALSE),
      tuneGrid = temp$bestTune
    )

    # Return caret-like object
    fit.xgbTree$results = temp$Results_folds
    fit.xgbTree$pred = temp$Prediction_folds
    fit.xgbTree$resample = temp$Resample_matrix
    fit.xgbTree$bestTune = temp$bestTune

    # Predictions in trained model
    predictions.xgboost = predict(fit.xgbTree, newdata = training_sets$XGboost, type = "prob") %>%
      data.frame() %>%
      dplyr::select(yes) %>%
      dplyr::rename(XGboost = yes) #Predictions of model (already ordered)

    parallel::stopCluster(cl)  # stop the cluster after parallel execution
    unregister_dopar() #Stop Dopar from running in the background

  }
  ############################################################## Save models

  ensembleResults <- list(BAG = fit.treebag,
                          RF = fit.rf,
                          C50 = fit.c50,
                          GLMNET = fit.glmnet,
                          KNN = fit.knn,
                          CART = fit.cart,
                          LASSO = fit.lasso,
                          RIDGE = fit.ridge,
                          SVM_radial = fit.svm_radial,
                          SVM_linear = fit.svm_linear,
                          XGboost = fit.xgbTree)

  ml_methods = list(BAG = "treebag",
                    RF = "rf",
                    C50 = "C5.0",
                    GLMNET = "glmnet",
                    KNN = "knn",
                    CART = "rpart",
                    LASSO = "glmnet",
                    RIDGE = "glmnet",
                    SVM_radial = "svmRadial",
                    SVM_linear = "svmLinear",
                    XGboost = "xgbTree")

  model_predictions = list(BAG = predictions.bag,
                           RF = predictions.rf,
                           C50 = predictions.c50,
                           GLMNET = predictions.glmnet,
                           KNN = predictions.knn,
                           CART = predictions.cart,
                           LASSO = predictions.lasso,
                           RIDGE = predictions.ridge,
                           SVM_radial = predictions.svm_radial,
                           SVM_linear = predictions.svm_linear,
                           XGboost = predictions.xgboost)

  #Remove models with same predictions across samples (not able to make distinction)
  model_predictions <- lapply(model_predictions, function(df) {
    df = df %>%
      dplyr::select(dplyr::where(~ dplyr::n_distinct(.) > 1))

    if(ncol(df) == 0){
      df = NULL
    }

    return(df)
  })

  model_predictions = Filter(Negate(is.null), model_predictions) #Discard not useful predictions
  ensembleResults = ensembleResults[names(model_predictions)] #Discard not useful models based on predictions

  model_predictions = do.call(cbind, model_predictions) #Join as data frame

  #Clean memory
  rm(fit.treebag, fit.rf, fit.c50, fit.knn, fit.cart, fit.glmnet, fit.lasso, fit.ridge, fit.svm_radial, fit.svm_linear)
  gc()


  if(stacking){
    features = colnames(train_data)[colnames(train_data) != "target"]

    #Base models using ML models with best accuracy or AUC from each family
    if(metric == "Accuracy"){
      base_models = compute_cv_accuracy(ensembleResults, base_models = T, file_name = file_name, return = return)
    }else if(metric == "AUROC" || metric == "AUPRC"){
      base_models = compute_cv_AUC(ensembleResults, base_models = T, file_name = file_name, AUC_type = metric, return = return)
    }

    cat("Meta-learners ML model based on GLM\n")

    features_predictions = model_predictions %>%
      t() %>%
      data.frame() %>%
      tibble::rownames_to_column("Models") %>%
      dplyr::filter(grepl(paste0("\\b(", paste(base_models$Base_models, collapse = "|"), ")\\b"), Models)) %>%
      tibble::column_to_rownames("Models") %>%
      t() %>%
      data.frame()

    meta_features = cbind(features_predictions, "true_label" = train_data$target)

    trainControl <- caret::trainControl(index = multifolds, method="repeatedcv", number=k_folds, repeats=n_rep, verboseIter = F, allowParallel = F, classProbs = TRUE, savePredictions=T)
    meta_learner <- caret::train(true_label ~ ., data = meta_features, method = "glmnet", trControl = trainControl) #Staking based on simple logistic regression

    output = list("Meta_learner" = meta_learner, "Base_models" = base_models$Base_models, "ML_models" = ensembleResults)

    ####################################################################### To be done, which output to retrieve when stacking is done? Multiple ML models used different cell groups depending on optimization
    # if(is.null(custom_output) == F){
    #   output[[length(output)+1]] = custom_outputs[[top_model]]
    #   names(output)[length(output)] = "Custom_output"
    # }

  }else{

    #Top model with best accuracy or AUC
    if(metric == "Accuracy"){
      metrics = compute_cv_accuracy(ensembleResults, file_name = file_name, return = return)
    }else if(metric == "AUROC" || metric == "AUPRC"){
      metrics = compute_cv_AUC(ensembleResults, file_name = file_name, AUC_type = metric, return = return)
    }

    top_model = metrics[["Top_model"]]
    AUROC_median = metrics[["AUROC"]]
    AUPRC_median = metrics[["AUPRC"]]

    model = ensembleResults[[top_model]]

    cat("Best ML model found: ", top_model, "\n")

    cat("Returning model trained\n")

    output = list("Model" = model, "ML_Models" = ensembleResults, "AUROC_median" = AUROC_median, "AUPRC_median" = AUPRC_median)

    if(!is.null(custom_outputs) && !any(sapply(custom_outputs, is.null))){ #Check whether custom_output exists or not
      output[[length(output)+1]] = custom_outputs[[top_model]]
      names(output)[length(output)] = "Custom_output"
    }

  }


  return(output)

}

#' Train and evaluate machine learning models on previously constructed k folds
#'
#' This function performs k-fold cross-validation using custom folds created from custom functions to be used for cohort-dependent algorithms (see vignette for more information about this).
#' It supports hyperparameter tuning over a grid and returns a model object that mimicks the caret's training output, including performance metrics and predictions.
#'
#' @param processed_folds A list of folds. Each fold contains processed training and test data with features.
#' @param ml_method A character string indicating the machine learning model to use, as supported by the `caret` package (e.g., `"rf"`, `"svmRadial"`, `"glmnet"`).
#' @param tuneGrid Optional. A data frame specifying the grid of hyperparameters to evaluate. If `NULL`, a default grid of length 3 is generated using caret's `getModelInfo()`.
#'
#' @return A list with the following components:
#' \itemize{
#'   \item \code{Results_folds}: A data frame summarizing average cross-validated Accuracy, Kappa, and their standard deviations for each hyperparameter combination.
#'   \item \code{Prediction_folds}: A data frame of predictions from each fold, including class probabilities, observed and predicted labels, and hyperparameter values.
#'   \item \code{Resample_matrix}: A data frame summarizing Accuracy and Kappa per fold for the best-tuned model.
#'   \item \code{Besttune}: List of optimized hyperparameters.
#' }
#'
#' @details
#' This function performs the following:
#' \enumerate{
#'   \item Trains models for each fold and hyperparameter combination.
#'   \item Predicts on the held-out test data of each fold.
#'   \item Aggregates prediction results and evaluates Accuracy and Kappa for each fold and hyperparameter set.
#'   \item Selects the best-performing hyperparameter set based on mean Accuracy across folds.
#'   \item Trains the final model on the full dataset using the selected hyperparameters.
#' }
#'
#' @importFrom dplyr tibble bind_cols group_by summarise rename ungroup select desc slice_max arrange all_of across
#' @importFrom tidyr unnest_wider
#' @importFrom caret train trainControl getModelInfo
#' @importFrom stats predict
#' @importFrom rlang .data
#'
#' @export
#'
compute_custom_k_fold_CV <- function(processed_folds, ml_method, tuneGrid) {

  train_data = processed_folds[["train_data"]]
  test_data = processed_folds[["test_data"]]

  all_preds = list()

  for (grid_row in seq(nrow(tuneGrid))) {
    # Extract hyperparameters
    hp <- tuneGrid[grid_row, , drop = FALSE]

    # Train model
    model <- suppressWarnings({caret::train(
      target ~ .,
      data = train_data,
      method = ml_method,
      trControl = caret::trainControl(method = "none", classProbs = TRUE, allowParallel = (ml_method != "xgbTree")),
      tuneGrid = hp,
      metric = "Accuracy"
    )})

    # Predict
    test_data <- test_data[, colnames(test_data) %in% model$coefnames]
    probs <- stats::predict(model, newdata = test_data, type = "prob")
    preds <- stats::predict(model, newdata = test_data)

    # Prepare results
    rownames(hp) <- NULL
    hp <- tibble::as_tibble(hp)

    pred_df <- dplyr::tibble(
      rowIndex = processed_folds$rowIndex,
      Resample = processed_folds$fold_name,
      obs = processed_folds$obs_test,
      pred = preds
    ) %>%
      dplyr::bind_cols(hp, probs)

    if (!is.null(processed_folds$params)) {
      pred_df <- dplyr::bind_cols(pred_df, processed_folds$params)
    }

    all_preds[[grid_row]] = pred_df
  }

  ## Combine predictions
  pred_df_all <- do.call(rbind, all_preds)

  return(list(pred_df_all, names(tuneGrid)))
}


#' Train machine learning or survival models with optional stacking and custom cross-validation
#'
#' This function trains one or more machine learning models using repeated k-fold cross-validation,
#' with optional model stacking, feature selection, and support for both classification and survival tasks.
#' It allows flexible cross-validation schemes, including:
#' \itemize{
#'   \item Standard stratified k-fold cross-validation
#'   \item Leave-One-Dataset-Out (LODO) stratified folds by cohort
#'   \item User-defined custom fold construction via a \code{fold_construction_fun}
#' }
#'
#' The function supports both classification and survival analysis pipelines by setting
#' \code{task_type = "classification"} or \code{task_type = "survival"}.
#'
#' @param features_train A data frame containing the features used for training (samples in rows, features in columns).
#' @param task_type Character. Specifies the type of prediction task. Either \code{"classification"} or \code{"survival"}.
#' @param target_var Vector. The target variable to predict (required for classification tasks).
#' @param trait.positive Value in \code{target_var} that represents the positive class (used for metrics like AUROC, AUPRC, confusion matrix, and SHAP values).
#'   This ensures that all performance metrics and interpretability analyses consistently treat the correct class as positive,
#'   which is important when datasets encode classes differently (e.g., "0"/"1", "R"/"NR", "pos"/"neg").
#' @param time_var Character. The name of the survival time variable (required for survival models).
#' @param event_var Character. The name of the event indicator variable (required for survival models; 1 = event occurred, 0 = censored).
#' @param metric Character. Performance metric used for model selection and tuning. Supported values are:
#'   \itemize{
#'     \item \code{"Accuracy"} — classification accuracy
#'     \item \code{"AUROC"} — area under the ROC curve
#'     \item \code{"AUPRC"} — area under the precision-recall curve
#'     \item \code{"C-index"} — concordance index (for survival tasks)
#'   }
#' @param stack Logical. Whether to perform model stacking (ensemble meta-learning). Default is \code{FALSE}.
#' @param k_folds Integer. Number of folds to use for cross-validation.
#' @param n_rep Integer. Number of repetitions for cross-validation (repeated CV).
#' @param LODO Logical. If \code{TRUE}, constructs cross-validation folds stratified by cohort (Leave-One-Dataset-Out scheme).
#' @param batch_var Character. Batch membership for each sample. Required if \code{LODO = TRUE}.
#' @param file_name Character. File name used for saving performance plots in the \code{"Results/"} directory.
#' @param ncores Integer. Number of CPU cores to use for parallelization. Defaults to \code{parallel::detectCores() - 1}.
#' @param return Logical. Whether to return and save the generated plots. Default is \code{FALSE}.
#'
#' @param fold_construction_fun Function. Optional user-defined function to construct cross-validation folds.
#'   This enables full control over how data splits and feature transformations are created.
#'   The function must accept a \code{bestune} argument:
#'   \itemize{
#'     \item If \code{bestune = NULL}, the function explores a parameter grid across folds
#'       (executed in parallel via \code{foreach}).
#'     \item If \code{bestune} is provided, optimized parameters are applied to the full dataset
#'       to rebuild features before final training.
#'   }
#'   The fold constructor should save individual folds as \code{"Results/fold_*.rds"} objects containing:
#'   \itemize{
#'     \item \code{train_data} — training data for that fold
#'     \item \code{test_data} — testing data for that fold
#'     \item \code{obs_test} — observed target or survival outcomes
#'     \item \code{params} — parameter combination used (if applicable)
#'   }
#'
#' @param fold_construction_args_fixed List. Arguments passed to \code{fold_construction_fun} that remain fixed
#'   during both cross-validation and final training (e.g., annotation files, normalization flags, etc.).
#' @param fold_construction_args_tunable List. Arguments passed to \code{fold_construction_fun} that define
#'   hyperparameters to tune during cross-validation. Each element should contain one or more candidate values.
#'
#' @details
#' The function supports:
#' \itemize{
#'   \item Automatic feature preprocessing (e.g., correlation filtering, low-variance removal).
#'   \item Parallelized cross-validation across folds and repetitions.
#'   \item Integration with custom model pipelines (e.g., CellTFusion, pathway-based deconvolution).
#'   \item Unified handling of both survival and classification models.
#' }
#'
#' When a custom fold constructor is provided via \code{fold_construction_fun}, the default stratified k-fold
#' logic is bypassed, and the function will instead iterate through all \code{Results/fold_*.rds} files generated
#' by the custom routine. This allows hybrid pipelines combining biological preprocessing (e.g., CellTFusion)
#' with downstream model fitting.
#'
#' @return A list containing:
#' \itemize{
#'   \item Trained model(s) or meta-learner (if \code{stack = TRUE})
#'   \item Feature set used for model training
#'   \item Cross-validation performance results and plots
#'   \item Best hyperparameter configuration (if applicable)
#' }
#'
#' @seealso [compute_k_fold_CV_survival()], [compute_ml_survival()], [compute_cv_CINDEX()]
#'
#' @export
compute_features.training.ML = function(features_train, task_type = c("classification", "survival"), target_var = NULL, trait.positive = NULL,
                                        time_var = NULL, event_var = NULL, metric = NULL, stack = FALSE, k_folds = 10, n_rep = 5, LODO = FALSE,
                                        batch_var = NULL, file_name = NULL, ncores = NULL, return = FALSE,
                                        fold_construction_fun = NULL, fold_construction_args_fixed = NULL, fold_construction_args_tunable = NULL){

  # ---------------------------------------------------------------------------
  # Validate task_type and required arguments
  # ---------------------------------------------------------------------------

  if (task_type == "classification") {
    if (is.null(target_var) || is.null(trait.positive)) {
      stop("For classification, both `target_var` and `trait.positive` must be provided.\n")
    }
  } else if (task_type == "survival") {
    if (is.null(time_var) || is.null(event_var)) {
      stop("For survival, both `time_var` and `event_var` must be provided.\n")
    }

    if (!is.null(metric)) {
      warning("`metric` was set, but the task is survival. The metric argument will be ignored.\n")
    }

  } else {
    stop("Invalid task_type. Must be either 'classification' or 'survival'.\n")
  }

  # ---------------------------------------------------------------------------
  # === CASE 1: CLASSIFICATION TASK ==========================================
  # ---------------------------------------------------------------------------

  if (task_type == "classification") {

    #Set training set for classification
    train_data = features_train %>%
      data.frame() %>%
      dplyr::mutate(Trait = target_var,
                    target = as.factor(ifelse(Trait == trait.positive, 'yes', 'no'))) %>%
      dplyr::select(-Trait)

    train_data$target <- factor(train_data$target, levels = c("no", "yes"))  # Order, just in case to ensure positive class is not well defined

    if(LODO == T){
      train_data = train_data %>%
        dplyr::mutate(dataset = batch_var)
    }

    #Cross-validation training
    training = compute_k_fold_CV(train_data, k_folds = k_folds, n_rep = n_rep, metric = metric, stacking = stack,
                                 file_name = file_name, LODO = LODO, ncores = ncores, return= return,
                                 fold_construction_fun = fold_construction_fun, fold_construction_args_fixed = fold_construction_args_fixed,
                                 fold_construction_args_tunable = fold_construction_args_tunable)

  }

  # ---------------------------------------------------------------------------
  # === CASE 2: SURVIVAL ANALYSIS TASK =======================================
  # ---------------------------------------------------------------------------
  else if (task_type == "survival"){

    #Prepare training data
    train_data <- features_train %>%
      data.frame() %>%
      dplyr::mutate(
        time  = time_var, ## standarize names
        event = as.numeric(event_var) ## standarize names
      )

    if (LODO == TRUE) {
      train_data = train_data %>%
        dplyr::mutate(dataset = batch_var)
    }

    #Split into features + outcomes
    df_outcome  <- train_data %>% dplyr::select(time, event)
    df_features <- train_data %>% dplyr::select(-time, -event)

    #Run survival cross-validation and tuning
    training = compute_k_fold_CV_survival(
      df_features = df_features,
      df_outcome  = df_outcome,
      outcome_col = "time",
      event_col   = "event",
      k_folds = k_folds,
      n_rep = n_rep,
      ncores = ncores,
      file_name   = file_name,
      fold_construction_fun = fold_construction_fun,
      fold_construction_args_fixed = fold_construction_args_fixed,
      fold_construction_args_tunable = fold_construction_args_tunable
    )

  }

  ####################################################Predicting
  if(length(training)!=0){
    return(training)
  }else{  #No features are selected as predictive
    message("No features selected as predictive after Boruta runs. No model returned.")
    return(NULL)
  }

}

#' Train and evaluate machine learning models for classification or survival analysis
#'
#' This function trains and evaluates machine learning models using cross-validation on training data
#' and then tests performance on independent test data. It supports both **classification** and
#' **survival analysis** tasks, including hyperparameter tuning, model stacking, and cohort-based
#' (Leave-One-Dataset-Out, LODO) validation. For survival models, it computes the **C-index**
#' and generates Kaplan–Meier plots stratified by predicted risk.
#'
#' @param features_train A data frame or matrix of predictor variables used for training
#'   (rows as samples, columns as features).
#' @param features_test A data frame or matrix of predictor variables used for testing.
#' @param clinical A data frame containing clinical or outcome information.
#'   Row names must match those of \code{features_train} and \code{features_test}.
#' @param task_type Character. Type of task: either \code{"classification"} or \code{"survival"}.
#' @param trait Character. Name of the column in \code{clinical} used as the target variable
#'   (required for classification tasks).
#' @param trait.positive Value in \code{target_var} that represents the positive class (used for metrics like AUROC, AUPRC, confusion matrix, and SHAP values).
#'   This ensures that all performance metrics and interpretability analyses consistently treat the correct class as positive,
#' @param time_var Character. Name of the column in \code{clinical} containing survival or follow-up time
#'   (required for survival tasks).
#' @param event_var Character. Name of the column in \code{clinical} indicating event occurrence
#'   (1 = event occurred, 0 = censored; required for survival tasks).
#' @param metric Character. Performance metric for model tuning and selection.
#'   Supported options for classification: \code{"Accuracy"}, \code{"AUROC"}, \code{"AUPRC"}.
#'   For survival models, performance is evaluated using the concordance index (C-index).
#' @param stack Logical. Whether to perform model stacking (default = \code{FALSE}).
#' @param k_folds Integer. Number of folds for cross-validation (default = 10).
#' @param n_rep Integer. Number of repetitions for cross-validation (default = 5).
#' @param LODO Logical. If \code{TRUE}, performs Leave-One-Dataset-Out (LODO) cross-validation
#'   based on cohort identifiers.
#' @param batch_id Column name indicating where the cohort or batch membership for each sample is.
#'   Required if \code{LODO = TRUE}.
#' @param file_name Character. Base name used to save plots and results under the \code{Results/} directory.
#'   For survival tasks, this will be used to create a Kaplan–Meier plot named
#'   \code{"Results/Survival_KM_<file_name>.pdf"}.
#' @param ncores Integer. Number of CPU cores to use for parallelization. If not specified,
#'   defaults to \code{parallel::detectCores() - 1}.
#' @param maximize Character. Metric to maximize when selecting the optimal classification threshold
#'   (options: \code{"Accuracy"}, \code{"Precision"}, \code{"Recall"}, \code{"Specificity"},
#'   \code{"Sensitivity"}, \code{"F1"}, or \code{"MCC"}). Default = \code{"Accuracy"}.
#' @param fold_construction_fun Function. Optional custom function to construct cross-validation folds.
#'   Must accept a \code{bestune} argument internally for optimized hyperparameter injection.
#' @param fold_construction_args_fixed List. Fixed arguments passed to \code{fold_construction_fun},
#'   used in both cross-validation and final model training.
#' @param fold_construction_args_tunable List. Tunable arguments passed to \code{fold_construction_fun},
#'   defining hyperparameters to explore during cross-validation.
#' @param return Logical. Whether to return and save plots (default = \code{FALSE}).
#'
#' @details
#' For **classification tasks**, this function performs cross-validation tuning based on the chosen
#' performance metric (e.g., Accuracy, AUROC, or AUPRC), followed by test-set evaluation and ROC/PR curve plotting.
#'
#' For **survival analysis tasks**, it performs model selection using the C-index, refits the best
#' model on the full training data, evaluates the test-set C-index, and plots Kaplan–Meier survival
#' curves across quantile-based risk strata (Low/Medium/High risk). The C-index and log-rank test
#' p-value are displayed on the plot.
#'
#' @return A named list containing:
#' \describe{
#'   \item{Model}{The trained model or workflow (classification) or refitted best model (survival).}
#'   \item{Metrics}{Performance metrics computed on the test data.}
#'   \item{AUC}{For classification tasks, a list containing AUROC and AUPRC values.}
#'   \item{Prediction}{Predicted class probabilities (classification) or risk scores (survival).}
#'   \item{CV_Results}{Cross-validation results, including median and MAD of C-index (survival).}
#'   \item{Test_CINDEX}{Concordance index for the test data (survival only).}
#'   \item{KM_Plot}{Kaplan–Meier plot object (if \code{return = TRUE}).}
#' }
#'
#' @examples
#' \dontrun{
#' # --- Classification Example ---
#' results_classif <- compute_features.ML(
#'   features_train = X_train,
#'   features_test  = X_test,
#'   clinical       = clin_df,
#'   task_type      = "classification",
#'   trait          = "Response",
#'   trait.positive = "Responder",
#'   k_folds        = 5,
#'   n_rep          = 1,
#'   file_name      = "classification_example",
#'   return         = TRUE
#' )
#'
#' # --- Survival Example ---
#' results_surv <- compute_features.ML(
#'   features_train = X_train,
#'   features_test  = X_test,
#'   clinical       = clin_df,
#'   task_type      = "survival",
#'   time_var       = "time",
#'   event_var      = "status",
#'   k_folds        = 5,
#'   n_rep          = 1,
#'   file_name      = "cox_survival_example",
#'   return         = TRUE
#' )
#' }
#'
#' @export
compute_features.ML <- function(features_train, features_test, clinical,
                                task_type = c("classification", "survival"),
                                trait = NULL, trait.positive = NULL,
                                time_var = NULL, event_var = NULL,
                                metric = "Accuracy", stack = FALSE,
                                k_folds = 10, n_rep = 5, LODO = FALSE,
                                batch_id = NULL, file_name = NULL, ncores = NULL,
                                return = FALSE,
                                fold_construction_fun = NULL,
                                fold_construction_args_fixed = NULL,
                                fold_construction_args_tunable = NULL){

  # ---------------------------------------------------------------------------
  # === CASE 1: CLASSIFICATION TASK ==========================================
  # ---------------------------------------------------------------------------
  if (task_type == "classification") {

    # Train cohort
    traitData_train = clinical[rownames(clinical)%in%rownames(features_train), ]

    # Test cohort
    traitData_test = clinical[rownames(clinical)%in%rownames(features_test), ]

    ####################################################Training

    #Set training set
    train_data = features_train %>%
      data.frame() %>%
      dplyr::mutate(Trait = traitData_train[,trait],
                    target = as.factor(ifelse(Trait == trait.positive, 'yes', 'no'))) %>%
      dplyr::select(-Trait)

    train_data$target <- factor(train_data$target, levels = c("no", "yes"))  # Order, just in case to ensure positive class is not well defined

    if(LODO == T){
      train_data = train_data %>%
        dplyr::mutate(dataset = traitData_train[,batch_id])
    }

    #Cross-validation training
    training = compute_k_fold_CV(train_data, k_folds = k_folds, n_rep = n_rep, metric = metric, stacking = stack,
                                 file_name = file_name, LODO = LODO, ncores = ncores, return= return,
                                 fold_construction_fun = fold_construction_fun, fold_construction_args_fixed = fold_construction_args_fixed,
                                 fold_construction_args_tunable = fold_construction_args_tunable)

    ####################################################Predicting
    if(length(training)!=0){
      ####################### Testing set

      if(stack){
        prediction = compute_prediction(training, features_test, traitData_test[,trait], trait.positive, stack = TRUE, file.name = file_name, return = return)
      }else{
        prediction = compute_prediction(training, features_test, traitData_test[,trait], trait.positive, stack = FALSE, file.name = file_name, return = return)
      }

      auc_roc_score = prediction[["AUC"]][["AUROC"]]
      auc_prc_score = prediction[["AUC"]][["AUPRC"]]

      metrics = prediction[["Metrics"]]
      predictions = prediction[["Predictions"]]

      return(list(Model = training, Metrics = metrics, AUC = list(AUROC = auc_roc_score, AUPRC = auc_prc_score), Prediction = predictions))
    }else{  #No features are selected as predictive

      message("No features selected as predictive after Boruta runs. No model returned.")

      return(NULL)
    }

  }

  # ---------------------------------------------------------------------------
  # === CASE 2: SURVIVAL ANALYSIS TASK =======================================
  # ---------------------------------------------------------------------------
  if (task_type == "survival") {

    # ---------------------------- Data prep ----------------------------------
    train_data <- features_train %>%
      as.data.frame() %>%
      dplyr::mutate(
        time = clinical[rownames(features_train), time_var],
        event = clinical[rownames(features_train), event_var]
      )

    test_data <- features_test %>%
      as.data.frame() %>%
      dplyr::mutate(
        time = clinical[rownames(features_test), time_var],
        event = clinical[rownames(features_test), event_var]
      )

    if (LODO) {
      train_data = train_data %>%
        dplyr::mutate(dataset = clinical[,batch_id])
    }

    df_outcome  <- train_data %>% dplyr::select(time, event)
    df_features <- train_data %>% dplyr::select(-time, -event)

    # ---------------------------- Cross-validation ---------------------------
    training <- compute_k_fold_CV_survival(
      df_features = df_features,
      df_outcome  = df_outcome,
      outcome_col = "time",
      event_col   = "event",
      k_folds = k_folds,
      n_rep = n_rep,
      ncores = ncores,
      file_name   = file_name
    )

    # ---------------------------- Refit best model ---------------------------

    ### Prediction

    outcome_test  <- test_data %>% dplyr::select(time, event)
    df_test <- test_data %>% dplyr::select(-time, -event)

    pred <- compute_prediction(model = training$Model,
                               test_data = df_test,
                               task_type = task_type,
                               time_var = outcome_test$time,
                               event_var = outcome_test$event,
                               file.name = file_name,
                               return = return)

    return(list(Model = training, C_index = pred$C_index, Prediction = pred$Predictions))

  }

}

#' Plot Pooled AUROC and AUPRC Performance Curves
#'
#' This function reads multiple `.rds` files containing machine learning results, pools the AUROC and AUPRC metrics,
#' and generates boxplots summarizing performance across iterations. Median values are annotated on the plots.
#'
#' @param file.name Character. Name to use when saving the plots (used as a prefix in the output file names).
#' @param folder_path Character. Path to the directory containing the `.rds` files with ML model results.
#'
#' @return
#' Saves two PDF files in the \code{Results/} directory:
#' \itemize{
#'   \item Boxplot of AUROC values with median annotation
#'   \item Boxplot of AUPRC values with median annotation
#' }
#' No value is returned to the R environment.
#'
#' @details
#' Each `.rds` file is expected to contain a list with a \code{result$AUC} element that includes
#' both \code{AUROC} and \code{AUPRC} values.
#'
#'
get_pooled_roc_curves = function(file.name, folder_path){

  # Get a list of all RDS files in the folder
  res <- list.files(folder_path, pattern = "\\.rds$", full.names = TRUE)

  # Initialize cumulative data frame
  cumulative_data <- data.frame(AUC_roc = numeric(),
                                AUC_prc = numeric(),
                                Cohort = character(),
                                stringsAsFactors = FALSE)

  for (file in res) {
    model <- readRDS(file)

    auc_roc <- model[["result"]][["AUC"]][["AUROC"]]
    auc_prc <- model[["result"]][["AUC"]][["AUPRC"]]

    # Append metrics to cumulative data
    cumulative_data <- rbind(cumulative_data,
                             data.frame(AUC_roc = auc_roc,
                                        AUC_prc = auc_prc,
                                        Cohort = file.name))
  }

  #########Boxplot
  iterations = nrow(cumulative_data)

  median_auc_roc = cumulative_data %>%
    dplyr::group_by(Cohort) %>%
    dplyr::summarize(medianAUROC = median(AUC_roc))

  median_auc_prc = cumulative_data %>%
    dplyr::group_by(Cohort) %>%
    dplyr::summarize(medianAUPRC = median(AUC_prc))

  # Plot boxplot with median AUC annotations
  plot_roc = ggplot2::ggplot(cumulative_data, ggplot2::aes(x = Cohort, y = AUC_roc, fill = Cohort)) +
    ggplot2::geom_boxplot() +
    ggplot2::labs(title = paste0("Distribution of AUROC values across ", iterations, " splits"),
         x = "Model",
         y = "AUROC") +
    ggplot2::theme_minimal() +
    ggplot2::theme(legend.position = "right") +
    ggplot2::geom_text(data = median_auc_roc, ggplot2::aes(x = Cohort, y = max(cumulative_data$AUC_roc),
                                                           label = paste("Median AUROC:", round(medianAUROC, 3))),
                       size = 4, color = "black", vjust = -0.5)

  grDevices::pdf(paste0("Results/Boxplot_AUROC_performance_", file.name, ".pdf"))
  print(plot_roc)
  grDevices::dev.off()

  plot_prc = ggplot2::ggplot(cumulative_data, ggplot2::aes(x = Cohort, y = AUC_prc, fill = Cohort)) +
    ggplot2::geom_boxplot() +
    ggplot2::labs(title = paste0("Distribution of AUPRC values across ", iterations, " splits"),
         x = "Model",
         y = "AUPRC") +
    ggplot2::theme_minimal() +
    ggplot2::theme(legend.position = "right") +
    ggplot2::geom_text(data = median_auc_prc, ggplot2::aes(x = Cohort, y = max(cumulative_data$AUC_prc),
                                                           label = paste("Median AUPRC:", round(medianAUPRC, 3))),
              size = 4, color = "black", vjust = -0.5)

  grDevices::pdf(paste0("Results/Boxplot_AUPRC_performance_", file.name, ".pdf"))
  print(plot_prc)
  grDevices::dev.off()

}

#' Plot Pooled AUROC and AUPRC Boxplots Across Multiple Folders
#'
#' This function aggregates AUROC and AUPRC metrics from multiple folders (typically corresponding to different cohorts or models),
#' and generates boxplots comparing model performance across groups.
#'
#' @param folder_paths Character vector. Paths to folders containing `.rds` files with ML model performance results.
#' @param file_name Character. Prefix for the saved PDF files containing the plots.
#' @param width Numeric. Width of the saved plots in inches. Default is 12.
#' @param height Numeric. Height of the saved plots in inches. Default is 8.
#'
#' @return
#' Saves two PDF files in the \code{Results/} directory:
#' \itemize{
#'   \item \code{Boxplots_AUROC_performance_<file_name>.pdf}
#'   \item \code{Boxplots_AUPRC_performance_<file_name>.pdf}
#' }
#' No object is returned to the R environment.
#'
#' @details
#' Each `.rds` file is expected to contain a list object with a \code{result$AUC} element,
#' which includes numeric values for both \code{AUROC} and \code{AUPRC}.
#' Folder names are used as grouping labels in the plots.
#'
#' Red dashed lines are drawn at a fixed reference value (e.g., 0.7) for visual interpretation.
#'
#'
get_pooled_boxplots = function(folder_paths, file_name, width = 12, height = 8) {

  # Initialize cumulative data frame
  cumulative_data <- data.frame(AUC_roc = numeric(),
                                AUC_prc = numeric(),
                                Cohort = character(),
                                Folder = character(),
                                stringsAsFactors = FALSE)

  for (folder_path in folder_paths) {

    # Extract folder name for labeling
    folder_name <- basename(folder_path)

    # Get a list of all RDS files in the folder
    res <- list.files(folder_path, pattern = "\\.rds$", full.names = TRUE)

    for (file in res) {
      model <- readRDS(file)

      auc_roc <- model[["result"]][["AUC"]][["AUROC"]]
      auc_prc <- model[["result"]][["AUC"]][["AUPRC"]]

      # Append metrics to cumulative data
      cumulative_data <- rbind(cumulative_data,
                               data.frame(AUC_roc = auc_roc,
                                          AUC_prc = auc_prc,
                                          Cohort = basename(file),
                                          Folder = folder_name))
    }
  }

  ######### AUROC Boxplots #########
  grDevices::pdf(paste0("Results/Boxplots_AUROC_performance_", file_name, ".pdf"), width = width, height = height)

  folder_data <- cumulative_data %>% dplyr::filter(Folder == folder_name)
  iterations <- nrow(folder_data)

  plot_roc <- ggplot2::ggplot(cumulative_data, ggplot2::aes(x = Folder, y = AUC_roc, fill = Folder)) +
    ggplot2::geom_boxplot() +
    ggplot2::geom_hline(yintercept = 0.7, linetype = "dashed", color = "red", linewidth = 1) +  # Red horizontal line
    ggplot2::coord_cartesian(ylim = c(0.2, 0.9)) +  # Set y-axis limits
    ggplot2::labs(title = paste0("ML models using TME features across ", iterations, " iterations"),
                  x = "Features",
                  y = "AUROC") +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      text = element_text(size = 16),       # Increase overall text size
      axis.text = element_text(size = 14),  # Increase axis tick labels
      axis.text.x = element_text(angle = 45, hjust = 1),  # Rotate x-axis labels diagonally
      axis.title = element_text(size = 16), # Increase axis titles
      plot.title = element_text(size = 18, face = "bold"), # Increase title size
      legend.text = element_text(size = 14),
      legend.title = element_text(size = 16)
    )

  print(plot_roc)

  grDevices::dev.off()

  ######### AUPRC Boxplots #########
  grDevices::pdf(paste0("Results/Boxplots_AUPRC_performance_", file_name, ".pdf"), width = width, height = height)

  plot_prc <- ggplot2::ggplot(cumulative_data, ggplot2::aes(x = Folder, y = AUC_prc, fill = Folder)) +
    ggplot2::geom_boxplot() +
    ggplot2::geom_hline(yintercept = 0.7, linetype = "dashed", color = "red", linewidth = 1) +  # Red horizontal line
    ggplot2::coord_cartesian(ylim = c(0, 1)) +  # Set y-axis limits
    ggplot2::labs(title = paste0("LODO analysis - AUPRC Distribution (", iterations, " splits)"),
         x = "Model",
         y = "AUPRC") +
    ggplot2::theme_minimal()

  print(plot_prc)

  grDevices::dev.off()

}

#' Compute Cross-Validation Accuracy for ML Models
#'
#' This function extracts cross-validated accuracy values from a list of trained machine learning models,
#' summarizes their median and standard deviation, and optionally plots a bar chart or selects base models for stacking.
#'
#' @param models A named list of trained ML models, each with a \code{resample} element containing cross-validated accuracy.
#' @param file_name (Optional) Character string specifying the filename prefix for the saved accuracy plot (PDF format).
#' @param base_models Logical. If \code{TRUE}, the function selects and returns base models using \code{choose_base_models()} for stacking.
#' @param return Logical. If \code{TRUE}, the function saves a barplot of the model accuracy values in the \code{Results/} directory.
#'
#' @return A list containing:
#' \itemize{
#'   \item \code{Accuracy}: A data frame with the median and standard deviation of accuracy for each model.
#'   \item \code{Top_model}: A character string naming the model with the highest median accuracy.
#'   \item \code{Base_models} (optional): A character vector of selected base models if \code{base_models = TRUE}.
#' }
#'
#' @details
#' This function assumes that each model in the list has a \code{$resample} component containing
#' a column named \code{Accuracy}. It calculates the median and standard deviation of accuracy
#' for each model and creates a barplot (if \code{return = TRUE}) with error bars.
#'
#' If \code{base_models = TRUE}, it calls a helper function \code{choose_base_models()} to select
#' models for use in stacking.
#'
#'
compute_cv_accuracy = function(models, file_name = NULL, base_models = FALSE, return = TRUE){

  # Bind accuracy values from each model
  accuracy = list()
  for (i in 1:length(models)) {
    accuracy[[i]] = models[[i]]$resample %>%
      dplyr::mutate(model = names(models)[i])
    names(accuracy)[i] = names(models)[i]
  }
  accuracy_data = do.call(rbind, accuracy)

  # Retrieve top model based on accuracy
  res_accuracy <- accuracy_data %>%
    dplyr::group_by(model) %>%
    dplyr::summarise(
      Mean_Accuracy = median(.data$Accuracy),
      MAD_Accuracy = stats::mad(.data$Accuracy)
    ) %>%
    dplyr::arrange(desc(Mean_Accuracy))

  top_model = res_accuracy %>%
    dplyr::slice(1) %>%
    dplyr::pull(model)

  if(return){
    grDevices::pdf(paste0("Results/Accuracy_CV_methods_", file_name, ".pdf"), width = 10)

    plot(ggplot2::ggplot(res_accuracy, ggplot2::aes(x = model, y = Mean_Accuracy)) +
           ggplot2::geom_bar(stat = "identity", position = ggplot2::position_dodge(), width = 0.6, fill = "#1f78b4") + # professional blue
           ggplot2::geom_errorbar(ggplot2::aes(ymin = Mean_Accuracy - MAD_Accuracy, ymax = Mean_Accuracy + MAD_Accuracy),
                                  width = 0.2, position = position_dodge(0.6)) +
           ggplot2::labs(title = "Performance of Models",
                         x = "Model",
                         y = "Median Accuracy") +
           ggplot2::theme_minimal() +
           ggplot2::theme(legend.position = "none") +
           ggplot2::scale_y_continuous(breaks = seq(0, 1, by = 0.05)))

    grDevices::dev.off()
  }

  if(base_models == T){
    cat("Choosing base models for stacking.......................................\n\n")
    base_models = choose_base_models(models, metric = "Accuracy")
    cat("Models chosen are:", paste0(base_models, collapse = ", "), "\n\n")
    return(list("Accuracy" = res_accuracy, "Top_model" = top_model, "Base_models" = base_models))
  }else{
    return(list("Accuracy" = res_accuracy, "Top_model" = top_model))
  }


}

#' Compute Cross-Validated AUC Values for Machine Learning Models
#'
#' This function computes cross-validated AUROC and AUPRC scores for a list of trained machine learning models.
#' It can also save performance barplots and optionally select base models for stacking.
#'
#' @param models A named list of trained machine learning models. Each model should contain a \code{resample}
#'   data frame with AUROC and AUPRC values from cross-validation.
#' @param file_name (Optional) Character string. Used as the prefix for the plot filenames if \code{save_plot = TRUE}.
#' @param base_models Logical. If \code{TRUE}, selects a subset of models as base learners for stacking using the \code{choose_base_models()} function.
#' @param AUC_type Character. Either \code{"AUROC"} or \code{"AUPRC"}; determines which metric is used to select the top-performing model.
#' @param return Logical. Whether to return the results and generated plots.
#'
#' @return A list containing:
#' \describe{
#'   \item{\code{AUROC}}{A data frame with median and standard deviation of AUROC values for each model.}
#'   \item{\code{AUPRC}}{A data frame with median and standard deviation of AUPRC values for each model.}
#'   \item{\code{Top_model}}{The name of the model with the highest median value for the selected metric (\code{AUC_type}).}
#'   \item{\code{Base_models}}{(Optional) A character vector of selected base models for stacking, returned if \code{base_models = TRUE}.}
#' }
#'
#'
#' @examples
#' \dontrun{
#' res <- compute_cv_AUC(
#'   models = ml_models,
#'   file_name = "Model_Performance",
#'   base_models = TRUE,
#'   AUC_type = "AUROC",
#'   save_plot = TRUE
#' )
#' }
compute_cv_AUC = function(models, file_name = NULL, base_models = FALSE, AUC_type = "AUROC", return = TRUE){

  if(!(AUC_type %in% c("AUROC", "AUPRC"))){
    stop("AUC type provided don't correspond neither to ROC or PRC")
  }

  #Bind AUROC values from each model
  auroc = list()
  for (i in 1:length(models)){
    auroc[[i]] = models[[i]]$resample %>% #we use the resample matrix and not directly the results matrix as some have hyperparameters so we will need to define best on the tuned parameter (=more code) - resample matrix is made based on the best tuning
      dplyr::mutate(model = names(models)[i])
    names(auroc)[i] = names(models)[i]
  }
  auroc_data = do.call(rbind, unname(auroc))

  #Bind AUPRC values from each model
  auprc = list()
  for (i in 1:length(models)){
    auprc[[i]] = models[[i]]$resample %>% #we use the resample matrix and not directly the results matrix as some have hyperparameters so we will need to define best on the tuned parameter (=more code) - resample matrix is made based on the best tuning
      dplyr::mutate(model = names(models)[i])
    names(auprc)[i] = names(models)[i]
  }
  auprc_data = do.call(rbind, auprc)

  ### Calculate AUROC and AUPRC across resamples for each model
  res_auroc <- auroc_data %>%
    dplyr::group_by(model) %>%
    dplyr::summarise(
      Mean_AUROC = median(.data$AUROC),
      MAD_AUROC = stats::mad(.data$AUROC)
    ) %>%
    dplyr::arrange(desc(Mean_AUROC))

  res_auprc <- auprc_data %>%
    dplyr::group_by(model) %>%
    dplyr::summarise(
      Mean_AUPRC = median(.data$AUPRC),
      MAD_AUPRC = stats::mad(.data$AUPRC)
    ) %>%
    dplyr::arrange(desc(Mean_AUPRC))


  if(return){
    grDevices::pdf(paste0("Results/AUROC_CV_methods_", file_name, ".pdf"), width = 10)
    plot(ggplot2::ggplot(res_auroc, ggplot2::aes(x = model, y = Mean_AUROC)) +
           ggplot2::geom_bar(stat = "identity", position = ggplot2::position_dodge(), width = 0.6, fill = "#1f78b4") + # professional blue
           ggplot2::geom_errorbar(aes(ymin = Mean_AUROC - MAD_AUROC, ymax = Mean_AUROC + MAD_AUROC),
                                  width = 0.2, position = position_dodge(0.6)) +
           ggplot2::labs(title = "Performance of Models",
                         x = "ML Model",
                         y = "Median AUROC across resamples") +
           ggplot2::theme_minimal() +
           ggplot2::theme(legend.position = "none") +
           ggplot2::scale_y_continuous(breaks = seq(0, 1, by = 0.05)))
    grDevices::dev.off()

    grDevices::pdf(paste0("Results/AUPRC_CV_methods_", file_name, ".pdf"), width = 10)
    plot(ggplot2::ggplot(res_auprc, ggplot2::aes(x = model, y = Mean_AUPRC)) +
           ggplot2::geom_bar(stat = "identity", position = ggplot2::position_dodge(), width = 0.6, fill = "#1f78b4") + # same professional blue
           ggplot2::geom_errorbar(aes(ymin = Mean_AUPRC - MAD_AUPRC, ymax = Mean_AUPRC + MAD_AUPRC),
                                  width = 0.2, position = ggplot2::position_dodge(0.6)) +
           ggplot2::labs(title = "Performance of Models",
                         x = "ML Model",
                         y = "Median AUPRC across resamples") +
           ggplot2::theme_minimal() +
           ggplot2::theme(legend.position = "none") +
           ggplot2::scale_y_continuous(breaks = seq(0, 1, by = 0.05)))
    grDevices::dev.off()

  }

  #Retrieve top model based on AUROC or AUPRC
  if(AUC_type == "AUROC"){
    top = res_auroc
  }else{
    top = res_auprc
  }

  top_model = top %>%
    dplyr::slice(1) %>%
    dplyr::pull(model)

  if(base_models == TRUE){ ### For stacking
    cat("Choosing base models for stacking.......................................\n\n")
    base_models = choose_base_models(models, metric = AUC_type)
    cat("Models chosen are:", paste0(base_models, collapse = ", "), "\n\n")
    return(list("AUROC" = res_auroc, "AUPRC" = res_auprc, "Top_model" = top_model, "Base_models" = base_models))
  }else{
    return(list("AUROC" = res_auroc, "AUPRC" = res_auprc, "Top_model" = top_model))
  }

}

#' Choose Top Base Models for Stacking Based on Accuracy or AUC Scores
#'
#' This function selects three base models for stacking based on either Accuracy or AUC metrics. It chooses the top models from
#' different categories (e.g., tree-based methods, linear models, instance-based methods) according to the specified metric.
#'
#' @param models A list of trained machine learning models. Each model must contain a \code{resample} data frame with
#'   performance metrics (Accuracy, AUROC, AUPRC) from cross-validation.
#' @param metric A character string specifying the metric to use for model selection. Can be either "Accuracy", "AUROC", or "AUPRC".
#'   Default is "Accuracy".
#'
#' @return A character vector containing the names of the top models selected based on the specified metric.
#'
#'
#' @examples
#' \dontrun{
#' base_models = choose_base_models(models = ml_models, metric = "AUROC")
#' }
choose_base_models = function(models, metric = "Accuracy"){

  #Bind accuracy values from each model
  resample_df = list()
  for (i in 1:length(models)){
    resample_df[[i]] = models[[i]]$resample %>%
      dplyr::mutate(model = names(models)[i])
    names(resample_df)[i] = names(models)[i]
  }
  resample_df = do.call(rbind, resample_df)

  if(metric == "Accuracy"){
    #Prepare data frame for ploting
    resample_df <- resample_df %>%
      dplyr::group_by(model) %>%
      dplyr::summarise(Accuracy = median(Accuracy))
  }else if(metric == "AUROC"){
    #Prepare data frame for ploting
    resample_df <- resample_df %>%
      dplyr::group_by(model) %>%
      dplyr::summarise(AUROC = median(AUROC))
  }else if(metric == "AUPRC"){
    #Prepare data frame for ploting
    resample_df <- resample_df %>%
      dplyr::group_by(model) %>%
      dplyr::summarise(AUPRC = median(AUPRC))
  }

  resample_df <- resample_df %>%
    dplyr::mutate(Category = dplyr::case_when(
      model %in% c("BAG", "C50", "CART", "RF", "XGboost") ~ "Tree-based Methods",
      model %in% c("GLM", "LDA", "GLMNET", "LASSO", "RIDGE") ~ "Linear Models",
      model %in% c("KNN", "SVM_linear", "SVM_radial") ~ "Instance-based Methods",
      TRUE ~ "Other"  # In case there are models not in the above lists
    ))

  if(metric == "Accuracy"){
    groupped_df <- resample_df %>%
      dplyr::group_by(Category) %>%
      dplyr::filter(Accuracy == max(Accuracy)) %>%
      dplyr::ungroup()
  }else if(metric == "AUROC"){
    groupped_df <- resample_df %>%
      dplyr::group_by(Category) %>%
      dplyr::filter(AUROC == max(AUROC)) %>%
      dplyr::ungroup()
  }else if(metric == "AUPRC"){
    groupped_df <- resample_df %>%
      dplyr::group_by(Category) %>%
      dplyr::filter(AUPRC == max(AUPRC)) %>%
      dplyr::ungroup()
  }else{
    stop("No metric defined")
  }

  #Retrieve top model based on accuracy/auc
  base_models <- groupped_df %>%
    dplyr::pull(model)

  return(base_models)
}

calculate_auc_roc_resample = function(obs, pred){

  prob_obs = data.frame("yes" = pred, "obs" = obs)

  prob_obs = prob_obs %>%
    dplyr::arrange(dplyr::desc(pred)) %>% #need to be arrange for apply cumulative sum
    dplyr::mutate(is_yes = (obs == "yes"),
            tp = cumsum(is_yes), #true positive above the threshold - cumulative sum to refer to the threshold
            fp = cumsum(!is_yes), #false positive above the threshold - cumulative sum to refer to the threshold
            fpr = fp/sum(obs == 'no'),
            tpr = tp/sum(obs == 'yes'))

  auc_value = calculate_auroc(prob_obs$fpr, prob_obs$tpr)

  return(auc_value)
}

calculate_auc_prc_resample = function(obs, pred) {

  prob_obs = data.frame("yes" = pred, "obs" = obs)

  prob_obs = prob_obs %>%
    dplyr::arrange(dplyr::desc(pred)) %>% # Sort by predicted probability
    dplyr::mutate(is_yes = (obs == "yes"),
                  tp = cumsum(is_yes), # True positives
                  fp = cumsum(!is_yes), # False positives
                  fn = sum(obs == "yes") - tp, # False negatives
                  precision = tp / (tp + fp), # Precision
                  recall = tp / (tp + fn)) # Recall

  # Now use your existing calculate_auprc function to calculate AUC-PRC
  auc_prc_value = calculate_auprc(prob_obs$recall, prob_obs$precision)

  return(auc_prc_value)
}

calculate_accuracy_kappa_resample <- function(obs, pred) {
  # Ensure input is factor with the same levels
  obs <- factor(obs)
  pred <- factor(pred, levels = levels(obs))

  # Confusion matrix components
  cm <- table(pred, obs)

  # Total observations
  total <- sum(cm)

  # Accuracy = (TP + TN) / Total
  accuracy <- sum(diag(cm)) / total

  # Row and column totals
  row_totals <- rowSums(cm)
  col_totals <- colSums(cm)

  # Expected accuracy (by chance)
  expected_accuracy <- sum(row_totals * col_totals) / (total ^ 2)

  # Kappa = (observed - expected) / (1 - expected)
  kappa <- (accuracy - expected_accuracy) / (1 - expected_accuracy)

  return(c(Accuracy_resample = accuracy, Kappa_resample = kappa))
}

#' Compute F1 Score
#'
#' The F1 score is the harmonic mean of precision and recall, and is used to evaluate the balance between the two metrics.
#' It is particularly useful when the class distribution is imbalanced.
#'
#' @param metrics A vector of predicted class labels or probabilities.
#' @param target A vector of true class labels.
#'
#' @return The F1 score, a numeric value between 0 and 1.
#' @export
#'
calculate_f1 = function(metrics, target) {

  confusion_values <- calculate_confusion_values(metrics, target)
  TP <- confusion_values$TP
  FN <- confusion_values$FN
  FP <- confusion_values$FP

  precision <- TP / (TP + FP)
  recall <- TP / (TP + FN)

  f1_score = ifelse((precision + recall) > 0, 2 * (precision * recall) / (precision + recall), 0) # F1 Score

  return(f1_score)
}

#' Compute Matthews Correlation Coefficient (MCC) Score
#'
#' The Matthews correlation coefficient (MCC) is a metric used for binary classification problems.
#' It takes into account true and false positives and negatives, and is considered a balanced metric.
#'
#' @param metrics A vector of predicted class labels or probabilities.
#' @param target A vector of true class labels.
#'
#' @return The MCC score, a numeric value between -1 and 1.
#' @export
#'
calculate_mcc = function(metrics, target) {

  confusion_values <- calculate_confusion_values(metrics, target)
  TP <- confusion_values$TP
  FN <- confusion_values$FN
  FP <- confusion_values$FP
  TN <- confusion_values$TN

  numerator = (TP * TN) - (FP * FN)
  denominator = sqrt((TP + FP) * (TP + FN) * (TN + FP) * (TN + FN))

  mcc = ifelse(denominator > 0, numerator / denominator, 0) # MCC Calculation

  return(mcc)
}

#' Calculate Sensitivity and Specificity Values
#'
#' This function calculates sensitivity (recall), specificity, and other related metrics (accuracy, precision, recall, F1 score, MCC)
#' from predicted and true class labels.
#'
#' @param predictions A vector of predicted class labels or probabilities.
#' @param observed A vector of true class labels.
#' @param ml.model The trained machine learning model used to generate predictions.
#'
#' @return A data frame containing sensitivity, specificity, precision, recall, F1 score, MCC, and other metrics.
#' @export
#'
get_sensitivity_specificity = function(predictions, observed, ml.model){
  prob_obs = dplyr::bind_cols(predictions, observed = observed)

  prob_obs = prob_obs %>%
    dplyr::arrange(dplyr::desc(yes)) %>% #need to be arrange for apply cumulative sum
    dplyr::mutate(is_yes = (observed == "yes"),
                  tp = cumsum(is_yes), #true positive above the threshold - cumulative sum to refer to the threshold
                  fp = cumsum(!is_yes), #false positive above the threshold - cumulative sum to refer to the threshold
                  Sensitivity = tp/sum(observed == 'yes'),
                  fpr = fp/sum(observed == 'no'),
                  Specificity = 1 - fpr) %>%
    dplyr::select(yes, Sensitivity, Specificity, fpr) %>%
    dplyr::mutate(model = ml.model)

  prob_obs = prob_obs %>%
    dplyr::mutate(Accuracy = calculate_accuracy(., observed),
                  Precision = calculate_precision(., observed),
                  Recall = calculate_recall(., observed),
                  F1 = calculate_f1(., observed),
                  MCC = calculate_mcc(., observed)) %>%
    dplyr::select(yes, model, dplyr::everything())

  return(prob_obs)

}

#' Calculate AUC from ROC Curve
#'
#' This function calculates the Area Under the Receiver Operating Characteristic (ROC) curve.
#' It uses the trapezoidal rule to compute the AUC from the false positive rate (FPR) and sensitivity (true positive rate).
#'
#' @param fpr A numeric vector of false positive rates from the ROC curve.
#' @param sensitivity A numeric vector of sensitivities (true positive rates) from the ROC curve.
#'
#' @return The AUC score, a numeric value between 0 and 1.
#' @export
#'
calculate_auroc <- function(fpr, sensitivity) {
  #tpr is sensitivity

  # Sort by FPR to ensure trapezoidal rule is correctly applied (Already ordered)
  ordered <- order(fpr)
  fpr <- fpr[ordered]
  sensitivity <- sensitivity[ordered]

  auc <- 0
  for (i in 1:(length(fpr) - 1)) { #-1 to avoid NA cause last terms are TPR = 1 and FPR = 1
    # Trapezoidal rule: (TPR_i + TPR_{i+1}) / 2 * (FPR_{i+1} - FPR_i)
    auc <- auc + ((sensitivity[i+1] + sensitivity[i]) / 2) * (fpr[i+1] - fpr[i])
  }
  return(auc)

}

#' Calculate AUC from Precision-Recall Curve
#'
#' This function calculates the Area Under the Precision-Recall Curve (AUPRC).
#' It uses the trapezoidal rule to compute the AUPRC from the recall and precision values.
#'
#' @param recall A numeric vector of recall values (sensitivity) from the precision-recall curve.
#' @param precision A numeric vector of precision values from the precision-recall curve.
#'
#' @return The AUPRC score, a numeric value between 0 and 1.
#' @export
#'
calculate_auprc <- function(recall, precision) {
  # Sort by Recall to ensure trapezoidal rule is correctly applied
  ordered <- order(recall)
  recall <- recall[ordered]
  precision <- precision[ordered]

  auprc <- 0
  for (i in 1:(length(recall) - 1)) { # -1 to avoid NA from the last terms
    # Trapezoidal rule: (precision[i] + precision[i+1]) / 2 * (recall[i+1] - recall[i])
    auprc <- auprc + ((precision[i+1] + precision[i]) / 2) * (recall[i+1] - recall[i])
  }
  return(auprc)
}

#' Compute Weighted Feature Importance from Base Models and Meta-Learner for Stacking Models
#'
#' This function computes the feature importance by weighing the feature importances from
#' multiple base models in a stacking ensemble, combined with the meta-learner's model importance.
#' The final importance score for each feature is calculated by multiplying the base model's
#' feature importance with the meta-learner's weight for each base model.
#'
#' @param base_importance A list where each element corresponds to a base model and contains a data frame
#'                        with feature importances. Each data frame should have a column called `importance`
#'                        (either for the positive class or overall, depending on the type of model).
#' @param base_models A character vector with the names of the base models whose feature importances are
#'                    provided in `base_importance`.
#' @param meta_learner A `caret` object representing the trained meta-learner model. This model is used to
#'                     obtain weights for each base model, based on their performance in the ensemble.
#'
#' @return A data frame with two columns:
#'   \item{features}{The feature names.}
#'   \item{final_importance}{The final weighted importance score for each feature, calculated by summing
#'                           the weighted importances across all base models. Features are sorted in descending
#'                           order of their final importance score.}
#'
#' @details
#' The function extracts feature importance values from the base models, then computes the weighted importance
#' for each feature based on the meta-learner's performance. The meta-learner's feature importance is normalized
#' to ensure the sum of the importances across all models is 1, and it is used as the weight for each base model.
#' The feature importances from all base models are then aggregated and weighted by their respective meta-learner
#' importance scores.
#'
#' @seealso \code{\link[caret]{varImp}}
#'
#' @import caret
#' @import dplyr
#' @import tibble
#' @export
calculate_feature_importance_stacking = function(base_importance, base_models, meta_learner){

  #Extract features importance values within each base model for the meta-learner
  base_importance_list = list()
  for (i in 1:length(base_models)) {
    check = ncol(base_importance[[base_models[i]]][["importance"]])
    if(check > 1){ #Means importance is given for each class
      base_importance_list[[i]] = base_importance[[base_models[i]]][["importance"]] %>%
        tibble::rownames_to_column("features") %>%
        dplyr::select(features, yes) %>% #Take only importance for positive class
        dplyr::rename(importance = yes)
    }else{
      base_importance_list[[i]] = base_importance[[base_models[i]]][["importance"]] %>%
        tibble::rownames_to_column("features") %>%
        dplyr::rename(importance = Overall)
    }
    names(base_importance_list)[i] = base_models[i]
  }

  #Combine all base model importances in one data frame and add the model name
  combined_importance <- dplyr::bind_rows(
    lapply(names(base_importance_list), function(model) {
      base_importance_list[[model]] %>%
        data.frame() %>%
        dplyr::mutate(model = model)
    })
  )

  #Calculate base-models importance for the meta-learner
  meta_importance = caret::varImp(meta_learner, scale = FALSE)$importance %>%
    tibble::rownames_to_column("model")

  #Normalize the meta-learner's importance scores so they sum to 1
  meta_importance$Overall <- meta_importance$Overall / sum(meta_importance$Overall)

  #Combine features importance within base models with the overall importance for meta-learner
  combined_importance <- combined_importance %>%
    dplyr::left_join(meta_importance, by = "model") %>%
    dplyr::mutate(weighted_importance = importance * Overall) # importance is from base, Overall is from meta

  #Sum the weighted importance by feature across all models
  final_importance <- combined_importance %>%
    dplyr::group_by(features) %>%
    dplyr::summarise(final_importance = sum(weighted_importance, na.rm = TRUE)) %>%
    dplyr::arrange(desc(final_importance))

  return(final_importance)

}

#' Compute Prediction Metrics for a Trained Machine Learning Model
#'
#' This function computes prediction metrics for a given machine learning model, including
#' the confusion matrix, AUROC, AUPRC, and other performance metrics such as Accuracy, Sensitivity,
#' Specificity, Precision, Recall, F1 score, and MCC. The function also determines the optimal
#' classification threshold based on a chosen metric (e.g., Accuracy, F1, or AUROC) and generates
#' a confusion matrix plot.
#'
#' @param model The trained machine learning model returned from `compute_features.training.ML()` or `compute_features.ML()`.
#' @param test_data A matrix or data frame containing the testing dataset (features only).
#' @param target_var A character vector of true target values for the test data (the observed labels).
#' @param trait.positive Value in \code{target_var} that represents the positive class (used for metrics like AUROC, AUPRC, confusion matrix, and SHAP values).
#'   This ensures that all performance metrics and interpretability analyses consistently treat the correct class as positive
#' @param stack Logical. If stacking was used during model training, this parameter should be set to TRUE in order to use the meta-learner for prediction. Default is FALSE.
#' @param file.name A character string to specify the filename for saving the confusion matrix plot
#'                  (optional). If `NULL`, the plot is not saved.
#' @param maximize A character string indicating which metric to maximize when selecting the best
#'                 threshold for the confusion matrix. Options include "Accuracy", "Precision",
#'                 "Recall", "Specificity", "Sensitivity", "F1", or "MCC". Default is "Accuracy".
#' @param return Logical. Whether to return the results and generated plots.
#'
#' @return A list containing:
#' \item{Metrics}{A data frame with various performance metrics (Accuracy, Sensitivity, Specificity,
#'                Precision, Recall, F1 score, MCC) for each threshold.}
#' \item{AUC}{A list containing the AUROC and AUPRC values.}
#' \item{Predictions}{A data frame with the predicted probabilities for each class (e.g., `yes` or `no`).}
#'
#' @details
#' This function first generates predictions for the test dataset using the trained machine learning model.
#' It then calculates performance metrics for a range of threshold values and selects the threshold that maximizes
#' the chosen metric (e.g., Accuracy, F1 score, etc.). The function returns the metrics for the best threshold,
#' including AUROC and AUPRC, and produces a confusion matrix plot that compares predicted versus actual labels.
#'
#' The confusion matrix plot is saved as a PDF with the name `Confusion_Matrix_<file.name>.pdf` if a valid
#' `file.name` is provided.
#'
#'
#' @seealso \code{\link[caret]{confusionMatrix}}, \code{\link[caret]{varImp}}, \code{\link[ggplot2]{ggplot}}
#'
#' @import caret
#' @import dplyr
#' @import ggplot2
#' @import reshape2
#' @import grDevices
#' @export
compute_prediction = function(model, test_data, target_var = NULL, trait.positive = NULL, task_type = "classification", time_var = NULL, event_var = NULL, stack = FALSE, file.name = NULL, return = FALSE){

  # --------------------------------------
  # Classification task
  # --------------------------------------
  if(task_type == "classification"){
    if(is.null(target_var) | is.null(trait.positive)) stop("target_var and trait.positive must be provided for classification")
    target = as.factor(ifelse(target_var == trait.positive, 'yes', 'no'))
    target <- factor(target, levels = c("no", "yes"))  # Order (just in case) to ensure positive class is not well defined

    cat("Predicting target variable using provided ML model.................................................\n")

    if(stack == FALSE){
      test_data = test_data[,colnames(test_data)%in%model[["coefnames"]]] #Only use features defined in the model
      features <- colnames(test_data)
      method = model$method
      are_equal = dplyr::setequal(model[["coefnames"]], features)
      if(are_equal == F){
        stop("Testing set does not count with the same features as model")
      }
      #Predict target variable
      predict <- data.frame(stats::predict(model, test_data, type = "prob"))
    }else{
      super.learner = model$Meta_learner
      ml.models = model$ML_models
      base.models = model$Base_models
      method = "Meta-learner"

      #Learning from simple meta-learner
      base_predictions = list()
      for (i in 1:length(base.models)) {
        base_predictions[[i]] = stats::predict(ml.models[[base.models[i]]], test_data, type = "prob")$yes
        names(base_predictions)[i] = base.models[i]
      }

      base_predictions = do.call(cbind, base_predictions)

      predict = data.frame(stats::predict(super.learner, base_predictions, type = "prob"))
    }

    #Get metrics
    sens_spec = get_sensitivity_specificity(predict, target, method)

    auroc <- list(
      mean  = calculate_auroc(sens_spec$fpr, sens_spec$Sensitivity),
      lower = NA,
      upper = NA
    )

    auprc <- list(
      mean  = calculate_auprc(sens_spec$Recall, sens_spec$Precision),
      lower = NA,
      upper = NA
    )

    #Bootstrap confidence intervals
    boot_auc <- bootstrap_auc(
      predict = predict,
      target = target,
      method = method,
      B = 1000
    )

    auroc$lower <- boot_auc$AUROC$lower
    auroc$upper <- boot_auc$AUROC$upper

    auprc$lower <- boot_auc$AUPRC$lower
    auprc$upper <- boot_auc$AUPRC$upper

    if(return == TRUE){
      # Plot ROC and PRC
      get_curves(data = sens_spec,
                 spec = "Specificity",
                 sens = "Sensitivity",
                 reca = "Recall",
                 prec = "Precision",
                 color = "model",
                 auc_roc = auroc,
                 auc_prc = auprc,
                 file.name = file.name)
    }

    cat("Prediction finished!.................................................\n")

    return(list(Metrics = sens_spec, AUC = list("AUROC" = auroc, "AUPRC" = auprc), Predictions = predict))
  }

  # --------------------------------------
  # Survival task
  # --------------------------------------
  if(task_type == "survival"){
    if(is.null(time_var) | is.null(event_var)) stop("time_var and event_var must be provided for survival")

    cat("Predicting survival outcomes using provided ML model...\n")

    # Predict survival using pre-trained model
    if(!stack){
      test_data = test_data %>%
        dplyr::mutate("event" = event_var,
                      "time" = time_var)

      # Standard survival model
      prediction = predict_and_evaluate_survival(model$Model_object, test_data, "time", "event")

      # Kaplan-Meier plot
      if(return == TRUE){
        plot_survival_performance(df_test = test_data, prediction = prediction, n_groups = 2, file_name = file.name)
      }

      test_data$.pred = as.numeric(unlist(prediction$preds))

    } else {
      stop("Stacked survival models not yet implemented")
    }

    cat("Prediction finished!.................................................\n")

    return(prediction)
  }

}

#' Calculates accuracy values from prediction
#'
#' This function calculates the accuracy of a model based on the provided metrics and the true target values.
#' The accuracy is computed as the ratio of correct predictions (both true positives and true negatives)
#' to the total number of predictions.
#'
#' @param metrics A data frame with metrics obtained using `get_sensitivity_specificity()`, containing
#'   at least two columns: "Sensitivity" and "Specificity".
#' @param target A character vector containing the true values from the target variable.
#'   It should have the same length as the predictions.
#'
#' @return A numeric vector representing the accuracy values.
#'   The result is the fraction of correct predictions out of all predictions.
#'
calculate_accuracy <- function(metrics, target) {
  sensitivity = metrics[,"Sensitivity", drop = T]
  specificity = metrics[,"Specificity", drop = T]
  total_positives = sum(target == "yes")
  total_negatives = sum(target == "no")
  TP <- sensitivity * total_positives
  FN <- total_positives - TP
  TN <- specificity * total_negatives
  FP <- total_negatives - TN

  accuracy <- (TP + TN) / (TP + TN + FP + FN)

  return(accuracy)
}

#' Calculate confusion values
#'
#' This function computes the confusion matrix values (True Positives, False Negatives,
#' True Negatives, and False Positives) based on the given metrics and the true target values.
#'
#' @param metrics A data frame with metrics obtained using `get_sensitivity_specificity()`, containing
#'   at least two columns: "Sensitivity" and "Specificity".
#' @param target A character vector containing the true values from the target variable.
#'   It should have the same length as the predictions.
#'
#' @return A list containing four elements:
#'   - `TP`: True Positives
#'   - `FN`: False Negatives
#'   - `TN`: True Negatives
#'   - `FP`: False Positives
#'
calculate_confusion_values <- function(metrics, target) {
  sensitivity <- metrics[,"Sensitivity", drop = T]
  specificity <- metrics[,"Specificity", drop = T]

  # Count total positives and negatives
  total_positives <- sum(target == "yes")
  total_negatives <- sum(target == "no")

  # Calculate confusion matrix values
  TP <- sensitivity * total_positives
  FN <- total_positives - TP
  TN <- specificity * total_negatives
  FP <- total_negatives - TN

  return(list(TP = TP, FN = FN, TN = TN, FP = FP))
}

#' Calculate precision values
#'
#' This function calculates the precision of a model based on the provided metrics and the true target values.
#' Precision is defined as the ratio of true positive predictions to all positive predictions (true positives + false positives).
#'
#' @param metrics A data frame with metrics obtained using `get_sensitivity_specificity()`, containing
#'   at least two columns: "Sensitivity" and "Specificity".
#' @param target A character vector containing the true values from the target variable.
#'   It should have the same length as the predictions.
#'
#' @return A numeric vector representing the precision values.
#'   Precision is the fraction of true positive predictions among all positive predictions.
#'
calculate_precision <- function(metrics, target) {
  confusion_values <- calculate_confusion_values(metrics, target)
  TP <- confusion_values$TP
  FP <- confusion_values$FP

  precision <- TP / (TP + FP)

  return(precision)
}

#' Calculate recall values
#'
#' This function calculates the recall (also known as sensitivity) of a model based on the provided metrics and the true target values.
#' Recall is defined as the ratio of true positive predictions to all actual positive instances (true positives + false negatives).
#'
#' @param metrics A data frame with metrics obtained using `get_sensitivity_specificity()`, containing
#'   at least two columns: "Sensitivity" and "Specificity".
#' @param target A character vector containing the true values from the target variable.
#'   It should have the same length as the predictions.
#'
#' @return A numeric vector representing the recall values.
#'   Recall is the fraction of actual positive instances that were correctly predicted.
#'
calculate_recall <- function(metrics, target) {
  confusion_values <- calculate_confusion_values(metrics, target)
  TP <- confusion_values$TP
  FN <- confusion_values$FN

  # Calculate recall (sensitivity)
  recall <- TP / (TP + FN)

  return(recall)
}

#' Get performance curves
#'
#' This function generates and saves the Receiver Operating Characteristic (ROC) curve
#' and Precision-Recall curve based on the provided metrics. It also includes the AUC values
#' for both curves in the plot legends.
#'
#' @param data A data frame containing the prediction metrics.
#' @param spec The name of the column containing the specificity values.
#' @param sens The name of the column containing the sensitivity values.
#' @param reca The name of the column containing the recall values.
#' @param prec The name of the column containing the precision values.
#' @param color The name of the column containing the cohort names. Each cohort will have a
#'        corresponding color in the plot. Multiple cohorts will result in different curves.
#' @param auc_roc A numeric value representing the AUC for the ROC curve.
#' @param auc_prc A numeric value representing the AUC for the Precision-Recall curve.
#' @param file.name A character string used as the file name prefix for saving the plots.
#'
#' @return Saves two PDF plots: one for the ROC curve and one for the Precision-Recall curve
#'         in the "Results/" directory.
#' @export
#'
get_curves = function(data, spec, sens, reca, prec, color, auc_roc, auc_prc, file.name){

  data = data %>%
    dplyr::mutate(specificity = data[,spec],
                  sensitivity = data[,sens],
                  recall = data[,reca],
                  precision = data[,prec],
                  color = data[,color])

  ## ---- ROC Curve ----
  data$sens_lower <- auc_roc$lower
  data$sens_upper <- auc_roc$upper

  data$color_label <- paste0(data$model_color,
                             "\n(AUC-ROC = ", round(auc_roc$mean, 2),
                             " [", round(auc_roc$lower[1],2), "-", round(auc_roc$upper[1],2), "])")

  roc <- ggplot(data, aes(x = 1 - specificity, y = sensitivity, color = color_label)) +
    geom_line(size = 1.2) +
    geom_ribbon(aes(x = 1 - specificity, ymin = sens_lower, ymax = sens_upper, fill = color_label),
                alpha = 0.2, color = NA) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey") +
    scale_color_brewer(palette = "Set1") +
    scale_fill_brewer(palette = "Set1") +
    labs(title = "ROC Curve", x = "1 - Specificity", y = "Sensitivity") +
    theme_minimal(base_size = 14) +
    theme(legend.title = element_blank(),
          legend.position = "bottom",
          axis.title = element_text(face = "bold"))

  ## ---- PRC Curve ----
  data$prec_lower <- auc_prc$lower
  data$prec_upper <- auc_prc$upper

  data$color_label_prc <- paste0(data$model_color,
                                 "\n(AUC-PRC = ", round(auc_prc$mean, 2),
                                 " [", round(auc_prc$lower[1],2), "-", round(auc_prc$upper[1],2), "])")

  prc <- ggplot(data, aes(x = recall, y = precision, color = color_label_prc)) +
    geom_line(size = 1.2) +
    geom_ribbon(aes(x = recall, ymin = prec_lower, ymax = prec_upper, fill = color_label_prc),
                alpha = 0.2, color = NA) +
    scale_color_brewer(palette = "Set1") +
    scale_fill_brewer(palette = "Set1") +
    labs(title = "Precision-Recall Curve", x = "Recall", y = "Precision") +
    ylim(0,1) +
    theme_minimal(base_size = 14) +
    theme(legend.title = element_blank(),
          legend.position = "bottom",
          axis.title = element_text(face = "bold"))

  ## ---- Save PDFs ----
  grDevices::pdf(paste0("Results/ROC_curve_", file.name, ".pdf"), width = 6, height = 6)
  print(roc)
  grDevices::dev.off()

  grDevices::pdf(paste0("Results/PRC_curve_", file.name, ".pdf"), width = 6, height = 6)
  print(prc)
  grDevices::dev.off()

}

#' Extract ML models from a directory based on specific AUC score
#'
#' This function searches a directory for machine learning models and filters them based on a
#' specified AUC threshold for either the ROC or Precision-Recall curves. It returns a list of
#' model file names that meet the specified AUC criteria.
#'
#' @param folder_path A character string specifying the directory path where the machine learning
#'                    models are stored.
#' @param metric A character string indicating which AUC metric to use. Choose either "ROC" or "PRC".
#' @param AUC A numeric value representing the minimum acceptable AUC score for the models.
#'
#' @return A character vector with the file paths of the ML models that meet the AUC criteria.
#' @export
#'
find.ML.models = function(folder_path, metric, AUC){

  #Read ML models
  res <- list.files(folder_path, pattern = "\\.rds$", full.names = TRUE)

  # Remove ML models with result NULL or AUC scores below threshold (parameter = "AUC")
  valid_indices <- logical(length(res))
  for (i in seq_along(res)) {
    model <- tryCatch(readRDS(res[i]), error = function(e) NULL)
    if (!is.null(model) && !is.null(model[["result"]])) {
      auc_value <- if (metric == "ROC") {
        model[["result"]][["AUC"]][["AUROC"]]
      } else if (metric == "PRC") {
        model[["result"]][["AUC"]][["AUPRC"]]
      } else {
        stop("Metric selected for AUC doesn't exist. Choose between ROC or PRC")
      }
      valid_indices[i] <- auc_value <= AUC
    }
  }

  # Filter files based on valid indices
  res <- res[valid_indices]

  if (length(res) == 0) {
    stop("No ML models with AUC below ", AUC, " were found. Try with another value.")
  }

  return(res)

}

feature.importance.alignment = function(model){

  #positive and negative class are defined based on the factor() from trait (this have been defined in compute_ML() function already)
  importance = model[["result"]][["Variable_importance"]]
  features_values = model[["result"]][["Cell_groups"]][[1]]
  trait = model[["result"]][["Model"]][["trainingData"]]

  for (i in 1:ncol(features_values)) {
    logreg = stats::glm(trait$.outcome ~ features_values[,i], family = binomial) #Calculate logistic regression using cell groups values and trait outcome from training
    beta = logreg$coefficients[[2]] #Extract beta coefficient from regression
    idx = which(importance$features == colnames(features_values)[i])
    if(beta >= 0){ #If beta is positive means positive association to positive target variable
      importance$final_importance[idx] = importance$final_importance[idx]*1
    }else{ #If beta is negative means there is a reduced likelihood of being in the positive class
      importance$final_importance[idx] = importance$final_importance[idx]*-1
    }
  }

  return(importance)
}

compute_platt.scaling = function(obs, yes){ ####DEPRECATED
  data = data.frame(obs = obs, yes = yes) #Create df from obs and yes to avoid nested problems using dplyr() when grouping by resamples
  # Fit a logistic regression model
  glm_model = stats::glm(obs ~ yes, family = binomial, data = data)
  # Predict calibrated probabilities
  calibrated_prob = as.numeric(predict(glm_model, type = "response")) # "response" returns only probabilities for 'yes', we dont specify new_data argument cause we are predicting on the same set where the training was done

  return(as.numeric(calibrated_prob))
}

# This function creates stratified folds while preserving dataset proportions, useful when doing Leaving-one-dataset-out (LODO) approach
construct_stratified_cohort_folds = function(train_data, batch_id, target_id, k_folds, n_rep){

  ## Named cohort and target variable
  train_data = train_data %>%
    dplyr::mutate(dataset = as.factor(train_data[,batch_id]),
                  target = as.factor(train_data[,target_id]))

  # Create stratified folds across different repeats
  folds_list <- purrr::map(1:n_rep, ~{
    # Split within each dataset
    dataset_folds <- train_data %>%
      dplyr::group_by(dataset) %>% #Groups data by dataset
      dplyr::group_split() %>% #Splits into a list where each dataset is a separate dataframe
      purrr::map(~{
        # get global indices for this dataset subset
        global_idx <- which(train_data$dataset == unique(.x$dataset))

        #Inside each dataset, stratified k-fold splits are created ensuring class balance
        local_folds <- caret::createFolds(.x$target, k = k_folds, returnTrain = TRUE)

        # convert local fold indices to global indices
        purrr::map(local_folds, ~ global_idx[.x])
      })

    fold_indices <- vector("list", k_folds) #Creates an empty list for k folds

    for (i in seq_len(k_folds)) {
      # Extracts the i-th fold from each dataset
      # unlist() now merges globally-corrected indices from all datasets
      fold_indices[[i]] <- unlist(purrr::map(dataset_folds, ~ .x[[i]])) #Merges fold indices from all datasets into a single vector
    }

    return(fold_indices)
  })

  multifolds <- unlist(folds_list, recursive = FALSE)

  # Name folds like "Fold1.Rep1", "Fold2.Rep1", ..., "FoldK.RepN"
  fold_names <- unlist(lapply(1:n_rep, function(rep) {
    paste0("Fold", 1:k_folds, ".Rep", rep)
  }))
  names(multifolds) <- fold_names

  return(multifolds)
}

calculate_cv_metrics = function(ml_model, metric, hyperparameters = NULL){

  ## List hyperparameters
  if(is.null(hyperparameters)==F){
    group_vars = c("Resample", hyperparameters)
  }else{
    group_vars = "Resample"
  }

  # -------------------------------------------------------------------------

  ############# Prediction_folds

  # Stores all OUT-OF-FOLD predictions generated during CV.
  # For each resample (FoldX.RepY), the model is trained on the training folds and predictions are produced for the held-out fold.
  # Therefore, every row corresponds to the prediction for one sample when that sample was part of the test fold.
  # Columns:
  # - pred: predicted class label
  # - obs: true class label
  # - no / yes: predicted probabilities for each class
  # - rowIndex: index of the sample in the dataset
  # - parameter: hyperparameter configuration used
  # - Resample: identifier of the CV split (fold and repetition)
  # - AUROC / AUPRC: performance metrics computed using all predictions
  # within that resample.
  # NOTE: AUROC and AUPRC are constant within each Resample because they are calculated once per fold using all samples in that fold.

  ml_model$Prediction_folds = ml_model$Prediction_folds %>%
    dplyr::group_by(dplyr::across(dplyr::all_of(group_vars))) %>% ## Group by resample and parameters
    dplyr::mutate(AUROC = calculate_auc_roc_resample(obs, yes), # Calculate AUC-ROC if metric is "AUROC"
                  AUPRC = calculate_auc_prc_resample(obs, yes) # Calculate AUC-PRC if metric is "AUPRC"
    ) %>%
    dplyr::ungroup() %>%
    data.frame()

  # -------------------------------------------------------------------------

  ############# Prediction metrics

  # Compute classification metrics for each resample (fold × repeat)
  # Using the out-of-fold predictions stored in `Prediction_folds`, we compute threshold-based classification metrics such as:
  # Accuracy, Sensitivity (Recall), Specificity, Precision, F1-score and MCC.
  # Steps:
  #  1. Group predictions by the resampling identifiers (`group_vars`), which typically include the CV split (Resample) and the model parameter set.
  #  2. For each group (i.e., each fold/repeat), compute the confusion-matrix derived metrics using `get_sensitivity_specificity()`.
  #  3. Combine the resulting metrics back into a single data frame.

  metrics <- ml_model$Prediction_folds %>%
    dplyr::group_by(dplyr::across(dplyr::all_of(group_vars))) %>%
    dplyr::group_split() %>%
    purrr::map(~ get_sensitivity_specificity(.x, .x$obs, "test")) %>%
    dplyr::bind_rows() %>%
    dplyr::select(-model)

  ml_model$Prediction_folds = ml_model$Prediction_folds %>%
    dplyr::select(-yes) %>% #remove yes probabilities from resamples (only keep those ordered for calculated the metrics)
    dplyr::bind_cols(metrics) %>%
    dplyr::select(-pred, -obs, -no) %>%
    dplyr::select(Resample, yes, dplyr::everything())

  if(is.null(hyperparameters) == F){
    # -------------------------------------------------------------------------

    ############# Results_folds (when there are hyperparameters)

    # If the model includes hyperparameters, we first evaluate the performance of each hyperparameter combination using the out-of-fold predictions stored in `Prediction_folds`.
    # Step 1 — Evaluate hyperparameter combinations
    #    Predictions are grouped by the hyperparameter values (not by Resample).
    #    This allows aggregating performance across all CV resamples for each hyperparameter configuration.
    #    The median AUROC, AUPRC, and Accuracy are computed to obtain robust estimates of performance.
    # Step 2 — Select the best hyperparameter configuration
    #   The optimal parameter set is selected based on the chosen optimization metric (`metric`, typically AUROC or AUPRC).
    #   The row with the highest value is identified and stored in `ml_model$bestTune`.
    # Step 3 — Filter predictions for the best configuration
    #   The full prediction table is filtered to keep only the rows corresponding to the selected hyperparameter combination.
    #   This ensures that subsequent evaluation is performed only on predictions generated by the tuned model.
    # Step 4 — Compute performance per resample
    #   After selecting the best hyperparameters, performance metrics are recalculated per resample (Fold × Repeat).
    #   This produces a table where each row represents the performance of the tuned model on a specific cross-validation split.

    # Result:

    # - `ml_model$Results_folds` contains the average across all resamples per hyperparameter configuration.
    # - `ml_model$bestTune` contains the optimal hyperparameter configuration.
    # - `df_avg` contains AUROC, AUPRC, and Accuracy for each resample using the tuned model

    ## Average across all resamples and per hyperparameter configuration
    ml_model$Results_folds <- ml_model$Prediction_folds %>%
      dplyr::group_by(dplyr::across(dplyr::all_of(hyperparameters))) %>% ## here we only group by hyperparameter not resample (to choose best hyperparams combination)
      dplyr::summarise(
        AUROC = median(AUROC, na.rm = TRUE),
        AUPRC = median(AUPRC, na.rm = TRUE),
        Accuracy = median(Accuracy, na.rm = TRUE),
        .groups = "keep"  # Ensures all groups are retained
      ) %>%
      dplyr::ungroup()

    ##### Hyperparameter tuning
    tune = which.max(data.frame(ml_model$Results_folds)[,metric]) # Tuning parameter (select combination with top AUROC or AUPRC)

    ml_model$bestTune = ml_model$Results_folds %>%
      dplyr::slice(tune) %>%  # Extract the top row
      dplyr::select(dplyr::all_of(hyperparameters)) %>%
      dplyr::mutate(dplyr::across(where(is.numeric), as.numeric))%>%
      as.data.frame()

    filter_conditions <- ml_model$bestTune[1, , drop = FALSE] #Take tuned parameters

    ##### Filter predictions to only those with tuned hyperparams
    df_filtered <- ml_model$Prediction_folds # All predictions with all hyperparams combinations
    for (col in names(filter_conditions)) {
      df_filtered <- df_filtered[df_filtered[[col]] == filter_conditions[[col]], ] # Filter by keeping only rows where the column matches the corresponding value in filter_conditions (Continue refining until all conditions are applied)
    }

    ##### Average across resamples with tuned hyperparams
    df_avg = df_filtered %>%
      dplyr::group_by(Resample) %>%
      dplyr::summarise(
        AUROC = median(AUROC, na.rm = TRUE),
        AUPRC = median(AUPRC, na.rm = TRUE),
        Accuracy = median(Accuracy, na.rm = TRUE)
      ) %>%
      dplyr::ungroup()

  }else{

    # -------------------------------------------------------------------------

    ############# Results_folds (when there are no hyperparameters)

    # Aggregate cross-validation performance across resamples
    # Here we summarise the performance metrics obtained from the out-of-fold predictions stored in `Prediction_folds`.
    # Since metrics such as AUROC, AUPRC, and Accuracy are computed per resample (fold × repeat), we aggregate them to obtain a single
    # representative value for the model. The median is used instead of the mean to provide a more robust estimate that is less sensitive
    # to extreme values across folds (in other words, i am taking the mediam across resamples).
    # Steps:
    # 1. Compute the median AUROC, AUPRC, and Accuracy across all resamples.
    # 2. Remove hyperparameter-related columns and redundant performance summaries from `Results_folds`.
    # 3. Append the aggregated metrics (`df_avg`) to `Results_folds`.

    # Finally, `bestTune` is set to `"none"` because this workflow does not involve hyperparameter tuning; the model is trained with fixed settings.
    # Result:
    # - `Results_folds` now contains the aggregated cross-validation performance metrics for the model.

    # -------------------------------------------------------------------------

    df_avg <- ml_model$Prediction_folds %>%
      dplyr::summarise(
        AUROC = median(AUROC, na.rm = TRUE),
        AUPRC = median(AUPRC, na.rm = TRUE),
        Accuracy = median(Accuracy, na.rm = TRUE),
        .groups = "keep"  # Ensures all groups are retained
      ) %>%
      dplyr::ungroup()

    ml_model$Results_folds <- ml_model$Results_folds %>%
      dplyr::select(-dplyr::all_of(hyperparameters), -Accuracy, -Kappa, -AccuracySD, -KappaSD) %>%
      dplyr::bind_cols(df_avg)

    ml_model$bestTune = tibble::tibble(parameter = "none")
  }

  # -------------------------------------------------------------------------

  ############# Resample_matrix

  # This step integrates the per-resample metrics from `df_avg` into the model's `Resample_matrix`.
  # Result:
  # - `ml_model$Resample_matrix` now contains resample-level performance metrics for the tuned model

  ml_model$Resample_matrix = ml_model$Resample_matrix %>%
    dplyr::select(-Accuracy, -Kappa) %>%
    dplyr::arrange(match(Resample, df_avg$Resample)) %>%
    dplyr::select(-Resample) %>%
    dplyr::bind_cols(df_avg)  %>%
    { if ("Resample" %in% colnames(.)) dplyr::select(., -Resample) else . }


  return(list(
    Prediction_folds = ml_model$Prediction_folds,   # Contains all out-of-fold predictions per resample and hyperparameter combination, with computed metrics (AUROC, AUPRC, Accuracy, etc.).
    Resample_matrix = ml_model$Resample_matrix,     # Contains per-resample performance metrics (AUROC, AUPRC, Accuracy) for the tuned model.
    Results_folds = ml_model$Results_folds,         # Aggregated performance metrics (median AUROC, AUPRC, Accuracy) per hyperparameter combination.
    bestTune = ml_model$bestTune                    # The optimal hyperparameter configuration selected based on the chosen metric.
  ))

}


compute_shap_values <- function(model_trained, data_train, task_type = "classification", target_col = NULL,
                                trait.positive, time_col = NULL, event_col = NULL, n_cores = 2, file.name = NULL) {

  if(task_type == "classification"){
    if(is.null(target_col)) stop("target_col must be provided for classification tasks")

    data_train <- data_train %>%
      dplyr::mutate(
        target = factor(
          ifelse(data_train[[target_col]] == trait.positive, "yes", "no"),
          levels = c("no", "yes")
        )
      )

    X = data_train %>%
      dplyr::select(-dplyr::all_of(target_col))

    pred_fun = function(object, newdata){
      predict(object, newdata, type = "prob")[,yes]
    }

    # Determine resamples from the model object
    resamples <- unique(model_trained$pred$Resample)

    # Extract ML model
    method = model_trained$method

  }else if(task_type == "survival"){
    if(is.null(time_col) || is.null(event_col)) stop("time_col and event_col must be provided for survival tasks")
    X = data_train %>%
      dplyr::select(-dplyr::all_of(c(time_col, event_col)))

    # Determine resamples from the model object
    resamples <- unique(model_trained$Resample_matrix$Resample)

    # Extract ML model
    method = unique(model_trained$Resample_matrix$model)

    pred_fun <- function(object, newdata){

      prediction = predict_and_evaluate_survival(object, newdata, "time", "event")
      preds = as.numeric(unlist(prediction$preds))

      return(preds)
    }

  }

  gc() #clean memory before start

  cat("Computing SHAP values...\n\n")

  # Get tuned hyperparameters
  filter_conditions <- model_trained$bestTune[1, , drop = FALSE]

  # Filter predictions for best tune
  # if (any(filter_conditions != "none")) {
  #   for (col in names(filter_conditions)) {
  #     model_trained$pred <- model_trained$pred[model_trained$pred[[col]] == filter_conditions[[col]], ]
  #   }
  # }

  # Register parallel backend
  cl <- parallel::makeCluster(n_cores)
  doParallel::registerDoParallel(cl)

  importance_list <- foreach::foreach(resample = resamples, .packages = c("dplyr", "caret", "censored")) %dopar% {

      source("~/Documents/pipeML/R/machine_learning.R")

      if(task_type == "classification"){

        test_index <- model_trained$pred %>%
          dplyr::filter(Resample == resample) %>%
          dplyr::pull(rowIndex)

        train_data_fold <- data_train[-test_index, ]
        test_data_fold  <- data_train[test_index, ]

        fit <- if (any(filter_conditions != "none")) {

          caret::train(
            target ~ .,
            data = train_data_fold,
            method = method,
            trControl = caret::trainControl(method = "none", classProbs = TRUE),
            tuneGrid = filter_conditions,
            metric = "Accuracy"
          )

        } else {
          caret::train(
            target ~ .,
            data = train_data_fold,
            method = method,
            trControl = caret::trainControl(method = "none", classProbs = TRUE),
            metric = "Accuracy"
          )
        }

        # Separating features from target
        X_train <- train_data_fold[, setdiff(names(train_data_fold), "target")]
        X_test  <- test_data_fold[, setdiff(names(test_data_fold), "target")]

        pred_probs <- pred_fun(fit, X_test)

      }else if(task_type == "survival"){

        test_index <- model_trained$Resample_matrix %>%
          dplyr::filter(Resample == resample) %>%
          dplyr::pull(rowIndex)

        train_data_fold <- data_train[-test_index, ]
        test_data_fold  <- data_train[test_index, ]

        res = compute_ml_survival(train_data_fold, test_data_fold, outcome_col = time_col,
                                  event_col = event_col, model = method,
                                  models_hyperparameters = if (is.null(filter_conditions)) NULL else
                                    list(filter_conditions),
                                  return_model = T)

        fit = res$Model

        # Separating features from target
        X_test <- test_data_fold %>% select(-all_of(c(time_col, event_col)))
        X_train <- train_data_fold %>% select(-all_of(c(time_col, event_col)))

        pred_probs = res$Metrics$predictions
      }



      # Check for trivial predictions (same for all "yes"/"no")
      if (length(unique(round(pred_probs, 5))) > 1) { ## equal prob in all the samples
        # If predictions are not trivial, compute SHAP values
        shap_values <- fastshap::explain(
          object = fit,
          X = X_train,
          pred_wrapper = pred_fun,
          newdata = X_test,
          nsim = 1000,
          adjust = TRUE
        )

        shap_values <- as.data.frame(shap_values)
        shap_values$Samples <- rownames(X_test)

        shap_values

      } else {
        message("Skipping resample ", resample, " (constant predictions)")
        NULL
      }

  }

  # Stop the cluster after parallel execution
  parallel::stopCluster(cl)
  unregister_dopar() # Stop Dopar from running in the background

  gc()

  importance_list = importance_list[!sapply(importance_list, is.null)]

  # Check if trivial predictions were found
  if (length(importance_list)==0) {
    warning("Trivial predictions were found in some resamples. SHAP values cannot be calculated")
    return(NULL)
  }

  # Combine and summarize importance results
  shap_df <- dplyr::bind_rows(importance_list)

  plot_shap_stability(shap_df, file.name = file.name)

  shap_df <- shap_df %>%
    dplyr::group_by(Samples) %>%
    dplyr::summarise(
      dplyr::across(where(is.numeric), \(x) median(x, na.rm = TRUE)),
      .groups = "drop"
    ) %>%
    dplyr::arrange(match(Samples, rownames(data_train))) %>%
    data.frame()

  rownames(shap_df) <- shap_df$Samples
  shap_df$Samples <- NULL

  cat("SHAP analisis finished! \n\n")

  return(shap_df)
}

plot_shap_stability <- function(shap_df, file.name){

  shap_long <- shap_df %>%
    tidyr::pivot_longer(
      cols = -Samples,
      names_to = "feature",
      values_to = "shap_value"
    )

  # mean |SHAP| per feature per sample occurrence
  shap_summary <- shap_long %>%
    dplyr::mutate(abs_shap = abs(shap_value)) %>%
    dplyr::group_by(feature) %>%
    dplyr::summarise(
      mean_importance = mean(abs_shap, na.rm = TRUE),
      sd_importance   = stats::sd(abs_shap, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    dplyr::arrange(mean_importance)

  p <- ggplot2::ggplot(shap_summary,
              ggplot2::aes(x = stats::reorder(feature, mean_importance),
                           y = mean_importance)) +
    ggplot2::geom_col(fill = "steelblue") +
    ggplot2::geom_errorbar(
      ggplot2::aes(
        ymin = mean_importance - sd_importance,
        ymax = mean_importance + sd_importance
      ),
      width = 0.2
    ) +
    ggplot2::coord_flip() +
    ggplot2::theme_minimal(base_size = 14) +
    ggplot2::labs(
      title = "SHAP feature importance stability across resamples",
      x = "Feature",
      y = "Mean |SHAP| importance"
    )

  grDevices::pdf(paste0("Results/SHAP_stability_resample_", file.name, ".pdf"))
  print(p)
  dev.off()

}

#' Compute variable importance using SHAP values
#'
#' Computes the variable importance for a machine learning model using SHAP (SHapley Additive exPlanations) values.
#'
#' @param model A trained machine learning model.
#' @param stacking A logical value indicating whether the model was trained using stacking (default is FALSE).
#' @param n_cores An integer specifying the number of workers to use for parallel computation (default is 2).
#'
#' @return A data frame with SHAP values representing the variable importance for each feature.
#'
#' @details If `stacking` is TRUE, the function computes the SHAP values for each base model in the stacked ensemble model and
#' averages them. If `stacking` is FALSE, the function computes the SHAP values for the provided single machine learning model.
#' The computed SHAP values are returned as a data frame with features as rows and samples as columns.
#'
#' @export
#'
compute_variable.importance = function(model, stacking = FALSE, n_cores = 2){

  if(stacking == TRUE){
    base_models = model$Base_models
    ml_models = model$ML_models
    train_data = model$ML_models$BAG$trainingData %>% #All training datasets are the same so BAG it's choosing randomly
      dplyr::rename(target = .outcome)

    importance = list() #Save variable importance of each base model
    for (i in 1:length(base_models)) {
      importance[[i]] = compute_shap_values(ml_models[[base_models[i]]], train_data, ml_models[[base_models[i]]]$method, n_cores) ## Compute SHAP values
    }
    non_null_importance = Filter(Negate(is.null), importance)

    if (length(non_null_importance) == 0) {
      return(NULL)
    } else {
      importance_df <- Reduce(`+`, non_null_importance) / length(non_null_importance)
    }
  }else{
    train_data = model$Model$trainingData %>%
      dplyr::rename(target = .outcome)
    ml_method = model$Model$method

    importance_df = compute_shap_values(model$Model, train_data, ml_method, n_cores) ## Compute SHAP values for variable importance
  }

  return(importance_df)
}


unregister_dopar <- function() {
  if (!is.null(foreach::getDoParRegistered())) {
    # switch back to sequential backend
    foreach::registerDoSEQ()
    gc()
  }
}

model_boruta_selection <- function(model,
                                   boruta_iterations = 100,
                                   fix_boruta = TRUE,
                                   file_name = NULL,
                                   boruta_threshold = 0.8,
                                   tentative = FALSE) {
  cat("Feature selection using Boruta...............................................................\n\n")

  # Run Boruta
  res_boruta <- feature.selection.boruta(
    data = model,
    iterations = boruta_iterations,
    fix = fix_boruta,
    file_name = file_name,
    doParallel = FALSE,
    workers = NULL,
    threshold = boruta_threshold,
    return = FALSE
  )

  # Decide which features to keep
  if (!tentative) {
    if (length(res_boruta$Confirmed) <= 1) {
      message("Not enough features selected for training a model")
      return(NULL)
    } else {
      cat("\nKeeping only features confirmed in more than 80% of the times for training...............................................................\n\n")
      cat("If you want to consider also tentative features, please specify tentative = TRUE in the parameters.\n\n")
      selected_model <- model[, colnames(model) %in% res_boruta$Confirmed, drop = FALSE] %>%
        dplyr::mutate(target = model$target)
    }
  } else {
    total_features <- length(res_boruta$Confirmed) + length(res_boruta$Tentative)
    if (total_features <= 1) {
      message("Not enough features selected for training a model")
      return(NULL)
    } else {
      cat("Keeping features confirmed and tentative in more than 80% of the times for training...............................................................\n\n")
      selected_model <- model[, colnames(model) %in% c(res_boruta$Confirmed, res_boruta$Tentative), drop = FALSE] %>%
        dplyr::mutate(target = model$target)
    }
  }

  return(selected_model)
}

#' Preprocess Features for Machine Learning
#'
#' This function preprocesses a dataset by removing features with
#' near-zero variance, highly correlated features, and features that are
#' constant within any class of the target variable. It ensures that the
#' resulting feature set is more suitable for machine learning models.
#'
#' @param data A data frame containing predictor features and the target variable.
#' @param target_col A character string specifying the name of the target column
#'   in \code{data}. Default is \code{"target"}.
#' @param cor_thresh A numeric value between 0 and 1 specifying the correlation
#'   threshold for removing highly correlated features. Default is \code{0.9}.
#' @param time_var A character string specifying the name of the time-to-event
#'   column in \code{data}. Used only for survival analysis.
#' @param event_var A character string specifying the name of the event indicator
#'   column in \code{data} (e.g., 0/1). Used only for survival analysis.
#'
#' @details
#' The preprocessing steps include:
#' \enumerate{
#'   \item Removing near-zero variance features (using \code{caret::nearZeroVar}).
#'   \item Removing highly correlated features above the specified threshold
#'         (using \code{caret::findCorrelation}).
#'   \item Removing features that are constant within any class of the target
#'         variable (i.e., provide no discriminatory power across classes).
#' }
#'
#' After preprocessing, the target column is re-attached to the dataset.
#'
#' @return A data frame with the preprocessed features and the target column.
#'
#' @examples
#' \dontrun{
#' library(caret)
#' library(dplyr)
#'
#' set.seed(123)
#' df <- data.frame(
#'   feature1 = c(1, 1, 1, 1, 1),             # constant
#'   feature2 = c(1, 2, 3, 4, 5),             # numeric
#'   feature3 = c(1, 2, 3, 4, 5) * 2,         # highly correlated with feature2
#'   target   = c("A", "A", "B", "B", "B")
#' )
#'
#' clean_df <- preprocess_features(df, target_col = "target", cor_thresh = 0.9)
#' }
#'
#' @importFrom dplyr select all_of
#' @importFrom caret nearZeroVar findCorrelation
#'
#' @export
preprocess_features <- function(data,
                                target_col = NULL,
                                time_var = NULL,
                                event_var = NULL,
                                cor_thresh = 0.9) {
  # ---- Detect outcome type ----
  if (!is.null(target_col) && target_col %in% names(data)) {
    target <- data[[target_col]]
    outcome_type <- "classification"
    data <- data %>% dplyr::select(-dplyr::all_of(target_col))
  } else if (!is.null(time_var) && !is.null(event_var)) {
    if (!time_var %in% names(data) || !event_var %in% names(data)) {
      stop("time_var and event_var must be column names in 'data'.")
    }
    time_col <- data[[time_var]]
    event_col <- data[[event_var]]
    outcome_type <- "survival"
    data <- data %>% dplyr::select(-dplyr::all_of(c(time_var, event_var)))
  } else {
    stop("You must provide either 'target_col' or both 'time_var' and 'event_var'.")
  }

  # ---- Step 1: Remove near-zero variance features ----
  if (ncol(data) > 0) {
    nzv <- caret::nearZeroVar(data, saveMetrics = TRUE)
    data <- data[, !nzv$nzv, drop = FALSE]
  }

  # ---- Step 2: Remove highly correlated features ----
  if (ncol(data) > 1) {
    cor_matrix <- stats::cor(data, use = "pairwise.complete.obs")
    high_cor <- caret::findCorrelation(cor_matrix, cutoff = cor_thresh)
    if (length(high_cor) > 0) {
      data <- data[, -high_cor, drop = FALSE]
    }
  }

  # ---- Step 3: Remove features constant within classes (if classification) ----
  if (outcome_type == "classification" && ncol(data) > 0) {
    constant_features <- c()
    for (i in seq_along(colnames(data))) {
      col_values <- data[[i]]
      for (cl in unique(target)) {
        class_vals <- col_values[target == cl]
        nzv_info <- caret::nearZeroVar(as.data.frame(class_vals), saveMetrics = TRUE)
        if (nzv_info$nzv) {
          constant_features <- c(constant_features, colnames(data)[i])
          break
        }
      }
    }
    if (length(constant_features) > 0) {
      data <- data[, !colnames(data) %in% constant_features, drop = FALSE]
    }
  }

  # ---- Step 4: Add back outcome columns ----
  if (outcome_type == "classification") {
    data[[target_col]] <- target
  } else if (outcome_type == "survival") {
    data[[time_var]] <- time_col
    data[[event_var]] <- event_col
  }

  return(data)
}

#' Train model with optimized hyperparameters for classification tasks
#'
#' This function wraps cross-validation, hyperparameter optimization,
#' and final training into a single workflow. It identifies the best
#' hyperparameters using a custom cross-validation function, reconstructs
#' the training set, preprocesses features, and retrains the model on the
#' complete training data with the selected hyperparameters.
#'
#' @param train_data A data frame containing the full training dataset,
#'   including predictors and the target variable.
#' @param optimized An object returned by \code{compute_custom_k_fold_CV}
#'   containing the optimized hyperparameters and cross-validation results.
#' @param ml_method A character string specifying the machine learning method
#'   to be passed to \code{caret::train}.
#' @param fold_construction_fun A function used to (re)construct training
#'   data partitions given the best hyperparameters.
#' @param fold_construction_args_fixed A named list of additional fixed arguments
#'   to pass to \code{fold_construction_fun}.
#'
#' @details
#' The workflow proceeds in the following steps:
#' \enumerate{
#'   \item Runs cross-validation using \code{compute_custom_k_fold_CV} to
#'         identify the best hyperparameter set.
#'   \item Reconstructs the training set using \code{fold_construction_fun} and
#'         the selected hyperparameters.
#'   \item Preprocesses the training features by removing near-zero variance,
#'         highly correlated, and constant-within-class features
#'         (via \code{\link{preprocess_features}}).
#'   \item Retrains the model on the complete training data with the optimized
#'         hyperparameters.
#' }
#'
#' The returned object mimics a \code{caret} model object but includes
#' additional elements derived from the custom cross-validation.
#'
#' @return A list with the following components:
#' \describe{
#'   \item{\code{Model}}{A trained \code{caret} model object with
#'         results, predictions, and resampling info attached.}
#'   \item{\code{training_set}}{The final preprocessed training dataset.}
#'   \item{\code{custom_output}}{Additional output from \code{fold_construction_fun}.}
#' }
#'
#' @examples
#' \dontrun{
#' library(caret)
#'
#' # Example placeholders
#' train_data <- your_training_data
#' fold_data  <- your_prepared_folds
#'
#' result <- wrapper_train_best_hyperparams_classification(
#'   train_data = train_data,
#'   fold_data  = fold_data,
#'   ml_method  = "rf",
#'   fold_construction_fun = your_fold_fun,
#'   fold_construction_args_fixed = list(arg1 = "value"),
#'   tuneGrid = expand.grid(mtry = 2:4),
#'   ncores = 4
#' )
#'
#' result$Model
#' }
#'
#' @seealso \code{\link{compute_custom_k_fold_CV}},
#'   \code{\link{preprocess_features}}
#'
#' @importFrom caret train trainControl
#' @export
wrapper_train_best_hyperparams_classification <- function(train_data, optimized, ml_method, fold_construction_fun, fold_construction_args_fixed) {

  # Extract optimized hyperparams
  besttune <- optimized$bestTune

  # Create complete training set using tuned hyperparams from custom function
  training_all <- do.call(fold_construction_fun,
                          c(list(data = train_data, bestune = optimized$bestTune), fold_construction_args_fixed))

  # Preprocess features
  training_set <- preprocess_features(training_all[[1]], cor_thresh = 0.9, target_col = "target")

  # Wrap correct model type for lasso and ridge
  if(ml_method == "lasso"){
    ml_method = "glmnet"
  }else if(ml_method == "ridge"){
    ml_method = "glmnet"
  }

  # Retrain ML model with tuned ML hyperparams
  fit <- caret::train(
    target ~ .,
    data = training_set,
    method = ml_method,
    trControl = caret::trainControl(method = "none", allowParallel = F, classProbs = TRUE),
    tuneGrid = besttune %>% select(-colnames(training_all[[3]]))
  )

  # Return caret-like object
  fit$Results_folds  <- optimized$Results_folds
  fit$Prediction_folds     <- optimized$Prediction_folds
  fit$Resample_matrix <- optimized$Resample_matrix ## Resample matrix contains the performance per resample with tuned param conf

  ##### training_all[[2]] needs to be filter to return only features from training_set (not possible because we need to generalize custom_fold function so sometimes the structure will be different)

  return(list(Model = fit, training_set = training_set, custom_output = c(training_all[[2]], list("Parameters" = training_all[[3]]))))
}

#' Train the Best Survival Model Using Optimized Hyperparameters
#'
#' Fits a survival model on the full training data using the optimal
#' hyperparameters obtained from nested cross-validation.
#' This wrapper ensures consistent retraining for different survival
#' model types (Cox, penalized Cox, AFT, tree-based, or ensemble models),
#' and supports preprocessing pipelines such as CellTFusion through a
#' user-provided fold construction function.
#'
#' @param train_data A data frame containing the original training data
#'   used for cross-validation.
#' @param optimized A list output from [`aggregate_results_survival()`] or
#'   [`compute_k_fold_CV_survival()`], containing the best-tuned parameters
#'   (`bestTune`) and model performance summaries.
#' @param ml_method Character string specifying the survival model to train.
#'   Must be one of:
#'   \itemize{
#'     \item `"cox_ph_survival"` — Cox proportional hazards model.
#'     \item `"proportional_hazards_glmnet"` — Penalized Cox (elastic net).
#'     \item `"survreg_flexsurv"` — Parametric AFT model.
#'     \item `"rand_forest_partykit"` — Random survival forest via `partykit`.
#'     \item `"rand_forest_aorsf"` — Oblique random survival forest.
#'     \item `"decision_tree_partykit"` — Single survival tree.
#'     \item `"bag_tree_rpart"` — Bagged CART-based survival trees.
#'     \item `"boost_tree_mboost"` — Gradient boosting for censored data.
#'   }
#' @param fold_construction_fun A custom function used to construct folds and
#'   preprocessed data (e.g., `prepare_CellTFusion_folds()`).
#'   Must accept arguments `data` and optionally `bestune`.
#' @param fold_construction_args_fixed A named list of fixed arguments to pass
#'   to `fold_construction_fun()` (e.g., paths, deconvolution matrices, etc.).
#' @param outcome_col Character string naming the survival time column
#'   (default = `"time"`).
#' @param event_col Character string naming the event indicator column
#'   (default = `"event"`).
#'
#' @details
#' This function performs the following steps:
#' \enumerate{
#'   \item Extracts the optimal hyperparameters from the `optimized` object.
#'   \item Reconstructs the training dataset using the provided
#'         `fold_construction_fun()`, including any custom preprocessing
#'         or feature generation.
#'   \item Applies the optimal hyperparameters to the model specification.
#'   \item Fits the final model using the full training data.
#' }
#'
#' If the selected model type has no tunable hyperparameters,
#' the function automatically detects this and proceeds with
#' the default model configuration.
#'
#' @return A named list containing:
#' \describe{
#'   \item{`Model`}{A list containing the fitted parsnip model object,
#'   resampling results, and tuning information.}
#'   \item{`training_set`}{The final preprocessed training dataset used for fitting.}
#'   \item{`custom_output`}{Additional data returned by the custom fold construction
#'   function (e.g., CellTFusion outputs or parameter tables).}
#' }
#'
#' @export
#'
wrapper_train_best_hyperparams_survival <- function(train_data,
                                                    optimized,
                                                    ml_method,
                                                    fold_construction_fun,
                                                    fold_construction_args_fixed,
                                                    outcome_col = "time",
                                                    event_col = "event") {

  # Extract optimized hyperparams
  besttune <- optimized$bestTune

  # Build full training data using tuned hyperparams
  training_all <- do.call(
    fold_construction_fun,
    c(list(data = train_data, bestune = besttune), fold_construction_args_fixed)
  )

  # training_all[[1]]: features
  # training_all[[2]]: custom CellTFusion output, etc.
  # training_all[[3]]: parameter table

  # Preprocess features
  training_set <- preprocess_features(training_all[[1]], cor_thresh = 0.9,
                                      time_var = "time", event_var = "event")

  # Define the model specification based on ml_method
  model <- ml_method
  model_spec <- NULL

  if (model == "cox_ph_survival") {
    model_spec <- parsnip::proportional_hazards(mode = "censored regression", engine = "survival")

  } else if (model == "proportional_hazards_glmnet") {
    model_spec <- parsnip::proportional_hazards(
      penalty = tune(), mixture = tune(),
      mode = "censored regression", engine = "glmnet"
    )

  } else if (model == "survreg_flexsurv") {
    model_spec <- parsnip::survival_reg(mode = "censored regression", engine = "flexsurv")

  } else if (model == "rand_forest_partykit") {
    model_spec <- parsnip::rand_forest(trees = tune(), mode = "censored regression", engine = "partykit")

  } else if (model == "rand_forest_aorsf") {
    model_spec <- parsnip::rand_forest(trees = tune(), mode = "censored regression", engine = "aorsf")

  } else if (model == "decision_tree_partykit") {
    model_spec <- parsnip::decision_tree(mode = "censored regression", engine = "partykit")

  } else if (model == "bag_tree_rpart") {
    model_spec <- parsnip::bag_tree(mode = "censored regression", engine = "rpart")

  } else if (model == "boost_tree_mboost") {
    model_spec <- parsnip::boost_tree(trees = tune(), mode = "censored regression", engine = "mboost")

  } else {
    stop("Unsupported survival model: ", model)
  }

  model_hyperparams <- if (is.data.frame(besttune)) {
    # Case 1: besttune is a data.frame --> means it has tunable params
    besttune %>%
      dplyr::select(-dplyr::all_of(colnames(training_all[[3]])))

  } else if (is.list(besttune)) {
    # Case 2: besttune is a list --> means came from fixed params
    besttune[setdiff(names(besttune), names(fold_construction_args_fixed))]

  } else {
    stop("`besttune` must be either a data.frame or a list.")
  }

  # If ML model doesnt have hyperparams
  if (is.null(model_hyperparams) || (is.data.frame(model_hyperparams) && ncol(model_hyperparams) == 0) || (is.list(model_hyperparams) && length(model_hyperparams) == 0)){
    model_hyperparams <- NULL
  }

  # Apply hyperparameters to model spec if available
  if (!is.null(model_hyperparams)) {
    model_spec <- do.call(parsnip::set_args, c(list(object = model_spec), model_hyperparams))
  }

  # Build survival formula
  formula_model <- as.formula(paste0("Surv(", outcome_col, ", ", event_col, ") ~ ."))

  wf <- workflows::workflow() %>%
    workflows::add_model(model_spec) %>%
    workflows::add_formula(formula_model)

  # Fit on full training set
  fitted <- parsnip::fit(wf, data = training_set)

  # Create a caret-like output
  fit <- list()
  fit$model <- model
  fit$fitted <- fitted
  fit$Results_folds <- optimized$Results_folds
  fit$Prediction_folds <- optimized$Prediction_folds
  fit$Resample_matrix <- optimized$Resample_matrix
  fit$bestTune <- besttune

  # Return everything consistent with the classification wrapper
  return(list(
    Model = fit,
    training_set = training_set,
    custom_output = c(training_all[[2]], list("Parameters" = training_all[[3]]))
  ))
}

#' Aggregate Nested Cross-Validation Results for Classification or Survival Tasks
#'
#' Aggregates performance metrics from nested cross-validation experiments for
#' either **classification** or **survival** models.
#' This function reads the resampled predictions (often saved per fold and
#' hyperparameter combination) and computes overall performance summaries,
#' identifies the best hyperparameter configuration, and collates
#' per-resample metrics for detailed inspection.
#'
#' @param all_loaded A nested list containing the results of cross-validation runs.
#'   Each element typically corresponds to a fold and may include one or more
#'   hyperparameter configurations and model predictions:
#'   \itemize{
#'     \item For classification: `all_loaded[[fold]][[param]][[model]][[1]]`
#'       contains predictions, and `[[2]]` the corresponding hyperparameters.
#'     \item For survival: `all_loaded[[fold]][[param]][[model]]` contains
#'       the C-index and hyperparameter columns.
#'   }
#' @param task Character string specifying the task type.
#'   Must be one of: `"classification"` or `"survival"`.
#'
#' @details
#' The function automatically detects whether the input structure includes
#' tunable hyperparameters (`has_params`) by inspecting the nested list depth.
#'
#' For **classification** tasks:
#' \itemize{
#'   \item Aggregates fold-level predictions and computes Accuracy and Kappa.
#'   \item Computes the median and MAD (robust SD) across resamples.
#'   \item Selects the hyperparameter set with the highest median Accuracy.
#' }
#'
#' For **survival** tasks:
#' \itemize{
#'   \item Aggregates per-fold C-index values across hyperparameter configurations.
#'   \item Computes the median and MAD of the C-index per configuration.
#'   \item Identifies and returns the best-performing parameter combination.
#' }
#'
#' The function is compatible with results produced by
#' [`compute_k_fold_CV_survival()`] and analogous classification CV pipelines.
#'
#' @return A list of length equal to the number of models evaluated.
#' Each element contains:
#' \describe{
#'   \item{`Prediction_folds`}{All predictions or C-index values per fold.}
#'   \item{`Results_folds`}{Aggregated performance summaries across folds.}
#'   \item{`bestTune`}{Best-performing hyperparameter combination.}
#'   \item{`Resample_matrix`}{Fold-level metrics for the best configuration.}
#' }
#'
#'
#' @export
#'
aggregate_results <- function(all_loaded, task = c("classification", "survival")) {

  # Detect structure (whether has tunable parameters or not)
  if(task == "classification"){
    has_params <-  length(all_loaded[[1]][[1]]) > 3
  }else{
    has_params <-  length(all_loaded[[1]][[1]]) > 5
  }

  n_folds  <- length(all_loaded)
  n_models <- if (has_params) length(all_loaded[[1]][[1]]) else length(all_loaded[[1]])
  n_params <- if (has_params) length(all_loaded[[1]]) else 1

  results <- vector("list", n_models)

  # ============================================================
  # === CLASSIFICATION =========================================
  # ============================================================
  if (task == "classification") {

    for (m in seq_len(n_models)) {

      all_preds   <- NULL
      hp_cols_all <- character()

      # ---- Gather predictions across folds ----
      for (f in seq_len(n_folds)) {
        if (has_params) {
          for (p in seq_len(n_params)) {
            preds <- all_loaded[[f]][[p]][[m]][[1]]
            hp    <- all_loaded[[f]][[p]][[m]][[2]]

            if (is.null(all_preds)) {
              all_preds <- preds
            } else {
              all_preds <- dplyr::bind_rows(all_preds, preds)
            }

            hp_cols_all <- union(hp_cols_all, hp)
          }
        } else {
          preds <- all_loaded[[f]][[m]][[1]]
          hp    <- all_loaded[[f]][[m]][[2]]

          if (is.null(all_preds)) {
            all_preds <- preds
          } else {
            all_preds <- dplyr::bind_rows(all_preds, preds)
          }

          hp_cols_all <- union(hp_cols_all, hp)
        }
      }

      rownames(all_preds) <- NULL

      if(has_params){
        # Add any extra columns if present
        extra_cols <- setdiff(
          names(all_preds),
          c("rowIndex", "Resample", "obs", "pred", "no", "yes", hp_cols_all)
        )
        if (length(extra_cols) > 0) {
          hp_cols_all <- c(hp_cols_all, extra_cols)
        }
      }

      # ---- Compute metrics per resample ----
      results_matrix <- all_preds %>%
        dplyr::group_by(dplyr::across(dplyr::all_of(hp_cols_all)), Resample) %>%
        dplyr::summarise(
          metrics = list(calculate_accuracy_kappa_resample(obs, pred)),
          .groups = "drop"
        ) %>%
        tidyr::unnest_wider(metrics) %>%
        dplyr::group_by(dplyr::across(dplyr::all_of(hp_cols_all))) %>%
        dplyr::summarise(
          Accuracy   = median(.data$Accuracy_resample),
          Kappa      = median(.data$Kappa_resample),
          AccuracySD = stats::mad(.data$Accuracy_resample),
          KappaSD    = stats::mad(.data$Kappa_resample),
          .groups = "keep"
        )

      # ---- Select best hyperparameters ----
      best_row <- results_matrix %>%
        dplyr::ungroup() %>%
        dplyr::arrange(dplyr::desc(Accuracy)) %>%
        dplyr::slice_max(Accuracy, n = 1, with_ties = FALSE)

      besttune <- best_row %>% dplyr::select(dplyr::all_of(hp_cols_all))

      # ---- Compute resample summaries for besttune ----
      resample_df <- all_preds %>%
        dplyr::inner_join(besttune, by = hp_cols_all) %>%
        dplyr::group_by(Resample) %>%
        dplyr::summarise(
          metrics = list(calculate_accuracy_kappa_resample(obs, pred)),
          .groups = "drop"
        ) %>%
        tidyr::unnest_wider(metrics) %>%
        dplyr::rename(
          Accuracy = Accuracy_resample,
          Kappa = Kappa_resample
        ) %>%
        dplyr::select(Accuracy, Kappa, Resample) %>%
        dplyr::arrange(Resample)

      results[[m]] <- list(
        Prediction_folds = all_preds,
        Results_folds    = results_matrix,
        bestTune         = besttune,
        Resample_matrix  = resample_df
      )
    }

    return(results)
  }

  # ============================================================
  # === SURVIVAL ===============================================
  # ============================================================
  if (task == "survival") {

    for (m in seq_len(n_models)) {

      all_preds   <- NULL
      hp_cols_all <- character()

      # ---- Gather predictions ----
      for (f in seq_len(n_folds)) {
        if (has_params) {
          for (p in seq_len(n_params)) {
            preds <- all_loaded[[f]][[p]][[m]]
            hp_cols <- setdiff(names(preds), c("predictions", "Resample", "rowIndex", "model", "c_index"))
            all_preds <- dplyr::bind_rows(all_preds, preds)
            hp_cols_all <- union(hp_cols_all, hp_cols)
          }
        } else {
          preds <- all_loaded[[f]][[m]]
          hp_cols <- setdiff(names(preds), c("predictions", "Resample",  "rowIndex", "model", "c_index"))
          all_preds <- dplyr::bind_rows(all_preds, preds)
          hp_cols_all <- union(hp_cols_all, hp_cols)
        }
      }

      rownames(all_preds) <- NULL

      if(ncol(all_preds)!=0){
        # ---- Summarize performance per hyperparameter configuration ----
        results_matrix <- all_preds %>%
          dplyr::group_by(dplyr::across(dplyr::all_of(hp_cols_all))) %>%
          dplyr::summarise(
            c_index_median = stats::median(c_index, na.rm = TRUE),
            c_index_mad    = stats::mad(c_index, constant = 1, na.rm = TRUE),
            .groups = "drop" ## keep only c_index and mad
          ) %>%
          dplyr::arrange(dplyr::desc(c_index_median))

        # ---- Select best configuration ----
        best_row <- results_matrix %>%
          dplyr::slice_max(c_index_median, n = 1, with_ties = FALSE)

        besttune <- best_row %>% dplyr::select(dplyr::all_of(hp_cols_all))

        # ---- Compute per-fold for best config ----
        if(length(hp_cols_all)>0){
          resample_df <- all_preds %>%
            dplyr::inner_join(besttune, by = hp_cols_all) %>%
            dplyr::arrange(Resample)
        }else{
          resample_df = all_preds %>%
            dplyr::arrange(Resample)
        }

        # ---- Store results ----
        results[[m]] <- list(
          Prediction_folds = all_preds,
          Results_folds    = results_matrix,
          bestTune         = besttune,
          Resample_matrix  = resample_df
        )
      }else{
        results[[m]] <- NULL
      }
    }

    return(results)
  }
}


### Helper function: get tuneGrid for models --->TO DO: might need to re-think how to choose how many values are gonna be evaluated
get_tune_grid = function(method, train_data){
  set.seed(123)

  if(method == "glmnet"){
    return(expand.grid(alpha = c(0,1), lambda = seq(0.001, 1, length = 20)))
  }
  if(method == "lasso"){
    return(expand.grid(alpha = 1, lambda = seq(0.001, 1, length = 20)))
  }
  if(method == "ridge"){
    return(expand.grid(alpha = 0, lambda = seq(0.001, 1, length = 20)))
  }
  if(method == "rf"){
    n_features <- ncol(train_data) - 1
    return(data.frame(mtry = unique(round(seq(n_features * 0.2, n_features * 0.9, length.out = 3)))))
  }
  if(method == "svmRadial"){
    # Typical small-to-moderate RBF widths + modest C range
    return(expand.grid(
      sigma = c(0.01, 0.05, 0.1),
      C     = c(0.5, 1, 2)
    ))
  }
  if(method == "treebag"){
    return(data.frame(parameter = "none"))
  }
  if(method == "C5.0"){
    return(expand.grid(
      trials = c(1, 5, 10),
      model  = "tree",
      winnow = c(TRUE, FALSE)
    ))
  }
  if(method == "knn"){
    # Odd ks to avoid ties; small-to-moderate neighborhood sizes
    return(expand.grid(k = c(3, 5, 7, 9, 11)))
  }
  if(method == "rpart"){
    # Coarse cp sweep across low/med/high regularization
    return(expand.grid(cp = seq(0.001, 0.1, length = 10)))
  }
  if(method == "svmLinear"){
    return(expand.grid(C = c(0.25, 0.5, 1, 2, 4)))
  }
  if(method == "xgbTree"){
    return(expand.grid(
      nrounds          = c(100, 300, 500),
      max_depth        = c(3, 6, 9),
      eta              = c(0.01, 0.1, 0.3),
      gamma            = 0,            # fixed default
      colsample_bytree = 0.8,          # fixed default
      min_child_weight = 1,            # fixed default
      subsample        = 0.8           # fixed default
    ))
  }

  stop(paste("No grid defined for method:", method))
}

# get_tune_grid = function(method, train_data){ #----> DEPRECATE cause it calls the training recursive times (with big train data) and it fall the stack protection from R
#   set.seed(123)
#
#   if(method == "glmnet"){
#     return(expand.grid(alpha = c(0,1), lambda = seq(0.001, 1, length = 20)))
#   }
#   if(method == "lasso"){
#     return(expand.grid(alpha = 1, lambda = seq(0.001, 1, length = 20)))
#   }
#   if(method == "ridge"){
#     return(expand.grid(alpha = 0, lambda = seq(0.001, 1, length = 20)))
#   }
#   if(method == "rf"){
#     grid_func <- caret::getModelInfo("rf")[["rf"]]$grid
#
#     grid_rf <- grid_func(
#       x = train_data[, -which(names(train_data) == "target")],
#       y = train_data$target,
#       len = 5 ### Think about module this as a parameter
#     )
#
#     # Filter and adjust numeric hyperparameters (avoid higher hyper than number of features)
#     n_features <- ncol(train_data[, -which(names(train_data) == "target")])
#
#     # Replace only values greater than n_features
#     for (param in names(grid_rf)) {
#       if (is.numeric(grid_rf[[param]])) {
#         invalid_idx <- which(unique(grid_rf[[param]]) >= n_features)
#         if (length(invalid_idx) > 0) {
#           # Identify the pattern in the original values
#           original_values <- unique(grid_rf[[param]])
#           pattern_length <- length(original_values)
#
#           # Build a replacement range within valid limits (20% - 90%)
#           replacements <- unique(round(seq(n_features * 0.2, n_features * 0.9, length.out = length(invalid_idx))))
#
#           # Replace values in hyperparam df
#           for(k in seq_along(invalid_idx)){ # If there are more than one value to replace
#             old_value <- grid_rf[[param]][invalid_idx[k]]
#             new_value <- replacements[k]
#             grid_rf[[param]][grid_rf[[param]] == old_value] <- new_value
#           }
#         }
#       }
#     }
#     return(grid_rf)
#   }
#   if(method == "svmRadial"){
#     grid_func <- caret::getModelInfo("svmRadial")[["svmRadial"]]$grid
#
#     ### Remove zero variance column
#     nzv <- caret::nearZeroVar(train_data, saveMetrics = TRUE)
#     train_data <- train_data[, !nzv$zeroVar, drop = FALSE]
#
#     tuneGrid <- grid_func(
#       x = train_data[, -which(names(train_data) == "target")],
#       y = train_data$target,
#       len = 5 ### Think about module this as a parameter
#     )
#     return(tuneGrid)
#   }else{
#     grid_func <- caret::getModelInfo(method)[[method]]$grid
#
#     tuneGrid <- grid_func(
#       x = train_data[, -which(names(train_data) == "target")],
#       y = train_data$target,
#       len = 5 ### Think about module this as a parameter
#     )
#     return(tuneGrid)
#   }
#
# }

#' Get default hyperparameter grids for supported survival models
#'
#' This function generates default hyperparameter grids for various survival models
#' compatible with the tidymodels framework. It uses the default parameter ranges
#' from the {dials} package and produces a sequence of evenly spaced values
#' within those ranges for each tunable hyperparameter.
#'
#' The function supports models such as Cox proportional hazards (regular and penalized),
#' parametric survival regression, decision trees, bagging, random forests, and gradient boosting.
#' For models without tunable parameters (e.g., classic Cox or AFT models),
#' the function returns \code{NULL}.
#'
#' @param model_name Character string specifying the model name.
#'   Supported options include:
#'   \itemize{
#'     \item \code{"cox_ph_survival"} – Classic Cox proportional hazards model
#'     \item \code{"proportional_hazards_glmnet"} – Penalized Cox (LASSO / Elastic Net)
#'     \item \code{"survreg_flexsurv"} – Parametric accelerated failure time (AFT)
#'     \item \code{"decision_tree_partykit"} – Single survival tree
#'     \item \code{"bag_tree_rpart"} – Bagged CART survival trees
#'     \item \code{"rand_forest_partykit"} – Random survival forest (ctree-based)
#'     \item \code{"rand_forest_aorsf"} – Oblique random survival forest
#'     \item \code{"boost_tree_mboost"} – Gradient boosting for survival
#'   }
#' @param train_x Optional data frame or matrix of training predictors.
#'   This is required for parameters that depend on the number of features,
#'   such as \code{mtry}, which determines the number of variables randomly
#'   sampled at each split in tree-based models.
#' @param levels Integer specifying how many values to generate per hyperparameter.
#'   Defaults to \code{5}. Must be at least 2.
#' @param v Integer. Number of folds for K-fold cross-validation (default = 5).
#'
#' @return A named list of hyperparameter grids.
#'   Each element is a numeric vector of sampled values for that parameter.
#'   Returns \code{NULL} for models without tunable hyperparameters.
#'
#' @details
#' The helper function \code{vs()} internally calls \code{dials::value_seq()}
#' to generate evenly spaced sequences of parameter values across their default ranges.
#' For data-dependent parameters (like \code{mtry}), \code{dials::finalize()}
#' is used to compute appropriate limits based on \code{train_x}.
#'
#' @examples
#' # Example 1: Random forest with feature-dependent hyperparameter
#' get_default_hyperparams("rand_forest_partykit", train_x = iris[, 1:4])
#'
#' # Example 2: Penalized Cox model with default parameter ranges
#' get_default_hyperparams("proportional_hazards_glmnet")
#' @export
get_default_hyperparams <- function(model_name, train_x = NULL, levels = 5, v = 5) {
  stopifnot(levels >= 2)

  # --- helpers -----------------------------------------------------------
  vs <- function(param) dials::value_seq(param, n = levels)  # evenly spaced across default range
  p  <- if (!is.null(train_x)) ncol(train_x) else NA_integer_ #number of columns (features)
  n_total <- if (!is.null(train_x)) nrow(train_x) else NA_integer_ #number of rows (samples)
  max_leaf <- if (!is.na(n_total)) floor(floor(n_total * (v - 1)/v) / 2) - 1 else NA_integer_

  # --- model-specific parameter grids -----------------------------------
  if (model_name == "cox_ph_survival") {
    return(NULL)

  } else if (model_name == "proportional_hazards_glmnet") {
    return(list(
      penalty = vs(dials::penalty()),
      mixture = vs(dials::mixture())
    ))

  } else if (model_name == "survreg_flexsurv") {
    return(NULL)

  } else if (model_name == "decision_tree_partykit") {
    vals <- list(
      cost_complexity = vs(dials::cost_complexity()),
      tree_depth      = vs(dials::tree_depth()),
      min_n           = vs(dials::min_n())
    )
    if (!is.na(max_leaf)) vals$min_n <- pmin(vals$min_n, max_leaf)
    return(vals)

  } else if (model_name == "bag_tree_rpart") {
    return(list(
      trees = vs(dials::trees())
    ))

  } else if (model_name %in% c("rand_forest_partykit", "rand_forest_aorsf")) {
    min_n_vals <- as.numeric(vs(dials::min_n()))
    # Cap at fold maximum
    if (!is.na(max_leaf)) min_n_vals <- pmin(min_n_vals, max_leaf)

    # mtry sequence
    if (!is.na(p)) mtry_vals <- as.numeric(vs(dials::finalize(dials::mtry(), train_x)))
    if (!is.na(p)) mtry_vals <- pmin(mtry_vals, p)

    vals <- list(
      trees = as.numeric(vs(dials::trees())),
      min_n = min_n_vals,
      mtry = mtry_vals
    )
    return(vals)

  } else if (model_name == "boost_tree_mboost") {
    return(list(
      trees          = vs(dials::trees()),
      min_n          = vs(dials::min_n()),
      tree_depth     = vs(dials::tree_depth()),
      learn_rate     = vs(dials::learn_rate()),
      loss_reduction = vs(dials::loss_reduction()),
      sample_size    = vs(dials::sample_size(range = c(0, 1))),
      stop_iter      = vs(dials::stop_iter())
    ))

  } else {
    stop("Unsupported model name: ", model_name)
  }
}



unregister_dopar <- function() {
  if (!is.null(foreach::getDoParRegistered())) {
    # switch back to sequential backend
    foreach::registerDoSEQ()
    gc()
  }
}

#' Train and Evaluate a Survival Model
#'
#' This function trains and evaluates a single survival model using the
#' **tidymodels** framework. It supports a range of survival model types,
#' including Cox, penalized Cox, AFT, random forests, bagged trees, and
#' gradient boosting models.
#'
#' The function automatically applies user-specified hyperparameters,
#' standardizes predictions from different engines, and evaluates model
#' performance using the **Concordance Index (C-index)** as the primary metric.
#' It ensures a consistent risk-score direction across models (higher values
#' indicate higher risk).
#'
#' @param df_train A data frame containing the training data, including
#'   survival time, event indicator, and predictor variables.
#' @param df_test A data frame containing the test data with the same structure
#'   and columns as `df_train`.
#' @param outcome_col Character string specifying the name of the survival
#'   time column.
#' @param event_col Character string specifying the name of the event indicator
#'   column (`1 = event occurred`, `0 = censored`).
#' @param model Character string specifying which survival model to train.
#'   Supported values include:
#'   \itemize{
#'     \item `"cox_ph_survival"` – Cox proportional hazards model ({survival})
#'     \item `"proportional_hazards_glmnet"` – Penalized Cox model (LASSO/Elastic Net)
#'     \item `"survreg_flexsurv"` – Parametric AFT model ({flexsurv})
#'     \item `"rand_forest_partykit"` – Random survival forest (ctree engine)
#'     \item `"rand_forest_aorsf"` – Oblique random survival forest
#'     \item `"decision_tree_partykit"` – Single survival tree
#'     \item `"bag_tree_rpart"` – Bagged survival trees
#'     \item `"boost_tree_mboost"` – Gradient boosting for survival data
#'   }
#' @param models_hyperparameters A list containing model hyperparameter values
#'   to apply. Typically created from a tuning grid or optimization step, e.g.:
#'   `list(list(trees = 500, min_n = 10))`. If `NULL`, default parameters are used.
#'
#' @return A tibble with two columns:
#' \describe{
#'   \item{`model`}{The model name as a character string.}
#'   \item{`c_index`}{Numeric value of the computed C-index on the test data.}
#' }
#'
#' @details
#' The function uses the **parsnip** interface from tidymodels to define,
#' train, and evaluate models.
#'
#' Workflow:
#' \enumerate{
#'   \item Defines the model specification (`parsnip::model_spec`).
#'   \item Optionally applies user-specified hyperparameters via
#'         `parsnip::set_args()`.
#'   \item Constructs a survival formula of the form `Surv(time, event) ~ .`.
#'   \item Fits the model using `parsnip::fit()` on the training data.
#'   \item Evaluates model performance using
#'         [predict_and_evaluate_survival()], which computes the C-index.
#' }
#'
#' @seealso [predict_and_evaluate_survival()]
#' @export
#'
compute_ml_survival <- function(df_train, df_test = NULL,
                                outcome_col, event_col,
                                model, models_hyperparameters, return_model = F){
  # ---------------------------------------------------------------------------
  # Define the model specification based on the chosen model
  # ---------------------------------------------------------------------------
  # Each model is created using the {parsnip} unified syntax.
  # Mode is always set to "censored regression" for survival analysis.
  # The engine determines the computational backend:
  #   - "survival" → Classical Cox or AFT
  #   - "glmnet"   → Penalized regression (LASSO/Elastic Net)
  #   - "partykit" → Tree-based survival (ctree)
  #   - "aorsf"    → Oblique random survival forest
  #   - "rpart"    → Bagged CART trees
  #   - "mboost"   → Gradient boosting for censored data
  #
  if (model == "cox_ph_survival") {
    # Cox proportional hazards (no tunable parameters)
    model_spec <- parsnip::proportional_hazards(
      mode = "censored regression",
      engine = "survival"
    )

  } else if (model == "proportional_hazards_glmnet") {
    # Penalized Cox (elastic net): tune penalty (λ) and mixture (α)
    model_spec <- parsnip::proportional_hazards(
      penalty = tune(),
      mixture = tune(),
      mode = "censored regression",
      engine = "glmnet"
    )

  } else if (model == "survreg_flexsurv") {
    # Parametric AFT model
    model_spec <- parsnip::survival_reg(
      mode = "censored regression",
      engine = "flexsurv"
    )

  } else if (model == "rand_forest_partykit") {
    # Random survival forest via conditional inference framework
    model_spec <- parsnip::rand_forest(
      trees = tune(),
      mode = "censored regression",
      engine = "partykit"
    )

  } else if (model == "rand_forest_aorsf") {
    # Oblique random survival forest (efficient implementation)
    model_spec <- parsnip::rand_forest(
      trees = tune(),
      mode = "censored regression",
      engine = "aorsf"
    )

  } else if (model == "decision_tree_partykit") {
    # Single survival decision tree
    model_spec <- parsnip::decision_tree(
      mode = "censored regression",
      engine = "partykit"
    )

  } else if (model == "bag_tree_rpart") {
    # Bagged CART-based survival trees
    model_spec <- parsnip::bag_tree(
      mode = "censored regression",
      engine = "rpart"
    )

  } else if (model == "boost_tree_mboost") {
    # Gradient boosting for censored data
    model_spec <- parsnip::boost_tree(
      trees = tune(),
      mode = "censored regression",
      engine = "mboost"
    )

  } else {
    # Unsupported model name
    stop("Unsupported model: ", model)
  }


  # ---------------------------------------------------------------------------
  # Apply user-defined hyperparameters (if provided)
  # ---------------------------------------------------------------------------
  # Hyperparameters are applied dynamically using do.call() and set_args(),
  # allowing flexible argument passing from a parameter grid.
  #
  # Example:
  #   models_hyperparameters = list(list(trees = 500, min_n = 10))
  #
  if (!is.null(models_hyperparameters)) {
    model_spec <- do.call(
      parsnip::set_args,
      c(list(object = model_spec), models_hyperparameters[[1]])
    )
  }

  # ---------------------------------------------------------------------------
  # Define workflow and fit the model
  # ---------------------------------------------------------------------------
  # Combine the model and survival formula into a workflow.
  # Formula format: Surv(time, event) ~ predictors
  # Ensures compatibility with survival engines during fitting.
  #
  formula_model <- stats::as.formula(
    paste0("Surv(", outcome_col, ", ", event_col, ") ~ .")
  )

  wf <- workflows::workflow() %>%
    workflows::add_model(model_spec) %>%
    workflows::add_formula(formula_model)

  # Fit model on training data
  fitted <- tryCatch(
    parsnip::fit(wf, data = df_train),
    error = function(e) {
      warning(paste("Model fitting failed for", model, ":", e$message))
      return(NULL)
    }
  )

  # If fitting failed, return NULL
  if (is.null(fitted)) {
    return(NULL)
  }

  # ---------------------------------------------------------------------------
  # Evaluate performance on test data (if given)
  # ---------------------------------------------------------------------------
  # Uses predict_and_evaluate_survival() to compute the Concordance Index (C-index),
  # which measures the model’s ability to correctly rank survival outcomes.
  #
  if(!is.null(df_test)){ ## Return train model + predictions
    prediction = predict_and_evaluate_survival(fitted, df_test, outcome_col, event_col)
    preds = as.numeric(unlist(prediction$preds))
    metric_val = prediction$c_index
    if(return_model){
      res = list(Model = fitted, Metrics = tibble::tibble(predictions = preds, c_index = metric_val))
    }else{
      res = tibble::tibble(predictions = preds, c_index = metric_val)
    }
  }else{ ## Only return train model
    res = fitted
  }


  return(res)
}


#' Nested Cross-Validation for Survival Models with Optional Custom Fold Construction
#'
#' Performs *nested cross-validation* to evaluate and tune multiple survival
#' models using the **tidymodels** ecosystem. Supports both standard
#' event-stratified cross-validation and *Leave-One-Domain-Out (LODO)* setups,
#' enabling cohort-balanced model evaluation. Hyperparameter grids are
#' automatically constructed for each model type.
#'
#' Depending on the inputs, the function can:
#' \enumerate{
#'   \item Build folds internally or accept custom folds from a user-defined function.
#'   \item Train survival models with or without hyperparameter tuning.
#'   \item Compute and aggregate the Concordance Index (C-index) across folds.
#'   \item Identify and retrain the top-performing model using optimal parameters.
#' }
#'
#' @param df_features A data frame of predictor variables (features).
#' @param df_outcome A data frame containing survival outcomes — typically including
#'   survival time and event indicator columns.
#' @param outcome_col Character string giving the name of the survival time column.
#' @param event_col Character string giving the name of the event indicator column
#'   (`0 = censored`, `1 = event`).
#' @param k_folds Integer. Number of folds for K-fold cross-validation (default = 5).
#' @param n_rep Integer. Number of repeated CV iterations (default = 1).
#' @param ncores Integer. Number of CPU cores to use for parallelization.
#' @param LODO Logical; if `TRUE`, performs Leave-One-Domain-Out cross-validation
#'   using `batch_id` to stratify samples by cohort.
#' @param batch_id Optional character string naming the column representing cohort or
#'   batch identifiers. Required if `LODO = TRUE`.
#' @param file_name Optional string specifying the suffix for the generated C-index
#'   summary PDF saved in the `"Results/"` directory.
#' @param fold_construction_fun Optional custom function for constructing data folds.
#'   Used to interface with external preprocessing workflows (e.g., CellTFusion).
#' @param fold_construction_args_fixed Optional list of fixed arguments passed to
#'   `fold_construction_fun()`.
#' @param fold_construction_args_tunable Optional list of tunable arguments passed to
#'   `fold_construction_fun()` during hyperparameter tuning.
#'
#' @details
#' Internally, the function:
#' \itemize{
#'   \item Merges predictors and outcomes into a single dataset.
#'   \item Creates stratified folds using **rsample**, either by event rate or by
#'     cohort × event combinations (if `LODO = TRUE`).
#'   \item Evaluates a predefined set of survival models:
#'     Cox PH, penalized Cox (glmnet), AFT (flexsurv), decision trees, bagged trees,
#'     and random forests.
#'   \item Aggregates the median and MAD of the C-index across resamples.
#'   \item Retrains the best-performing model with its optimal hyperparameters.
#' }
#'
#' When a custom fold construction function is provided via `fold_construction_fun`,
#' the function handles folds in parallel, saves intermediate results under
#' `"Results/"`, and returns additional outputs for advanced integration.
#'
#' @return A named list with the following elements:
#' \describe{
#'   \item{`Model`}{The best-performing survival model retrained on the full dataset.}
#'   \item{`ML_Models`}{All evaluated survival models with aggregated C-index results.}
#'   \item{`C_index_median`}{Median C-index of the top-performing model.}
#'   \item{`Custom_output`}{Optional list of custom outputs from fold construction.}
#' }
#'
#' @export
#'
compute_k_fold_CV_survival <- function(df_features, df_outcome, outcome_col, event_col, k_folds, n_rep, ncores,
                                       LODO = FALSE, batch_id = NULL, file_name = NULL, fold_construction_fun = NULL,
                                       fold_construction_args_fixed = NULL,
                                       fold_construction_args_tunable = NULL){

  # ---------------------------------------------------------------------------
  # Step 2: Merge features and outcomes
  # ---------------------------------------------------------------------------
  # Combine predictor variables (df_features) with survival outcomes (df_outcome)
  # to ensure consistent ordering before creating cross-validation folds.
  #
  df_all <- df_features %>%
    dplyr::bind_cols(df_outcome) ## df with time and event columns

  # ---------------------------------------------------------------------------
  # Step 3: Define cross-validation folds
  # ---------------------------------------------------------------------------
  # Generate stratified K-fold CV partitions using {rsample}.
  #   - If LODO = FALSE → standard stratification by event rate.
  #   - If LODO = TRUE  → stratify by cohort × event to preserve domain balance.
  #
  if (isTRUE(LODO)) {
    if (is.null(batch_id) || !batch_id %in% names(df_all)) {
      stop("When LODO = TRUE, you must provide 'batch_id' as a column name in df_all.")
    }

    batch_col <- batch_id

    # Create composite stratification variable: cohort × event
    df_all <- df_all %>%
      dplyr::mutate(strata = interaction(.data[[batch_col]], .data[[event_col]], drop = TRUE))

    folds <- rsample::vfold_cv(
      df_all,
      v = k_folds,
      repeats = n_rep,
      strata = "strata"
    )

  } else {
    cat("Creating stratified v-fold CV (stratified by event rate only)\n")

    folds <- rsample::vfold_cv(
      df_all,
      v = k_folds,
      repeats = n_rep,
      strata = all_of(event_col) # stratify by event indicator only
    )
  }

  ### -------------------------------------------------------------------------
  ### Step 4: Reformat folds to caret-compatible list format
  ### -------------------------------------------------------------------------
  # Convert the rsample object into a list of training indices, enabling reuse
  # of legacy caret-compatible aggregation utilities.
  #
  multifolds <- purrr::map(
    seq_len(nrow(folds)),
    ~ folds$splits[[.x]]$in_id
  ) %>%
    purrr::set_names(folds$id)

  # ---------------------------------------------------------------------------
  # Step 5: Define candidate survival models
  # ---------------------------------------------------------------------------
  # Each model is a tidymodels-compatible specification supporting censored data.
  #
  model_list <- c(
    "cox_ph_survival",              # classical Cox proportional hazards model
    "proportional_hazards_glmnet",  # penalized Cox (LASSO / elastic net)
    "survreg_flexsurv",             # parametric AFT regression
    "decision_tree_partykit",       # single survival tree
    "bag_tree_rpart",               # bagged decision trees
    #"rand_forest_partykit",        # random survival forest (ctree engine)
    "rand_forest_aorsf"             # oblique random survival forest
    #"boost_tree_mboost"            # gradient boosting for survival (optional)
  )

  # ---------------------------------------------------------------------------
  # Step 6: Generate default hyperparameter grids
  # ---------------------------------------------------------------------------
  # Construct model-specific grids (via get_default_hyperparams()).
  # Each element corresponds to a model with tunable parameter ranges.
  #
  model_grids <- purrr::map(model_list, ~ get_default_hyperparams(.x, train_x = df_features, v = k_folds))
  names(model_grids) <- model_list

  # Initialize master list to store everything in memory
  models_all_folds <- vector("list", k_folds*n_rep)

  # ---------------------------------------------------------------------------
  # Step 7: Run nested cross-validation and model training
  # ---------------------------------------------------------------------------
  # Depending on whether fold_construction_fun is defined, this section either:
  #   - Performs standard nested CV (train/test splits within this function), or
  #   - Executes an external custom fold construction pipeline (e.g., CellTFusion).
  #
  # Each fold produces performance results stored for later aggregation.
  if(is.null(fold_construction_fun)){
    for (fold_i in seq_along(multifolds)) { ### number of folds (k_fold x n_rep)
      cat("Running fold_i", fold_i, "\n")
      train_idx <- multifolds[[fold_i]]
      test_idx <- setdiff(seq_len(nrow(df_all)), train_idx)

      train_data <- df_all[train_idx, , drop = FALSE]
      test_data <- df_all[test_idx, , drop = FALSE]

      #Preprocessing features (remove collinear variables and no-variance)
      train_data <- preprocess_features(train_data, cor_thresh = 0.9, time_var = "time", event_var = "event")

      common_features <- intersect(
        colnames(train_data),
        colnames(test_data)
      )

      test_data <- test_data[, common_features, drop = FALSE]

      # Run all ML methods for this parameter configuration
      models <- lapply(
        model_list,
        function(method) {
          # Retrieve hyperparameter grid for the current model
          hyperparams <- model_grids[[method]]
          cat("Running ", method, "\n")
          # Fallback: if the model has no tunable parameters, create a single configuration
          if (is.null(hyperparams)) {
            param_grid <- tibble::tibble(.config_id = 1)
          } else {
            param_grid <- tidyr::expand_grid(!!!hyperparams) %>%
              dplyr::mutate(.config_id = row_number())
          }

          model_results <- purrr::map_dfr(1:nrow(param_grid), function(g) {

            current_params <- param_grid[g, , drop = FALSE]

            # Train model and evaluate C-index on the held-out fold
            trained <- compute_ml_survival(df_train = train_data,
                                           df_test = test_data,
                                           outcome_col = outcome_col, event_col = event_col,
                                           model = method,
                                           models_hyperparameters = if (is.null(hyperparams)) NULL else list(
                                             current_params %>% dplyr::select(-.config_id)))

            if(!is.null(trained)){
              trained_df <- trained %>%
                data.frame() %>%
                dplyr::mutate(
                  model = method,
                  Resample = paste0("Fold", fold_i),
                  rowIndex = test_idx,
                ) %>%
                dplyr::bind_cols(current_params %>% dplyr::select(-.config_id))%>%
                dplyr::relocate(c_index, .after = dplyr::last_col())
            }

          })

        }
      )

      # Store all results for this fold
      models_all_folds[[fold_i]] <- models

    }

    models = aggregate_results(models_all_folds, task = 'survival')
    names(models) <- model_list

    ## Sanity check (each param conf has to be evaluated in all resamples)
    for(i in 1:length(models)){
      if(!is.null(models[[i]])){
        hp_cols_all = names(models[[i]][["bestTune"]]) ### Hyperparameter names
        x = models[[i]][["Prediction_folds"]] %>%
          dplyr::distinct(Resample, dplyr::across(all_of(hp_cols_all))) %>%
          dplyr::count(dplyr::across(all_of(hp_cols_all)), name = "n_resamples") %>%
          dplyr::arrange(desc(n_resamples))

        # Expected number of resamples (folds × repeats)
        expected_resamples <- length(unique(models[[i]][["Prediction_folds"]]$Resample))

        # Sanity check
        if (any(x$n_resamples != expected_resamples)) {
          stop("Inconsistent number of resamples detected for parameter configuration\n",
               paste0(hp_cols_all, collapse = " "))
        }
      }
    }

    custom_outputs = NULL

  }else{
    # Custom fold construction (is running in parallel)
    do.call(fold_construction_fun, c(list(data = df_features, folds = multifolds), fold_construction_args_fixed, fold_construction_args_tunable))

    ### Extract the file names of the folds
    result_files <- list.files("Results", pattern = "^fold_.*\\.rds$", full.names = TRUE)

    # Iterate across folds and inside each subfold corresponding to each param combination (if exist)
    for (fold_i in seq_along(result_files)) { ### number of folds (k_fold x n_rep)

      result = readRDS(result_files[[fold_i]]) ## per resample

      # Each fold contains multiple parameter sets (list of lists) --> fold_construction_args_tunable != NULL
      if (!is.null(fold_construction_args_tunable)) {
        models_all_params <- vector("list", length(result))

        cl <- parallel::makeCluster(ncores)
        doParallel::registerDoParallel(cl)

        models_all_params <- foreach::foreach(parameter_i = seq_along(result),
                                              .packages = c("dplyr", "caret", "censored")) %dopar% {

                                                source("~/Documents/pipeML/R/machine_learning.R")

                                                train_data_i <- result[[parameter_i]][["train_data"]]
                                                test_data_i  <- result[[parameter_i]][["test_data"]]

                                                # Preprocessing features (remove collinear variables and no-variance)
                                                train_data_i <- preprocess_features(train_data_i, cor_thresh = 0.9,
                                                                                    time_var = "time", event_var = "event")

                                                # Replace in original train/test datasets
                                                result[[parameter_i]][["train_data"]] <- train_data_i

                                                common_features <- intersect(
                                                  colnames(train_data_i),
                                                  colnames(test_data_i)
                                                )

                                                result[[parameter_i]][["test_data"]] <- test_data_i[, common_features, drop = FALSE]

                                                # Run all ML methods for this parameter configuration
                                                models <- lapply(
                                                  model_list,
                                                  function(method) {
                                                    cat("Running model", method, "with param configuration", parameter_i, "and fold", fold_i, "\n")

                                                    # Retrieve hyperparameter grid for the current model
                                                    hyperparams <- model_grids[[method]]

                                                    # Fallback: if the model has no tunable parameters, create a single configuration
                                                    if (is.null(hyperparams)) {
                                                      param_grid <- tibble::tibble(.config_id = 1)
                                                    } else {
                                                      param_grid <- tidyr::expand_grid(!!!hyperparams) %>%
                                                        dplyr::mutate(.config_id = row_number())
                                                    }

                                                    model_results <- purrr::map_dfr(1:nrow(param_grid), function(g) {

                                                      current_params <- param_grid[g, , drop = FALSE]

                                                      # Train model and evaluate C-index on the held-out fold
                                                      trained <- compute_ml_survival(df_train = result[[parameter_i]][["train_data"]],
                                                                                     df_test = result[[parameter_i]][["test_data"]],
                                                                                     outcome_col = outcome_col, event_col = event_col,
                                                                                     model = method,
                                                                                     models_hyperparameters = if (is.null(hyperparams)) NULL else list(
                                                                                       current_params %>% dplyr::select(-.config_id)))

                                                      trained_df <- trained %>%
                                                        data.frame() %>%
                                                        dplyr::mutate(
                                                          model = method,
                                                          Resample = paste0("Fold", fold_i)
                                                        ) %>%
                                                        dplyr::bind_cols(current_params %>% dplyr::select(-.config_id))%>%
                                                        dplyr::bind_cols(result[[parameter_i]][["params"]]) %>%
                                                        dplyr::relocate(c_index, .after = dplyr::last_col())

                                                    })

                                                  }
                                                )

                                                models
                                              }

        parallel::stopCluster(cl)  # stop the cluster after parallel execution
        unregister_dopar() #Stop Dopar from running in the background

        # Store all parameter results for this fold
        models_all_folds[[fold_i]] <- models_all_params

      }else { # Custom function does not have hyperparams to tune --> fold_construction_args_tunable != NULL
        train_data_i <- result[["train_data"]]
        test_data_i <- result[["test_data"]]

        # Preprocessing
        train_data_i <- preprocess_features(train_data_i, cor_thresh = 0.9,
                                            time_var = "time", event_var = "event")

        # Replace
        result[["train_data"]] = train_data_i
        result[["test_data"]] = test_data_i[, setdiff(colnames(train_data_i), "target")]

        # Run all ML methods for this parameter configuration
        models <- lapply(
          model_list,
          function(method) {

            # Retrieve hyperparameter grid for the current model
            hyperparams <- model_grids[[method]]

            # Fallback: if the model has no tunable parameters, create a single configuration
            if (is.null(hyperparams)) {
              param_grid <- tibble::tibble(.config_id = 1)
            } else {
              param_grid <- tidyr::expand_grid(!!!hyperparams) %>%
                dplyr::mutate(.config_id = row_number())
            }

            model_results <- purrr::map_dfr(1:nrow(param_grid), function(g) {

              current_params <- param_grid[g, , drop = FALSE]

              # Train model and evaluate C-index on the held-out fold
              trained <- compute_ml_survival(df_train = result[["train_data"]],
                                             df_test = result[["test_data"]],
                                             outcome_col = outcome_col, event_col = event_col,
                                             model = method,
                                             models_hyperparameters = if (is.null(hyperparams)) NULL else list(
                                               current_params %>% dplyr::select(-.config_id)))

              trained_df <- trained %>%
                data.frame() %>%
                dplyr::mutate(
                  model = method,
                  Resample = paste0("Fold", fold_i)
                ) %>%
                dplyr::bind_cols(current_params %>% dplyr::select(-.config_id))%>%
                dplyr::relocate(c_index, .after = dplyr::last_col())

            })

          }
        )

        models_all_folds[[fold_i]] <- models

      }
    }

    # Step 8: Aggregate results and validate fold consistency
    # ---------------------------------------------------------------------------
    # Combines all model results using aggregate_results().
    # Ensures that each hyperparameter configuration was evaluated across all folds.
    models = aggregate_results(models_all_folds, task = 'survival')
    names(models) <- model_list

    ## Sanity check (each param conf has to be evaluated in all resamples)
    for(i in 1:length(models)){
      hp_cols_all = names(models[[i]][["bestTune"]]) ### Hyperparameter names
      x = models[[i]][["Prediction_folds"]] %>%
        dplyr::distinct(Resample, dplyr::across(all_of(hp_cols_all))) %>%
        dplyr::count(dplyr::across(all_of(hp_cols_all)), name = "n_resamples") %>%
        dplyr::arrange(desc(n_resamples))

      # Expected number of resamples (folds × repeats)
      expected_resamples <- length(unique(models[[i]][["Prediction_folds"]]$Resample))

      # Sanity check
      if (any(x$n_resamples != expected_resamples)) {
        stop("Inconsistent number of resamples detected for parameter configuration\n",
             paste0(hp_cols_all, collapse = " "))
      }
    }
  }

  if(!is.null(fold_construction_fun)){
    if (!is.null(fold_construction_args_tunable)){

      ################################ Train model with optimized hyperparameters

      optimized_models <- lapply(seq_along(model_list), function(i) {
        cat("\nRunning model...", model_list[i], "\n")
        wrapper_train_best_hyperparams_survival(
          train_data = df_features,
          optimized = models[[i]],
          ml_method = model_list[i],
          fold_construction_fun,
          fold_construction_args_fixed
        )
      })

      # Split components across lists
      training_sets <- lapply(optimized_models, `[[`, "training_set")
      custom_outputs <- lapply(optimized_models, `[[`, "custom_output")
      models           <- lapply(optimized_models, `[[`, "Model")

      # Assign pretty names
      names(training_sets) <- model_list
      names(custom_outputs) <- model_list
      names(models) <- model_list
    }else{

      optimized_models <- lapply(seq_along(model_list), function(i) {
        cat("\nRunning model...", model_list[i], "\n")

        temp = models[[i]]$bestTune
        models[[i]]$bestTune = c(temp, fold_construction_args_fixed)

        p = wrapper_train_best_hyperparams_survival(
            train_data = df_features,
            optimized = models[[i]],
            ml_method = model_list[i],
            fold_construction_fun,
            fold_construction_args_fixed
          )

        models[[i]]$bestTune = temp

        p

      })

      # Split components across lists
      training_sets <- lapply(optimized_models, `[[`, "training_set")
      custom_outputs <- lapply(optimized_models, `[[`, "custom_output")
      models           <- lapply(optimized_models, `[[`, "Model")

      # Assign pretty names
      names(training_sets) <- model_list
      names(custom_outputs) <- model_list
      names(models) <- model_list
    }
  }

  # Step 9: Retrain best model on full data (if tuning performed)
  # ---------------------------------------------------------------------------
  # Once optimal hyperparameters are found per model, retrain the top model
  # using wrapper_train_best_hyperparams_survival().

  #Top model with best accuracy or AUC
  metrics <- compute_cv_CINDEX(models, plot_results = TRUE, file_name = file_name)

  top_model = metrics[["Top_model"]]
  c_index_median = metrics[["CINDEX_summary"]] %>%
    dplyr::filter(model == top_model) %>%
    dplyr::pull(Median_CINDEX)

  model_metrics = models[[top_model]]

  cat("Best ML model found: ", top_model, "\n")

  ############# Re-train only 'best model' to return ML model using tuned parameters
  model_metrics[["Model_object"]] <- compute_ml_survival(
                                        df_train = df_all,
                                        outcome_col = "time",
                                        event_col   = "event",
                                        model = top_model,
                                        models_hyperparameters = list(model_metrics$bestTune)
                                      )

  cat("Returning model trained\n")

  output = list("Model" = model_metrics, "ML_Models" = models, "C_index_median" = c_index_median)

  if(!is.null(custom_outputs) && !any(sapply(custom_outputs, is.null))){ #Check whether custom_output exists or not
    output[[length(output)+1]] = custom_outputs[[top_model]]
    names(output)[length(output)] = "Custom_output"
  }

  return(output)
}

#' Summarize and Visualize C-index Results from Survival Model Cross-Validation
#'
#' This function aggregates and visualizes the C-index (concordance index)
#' results obtained from cross-validation of multiple survival models.
#' Each model should contain a `Resample_matrix` element with per-fold C-index
#' values. The function computes the median and MAD (median absolute deviation)
#' of the C-index for each model, identifies the top-performing model, and
#' optionally generates a bar plot summarizing model performance.
#'
#' @param models A named list of survival model objects, where each element
#'   corresponds to one fitted model. Each model must contain a
#'   `Resample_matrix` data frame with columns:
#'   \describe{
#'     \item{`c_index`}{C-index value per resample (numeric).}
#'     \item{`Resample`}{Fold or resample identifier (e.g., "Fold1", "Fold2").}
#'   }
#' @param file_name Optional character string used to name the output PDF file
#'   saved under `"Results/CINDEX_CV_methods_<file_name>.pdf"`. If `NULL`, the
#'   file is not named explicitly.
#' @param plot_results Logical; if `TRUE` (default), generates and saves a PDF
#'   bar plot showing median C-index ± MAD per model.
#'
#' @details
#' For each model:
#' \itemize{
#'   \item The **median C-index** represents the typical discrimination
#'     performance across folds.
#'   \item The **MAD (Median Absolute Deviation)** measures variability
#'     in C-index values (robust equivalent of standard deviation).
#' }
#'
#' The plot helps visualize and compare models based on survival prediction
#' performance, with error bars representing ± MAD.
#'
#' @return
#' A list with the following elements:
#' \describe{
#'   \item{`CINDEX_summary`}{A tibble summarizing median and MAD per model.}
#'   \item{`All_folds`}{A tibble with raw C-index values from all folds and models.}
#'   \item{`Top_model`}{Character string naming the model with the highest median C-index.}
#' }
#'
#' @examples
#' \dontrun{
#' models <- list(
#'   cox_ph_survival = list(
#'     Resample_matrix = tibble::tibble(
#'       c_index = c(0.72, 0.75, 0.70),
#'       Resample = c("Fold1", "Fold2", "Fold3")
#'     )
#'   ),
#'   rand_forest_aorsf = list(
#'     Resample_matrix = tibble::tibble(
#'       c_index = c(0.80, 0.82, 0.78),
#'       Resample = c("Fold1", "Fold2", "Fold3")
#'     )
#'   )
#' )
#' compute_cv_CINDEX(models, file_name = "example")
#' }
#'
#' @seealso [aggregate_results()], [predict_and_evaluate_survival()]
#' @export
#'
compute_cv_CINDEX <- function(models, file_name = NULL, plot_results = TRUE){
  # Step 1: Extract C-index values per fold from all models
  # ---------------------------------------------------------------------------
  # Each element in `models` must contain $Resample_matrix with columns:
  #   - c_index: numeric value per resample (e.g., per fold)
  #   - Resample: fold identifier (e.g., "Fold1", "Fold2", ...)
  #
  res_cindex <- purrr::map_df(names(models), function(m) {
    model_obj <- models[[m]]
    if (!is.null(model_obj$Resample_matrix) && "c_index" %in% colnames(model_obj$Resample_matrix)) {
      tibble::tibble(
        model = m,
        c_index = model_obj$Resample_matrix$c_index,
        Resample = model_obj$Resample_matrix$Resample
      )
    } else {
      tibble::tibble(model = m, c_index = NA_real_, Resample = NA_character_)
    }
  })

  # Remove missing values (in case some models failed)
  res_cindex <- res_cindex %>% dplyr::filter(!is.na(c_index))

  # ---------------------------------------------------------------------------
  # Step 2: Summarize performance statistics per model
  # ---------------------------------------------------------------------------
  summary_cindex <- res_cindex %>%
    dplyr::group_by(model) %>%
    dplyr::summarise(
      Median_CINDEX = stats::median(c_index, na.rm = TRUE),
      MAD_CINDEX = stats::mad(c_index, constant = 1, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    dplyr::arrange(dplyr::desc(Median_CINDEX))

  # ---------------------------------------------------------------------------
  # Step 3: Optional plotting
  # ---------------------------------------------------------------------------
  if (plot_results) {

    grDevices::pdf(paste0("Results/CINDEX_CV_methods_", file_name, ".pdf"), width = 10)

    print(
      ggplot2::ggplot(summary_cindex,
                      ggplot2::aes(x = model, y = Median_CINDEX)) +
        ggplot2::geom_bar(
          stat = "identity",
          width = 0.9,
          fill = "#1f78b4"
        ) +
        ggplot2::geom_errorbar(
          ggplot2::aes(
            ymin = Median_CINDEX - MAD_CINDEX,
            ymax = Median_CINDEX + MAD_CINDEX
          ),
          width = 0.15
        ) +
        ggplot2::labs(
          title = "Cross-Validation performance (C-index)",
          x = "ML Model",
          y = "Median C-index ± MAD"
        ) +
        ggplot2::scale_x_discrete(expand = ggplot2::expansion(mult = c(0, 0))) +
        ggplot2::scale_y_continuous(breaks = seq(0, 1, by = 0.05)) +
        ggplot2::theme_minimal(base_size = 16) +
        ggplot2::theme(
          legend.position = "none",
          axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, size = 14),
          axis.text.y = ggplot2::element_text(size = 14),
          axis.title = ggplot2::element_text(size = 16),
          plot.title = ggplot2::element_text(size = 18, face = "bold")
        )
    )

    grDevices::dev.off()
  }

  # ---------------------------------------------------------------------------
  # Step 4: Identify top-performing model
  # ---------------------------------------------------------------------------
  top_model <- summary_cindex %>%
    dplyr::slice_max(Median_CINDEX, n = 1, with_ties = FALSE) %>%
    dplyr::pull(model)

  # ---------------------------------------------------------------------------
  # Step 5: Return structured summary
  # ---------------------------------------------------------------------------
  return(list(CINDEX_summary = summary_cindex, All_folds = res_cindex, Top_model = top_model))
}

#' Plot and save survival performance of a model on test data
#'
#' This function groups individuals into risk strata (e.g., Low/Medium/High)
#' based on their predicted risk scores from a fitted survival model.
#' It then plots Kaplan–Meier survival curves for each risk group,
#' including a log-rank test for group separation and displays the
#' concordance index (C-index) of the model on the test data.
#'
#' @param df_test A data frame containing at least the following columns:
#'   \describe{
#'     \item{time}{Observed survival or follow-up time (numeric).}
#'     \item{event}{Event indicator (1 = event occurred, 0 = censored).}
#'     \item{.pred}{Predicted risk score or linear predictor from the model
#'       (higher values indicate higher risk).}
#'   }
#' @param c_index Numeric. The concordance index (C-index) computed on the test data.
#' @param n_groups Integer. Number of risk groups to stratify by (default = 3).
#'   Typically 3 groups correspond to "Low", "Medium", and "High" risk strata.
#' @param file_name Character (optional). If provided, the Kaplan–Meier plot will be
#'   saved to \code{"Results/Survival_KM_<file_name>.pdf"}.
#'
#' @details
#' Risk groups are defined by quantile-based cut points on the predicted risk scores.
#' The log-rank test is used to assess whether survival curves differ significantly
#' between risk strata. The C-index is displayed for interpretability.
#'
#' @return Invisibly returns the \code{ggsurvplot} object for further customization,
#'   and saves a PDF of the plot in the "Results/" directory if \code{file_name} is provided.
#'
#' @examples
#' \dontrun{
#' plot_survival_performance(
#'   df_test = test_data,
#'   c_index = 0.74,
#'   n_groups = 3,
#'   file_name = "cox_model_test"
#' )
#' }
#'
#' @export
plot_survival_performance <- function(df_test,
                                      prediction,
                                      n_groups = 3,
                                      file_name = NULL) {

  # ---------------------------------------------------------------------------
  # Add predictions to test data
  # ---------------------------------------------------------------------------

  df_test$.pred <- prediction$preds$.pred

  c_index <- prediction$c_index
  c_low   <- prediction$c_index_lower
  c_high  <- prediction$c_index_upper

  # ---------------------------------------------------------------------------
  # Check required columns
  # ---------------------------------------------------------------------------

  required_cols <- c("time", "event", ".pred")

  if (!all(required_cols %in% names(df_test))) {
    stop("df_test must include columns: time, event, and .pred")
  }

  # ---------------------------------------------------------------------------
  # Risk group labels
  # ---------------------------------------------------------------------------

  labels <- if (n_groups == 2) {
    c("Low risk", "High risk")
  } else if (n_groups == 3) {
    c("Low risk", "Medium risk", "High risk")
  } else {
    paste0("Group ", seq_len(n_groups))
  }

  # ---------------------------------------------------------------------------
  # Create risk groups
  # ---------------------------------------------------------------------------

  df_test <- df_test %>%
    dplyr::mutate(
      risk_group = cut(
        .pred,
        breaks = stats::quantile(
          .pred,
          probs = seq(0, 1, length.out = n_groups + 1),
          na.rm = TRUE
        ),
        include.lowest = TRUE,
        labels = labels
      )
    )

  # ---------------------------------------------------------------------------
  # Kaplan–Meier model
  # ---------------------------------------------------------------------------

  fit_km <- survival::survfit(
    survival::Surv(time, event) ~ risk_group,
    data = df_test
  )

  # ---------------------------------------------------------------------------
  # Log-rank test
  # ---------------------------------------------------------------------------

  logrank <- survival::survdiff(
    survival::Surv(time, event) ~ risk_group,
    data = df_test
  )

  p_val <- 1 - stats::pchisq(logrank$chisq, df = length(logrank$n) - 1)

  # ---------------------------------------------------------------------------
  # Subtitle with C-index + CI
  # ---------------------------------------------------------------------------

  subtitle_text <- paste0(
    "C-index: ",
    round(c_index, 3),
    " (95% CI ",
    round(c_low, 3), "-",
    round(c_high, 3), ")",
    " | Log-rank p = ",
    format.pval(p_val, digits = 3, eps = .001)
  )

  # ---------------------------------------------------------------------------
  # Plot KM curves
  # ---------------------------------------------------------------------------

  # Create ggsurvplot without printing
  plot_obj <- survminer::ggsurvplot(
    fit_km,
    data = df_test,
    risk.table = TRUE,
    conf.int = TRUE,
    pval = FALSE,
    ggtheme = ggplot2::theme_minimal(),
    palette = c("#1B9E77", "#7570B3", "#D95F02")[1:n_groups],
    legend.title = "Risk Group",
    legend.labs = levels(df_test$risk_group),
    title = paste0("Test-set Survival Performance_", file_name),
    subtitle = subtitle_text,
    risk.table.height = 0.25,
    tables.theme = ggplot2::theme_minimal(),
    print = FALSE
  )

  # For a single plot, no need to use arrange_ggsurvplots()
  grDevices::pdf(paste0("Results/Survival_KM_", file_name, ".pdf"), width = 8, height = 6)
  print(plot_obj)   # prints main plot + risk table together
  grDevices::dev.off()

}

#' Predict and Evaluate Survival Model Performance
#'
#' This function generates predictions from a fitted survival model and evaluates
#' its performance using the Concordance Index (C-index). It automatically
#' handles different prediction output types supported by various survival
#' modeling engines and standardizes predictions into a comparable numeric format.
#'
#' @param model_fit A fitted survival model object (typically created via
#'   a `parsnip` or `workflow` model specification).
#' @param data A data frame containing the predictors and survival outcome
#'   variables used for prediction and evaluation.
#' @param outcome_col Character string specifying the name of the survival time
#'   column in `data`. Default is `"time"`.
#' @param event_col Character string specifying the name of the event indicator
#'   column in `data`. Default is `"event"`.
#'
#' @details
#' The function attempts to generate predictions using multiple possible
#' prediction types, depending on what the model supports:
#' \itemize{
#'   \item `"linear_pred"` — Linear predictor or log hazard (e.g., Cox models),
#'     where higher values indicate higher risk.
#'   \item `"time"` — Expected survival time (e.g., AFT models),
#'     where higher values imply longer survival (automatically reversed).
#'   \item `"survival"` — Survival probability at a specific evaluation time,
#'     where higher values indicate better survival (automatically reversed).
#' }
#'
#' It standardizes the prediction output into a tibble with one numeric column
#' `.pred`, ensuring consistency across engines.
#' The function then computes the Concordance Index (C-index) using
#' `yardstick::concordance_survival_vec()` to assess how well the model ranks
#' predicted survival relative to actual outcomes.
#'
#' @return
#' A list with two elements:
#' \describe{
#'   \item{`preds`}{A tibble containing standardized numeric predictions.}
#'   \item{`c_index`}{Numeric value representing the computed C-index.}
#' }
#'
#' @examples
#' \dontrun{
#' fitted_model <- parsnip::fit(
#'   parsnip::proportional_hazards(engine = "survival"),
#'   Surv(time, event) ~ .,
#'   data = lung
#' )
#' results <- predict_and_evaluate_survival(fitted_model, lung)
#' results$c_index
#' }
#'
#' @seealso [yardstick::concordance_survival_vec()], [aggregate_results()]
#' @export
#'
predict_and_evaluate_survival <- function(model_fit,
                                          data,
                                          outcome_col = NULL,
                                          event_col = NULL) {


  # ---------------------------------------------------------------------------
  # Generate predictions on test data
  # ---------------------------------------------------------------------------
  # The code tries different prediction "types" depending on what the fitted model supports.
  #
  # 1. "linear_pred"  → Risk score or log hazard (e.g., Cox model)
  #      - Higher values = higher risk (shorter survival)
  #      - Direction: higher = worse outcome
  #
  # 2. "time"         → Expected survival time (e.g., parametric AFT models)
  #      - Higher values = longer expected lifetime
  #      - Direction: higher = better outcome (will be reversed later)
  #
  # 3. "survival"     → Survival probability at a given evaluation time
  #      - Higher values = more likely to survive (less risk)
  #      - Direction: higher = better outcome (will be reversed later)
  #
  # The try–catch structure ensures the code works for any survival model:
  # - Try risk-based prediction first
  # - If not supported, try time-based prediction
  # - If that fails, fall back to survival probability at a fixed eval_time
  # ---------------------------------------------------------------------------
  pred_type = NULL

  preds <- tryCatch({
    pred_type = "linear_pred"
    predict(model_fit, new_data = data, type = "linear_pred")
  }, error = function(e1) {
    tryCatch({
      pred_type = "time"
      predict(model_fit, new_data = data, type = "time")
    }, error = function(e2) {
      # Only try "survival" prediction if outcome/event columns exist
      if (!is.null(outcome_col) && !is.null(event_col)) {
        eval_time <- stats::median(data[[outcome_col]], na.rm = TRUE)
        pred_type <- "survival"
        predict(model_fit, new_data = data, type = "survival", eval_time = eval_time)
      } else {
        message("Outcome or event column is NULL; skipping survival-type prediction")
        return(NULL)
      }
    })
  })

  # ---------------------------------------------------------------------------
  # Standardize prediction output into a numeric tibble
  # ---------------------------------------------------------------------------
  # Prediction outputs vary by engine:
  #   - Some return numeric vectors (risk scores)
  #   - Others return matrices (e.g., survival curves across time)
  #   - Some may return lists of probabilities
  #
  # This section standardizes all prediction outputs into a tibble with a
  # single numeric column `.pred`, computing the median when multiple
  # predictions per sample exist.
  #
  if ("matrix" %in% class(preds[[1]])) {
    preds <- tibble::tibble(.pred = apply(as.data.frame(preds[[1]]), 1, median, na.rm = TRUE))
  } else if (!is.numeric(preds[[1]])) {
    preds <- tibble::tibble(.pred = apply(as.matrix(preds), 1, median, na.rm = TRUE))
  } else {
    pred_col <- names(preds)[1]
    preds <- preds %>% dplyr::rename(.pred = dplyr::all_of(pred_col))
  }

  # ---------------------------------------------------------------------------
  # Ensure direction consistency (higher = higher risk)
  # ---------------------------------------------------------------------------

  # Reverse predictions if type implies "more = better survival"
  if (pred_type %in% c("time", "survival")) {
    preds$.pred <- -preds$.pred
  }

  # ---------------------------------------------------------------------------
  # Evaluate model performance using the Concordance Index (C-index) only if outcome/event are available
  # ---------------------------------------------------------------------------
  # The C-index measures how well the model ranks survival times relative to
  # true outcomes (i.e., discrimination ability). It ranges from 0.5 (random)
  # to 1 (perfect discrimination).
  #
  # The function `censored::concordance_survival_vec()` computes the metric.
  # If the computation fails, the value is safely returned as NA.
  #

  if (!is.null(outcome_col) && !is.null(event_col)) {

    df_eval <- data
    df_eval$.pred <- preds$.pred

    cindex_res <- compute_cindex_ci(
      data = df_eval,
      time_col = outcome_col,
      event_col = event_col,
      pred_col = ".pred"
    )

    metric_val <- cindex_res$c_index
    ci_lower <- cindex_res$CI_lower
    ci_upper <- cindex_res$CI_upper

  } else {

    metric_val <- NA_real_
    ci_lower <- NA_real_
    ci_upper <- NA_real_

  }


  # ---------------------------------------------------------------------------
  # Return the model name and computed C-index
  # ---------------------------------------------------------------------------
  # The function returns a tibble containing the model name and the computed
  # C-index, making it easy to aggregate and compare across models and folds.
  #
  return(list(preds = preds, c_index = metric_val,
              c_index_lower = ci_lower, c_index_upper = ci_upper))

}

bootstrap_auc <- function(predict, target, method, B = 1000, seed = 123){

  set.seed(seed)
  n <- length(target)

  auroc_vals <- numeric(B)
  auprc_vals <- numeric(B)

  for(b in seq_len(B)){

    idx <- sample(seq_len(n), replace = TRUE)

    predict_b <- predict[idx, , drop = FALSE]
    target_b  <- target[idx]

    sens_spec_b <- get_sensitivity_specificity(
      predictions = predict_b,
      observed = target_b,
      ml.model = method
    )

    auroc_vals[b] <- calculate_auroc(
      sens_spec_b$fpr,
      sens_spec_b$Sensitivity
    )

    auprc_vals[b] <- calculate_auprc(
      sens_spec_b$Recall,
      sens_spec_b$Precision
    )
  }

  list(
    AUROC = list(
      mean  = mean(auroc_vals),
      lower = quantile(auroc_vals, 0.025),
      upper = quantile(auroc_vals, 0.975),
      values = auroc_vals
    ),
    AUPRC = list(
      mean  = mean(auprc_vals),
      lower = quantile(auprc_vals, 0.025),
      upper = quantile(auprc_vals, 0.975),
      values = auprc_vals
    )
  )
}

compute_cindex_ci <- function(data, time_col = "time", event_col = "event",
                              pred_col = ".pred", n_boot = 1000, seed = 123){

  set.seed(seed)

  n <- nrow(data)

  # Observed C-index
  c_index_obs <- yardstick::concordance_survival_vec(
    truth = survival::Surv(data[[time_col]], data[[event_col]]),
    estimate = data[[pred_col]]
  )

  # Bootstrap distribution
  cindex_boot <- replicate(n_boot, {

    idx <- sample(seq_len(n), replace = TRUE)

    boot_data <- data[idx, ]

    yardstick::concordance_survival_vec(
      truth = survival::Surv(boot_data[[time_col]], boot_data[[event_col]]),
      estimate = boot_data[[pred_col]]
    )

  })

  # Confidence interval
  ci <- stats::quantile(cindex_boot, probs = c(0.025, 0.975), na.rm = TRUE)

  tibble::tibble(
    c_index = c_index_obs,
    CI_lower = ci[1],
    CI_upper = ci[2]
  )
}

