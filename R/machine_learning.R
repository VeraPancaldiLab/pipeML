
dir.create("Results", showWarnings = FALSE, recursive = TRUE)

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
compute.boruta <- function(data, seed, fix = TRUE) {

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
#' @param data A data frame with the column "target" as the response and other columns as features.
#' @param iterations Integer. The number of Boruta iterations to perform.
#' @param fix Logical. If TRUE, applies TentativeRoughFix() to resolve tentative features after each iteration.
#' @param doParallel Logical. Whether to use parallel processing.
#' @param workers Integer. Number of CPU cores to use for parallel execution. If NULL, uses all available cores minus one.
#' @param file_name A string for naming output plots and CSV files saved in the "Results/" directory.
#' @param threshold A numeric value between 0 and 1. A feature must be confirmed in more than \code{threshold * iterations} to be finally labeled as confirmed.
#' @param return Logical. Whether to save the resulting plots in the "Results/" directory.
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
#' @examples
#' res_boruta <- feature.selection.boruta(
#'   data = training_set,
#'   iterations = 10,
#'   fix = FALSE,
#'   doParallel = TRUE,
#'   workers = 4,
#'   threshold = 0.8,
#'   file_name = "Test",
#'   return = FALSE
#' )
feature.selection.boruta <- function(data, iterations = NULL, fix, doParallel = F, workers=NULL, file_name = NULL, threshold = NULL, return) {
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

      message("Running ", iterations, " iterations of the Boruta algorithm using ", num_cores, " cores")

      res <- foreach::foreach(seed = sample.int(100000, iterations)) %dopar% {

        library(pipeML)

        tryCatch({
          # If successful, return the result and the seed
          list(result = compute.boruta(data, seed, fix),
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
      message("Running ", iterations, " iterations of the Boruta algorithm")
      res = list()
      for (i in 1:iterations) {
        res[[i]] = compute.boruta(data, seed = sample.int(100000, 1), fix)
      }

      # Extract the first sublist of each element
      matrix_of_importance <- lapply(res, function(x) x[[1]])

      # Extract the second sublist of each element
      features_labels <- lapply(res, function(x) x[[2]])

      res = merge_boruta_results(matrix_of_importance, features_labels, file_name = file_name, iterations = iterations, threshold = threshold, return = return)
    }

  }

  return(res)

}

#' Perform repeated stratified k-fold cross-validation for model training and tuning
#'
#' This function performs repeated stratified k-fold cross-validation on a dataset to train and tune hyperparameters for 13 machine learning methods. Optionally, it can also perform model stacking and Boruta-based feature selection. Performance is evaluated using user-specified metrics such as Accuracy, AUROC, or AUPRC.
#'
#' @param model A data frame containing features and a target column named 'target' corresponding to the response variable to predict.
#' @param k_folds Integer. Number of folds for k-fold cross-validation. Default is 5.
#' @param n_rep Integer. Number of repetitions of the k-fold cross-validation. Default is 100.
#' @param stacking Logical. Whether to perform model stacking. Default is FALSE.
#' @param metric Character. Metric used for hyperparameter tuning and model evaluation. Supported values are "Accuracy", "AUROC", and "AUPRC".
#' @param boruta Logical. Whether to apply Boruta for feature selection before model training. Note that many ML models handle feature importance internally, so prior selection is optional unless multicollinearity is a concern. Default is FALSE.
#' @param boruta_iterations Integer. Number of iterations to run Boruta. Since Boruta involves randomness, repeated runs improve consistency. Default is 100.
#' @param fix_boruta Logical. Whether to fix Boruta’s internal parameters. See `compute.boruta()` for details.
#' @param tentative Logical. Whether to include tentative features as confirmed in the training dataset.
#' @param boruta_threshold Numeric. Threshold for confirming features after multiple Boruta iterations. For example, 0.8 means features must be confirmed in at least 80% of iterations. Default is 0.8.
#' @param file_name Character. File name used for saving output plots in the `Results/` directory.
#' @param LODO Logical. If TRUE, performs Leave-One-Dataset-Out (LODO) cross-validation by stratifying folds based on cohort membership.
#' @param return Logical. Whether to return the results and generated plots.
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
#' @export
#'
#' @examples
#' training <- compute.k_fold_CV(
#'   train_data,
#'   k_folds = 5,
#'   n_rep = 100,
#'   metric = "Accuracy",
#'   stacking = TRUE,
#'   boruta = FALSE,
#'   boruta_iterations = 100,
#'   fix_boruta = FALSE,
#'   boruta_threshold = 0.8,
#'   file_name = "Test",
#'   return = TRUE
#' )
compute.k_fold_CV = function(model, k_folds, n_rep, stacking = FALSE, metric = "Accuracy", boruta, boruta_iterations = NULL, fix_boruta = NULL, tentative = FALSE, boruta_threshold = NULL, file_name = NULL, LODO = FALSE, return = FALSE){

  if(!(metric %in% c("AUROC", "AUPRC","Accuracy"))){
    stop("The metric assigned is not supported. Choose either accuracy or AUC.")
  }

  ######### Feature selection across folds
  if(boruta == T){
    cat("Feature selection using Boruta...............................................................\n\n")
    # Feature selection using Boruta
    res_boruta = feature.selection.boruta(model, iterations = boruta_iterations, fix = fix_boruta, file_name = file_name, doParallel = F, workers=NULL, threshold = boruta_threshold, return = return)

    if(tentative == F){
      if(length(res_boruta$Confirmed) <= 1){ #No enough features selected for training model
        message("No enough features selected for training a model")
        results = list()
        return(results)
      }else{
        cat("\nKeeping only features confirmed in more than 80% of the times for training...............................................................\n\n")
        cat("If you want to consider also tentative features, please specify tentative = TRUEin the parameters.\n\n")
        train_data = model[,colnames(model)%in%res_boruta$Confirmed, drop = F] %>%
          dplyr::mutate(target = model$target)
      }
    }else{
      sum_features = length(res_boruta$Confirmed) + length(res_boruta$Tentative)
      if(sum_features <= 1){
        message("No enough features selected for training a model")
        results = list()
        return(results)
      }else{
        cat("Keeping features confirmed and tentative in more than 80% of the times for training...............................................................\n\n")
        train_data = model[,colnames(model)%in%c(res_boruta$Confirmed, res_boruta$Tentative), drop = F] %>%
          dplyr::mutate(target = model$target)
      }
    }

    rm(res_boruta) #Clean memory
    gc()

  }else{
    train_data = model
  }

  rm(model) #Clean memory
  gc()

  cat("Training machine learning model...............................................................\n\n")

  ######### Machine Learning models

  cl <- parallel::makeCluster(4)
  doParallel::registerDoParallel(cl)

  ######### Stratify K fold cross-validation
  if(LODO == T){
    multifolds = construct_stratified_cohort_folds(train_data, 'dataset', 'target', k_folds = k_folds, n_rep = n_rep)
    train_data = train_data %>% dplyr::select(-dataset) #After creating multifolds we remove this variable to be able to train
  }else{
    multifolds = caret::createMultiFolds(train_data[,'target'], k = k_folds, times = n_rep) #repeated folds
  }

  trainControl <- caret::trainControl(index = multifolds, method="repeatedcv", number=k_folds, repeats=n_rep, verboseIter = F, allowParallel = T, classProbs = TRUE, savePredictions=T)

  ##################################################### ML models
  #To do: Re-calculate accuracy values based on tuning parameters optimized by the cv AUC - now the values are based on accuracy! be careful

  ################## Bagged CART
  fit.treebag <- caret::train(target~., data = train_data, method = "treebag", metric = "Accuracy",trControl = trainControl)

  ################## RF
  require(randomForest)
  fit.rf <- caret::train(target~., data = train_data, method = "rf", metric = "Accuracy",trControl = trainControl)

  ################## C5.0
  require(C50)
  fit.c50 <- caret::train(target~., data = train_data, method = "C5.0", metric = "Accuracy",trControl = trainControl)

  ################## LG - Logistic Regression
  fit.glm <- caret::train(target~., data = train_data, method="glm", metric="Accuracy",trControl=trainControl)

  ################## LDA - Linear Discriminate Analysis
  fit.lda <- caret::train(target~., data = train_data, method="lda", metric="Accuracy",trControl=trainControl)

  ################## GLMNET - Regularized Logistic Regression (Elastic net)
  fit.glmnet <- caret::train(target~., data = train_data, method="glmnet", metric="Accuracy",trControl=trainControl)

  ################## KNN - k-Nearest Neighbors
  fit.knn <- caret::train(target~., data = train_data, method="knn", metric="Accuracy",trControl=trainControl)

  ################## CART - Classification and Regression Trees (CART),
  fit.cart <- caret::train(target~., data = train_data, method="rpart", metric="Accuracy",trControl=trainControl)

  # NB - Naive Bayes (NB)
  #Grid = expand.grid(usekernel=TRUE,adjust=1,fL=c(0,0.2,0.5,0.8,1))
  #fit.nb <- train(target~., data = train_data, method="nb", metric="Accuracy",trControl=trainControl, tuneGrid=Grid)

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

  parallel::stopCluster(cl)  # stop the cluster after parallel execution
  unregister_dopar() #Stop Dopar from running in the background

  #### Set allowParallel = F
  #Disable parallelization in xgbTree cause the xgboost algorithm has its own internal parallelization, controlled by the nthread parameter — it uses this to speed up tree construction within a single model.
  #The caret package, on the other hand, can parallelize between models — for example, it can train different cross-validation folds or hyperparameter combinations at the same time if a parallel backend is registered
  #https://stackoverflow.com/questions/39528392/parallel-processing-with-xgboost-and-caret
  #If both are ON it can slower performance (lead to over-parallelization and CPU contention)
  trainControl <- caret::trainControl(index = multifolds, method="repeatedcv", number=k_folds, repeats=n_rep, verboseIter = F, allowParallel = F, classProbs = TRUE, savePredictions=T)

  fit.xgbTree <- caret::train(target~., data=train_data, method="xgbTree", metric = "Accuracy", trControl=trainControl)

  ####### Optimized based on metric (only AUC or Accuracy available)
  if(metric == "AUROC" || metric == "AUPRC"){

    # Store models in a named list
    models <- list(
      BAG = fit.treebag,
      RF = fit.rf,
      C50 = fit.c50,
      GLM = fit.glm,
      LDA = fit.lda,
      GLMNET = fit.glmnet,
      KNN = fit.knn,
      CART = fit.cart,
      LASSO = fit.lasso,
      RIDGE = fit.ridge,
      SVM_radial = fit.svm_radial,
      SVM_linear = fit.svm_linear,
      XGboost = fit.xgbTree
    )

    # Define corresponding hyperparameters
    hyperparams <- list(
      BAG = NULL,
      RF = "mtry",
      C50 = c("trials", "model", "winnow"),
      GLM = NULL,
      LDA = NULL,
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
      res <- calculate_cv_metrics(model, metric, hyperparams[[name]])
      # Update model fields
      model$pred <- res$Prediction
      model$resample <- res$Resamples
      model$results <- res$Results

      return(model)
    })

    # Restore updated models
    names(models) = names(hyperparams)
    fit.treebag <- models$BAG
    fit.rf <- models$RF
    fit.c50 <- models$C50
    fit.glm <- models$GLM
    fit.lda <- models$LDA
    fit.glmnet <- models$GLMNET
    fit.knn <- models$KNN
    fit.cart <- models$CART
    fit.lasso <- models$LASSO
    fit.ridge <- models$RIDGE
    fit.svm_radial <- models$SVM_radial
    fit.svm_linear <- models$SVM_linear
    fit.xgbTree <- models$XGboost


  }

  ###Prediction with best tuned hyper-parameters (Missing to add platt scaling to calibrated probabilities (when tested it didnt converge, need to be checked)) See https://www.cs.cornell.edu/~alexn/papers/calibration.icml05.crc.rev3.pdf

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

  predictions.glm = stats::predict(fit.glm, newdata = train_data, type = "prob") %>%
    data.frame() %>%
    dplyr::select(yes) %>%
    dplyr::rename(GLM = yes)  #Predictions of model (already ordered)

  ### LDA

  predictions.lda = stats::predict(fit.lda, newdata = train_data, type = "prob") %>%
    data.frame() %>%
    dplyr::select(yes) %>%
    dplyr::rename(LDA = yes)  #Predictions of model (already ordered)

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

  ############################################################## Save models

  ensembleResults <- list(BAG = fit.treebag,
                          RF = fit.rf,
                          C50 = fit.c50,
                          GLM = fit.glm,
                          LDA = fit.lda,
                          KNN = fit.knn,
                          CART = fit.cart,
                          GLMNET = fit.glmnet,
                          LASSO = fit.lasso,
                          RIDGE = fit.ridge,
                          SVM_radial = fit.svm_radial,
                          SVM_linear = fit.svm_linear,
                          XGboost = fit.xgbTree)

  ml_methods = list(BAG = "treebag",
                    RF = "rf",
                    C50 = "C5.0",
                    GLM = "glm",
                    LDA = "lda",
                    KNN = "knn",
                    CART = "rpart",
                    GLMNET = "glmnet",
                    LASSO = "glmnet",
                    RIDGE = "glmnet",
                    SVM_radial = "svmRadial",
                    SVM_linear = "svmLinear",
                    XGboost = "xgbTree")

  ############################################# These predictions are use for the meta-learner because it needs the predictions from the models in the complete dataset (might change in the future)
  model_predictions = list(BAG = predictions.bag,
                           RF = predictions.rf,
                           C50 = predictions.c50,
                           GLM = predictions.glm,
                           LDA = predictions.lda,
                           KNN = predictions.knn,
                           CART = predictions.cart,
                           GLMNET = predictions.glmnet,
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
  rm(fit.treebag, fit.rf, fit.c50, fit.glm, fit.lda, fit.knn, fit.cart, fit.glmnet, fit.lasso, fit.ridge, fit.svm_radial, fit.svm_linear, multifolds)
  gc()


  if(stacking){
    features = colnames(train_data)[colnames(train_data) != "target"]

    #Base models using ML models with best accuracy or AUC from each family
    if(metric == "Accuracy"){
      base_models = compute_cv_accuracy(ensembleResults, base_models = T, file_name = file_name, return = return)
    }else if(metric == "AUROC" || metric == "AUPRC"){
      base_models = compute_cv_AUC(ensembleResults, base_models = T, file_name = file_name, AUC_type = metric, return = return)
      # if(return == T){
      #   plot_cv_metrics(ensembleResults, file_name = file_name)
      # }
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

    meta_learner <- caret::train(true_label ~ ., data = meta_features, method = "glmnet", trControl = trainControl) #Staking based on simple logistic regression

    output = list("Features" = features, "Meta_learner" = meta_learner, "Base_models" = base_models$Base_models, "ML_models" = ensembleResults)

  }else{
    features = colnames(train_data)[colnames(train_data) != "target"] #Extract features used for model training

    #Top model with best accuracy or AUC
    if(metric == "Accuracy"){
      metrics = compute_cv_accuracy(ensembleResults, file_name = file_name, return = return)
    }else if(metric == "AUROC" || metric == "AUPRC"){
      metrics = compute_cv_AUC(ensembleResults, file_name = file_name, AUC_type = metric, return = return)
    }

    top_model = metrics[["Top_model"]]

    model = ensembleResults[[top_model]]

    cat("Best ML model found: ", top_model, "\n")

    cat("Returning model trained\n")

    output = list("Features" = features, "Model" = model, "ML_Models" = ensembleResults)
  }


  return(output)

}


#' Train machine learning models with optional stacking and feature selection
#'
#' This function trains one or more machine learning models using repeated k-fold cross-validation, with optional model stacking and feature selection using Boruta. It supports stratified cross-validation, including Leave-One-Dataset-Out (LODO) validation when cohort information is available.
#'
#' @param features_train A data frame containing the features used for training.
#' @param target_var A vector containing the target variable to predict.
#' @param trait.positive Value in \code{target_var} to be considered as the positive class.
#' @param metric Character. Metric used for hyperparameter tuning and model selection. Supported values are \code{"Accuracy"}, \code{"AUROC"}, and \code{"AUPRC"}.
#' @param stack Logical. Whether to perform model stacking. Default is \code{FALSE}.
#' @param k_folds Integer. Number of folds to use in cross-validation.
#' @param n_rep Integer. Number of repetitions of the cross-validation.
#' @param feature.selection Logical. Whether to apply Boruta feature selection before model training. Default is \code{FALSE}.
#' @param seed Integer. Random seed for reproducibility.
#' @param LODO Logical. If \code{TRUE}, constructs folds stratified by cohorts (Leave-One-Dataset-Out CV).
#' @param batch_id A vector indicating the cohort or batch for each sample (required only if \code{LODO = TRUE}).
#' @param file_name Character. File name used to save plots in the \code{Results/} directory.
#' @param return Logical. Whether to return and save the plots generated by the function.
#'
#' @return A list containing:
#' \itemize{
#'   \item Trained model (or meta-learner if \code{stack = TRUE})
#'   \item Features used in model training (all features if \code{feature.selection = FALSE})
#' }
#'
#' @export
#'
#' @examples
#' res <- compute.features.training.ML(
#'   x, y,
#'   trait.positive = "R",
#'   metric = "AUROC",
#'   stack = FALSE,
#'   k_folds = 10,
#'   n_rep = 5,
#'   feature.selection = FALSE,
#'   seed = 123,
#'   LODO = FALSE,
#'   batch_id = NULL,
#'   file_name = "Test",
#'   return = FALSE
#' )
compute.features.training.ML = function(features_train, target_var, trait.positive, metric = "Accuracy", stack, k_folds = 10, n_rep = 5, feature.selection = FALSE, seed, LODO = FALSE, batch_id = NULL, file_name = NULL, return = FALSE){

  set.seed(seed)

  #Set training set
  train_data = features_train %>%
    data.frame() %>%
    dplyr::mutate(Trait = target_var,
                  target = as.factor(ifelse(Trait == trait.positive, 'yes', 'no'))) %>%
    dplyr::select(-Trait)

  train_data$target <- factor(train_data$target, levels = c("no", "yes"))  # Order, just in case to ensure positive class is not well defined

  if(LODO == T){
    train_data = train_data %>%
      dplyr::mutate(dataset = batch_id)
  }

  #Cross-validation training (5 k-folds and 100 repetitions)
  training = compute.k_fold_CV(train_data, k_folds = k_folds, n_rep = n_rep, metric = metric, stacking = stack, boruta = feature.selection, boruta_iterations = 100, fix_boruta = FALSE, boruta_threshold = 0.8, file_name = file_name, LODO = LODO, return= return)

  ####################################################Predicting
  if(length(training)!=0){
    features = training[["Features"]] #Save selected features per partition (only useful if Boruta = TRUE- needs to be improve it)
    if(stack){
      model = training[["Meta_learner"]]
    }else{
      model = training[["Model"]] #Save best ML model based on the Accuracy/AUC from CV per partition
    }

    return(list(Model = training, Features = features))
  }else{  #No features are selected as predictive

    message("No features selected as predictive after Boruta runs. No model returned.")

    return(NULL)
  }

}

#' Train and evaluate machine learning models with optional stacking and feature selection
#'
#' This function trains machine learning models using cross-validation on training data and evaluates them on test data. It supports feature selection with Boruta, model stacking, and cohort-based (LODO) validation.
#'
#' @param features_train A data frame of features used for training the models.
#' @param features_test A data frame of features used for testing the models.
#' @param clinical A data frame containing clinical information, including the target variable and optionally a batch ID. Row names must match the sample identifiers in \code{features_train} and \code{features_test}.
#' @param trait Character. The name of the column in \code{clinical} corresponding to the target variable.
#' @param trait.positive Value in \code{trait} to be considered as the positive class.
#' @param metric Character. Metric used for hyperparameter tuning and model selection. Supported values are \code{"Accuracy"}, \code{"AUROC"}, and \code{"AUPRC"}.
#' @param stack Logical. Whether to apply model stacking. Default is \code{FALSE}.
#' @param k_folds Integer. Number of folds for cross-validation.
#' @param n_rep Integer. Number of cross-validation repetitions.
#' @param feature.selection Logical. Whether to apply Boruta feature selection before training. Default is \code{FALSE}.
#' @param seed Integer. Random seed for reproducibility.
#' @param LODO Logical. If \code{TRUE}, folds are constructed in a Leave-One-Dataset-Out (LODO) manner based on cohorts.
#' @param batch_id A vector indicating the cohort/batch for each sample (only required if \code{LODO = TRUE}).
#' @param file_name Character. Base name used to save plots in the \code{Results/} directory.
#' @param return Logical. Whether to return and save plots generated by the function.
#'
#' @return A list containing:
#' \itemize{
#'   \item Trained model (or meta-learner if stacking was used)
#'   \item Features used in model training
#'   \item Prediction performance metrics
#'   \item AUC scores (AUROC and AUPRC)
#'   \item Predicted class probabilities on test data
#' }
#'
#' @export
#'
#' @examples
#' res <- compute.features.ML(
#'   x, y, clinical,
#'   trait = "Response",
#'   trait.positive = "R",
#'   metric = "AUROC",
#'   stack = FALSE,
#'   k_folds = 10,
#'   n_rep = 5,
#'   feature.selection = FALSE,
#'   seed = 123,
#'   LODO = FALSE,
#'   batch_id = NULL,
#'   file_name = "Test",
#'   return = FALSE
#' )
compute.features.ML = function(features_train, features_test, clinical, trait, trait.positive, metric = "Accuracy", stack, k_folds = 10, n_rep = 5, feature.selection = FALSE, seed, LODO = FALSE, batch_id = NULL, file_name = NULL, return = FALSE){

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

  set.seed(seed)

  #Cross-validation training (5 k-folds and 100 repetitions)
  training = compute.k_fold_CV(train_data, k_folds = k_folds, n_rep = n_rep, metric = metric, stacking = stack, boruta = feature.selection, boruta_iterations = 100, fix_boruta = F, boruta_threshold = 0.8, file_name = file_name, LODO = LODO, return= return)

  ####################################################Predicting
  if(length(training)!=0){
    features = training[["Features"]] #Save selected features per partition (only useful if Boruta = TRUE - needs to be improve it)
    ####################### Testing set

    #Extract target variable
    target = traitData_test %>%
      dplyr::mutate(target = ifelse(traitData_test[,trait] == trait.positive, "yes", "no")) %>%
      dplyr::pull(target)

    target = factor(target, levels = c("no", "yes"))

    if(stack){
      model = training[["Meta_learner"]]
      prediction = compute.prediction.stacked(model, features_test, target, training[["ML_models"]], training[["Base_models"]])

    }else{
      model = training[["Model"]] #Save best ML model based on the Accuracy/AUC from CV per partition
      prediction = compute.prediction(model, features_test, target, file_name, maximize = "Accuracy")
    }

    auc_roc_score = prediction[["AUC"]][["AUROC"]]
    auc_prc_score = prediction[["AUC"]][["AUPRC"]]

    metrics = prediction[["Metrics"]]
    predictions = prediction[["Predictions"]]

    if(return == T){
      get_curves(metrics, "Specificity", "Sensitivity", "Recall", "Precision", "model", auc_roc_score, auc_prc_score, file_name)
    }


    return(list(Model = model, Features = features, Prediction_metrics = metrics, AUC = list(AUROC = auc_roc_score, AUPRC = auc_prc_score), Prediction = predictions))
  }else{  #No features are selected as predictive

    message("No features selected as predictive after Boruta runs. No model returned.")

    return(NULL)
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
#' @export
#'
#' @examples
#' get_pooled_roc_curves(file.name = "Combined_Model", folder_path = "Results/Models/")
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
#' @export
#'
#' @examples
#' get_pooled_boxplots(folder_paths = c("Results/Cohort1", "Results/Cohort2"),
#'                     file_name = "TME_Comparison")
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
#' @export
#'
#' @examples
#' res <- compute_cv_accuracy(models = ml_models, file_name = "MyModels", base_models = TRUE, return = TRUE)
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
      Mean_Accuracy = median(Accuracy),
      SD_Accuracy = sd(Accuracy)
    ) %>%
    dplyr::arrange(desc(Mean_Accuracy))

  top_model = res_accuracy %>%
    dplyr::slice(1) %>%
    dplyr::pull(model)

  if(return){
    grDevices::pdf(paste0("Results/Accuracy_CV_methods_", file_name, ".pdf"), width = 10)

    plot(ggplot2::ggplot(res_accuracy, ggplot2::aes(x = model, y = Mean_Accuracy, fill = model)) +
      ggplot2::geom_bar(stat = "identity", position = ggplot2::position_dodge(), width = 0.6) +
      ggplot2::geom_errorbar(ggplot2::aes(ymin = Mean_Accuracy - SD_Accuracy, ymax = Mean_Accuracy + SD_Accuracy),
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
#' @param save_plot Logical. If \code{TRUE}, generates and saves barplots of median AUROC and AUPRC values (with error bars) to the "Results" directory.
#'
#' @return A list containing:
#' \describe{
#'   \item{\code{AUROC}}{A data frame with median and standard deviation of AUROC values for each model.}
#'   \item{\code{AUPRC}}{A data frame with median and standard deviation of AUPRC values for each model.}
#'   \item{\code{Top_model}}{The name of the model with the highest median value for the selected metric (\code{AUC_type}).}
#'   \item{\code{Base_models}}{(Optional) A character vector of selected base models for stacking, returned if \code{base_models = TRUE}.}
#' }
#'
#' @export
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
  auroc_data = do.call(rbind, auroc)

  #Bind AUPRC values from each model
  auprc = list()
  for (i in 1:length(models)){
    auprc[[i]] = models[[i]]$resample %>% #we use the resample matrix and not directly the results matrix as some have hyperparameters so we will need to define best on the tuned parameter (=more code) - resample matrix is made based on the best tuning
      dplyr::mutate(model = names(models)[i])
    names(auprc)[i] = names(models)[i]
  }
  auprc_data = do.call(rbind, auprc)

  res_auroc <- auroc_data %>%
    dplyr::group_by(model) %>%
    dplyr::summarise(
      Mean_AUROC = median(AUROC),
      SD_AUROC = sd(AUROC)
    ) %>%
    dplyr::arrange(desc(Mean_AUROC))

  res_auprc <- auprc_data %>%
    dplyr::group_by(model) %>%
    dplyr::summarise(
      Mean_AUPRC = median(AUPRC),
      SD_AUPRC = sd(AUPRC)
    ) %>%
    dplyr::arrange(desc(Mean_AUPRC))


  if(return){
    grDevices::pdf(paste0("Results/AUROC_CV_methods_", file_name, ".pdf"), width = 10)
    plot(ggplot2::ggplot(res_auroc, ggplot2::aes(x = model, y = Mean_AUROC, fill = model)) +
           ggplot2::geom_bar(stat = "identity", position = ggplot2::position_dodge(), width = 0.6) +
           ggplot2::geom_errorbar(aes(ymin = Mean_AUROC - SD_AUROC, ymax = Mean_AUROC + SD_AUROC),
                         width = 0.2, position = position_dodge(0.6)) +
           ggplot2::labs(title = "Performance of Models",
                x = "Model",
                y = "Median AUROC") +
           ggplot2::theme_minimal() +
           ggplot2::theme(legend.position = "none") +
           ggplot2::scale_y_continuous(breaks = seq(0, 1, by = 0.05)))
    grDevices::dev.off()

    grDevices::pdf(paste0("Results/AUPRC_CV_methods_", file_name, ".pdf"), width = 10)
    plot(ggplot2::ggplot(res_auprc, ggplot2::aes(x = model, y = Mean_AUPRC, fill = model)) +
           ggplot2::geom_bar(stat = "identity", position = ggplot2::position_dodge(), width = 0.6) +
           ggplot2::geom_errorbar(aes(ymin = Mean_AUPRC - SD_AUPRC, ymax = Mean_AUPRC + SD_AUPRC),
                         width = 0.2, position = ggplot2::position_dodge(0.6)) +
           ggplot2::labs(title = "Performance of Models",
                x = "Model",
                y = "Median AUPRC") +
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

  if(base_models == TRUE){
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
#' @export
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
#' @examples
#' f1_scores = calculate_f1(predictions, target)
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
#' @examples
#' mcc_scores = calculate_mcc(predictions, target)
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
#' @examples
#' prediction_metrics = get_sensitivity_specificity(predictions, target, model)
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
#' @examples
#' sens_spec = get_sensitivity_specificity(predictions, target, model)
#' auroc = calculate_auroc(sens_spec$fpr, sens_spec$sensitivity)
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
#' @examples
#' sens_spec = get_sensitivity_specificity(predictions, target, model)
#' auprc = calculate_auprc(sens_spec$recall, sens_spec$precision)
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
#' @examples
#' # Assuming `base_importance` is a list of feature importance matrices from base models,
#' # `base_models` is a vector of model names, and `meta_learner` is the trained meta-learner model:
#' var_importance = calculate_feature_importance_stacking(base_importance, base_models, meta_learner)
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
#' @param model A trained machine learning model (e.g., from `caret` or other model training functions).
#' @param test_data A matrix or data frame containing the testing dataset (features only).
#' @param target A character vector of true target values for the test data (the observed labels).
#' @param file.name A character string to specify the filename for saving the confusion matrix plot
#'                  (optional). If `NULL`, the plot is not saved.
#' @param maximize A character string indicating which metric to maximize when selecting the best
#'                 threshold for the confusion matrix. Options include "Accuracy", "Precision",
#'                 "Recall", "Specificity", "Sensitivity", "F1", or "MCC". Default is "Accuracy".
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
#' @examples
#' # Example of usage with a trained model, testing data, and true labels:
#' prediction_results = compute.prediction(model = trained_model,
#'                                         test_data = testing_set,
#'                                         target = true_labels)
#'
#' @seealso \code{\link[caret]{confusionMatrix}}, \code{\link[caret]{varImp}}, \code{\link[ggplot2]{ggplot}}
#'
#' @import caret
#' @import dplyr
#' @import ggplot2
#' @import reshape2
#' @import grDevices
#' @export
compute.prediction = function(model, test_data, target, file.name = NULL, maximize = "Accuracy"){
  # Maximize: parameter for choosing threshold for confusing matrix: maximize sensitivity, specificity, F1, AUROC, AUPRC
  cat("Predicting target variable using provided ML model.................................................\n")

  if(!maximize %in% c("Accuracy", "Precision", "Recall", "Specificity", "Sensitivity", "F1", "MCC")){
    stop("Metric to maximize score to calculate confusion matrix not supported!")
  }

  features <- colnames(test_data)
  are_equal = dplyr::setequal(model[["coefnames"]], features)
  if(are_equal == T){
    #Predict target variable
    predict <- data.frame(stats::predict(model, test_data, type = "prob"))
    #predict$yes = compute.platt.scaling(target, predict$yes)
    #Get metrics
    sens_spec = get_sensitivity_specificity(predict, target, model$method)
    auroc = calculate_auroc(sens_spec$fpr, sens_spec$Sensitivity)
    auprc = calculate_auprc(sens_spec$Recall, sens_spec$Precision)

    ## Calculate confusion matrix based on best threshold
    cat("Choosing the threshold that maximizes", maximize ,"for calculating the confusion matrix...................................................\n")
    max_ind = which.max(sens_spec[,maximize])
    best_threshold <- sens_spec$yes[max_ind] #Find the threshold that maximizes F1 score
    cat("Best threshold: ", best_threshold, "\n")
    cat("Accuracy: ", round(sens_spec$Accuracy[max_ind]*100,3), "\n")
    cat("Sensitivity: ", round(sens_spec$Sensitivity[max_ind]*100,3), "\n")
    cat("Specificity: ", round(sens_spec$Specificity[max_ind]*100,3), "\n")
    cat("F1 score: ", round(sens_spec$F1[max_ind]*100, 3), "\n")
    cat("MCC score: ", round(sens_spec$MCC[max_ind]*100, 3), "\n")
    cat("Recall: ", round(sens_spec$Recall[max_ind]*100,3), "\n")
    cat("Precision: ", round(sens_spec$Precision[max_ind]*100, 3), "\n")
    predicted_classes <- ifelse(predict$yes >= best_threshold, "yes", "no") #Classify predictions based on the best threshold
    conf_matrix <- caret::confusionMatrix(factor(predicted_classes, levels = c("no", "yes")), factor(target, levels = c("no", "yes"))) #Calculate confusion matrix values using the predicted_classes and true_labels
    confusion_matrix <- as.data.frame(as.table(conf_matrix$table))
    colnames(confusion_matrix) <- c("Prediction", "Actual", "Count")

    confusion_matrix_melted <- reshape2::melt(confusion_matrix, id.vars = c("Prediction", "Actual"))

    p = ggplot2::ggplot(confusion_matrix_melted, ggplot2::aes(x = Actual, y = Prediction, fill = value)) +
      ggplot2::geom_tile(color = "black") +  # Add black border around tiles
      ggplot2::geom_text(ggplot2::aes(label = value), color = "black", size = 6) +  # Numbers in black
      ggplot2::scale_fill_gradient(low = "white", high = "red", limits = c(0, max(confusion_matrix_melted$value))) +
      ggplot2::labs(title = "Confusion Matrix", x = "Actual", y = "Prediction") +
      ggplot2::theme_minimal() +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)) +
      ggplot2::theme(axis.text.y = ggplot2::element_text(angle = 45, hjust = 1))

    grDevices::pdf(paste0("Results/Confusion_Matrix_", file.name, ".pdf"))
    print(p)
    grDevices::dev.off()

    return(list(Metrics = sens_spec, AUC = list("AUROC" = auroc, "AUPRC" = auprc), Predictions = predict))
  }else{
    message("Testing set does not count with the same features as model")
  }
}

#' Compute prediction using stacking approach
#'
#' @param super.learner Meta-learner model
#' @param test_data A matrix with the testing dataset
#' @param target Character vector with the true values
#' @param ml.models Machine learning models
#' @param base.models A character vector with the base models for the meta-learner
#'
#' @return A list containing
#'
#' - Prediction metrics
#' - AUROC and AUPRC
#' - Prediction values
#'
#' @export
#'
#' @examples
#'
#' prediction = compute.prediction.stacked(model, testing_set, target, ML_models, base_models)
#'
compute.prediction.stacked = function(super.learner, test_data, target, ml.models, base.models){

  #Learning from simple meta-learner
  base_predictions = list()
  for (i in 1:length(base.models)) {
    base_predictions[[i]] = stats::predict(ml.models[[base.models[i]]], test_data, type = "prob")$yes
    names(base_predictions)[i] = base.models[i]
  }

  base_predictions = do.call(cbind, base_predictions)

  prediction_simple = data.frame(stats::predict(super.learner, base_predictions, type = "prob"))

  #Metrics

  #Meta-learner simple
  sens_spec_simple = get_sensitivity_specificity(prediction_simple, target, "Meta-learner_simple")
  auroc_simple = calculate_auroc(sens_spec_simple$fpr, sens_spec_simple$Sensitivity)
  auprc_simple = calculate_auprc(sens_spec_simple$Recall, sens_spec_simple$Precision)

  return(list(Metrics = sens_spec_simple, AUC = list("AUROC" = auroc_simple, "AUPRC" = auprc_simple), Predictions = prediction_simple))

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
#' @export
#'
#' @examples
#' observed_values = c("yes", "yes", "no", "yes")
#' accuracy = calculate_accuracy(metrics, observed_values)
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
#' @export
#'
#' @examples
#' target = c("yes", "no", "yes", "no", "no")
#' metrics = get_sensitivity_specificity(predict, target, model)
#' confusion_values = calculate_confusion_values(metrics, target)
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
#' @export
#'
#' @examples
#' observed = c("yes", "no", "yes", "no", "no")
#' metrics = get_sensitivity_specificity(predict, target, model)
#' precision = calculate_precision(metrics, observed)
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
#' @export
#'
#' @examples
#' observed = c("yes", "no", "yes", "no", "no")
#' metrics = get_sensitivity_specificity(predict, target, model)
#' recall = calculate_recall(metrics, observed)
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
#' @examples
#' get_curves(metrics, "specificity", "sensitivity", "recall", "precision", "model", auc_roc_score, auc_prc_score, "Test")
#'
get_curves = function(data, spec, sens, reca, prec, color, auc_roc, auc_prc, file.name){

  data = data %>%
    dplyr::mutate(specificity = data[,spec],
                  sensitivity = data[,sens],
                  recall = data[,reca],
                  precision = data[,prec],
                  color = data[,color])

  #Add AUC scores to data frame
  data$color.roc <- paste(data$color, "\n(AUC-ROC =", round(auc_roc, 2), ")\n")

  # Plot the ROC curves
  roc = ggplot2::ggplot(data = data, ggplot2::aes(x = 1- specificity, y = sensitivity, color = color.roc)) +
    ggplot2::geom_line() +
    ggplot2::geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey") +
    ggplot2::labs(title = "ROC Curve", x = "1 - Specificity", y = "Sensitivity") +
    ggplot2::theme_minimal() +
    ggplot2::theme(legend.title = ggplot2::element_blank())

  #Add AUC scores to data frame
  data$color.prc <- paste(data$color, "\n(AUC-PRC =", round(auc_prc, 2), ")\n")

  # Plot recall curves
  recall = ggplot2::ggplot(data = data, aes(x = recall, y = precision, color = color.prc)) +
    ggplot2::geom_line() +
    ggplot2::labs(title = "Precision-Recall Curve", x = "Recall", y = "Precision") +
    ggplot2::ylim(0, 1) +
    ggplot2::theme_minimal() +
    ggplot2::theme(legend.title = ggplot2::element_blank())

  grDevices::pdf(paste0("Results/ROC_curve_", file.name, ".pdf"))
  print(roc)
  grDevices::dev.off()

  grDevices::pdf(paste0("Results/Recall_curve_", file.name, ".pdf"))
  print(recall)
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
#' @examples
#' find.ML.models("Results/ML_models", "ROC", 0.7)
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

  #positive and negative class are defined based on the factor() from trait (this have been defined in compute.ML() function already)
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

compute.platt.scaling = function(obs, yes){
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
  folds_list <- map(1:n_rep, ~{
    # Split within each dataset
    dataset_folds <- train_data %>%
      dplyr::group_by(dataset) %>% #Groups data by dataset
      dplyr::group_split() %>% #Splits into a list where each dataset is a separate dataframe
      purrr::map(~ caret::createFolds(.x$target, k = k_folds, returnTrain = TRUE))  #Inside each dataset, stratified k-fold splits are created using ensuring class balance

    fold_indices <- vector("list", k_folds) #Creates an empty list for k folds

    for (i in seq_len(k_folds)) {
      # Extracts the i-th fold from each dataset
      fold_indices[[i]] <- unlist(purrr::map(dataset_folds, ~ .x[[i]])) #Merges fold indices from all datasets into a single vector
    }

    return(fold_indices)
  })

  multifolds <- unlist(folds_list, recursive = FALSE)

  return(multifolds)
}

calculate_cv_metrics = function(ml_model, metric, hyperparameters = NULL){

  ## List hyperparameters
  if(is.null(hyperparameters)==F){
    group_vars = c("Resample", hyperparameters)
  }else{
    group_vars = "Resample"
  }

  ## Calculate AUCs and integrate into prediction matrix
  ml_model$pred = ml_model$pred %>%
    dplyr::group_by(dplyr::across(dplyr::all_of(group_vars))) %>%
    dplyr::mutate(calibrated_yes = compute.platt.scaling(obs, yes), # Apply Plat scaling to calibrate probabilities
                  AUROC= calculate_auc_roc_resample(obs, yes), # Calculate AUC-ROC if metric is "AUROC"
                  AUPRC = calculate_auc_prc_resample(obs, yes) # Calculate AUC-PRC if metric is "AUPRC"
    ) %>%
    dplyr::ungroup() %>%
    data.frame() %>%
    dplyr::select(-yes) %>%
    dplyr::rename(yes = calibrated_yes)

  ## Calculate prediction metrics (Accuracy, Recall, Precision, F1, MCC)
  metrics <- ml_model$pred %>%
    dplyr::group_by(dplyr::across(dplyr::all_of(group_vars))) %>%
    dplyr::group_split() %>%
    purrr::map(~ get_sensitivity_specificity(.x, .x$obs, "test")) %>%
    dplyr::bind_rows() %>%
    dplyr::select(-model)

  ml_model$pred = ml_model$pred %>%
    dplyr::select(-yes) %>% #remove yes probabilities from resamples (only keep those ordered for calculated the metrics)
    dplyr::bind_cols(metrics) %>%
    dplyr::select(-pred, -obs, -no) %>%
    dplyr::select(Resample, yes, dplyr::everything())


  if(is.null(hyperparameters) == F){
    ## Integrate average CV metrics across repetitions into resample matrix
    df_avg <- ml_model$pred %>%
      dplyr::group_by(dplyr::across(dplyr::all_of(hyperparameters))) %>%
      dplyr::summarise(
        Sensitivity = mean(Sensitivity, na.rm = TRUE),
        Specificity = mean(Specificity, na.rm = TRUE),
        Accuracy = mean(Accuracy, na.rm = TRUE),
        AUROC = mean(AUROC, na.rm = TRUE),
        AUPRC = mean(AUPRC, na.rm = TRUE),
        Precision = mean(Precision, na.rm = TRUE),
        Recall = mean(Recall, na.rm = TRUE),
        F1 = mean(F1, na.rm = TRUE),
        MCC = mean(MCC, na.rm = TRUE),
        Sensitivity_SD = sd(Sensitivity, na.rm = TRUE),
        Specificity_SD = sd(Specificity, na.rm = TRUE),
        Accuracy_SD = sd(Accuracy, na.rm = TRUE),
        AUROC_SD = sd(AUROC, na.rm = TRUE),
        AUPRC_SD = sd(AUPRC, na.rm = TRUE),
        Precision_SD = sd(Precision, na.rm = TRUE),
        Recall_SD = sd(Recall, na.rm = TRUE),
        F1_SD = sd(F1, na.rm = TRUE),
        MCC_SD = sd(MCC, na.rm = TRUE),
        .groups = "keep"  # Ensures all groups are retained
      ) %>%
      dplyr::ungroup()

    ml_model$results <- ml_model$results %>%
      dplyr::select(-dplyr::all_of(hyperparameters), -Accuracy, -Kappa, -AccuracySD, -KappaSD) %>%
      dplyr::bind_cols(df_avg)

    tune = which.max(ml_model$results[,metric])  #Tuning parameter (select combination with top AUROC or AUPRC)

    ml_model$bestTune = ml_model$results %>%
      dplyr::slice(tune) %>%  # Extract only the row with the tuned value
      dplyr::select(dplyr::all_of(hyperparameters))

    filter_conditions <- ml_model$bestTune[1, , drop = FALSE] #Take tuned parameters

    ## Integrate average CV metrics across repetitions only in tuned parameters into resample matrix
    df_filtered <- ml_model$pred

    for (col in names(filter_conditions)) {
      df_filtered <- df_filtered[df_filtered[[col]] == filter_conditions[[col]], ] #Filter by keeping only rows where the column matches the corresponding value in filter_conditions (Continue refining until all conditions are applied)
    }

    df_avg = df_filtered %>%
      dplyr::group_by(Resample) %>%
      dplyr::summarise(
        Sensitivity = mean(Sensitivity, na.rm = TRUE),
        Specificity = mean(Specificity, na.rm = TRUE),
        Accuracy = mean(Accuracy, na.rm = TRUE),
        AUROC = mean(AUROC, na.rm = TRUE),
        AUPRC = mean(AUPRC, na.rm = TRUE),
        Precision = mean(Precision, na.rm = TRUE),
        Recall = mean(Recall, na.rm = TRUE),
        F1 = mean(F1, na.rm = TRUE),
        MCC = mean(MCC, na.rm = TRUE)
      ) %>%
      dplyr::ungroup()



  }else{
    ## Integrate average CV metrics across repetitions into resample matrix
    df_avg <- ml_model$pred %>%
      dplyr::group_by(Resample) %>%
      dplyr::summarise(
        Sensitivity = mean(Sensitivity, na.rm = TRUE),
        Specificity = mean(Specificity, na.rm = TRUE),
        Accuracy = mean(Accuracy, na.rm = TRUE),
        AUROC = mean(AUROC, na.rm = TRUE),
        AUPRC = mean(AUPRC, na.rm = TRUE),
        Precision = mean(Precision, na.rm = TRUE),
        Recall = mean(Recall, na.rm = TRUE),
        F1 = mean(F1, na.rm = TRUE),
        MCC = mean(MCC, na.rm = TRUE)
      ) %>%
      dplyr::ungroup()
  }

  ml_model$resample = ml_model$resample %>%
    dplyr::select(-Accuracy, -Kappa) %>%
    dplyr::arrange(match(Resample, df_avg$Resample)) %>%
    dplyr::select(-Resample) %>%
    dplyr::bind_cols(df_avg)

  if(is.null(hyperparameters) == T){
    ## Integrate average CV metrics into results
    ml_model$results = ml_model$results %>%
      dplyr::select(-Accuracy, -AccuracySD, -Kappa, -KappaSD) %>%
      dplyr::mutate(Sensitivity = mean(ml_model$resample$Sensitivity, na.rm = TRUE),
                    Specificity = mean(ml_model$resample$Specificity, na.rm = TRUE),
                    Accuracy = mean(ml_model$resample$Accuracy, na.rm = TRUE),
                    AUROC = mean(ml_model$resample$AUROC, na.rm = TRUE),
                    AUPRC = mean(ml_model$resample$AUPRC, na.rm = TRUE),
                    Precision = mean(ml_model$resample$Precision, na.rm = TRUE),
                    Recall = mean(ml_model$resample$Recall, na.rm = TRUE),
                    F1 = mean(ml_model$resample$F1, na.rm = TRUE),
                    MCC = mean(ml_model$resample$MCC, na.rm = TRUE),
                    Sensitivity_SD = sd(ml_model$resample$Sensitivity, na.rm = TRUE),
                    Specificity_SD = sd(ml_model$resample$Specificity, na.rm = TRUE),
                    Accuracy_SD = sd(ml_model$resample$Accuracy, na.rm = TRUE),
                    AUC_SD = sd(ml_model$resample$AUC, na.rm = TRUE),
                    Precision_SD = sd(ml_model$resample$Precision, na.rm = TRUE),
                    Recall_SD = sd(ml_model$resample$Recall, na.rm = TRUE),
                    F1_SD = sd(ml_model$resample$F1, na.rm = TRUE),
                    MCC_SD = sd(ml_model$resample$MCC, na.rm = TRUE))
  }


  return(list(Prediction = ml_model$pred, Resamples = ml_model$resample, Results = ml_model$results))


}


compute_shap_values <- function(model_trained, data_train, method, n_cores = 2) {
  gc() #clean memory before start

  cat("Computing SHAP values for feature importance...............................................................\n\n")

  # Get tuned hyperparameters
  filter_conditions <- model_trained$bestTune[1, , drop = FALSE]

  # Filter predictions for best tune
  if (any(filter_conditions != "none")) {
    for (col in names(filter_conditions)) {
      model_trained$pred <- model_trained$pred[model_trained$pred[[col]] == filter_conditions[[col]], ]
    }
  }

  resamples <- unique(model_trained$pred$Resample)

  # Register parallel backend
  cl <- parallel::makeCluster(n_cores)
  doParallel::registerDoParallel(cl)

  # Flag to track trivial predictions
  trivial_predictions_found <- FALSE

  importance_list <- foreach::foreach(resample = resamples, .export = c('trivial_predictions_found')) %dopar% {

    if (trivial_predictions_found) {
      return(NULL)  # Skip further iterations if trivial predictions were found in any previous resample
    }

    tryCatch({
      library(pipeML)
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

      X_test <- test_data_fold[, setdiff(colnames(test_data_fold), "target")]

      # Get predicted probabilities
      pred_probs <- stats::predict(fit, newdata = X_test, type = "prob")

      # Check for trivial predictions (same for all "yes"/"no")
      if (length(unique(round(pred_probs[, 2], 5))) > 1) {
        # If predictions are not trivial, compute SHAP values
        shap_values <- kernelshap::kernelshap(
          fit,
          X = X_test,
          bg_X = train_data_fold[, setdiff(colnames(train_data_fold), "target")],
          exact = F,
          type = "prob",
          parallel = F
        )

        # Return importance dataframe for positive class ("yes")
        shap_values$S$yes %>%
          data.frame() %>%
          dplyr::mutate(Samples = rownames(X_test))

      } else {
        # Set the flag for trivial predictions found
        trivial_predictions_found = TRUE
        stop("Error on SHAP: model predictions are constant!")  # Stop further computation for this worker
      }
    }, error = function(e) {
      # Handle any errors here and print them to help with debugging
      message("Error encountered in resample: ", resample, " - ", e$message)
      return(NULL)  # Return NULL in case of an error, so the loop continues
    })
  }

  # Stop the cluster after parallel execution
  parallel::stopCluster(cl)
  unregister_dopar() # Stop Dopar from running in the background

  gc()

  # Check if trivial predictions were found
  if (any(sapply(importance_list, is.null))) {
    stop("Trivial predictions were found in some resamples. SHAP values cannot be calculated")
  }

  # Combine and summarize importance results
  importance_df <- do.call(rbind, importance_list) %>%
    dplyr::group_by(Samples) %>%
    dplyr::summarise(dplyr::across(dplyr::where(is.numeric), mean), .groups = "drop") %>%
    dplyr::arrange(match(Samples, rownames(data_train))) %>%
    tibble::column_to_rownames("Samples")

  return(importance_df)
}

#' Plot SHAP values
#'
#' Plots the variable importance based on SHAP (SHapley Additive exPlanations) values
#' and saves the plot to the `Results/` directory.
#'
#' @param shap_df A data frame or matrix containing SHAP values where rows represent samples and columns represent features.
#' @param ml_model A character string representing the name of the machine learning model.
#' @param file_name A character string representing the name of the file to save the plot in the `Results/` directory.
#' @param width A numeric value specifying the width of the plot in inches (default is 10).
#' @param height A numeric value specifying the height of the plot in inches (default is 10).
#'
#' @return A plot saved as a PDF file in the `Results/` directory showing the variable importance of the machine learning model.
#'
#' @details This function generates a bar plot of the SHAP values, where the features are sorted by their mean SHAP value. The plot distinguishes
#' between features that increase the predicted outcome (colored in red) and those that decrease the predicted outcome (colored in blue).
#' The plot is saved as a PDF file in the `Results/` directory, with the filename specified by the user.
#'
#' @export
#'
#' @examples
#' shap_df = data.frame(feature1 = rnorm(100), feature2 = rnorm(100))
#' ml_model = "RandomForest"
#' plot_shap_values(shap_df, ml_model, "shap_plot")
#'
plot_shap_values = function(shap_df, ml_model, file_name, width = 10, height = 10){

  #shap_df : shap values (samples as rows and features as columns)

  shap_long <- shap_df %>%
    tidyr::pivot_longer(cols = dplyr::everything(), names_to = "feature", values_to = "shap_value")

  ### Average SHAP values per feature across samples
  shap_summary <- shap_long %>%
    dplyr::group_by(feature) %>%
    dplyr::summarise(mean_shap = mean(shap_value, na.rm = TRUE)) %>%
    dplyr::filter(abs(mean_shap) > summary(abs(mean_shap))[3]) %>%
    dplyr::mutate(direction = ifelse(mean_shap > 0, "Increase", "Decrease"))

  p = ggplot2::ggplot(shap_summary, ggplot2::aes(x = stats::reorder(feature, mean_shap), y = mean_shap, fill = direction)) +
    ggplot2::geom_col() +
    ggplot2::coord_flip() +
    ggplot2::scale_fill_manual(values = c("Increase" = "red", "Decrease" = "blue")) +
    ggplot2::theme_minimal(base_size = 14) +
    ggplot2::labs(
      title = paste0("Variable importance for ML model ", ml_model),
      x = NULL,
      y = NULL,
      fill = NULL
    ) +
    ggplot2::geom_hline(yintercept = 0, color = "black")

  grDevices::pdf(paste0("Results/Feature_importance_", ml_model, "_", file_name, ".pdf"), width = width, height = height)
  print(p)
  grDevices::dev.off()

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
#' @examples
#' importance = compute.variable.importance(ml_model, stacking = FALSE, n_cores = 2)
#'
compute.variable.importance = function(model, stacking = FALSE, n_cores = 2){

  if(stacking == TRUE){
    base_models = model$Model$Base_models
    ml_models = model$Model$ML_models
    train_data = model$Model$ML_models$BAG$trainingData %>% #All training datasets are the same so BAG it's choosing randomly
      dplyr::rename(target = .outcome)

    importance = list() #Save variable importance of each base model
    for (i in 1:length(base_models)) {
      importance[[i]] = compute_shap_values(ml_models[[base_models[i]]], train_data, ml_models[[base_models[i]]]$method, n_cores) ## Compute SHAP values
    }
    importance_df <- Reduce(function(x, y) (x + y) / length(importance), importance) #Take the mean importance
  }else{
    ml_model = model$Model$Model
    train_data = model$Model$Model$trainingData %>%
      dplyr::rename(target = .outcome)
    ml_method = model$Model$Model$method

    importance_df = compute_shap_values(ml_model, train_data, ml_method, n_cores) ## Compute SHAP values for variable importance
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
