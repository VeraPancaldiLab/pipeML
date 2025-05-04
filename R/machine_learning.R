
##Basic function to compute boruta algorithm
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

##Merge results from boruta iterations
merge_boruta_results = function(importance_values, decisions, file_name, iterations, threshold, return = T){

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

##Function for iteratively running boruta (parallelization available)
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
      # arg_list <- replicate(iterations, list(data, sample.int(100000, 1), fix), simplify = FALSE)
      # system.time({
      #   res <- mclapply(arg_list, function(x) {
      #     do.call(compute.boruta, x)
      #   }, mc.cores = num_cores)
      # })

      res <- foreach::foreach(seed = sample.int(100000, iterations)) %dopar% {

        source("src/environment_set.R")

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

# Create folds for each repetition with different seeds (Not used anymore, replace by Multifolds)
create_folds_for_repetitions <- function(data, k_folds, n_rep) {
  all_folds <- list()
  for (rep in 1:n_rep) {
    set.seed(sample.int(100000, 1)) # Change seed for each repetition
    folds <- caret::createFolds(data$target, k = k_folds, returnTrain = TRUE, list = TRUE)
    all_folds[[rep]] <- folds
  }
  return(all_folds)
}

##Compute boruta during each fold (DEPRECATED - might give overfitting)
computed.boruta.kfolds = function(folds, data_model, boruta_iterations, fix_boruta, tentative, threshold, file_name){

  folds_threshold = 0.8*length(folds)
  features_folds = list()

  for (i in 1:length(folds)) {
    message("Feature selection using Boruta...............................................................\n\n")
    training_set = data_model[folds[[i]],]
    res_boruta = feature.selection.boruta(training_set, iterations = boruta_iterations, fix = fix_boruta, thres = threshold, file_name = file_name, return = T)

    if(tentative == F){
      if(length(res_boruta$Confirmed) <= 1){
        features_folds[[i]] = list()
        message("\nNo features were confirmed in more than ", round(threshold*100) ,"% of the times for training in this specific fold.......................\n\n")
      }
      message("\nKeeping only features confirmed in more than", round(threshold*100) ,"% of the times for training in this specific fold......................\n\n")
      message("If you want to consider also tentative features, please specify tentative = T in the parameters.\n\n")
      features_folds[[i]] = res_boruta$Confirmed
    }else{
      sum_features = length(res_boruta$Confirmed) + length(res_boruta$Tentative)
      if(sum_features <= 1){
        features_folds[[i]] = list()
        message("\nNo features were confirmed in more than ", round(thresh*100) ,"% of the times for training in this specific fold.......................\n\n")
      }
      message("\nKeeping only features confirmed and tentative in more than", round(thresh*100) ,"% of the times for training in this specific fold............................\n\n")
      features_folds[[i]] = c(res_boruta$Confirmed, res_boruta$Tentative)
    }
  }

  all_features <- unlist(features_folds)
  feature_freq <- table(all_features)
  selected_features <- names(feature_freq[feature_freq >= folds_threshold])

  if(length(selected_features)<=1){
    message("No features selected meet the requirements. Try with different parameter values.")
  }else{
    return(selected_features)
  }

}

#Main function for CV training using 13 ML models
compute.k_fold_CV = function(model, k_folds, n_rep, stacking = F, metric = "Accuracy", boruta, boruta_iterations = NULL, fix_boruta = NULL, tentative = F, boruta_threshold = NULL, file_name = NULL, LODO = F, return){

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
        cat("If you want to consider also tentative features, please specify tentative = T in the parameters.\n\n")
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
  #folds <- createFolds(train_data[,'target'], k = k_folds, returnTrain = T, list = T) #this for single folds
  if(LODO == T){
    multifolds = construct_stratified_cohort_folds(train_data, 'dataset', 'target', k_folds = k_folds, n_rep = n_rep)
    train_data = train_data %>% dplyr::select(-dataset) #After creating multifolds we remove this variable to be able to train
  }else{
    multifolds = caret::createMultiFolds(train_data[,'target'], k = k_folds, times = n_rep) #repeated folds
  }

  trainControl <- caret::trainControl(index = multifolds, method="repeatedcv", number=k_folds, repeats=n_rep, verboseIter = F, allowParallel = T, classProbs = TRUE, savePredictions=T)

  ######### Feature selection across folds (DEPRECATED)
  # if(boruta == T){
  #   features = computed.boruta.kfolds(multifolds, model, boruta_iterations = boruta_iterations, fix_boruta = fix_boruta, tentative = tentative, threshold = boruta_threshold, file_name = file_name)
  #   if(is.null(features)==T){
  #     cat("No features selected across folds after Boruta. Try different parameter values")
  #     return(NULL)
  #   }else{
  #     train_data = model[,colnames(model)%in%features, drop = F] %>%
  #       mutate(target = model$target)
  #   }
  # }else{
  #   train_data = model
  # }


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
      # if(return == T){
      #   plot_cv_metrics(ensembleResults, file_name = file_name)
      # }
    }

    top_model = metrics[["Top_model"]]

    model = ensembleResults[[top_model]]

    cat("Best ML model found: ", top_model, "\n")

    cat("Returning model trained\n")

    output = list("Features" = features, "Model" = model, "ML_Models" = ensembleResults)
  }


  return(output)

}


#Main ML Pipeline for one partition
compute.ML = function(raw.counts, normalized = F, deconv, tf.universe, clinical, trait, trait.positive, partition, metric = "Accuracy", stack, feature.selection = F, seed, file_name = NULL, return = F){

  set.seed(seed)

  # Do stratified partition
  index = caret::createDataPartition(clinical[,trait], times = 1, p = partition, list = FALSE)

  # Normalize counts
  if(normalized == T){
    norm.counts = data.frame(ADImpute::NormalizeTPM(raw.counts, log = T))
  }else{
    norm.counts = raw.counts
  }

  # Train cohort
  traitData_train = clinical[index, ]
  raw.counts_train = raw.counts[,index]
  counts.normalized_train = norm.counts[,index]
  deconv_train = deconv[index,]
  tfs_train = compute.TFs.activity(counts.normalized_train, universe = tf.universe)

  # Test cohort
  traitData_test = clinical[-index, ]
  raw.counts_test = raw.counts[,-index]
  counts.normalized_test = norm.counts[,-index]
  deconv_test = deconv[-index,]
  tfs_test = compute.TFs.activity(counts.normalized_test, universe = tf.universe)

  ######################################################################################################################## Training
  network = compute.WTCNA(tfs_train, corr_mod = 0.7, clustering.method = "ward.D2", return = return)
  dt = compute.deconvolution.analysis(deconv_train, corr = 0.7, seed = 123, return = return, file_name = file_name)

  set.seed(seed)

  ##### Construct cell groups
  cell_groups = construct_cell_groups(counts.normalized_train, tfs_train, deconv_train, network, dt, traitData_train, pval = 0.05, high_corr_groups = 0.9, clustering.method = "ward.D2", trait = trait, positive = trait.positive)

  #Set training set
  train_data = cell_groups[[1]] %>%
    data.frame() %>%
    dplyr::mutate(Trait = traitData_train[,trait],
                  target = as.factor(ifelse(Trait == trait.positive, 'yes', 'no'))) %>%
    dplyr::select(-Trait)

  train_data$target <- factor(train_data$target, levels = c("no", "yes"))  # Order, just in case to ensure positive class is not well defined

  #Cross-validation training (5 k-folds and 100 repetitions)
  set.seed(seed)
  training = compute.k_fold_CV(train_data, k_folds = 5, n_rep = 5, metric = metric, stacking = stack, boruta = feature.selection, boruta_iterations = 100, fix_boruta = F, boruta_threshold = 0.8, file_name = file_name, LODO = LODO, batch_id = batch_id, return= return)

  ####################################################Predicting
  if(length(training)!=0){
    network = network #Save TF network per partition
    features = training[["Features"]] #Save selected features per partition (only useful if Boruta = T - needs to be improve it)
    deconvolution_subgroups = dt[["Deconvolution groups - Linear-based correlation"]] #Save cell subgroups
    ####################### Testing set
    testing_set = compute.test.set(dt, cell_groups, names(cell_groups[[2]]), deconv_test) #Compute features in testing set
    #testing_set = replicate_deconvolution_subgroups(dt, deconv_test)

    #Extract target variable
    target = traitData_test %>%
      dplyr::mutate(target = ifelse(traitData_test[,trait] == trait.positive, "yes", "no")) %>%
      dplyr::pull(target)

    target = factor(target, levels = c("no", "yes"))

    if(stack){
      model = training[["Meta_learner"]]
      prediction = compute.prediction.stacked(model, testing_set, target, training[["ML_models"]], training[["Base_models"]])
    }else{
      model = training[["Model"]] #Save best ML model based on the Accuracy/AUC from CV per partition
      prediction = compute.prediction(model, testing_set, target)
    }

    auc_roc_score = prediction[["AUC"]][["AUROC"]]
    auc_prc_score = prediction[["AUC"]][["AUPRC"]]

    metrics = prediction[["Metrics"]]
    predictions = prediction[["Predictions"]]

    if(return == T){
      get_curves(metrics, "Specificity", "Sensitivity", "Recall", "Precision", "model", auc_roc_score, auc_prc_score, file_name)
    }

    gc() #Clean garbage

    return(list(Model = model, Features = features, Cell_groups = cell_groups, Deconvolution_subgroups = deconvolution_subgroups, TF_network = network, Prediction_metrics = metrics, AUC = list(AUROC = auc_roc_score, AUPRC = auc_prc_score), Prediction = predictions))
  }else{  #No features are selected as predictive

    gc() #Clean garbage

    message("No features selected as predictive after Boruta runs. No model returned.")

    return(NULL)
  }

}

compute.TME.features.ML = function(raw.counts, normalized = F, TME, clinical, trait, trait.positive, partition, metric = "Accuracy", stack, feature.selection = F, seed, file_name = NULL, return = F){

  set.seed(seed)

  # Do stratified partition
  index = caret::createDataPartition(clinical[,trait], times = 1, p = partition, list = FALSE)

  # Normalize counts
  if(normalized == T){
    norm.counts = data.frame(ADImpute::NormalizeTPM(raw.counts, log = T))
  }else{
    norm.counts = raw.counts
  }

  # Train cohort
  traitData_train = clinical[index, ]
  raw.counts_train = raw.counts[,index]
  counts.normalized_train = norm.counts[,index]

  # Test cohort
  traitData_test = clinical[-index, ]
  raw.counts_test = raw.counts[,-index]
  counts.normalized_test = norm.counts[,-index]

  ### TME features
  hallmarks_of_immune_response <- c("CYT", "Roh_IS", "chemokines", "Davoli_IS", "IFNy", "Ayers_expIS", "RIR", "TLS")

  # Training
  immune_response_scores.training <- compute_scores_immune_response(RNA_tpm = counts.normalized_train, selected_scores = hallmarks_of_immune_response)
  cell_fractions.training <- compute_cell_fractions(RNA_tpm = counts.normalized_train)
  pathway_activities.training <- compute_pathway_activity(RNA_counts = counts.normalized_train, remove_sig_genes_immune_response = TRUE)
  tf_activities.training <- compute.TFs.activity(counts.normalized_train)

  colnames(cell_fractions.training) <- gsub(" |\\+", "_", colnames(cell_fractions.training))
  colnames(pathway_activities.training) = gsub(" |\\+|-", "_", colnames(pathway_activities.training))

  # Testing
  immune_response_scores.testing <- compute_scores_immune_response(RNA_tpm = counts.normalized_test, selected_scores = hallmarks_of_immune_response)
  cell_fractions.testing <- compute_cell_fractions(RNA_tpm = counts.normalized_test)
  pathway_activities.testing <- compute_pathway_activity(RNA_counts = counts.normalized_test, remove_sig_genes_immune_response = TRUE)
  tf_activities.testing <- compute.TFs.activity(counts.normalized_test)

  colnames(cell_fractions.testing) <- gsub(" |\\+", "_", colnames(cell_fractions.testing))
  colnames(pathway_activities.testing) = gsub(" |\\+|-", "_", colnames(pathway_activities.testing))

  #Set training and testing set
  if(TME == "Immunescores"){
    train_data = immune_response_scores.training %>%
      data.frame() %>%
      dplyr::mutate(Trait = traitData_train[,trait],
                    target = as.factor(ifelse(Trait == trait.positive, 'yes', 'no'))) %>%
      dplyr::select(-Trait)

    train_data$target <- factor(train_data$target, levels = c("no", "yes"))  # Order, just in case to ensure positive class is not well defined

    testing_set = immune_response_scores.testing
  }else if(TME == "TF"){
    train_data = tf_activities.training %>%
      data.frame() %>%
      dplyr::mutate(Trait = traitData_train[,trait],
                    target = as.factor(ifelse(Trait == trait.positive, 'yes', 'no'))) %>%
      dplyr::select(-Trait)

    train_data$target <- factor(train_data$target, levels = c("no", "yes"))  # Order, just in case to ensure positive class is not well defined

    testing_set = tf_activities.testing
  }else if(TME == "Pathways"){
    train_data = pathway_activities.training %>%
      data.frame() %>%
      dplyr::mutate(Trait = traitData_train[,trait],
                    target = as.factor(ifelse(Trait == trait.positive, 'yes', 'no'))) %>%
      dplyr::select(-Trait)

    train_data$target <- factor(train_data$target, levels = c("no", "yes"))  # Order, just in case to ensure positive class is not well defined

    testing_set = pathway_activities.testing
  }else if(TME == "Deconvolution"){
    train_data = cell_fractions.training %>%
      data.frame() %>%
      dplyr::mutate(Trait = traitData_train[,trait],
                    target = as.factor(ifelse(Trait == trait.positive, 'yes', 'no'))) %>%
      dplyr::select(-Trait)

    train_data$target <- factor(train_data$target, levels = c("no", "yes"))  # Order, just in case to ensure positive class is not well defined

    testing_set = cell_fractions.testing
  }

  #Cross-validation training (5 k-folds and 100 repetitions)
  set.seed(seed)
  training = compute.k_fold_CV(train_data, k_folds = 5, n_rep = 5, metric = metric, stacking = stack, boruta = feature.selection, boruta_iterations = 100, fix_boruta = F, boruta_threshold = 0.8, file_name = file_name, LODO = LODO, batch_id = batch_id, return= return)

  ####################################################Predicting
  if(length(training)!=0){
    features = training[["Features"]] #Save selected features per partition (only useful if Boruta = T - needs to be improve it)
    #Extract target variable
    target = traitData_test %>%
      dplyr::mutate(target = ifelse(traitData_test[,trait] == trait.positive, "yes", "no")) %>%
      dplyr::pull(target)

    target = factor(target, levels = c("no", "yes"))

    if(stack){
      model = training[["Meta_learner"]]
      prediction = compute.prediction.stacked(model, testing_set, target, training[["ML_models"]], training[["Base_models"]])
    }else{
      model = training[["Model"]] #Save best ML model based on the Accuracy/AUC from CV per partition
      prediction = compute.prediction(model, testing_set, target)
    }

    auc_roc_score = prediction[["AUC"]][["AUROC"]]
    auc_prc_score = prediction[["AUC"]][["AUPRC"]]

    metrics = prediction[["Metrics"]]
    predictions = prediction[["Predictions"]]

    if(return == T){
      get_curves(metrics, "Specificity", "Sensitivity", "Recall", "Precision", "model", auc_roc_score, auc_prc_score, file_name)
    }

    gc() #Clean garbage

    return(list(Model = model, Features = features, Prediction_metrics = metrics, AUC = list(AUROC = auc_roc_score, AUPRC = auc_prc_score), Prediction = predictions))
  }else{  #No features are selected as predictive

    gc() #Clean garbage

    message("No features selected as predictive after Boruta runs. No model returned.")

    return(NULL)
  }

}

compute.training.ML = function(raw.counts, normalized = F, deconv_train, traitData_train, trait, trait.positive, metric = "Accuracy", stack, k_folds = 10, n_rep = 5, feature.selection = F, seed, LODO = F, batch_id = NULL, file_name = NULL, return = F){

  # Normalize counts
  if(normalized == T){
    norm.counts = data.frame(ADImpute::NormalizeTPM(raw.counts, log = T))
  }else{ #If they are already normalized
    norm.counts = raw.counts
  }

  # Train cohort
  tfs_train = compute.TFs.activity(norm.counts)

  ######################################################################################################################## Training
  network = compute.WTCNA(tfs_train, corr_mod = 0.7, clustering.method = "ward.D2", return = return)
  dt = compute.deconvolution.analysis(deconv_train, corr = 0.7, seed = 123, return = return, file_name = file_name)

  set.seed(seed)

  ##### Construct cell groups
  cell_groups = construct_cell_groups(norm.counts, tfs_train, deconv_train, network, dt, traitData_train, pval = 0.05, high_corr_groups = 0.9, clustering.method = "ward.D2", trait = trait, positive = trait.positive)

  #Set training set
  train_data = cell_groups[[1]] %>%
    data.frame() %>%
    dplyr::mutate(Trait = traitData_train[,trait],
                  target = as.factor(ifelse(Trait == trait.positive, 'yes', 'no'))) %>%
    dplyr::select(-Trait)

  train_data$target <- factor(train_data$target, levels = c("no", "yes"))  # Order, just in case to ensure positive class is not well defined

  if(LODO == T){
    train_data = train_data %>%
      dplyr::mutate(dataset = traitData_train[,batch_id])
  }

  #Cross-validation training (5 k-folds and 100 repetitions)
  training = compute.k_fold_CV(train_data, k_folds = k_folds, n_rep = n_rep, metric = metric, stacking = stack, boruta = feature.selection, boruta_iterations = 100, fix_boruta = F, boruta_threshold = 0.8, file_name = file_name, LODO = LODO, return= return)

  ####################################################Predicting
  if(length(training)!=0){
    network = network #Save TF network per partition
    features = training[["Features"]] #Save selected features per partition (only useful if Boruta = T - needs to be improve it)
    deconvolution_subgroups = dt #Save cell subgroups
    if(stack){
      model = training[["Meta_learner"]]
    }else{
      model = training[["Model"]] #Save best ML model based on the Accuracy/AUC from CV per partition
    }

    return(list(Model = training, Features = features, Cell_groups = cell_groups, Deconvolution_subgroups = deconvolution_subgroups, TF_network = network))
  }else{  #No features are selected as predictive

    message("No features selected as predictive after Boruta runs. No model returned.")

    return(NULL)
  }

}

compute.features.training.ML = function(features_train, traitData_train, trait, trait.positive, metric = "Accuracy", stack, k_folds = 10, n_rep = 5, feature.selection = F, seed, LODO = F, batch_id = NULL, file_name = NULL, return = F){

  set.seed(seed)

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

  #Cross-validation training (5 k-folds and 100 repetitions)
  training = compute.k_fold_CV(train_data, k_folds = k_folds, n_rep = n_rep, metric = metric, stacking = stack, boruta = feature.selection, boruta_iterations = 100, fix_boruta = F, boruta_threshold = 0.8, file_name = file_name, LODO = LODO, return= return)

  ####################################################Predicting
  if(length(training)!=0){
    features = training[["Features"]] #Save selected features per partition (only useful if Boruta = T - needs to be improve it)
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

compute.features.ML = function(features_train, features_test, clinical, trait, trait.positive, metric = "Accuracy", stack, k_folds = 10, n_rep = 5, feature.selection = F, seed, LODO = F, batch_id = NULL, file_name = NULL, return = F){

  # Train cohort
  traitData_train = clinical[rownames(clinical)%in%rownames(features_train), ]

  # Test cohort
  traitData_test = clinical[rownames(clinical)%in%rownames(features_test), ]

  ###############################################################################################################################################################################

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
    features = training[["Features"]] #Save selected features per partition (only useful if Boruta = T - needs to be improve it)
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


compute.deconvolution.ML = function(deconv, clinical, trait, trait.positive, partition, metric = "Accuracy", stack, feature.selection = F, seed, file_name = NULL, return = F){

  set.seed(seed)

  # Do stratified partition
  index = caret::createDataPartition(clinical[,trait], times = 1, p = partition, list = FALSE)

  # Train cohort
  traitData_train = clinical[index, ]
  deconv_train = deconv[index,]
  dt = compute.deconvolution.analysis(deconv_train, corr = 0.7, seed = 123, return = return, file_name = file_name)

  set.seed(seed)

  # Test cohort
  traitData_test = clinical[-index, ]
  deconv_test = deconv[-index,]


  ###############################################################################################################################################################################

  ####################################################Training

  #Set training set
  train_data = dt[[1]] %>%
    data.frame() %>%
    dplyr::mutate(Trait = traitData_train[,trait],
                  target = as.factor(ifelse(Trait == trait.positive, 'yes', 'no'))) %>%
    dplyr::select(-Trait)

  train_data$target <- factor(train_data$target, levels = c("no", "yes"))  # Order, just in case to ensure positive class is not well defined

  #Cross-validation training (5 k-folds and 100 repetitions)
  training = compute.k_fold_CV(train_data, k_folds = 5, n_rep = 5, metric = metric, stacking = stack, boruta = feature.selection, boruta_iterations = 100, fix_boruta = F, boruta_threshold = 0.8, file_name = file_name, LODO = LODO, batch_id = batch_id, return= return)

  ####################################################Predicting
  if(length(training)!=0){
    features = training[["Features"]] #Save selected features per partition (only useful if Boruta = T - needs to be improve it)
    deconvolution_subgroups = dt[["Deconvolution groups - Linear-based correlation"]] #Save cell subgroups
    ####################### Testing set
    testing_set = replicate_deconvolution_subgroups(dt, deconv_test) #Replicate subgroups

    #Extract target variable
    target = traitData_test %>%
      dplyr::mutate(target = ifelse(traitData_test[,trait] == trait.positive, "yes", "no")) %>%
      dplyr::pull(target)

    target = factor(target, levels = c("no", "yes"))

    if(stack){
      model = training[["Meta_learner"]]
      prediction = compute.prediction.stacked(model, testing_set, target, training[["ML_models"]], training[["Base_models"]])

    }else{
      model = training[["Model"]] #Save best ML model based on the Accuracy/AUC from CV per partition
      prediction = compute.prediction(model, testing_set, target)
    }

    auc_roc_score = prediction[["AUC"]][["AUROC"]]
    auc_prc_score = prediction[["AUC"]][["AUPRC"]]

    metrics = prediction[["Metrics"]]
    predictions = prediction[["Predictions"]]

    if(return == T){
      get_curves(metrics, "Specificity", "Sensitivity", "Recall", "Precision", "model", auc_roc_score, auc_prc_score, file_name)
    }


    return(list(Model = model, Features = features, Deconvolution_subgroups = deconvolution_subgroups, Prediction_metrics = metrics, AUC = list(AUROC = auc_roc_score, AUPRC = auc_prc_score), Prediction = predictions))
  }else{  #No features are selected as predictive

    message("No features selected as predictive after Boruta runs. No model returned.")

    return(NULL)
  }

}


compute.raw.deconvolution.ML = function(deconv, clinical, trait, trait.positive, partition, metric = "Accuracy", stack, feature.selection = F, seed, file_name = NULL, return = F){

  set.seed(seed)

  # Do stratified partition
  index = caret::createDataPartition(clinical[,trait], times = 1, p = partition, list = FALSE)

  # Train cohort
  traitData_train = clinical[index, ]
  deconv_train = deconv[index,]

  # Test cohort
  traitData_test = clinical[-index, ]
  deconv_test = deconv[-index,]


  ###############################################################################################################################################################################

  ####################################################Training

  #Set training set
  train_data = deconv_train %>%
    data.frame() %>%
    dplyr::mutate(Trait = traitData_train[,trait],
                  target = as.factor(ifelse(Trait == trait.positive, 'yes', 'no'))) %>%
    dplyr::select(-Trait)

  train_data$target <- factor(train_data$target, levels = c("no", "yes"))  # Order, just in case to ensure positive class is not well defined

  #Cross-validation training (5 k-folds and 100 repetitions)
  training = compute.k_fold_CV(train_data, k_folds = 5, n_rep = 5, metric = metric, stacking = stack, boruta = feature.selection, boruta_iterations = 100, fix_boruta = F, boruta_threshold = 0.8, file_name = file_name, LODO = LODO, batch_id = batch_id, return= return)

  ####################################################Predicting
  if(length(training)!=0){
    features = training[["Features"]] #Save selected features per partition (only useful if Boruta = T - needs to be improve it)
    ####################### Testing set
    testing_set = deconv_test[,features] #Replicate subgroups

    #Extract target variable
    target = traitData_test %>%
      dplyr::mutate(target = ifelse(traitData_test[,trait] == trait.positive, "yes", "no")) %>%
      dplyr::pull(target)

    target = factor(target, levels = c("no", "yes"))

    if(stack){
      model = training[["Meta_learner"]]
      prediction = compute.prediction.stacked(model, testing_set, target, training[["ML_models"]], training[["Base_models"]])

    }else{
      model = training[["Model"]] #Save best ML model based on the Accuracy/AUC from CV per partition
      prediction = compute.prediction(model, testing_set, target)
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

compute.simple.ML = function(train, test, target, metric = "Accuracy", stack, feature.selection = F, file_name = NULL, return = F){

  #Cross-validation training (5 k-folds and 100 repetitions)
  training = compute.k_fold_CV(train, k_folds = 5, n_rep = 5, metric = metric, stacking = stack, boruta = feature.selection, boruta_iterations = 100, fix_boruta = F, boruta_threshold = 0.8, file_name = file_name, LODO = LODO, batch_id = batch_id, return= return)

  if(length(training)!=0){
    ####################### Testing
    if(stack){
      model = training[["Meta_learner"]]
      prediction = compute.prediction.stacked(model, test, target, training[["ML_models"]], training[["Base_models"]])

    }else{
      model = training[["Model"]] #Save best ML model based on the Accuracy/AUC from CV per partition
      prediction = compute.prediction(model, test, target)
    }

    auc_roc_score = prediction[["AUC"]][["AUROC"]]
    auc_prc_score = prediction[["AUC"]][["AUPRC"]]

    metrics = prediction[["Metrics"]]
    predictions = prediction[["Predictions"]]

    if(return == T){
      get_curves(metrics, "Specificity", "Sensitivity", "Recall", "Precision", "model", auc_roc_score, auc_prc_score, file_name)
    }

    return(list(Model = model, Prediction_metrics = metrics, AUC = list(AUROC = auc_roc_score, AUPRC = auc_prc_score), Prediction = predictions))
  }else{  #No features are selected as predictive

    message("No features selected as predictive after Boruta runs. No model returned.")

    return(NULL)
  }

}

compute.LODO.ML = function(raw.counts, normalized = F, deconv, clinical, trait, trait.positive, trait.out, out = NULL, metric = "Accuracy", maximize_threshold = "Accuracy", stack, feature.selection = T, file_name = NULL, seed, k_folds, n_rep, return = T){

  clinical = clinical %>%
    dplyr::mutate(trait.out = clinical[,trait.out]) ##Just a way to make work "filter" cause it does not allow "" variables (might change after)

  # Normalize counts
  if(normalized == T){
    norm.counts = data.frame(ADImpute::NormalizeTPM(raw.counts))
  }else{ #If they are already normalized
    norm.counts = raw.counts
  }

  # Test cohort
  traitData_test = clinical %>%
    dplyr::filter(trait.out == out)
  raw.counts_test = raw.counts[,colnames(raw.counts)%in%rownames(traitData_test)]
  counts.normalized_test = norm.counts[,colnames(norm.counts)%in%rownames(traitData_test)]
  deconv_test = deconv[rownames(deconv)%in%rownames(traitData_test),]
  tfs_test = compute.TFs.activity(counts.normalized_test)

  # Train cohort
  traitData_train = clinical %>%
    dplyr::filter(trait.out != out)
  raw.counts_train = raw.counts[,colnames(raw.counts)%in%rownames(traitData_train)]
  counts.normalized_train = norm.counts[,colnames(norm.counts)%in%rownames(traitData_train)]
  deconv_train = deconv[rownames(deconv)%in%rownames(traitData_train),]
  tfs_train = compute.TFs.activity(counts.normalized_train)

  ######################################################################################################################## Training
  network = compute.WTCNA(tfs_train, corr_mod = 0.8, clustering.method = "ward.D2", return = return, minMod = 3)
  dt = compute.deconvolution.analysis(deconv_train, corr = 0.7, seed = 123, return = return, file_name = file_name)

  set.seed(seed)

  ##### Construct cell groups
  cell_groups = construct_cell_groups(counts.normalized_train, tfs_train, deconv_train, network, dt, traitData_train, pval = 0.05, high_corr_groups = 0.9, clustering.method = "ward.D2", trait = trait, positive = trait.positive)

  #Set training set
  train_data = cell_groups[[1]] %>%
    data.frame() %>%
    dplyr::mutate(Trait = traitData_train[,trait],
                  target = as.factor(ifelse(Trait == trait.positive, 'yes', 'no'))) %>%
    dplyr::select(-Trait)

  train_data$target <- factor(train_data$target, levels = c("no", "yes"))  # Order, just in case to ensure positive class is not well defined

  #Cross-validation training (5 k-folds and 100 repetitions)
  set.seed(seed)
  training = compute.k_fold_CV(train_data, k_folds = k_folds, n_rep = n_rep, metric = metric, stacking = stack, boruta = feature.selection, boruta_iterations = 100, fix_boruta = F, boruta_threshold = 0.8, file_name = file_name, return= return)

  ####################################################Predicting
  if(length(training)!=0){
    network = network #Save TF network per partition
    features = training[["Features"]] #Save selected features per partition (only useful if Boruta = T - needs to be improve it)
    deconvolution_subgroups = dt[["Deconvolution groups - Linear-based correlation"]] #Save cell subgroups
    ####################### Testing set
    testing_set = compute.test.set(dt, cell_groups, names(cell_groups[[2]]), deconv_test) #Compute features in testing set
    #testing_set = replicate_deconvolution_subgroups(dt, deconv_test)

    #Extract target variable
    target = traitData_test %>%
      dplyr::mutate(target = ifelse(traitData_test[,trait] == trait.positive, "yes", "no")) %>%
      dplyr::pull(target)

    target = factor(target, levels = c("no", "yes"))

    if(stack){
      model = training[["Meta_learner"]]
      prediction = compute.prediction.stacked(model, testing_set, target, training[["ML_models"]], training[["Base_models"]])

    }else{
      model = training[["Model"]] #Save best ML model based on the Accuracy/AUC from CV per partition
      prediction = compute.prediction(model, testing_set, target, file_name, maximize = maximize_threshold)
    }

    auc_roc_score = prediction[["AUC"]][["AUROC"]]
    auc_prc_score = prediction[["AUC"]][["AUPRC"]]

    metrics = prediction[["Metrics"]]
    predictions = prediction[["Predictions"]]

    if(return == T){
      get_curves(metrics, "Specificity", "Sensitivity", "Recall", "Precision", "model", auc_roc_score, auc_prc_score, file_name)
    }

    gc() #Clean garbage

    return(list(Model = model, Features = features, Cell_groups = cell_groups, Deconvolution_subgroups = deconvolution_subgroups, TF_network = network, Prediction_metrics = metrics, AUC = list(AUROC = auc_roc_score, AUPRC = auc_prc_score), Prediction = predictions))
  }else{  #No features are selected as predictive

    gc() #Clean garbage

    message("No features selected as predictive after Boruta runs. No model returned.")

    return(NULL)
  }


}

unregister_dopar <- function() {
  env <- foreach:::.foreachGlobals
  rm(list=ls(name=env), pos=env)
  gc()
}

#Main ML pipeline for several partitions
compute.bootstrap.ML = function(raw.counts, normalized = F, deconv, clinical, trait, trait.positive, partition = 0.8, metric = "Accuracy", iterations, feature.selection = F, stack, workers = NULL, file.name = NULL, return = F){

  folder = paste0("Results/", file.name)
  dir.create(file.path(getwd(), folder))

  if(is.null(iterations) == T){
    stop("No iterations specified, please set a number")
  }else{
    if(is.null(workers)==T){
      num_cores <- parallel::detectCores() - 1
    }else{
      num_cores <- workers
    }

    cl = parallel::makeCluster(num_cores) #Forking just copy the R session in its current state. - makeCluster() all must be exported (copied) to the clusters, which can add some overhead
    doParallel::registerDoParallel(cl)

    message("\nRunning ", iterations, " splits for training and test using ", num_cores, " cores")

    tfs.universe = decoupleR::get_collectri(organism = 'human', split_complexes = F)

    # Run foreach loop using each random seed directly
    foreach::foreach(iteration = seq_len(iterations), random.seed = sample.int(100000, iterations)) %dopar% {

      # Use absolute path for the source file to avoid path issues
      source("src/environment_set.R")

      # Run foreach loop using each random seed directly
      tryCatch({
        # Compute the result with the current random seed
        result <- compute.ML(
          raw.counts, normalized, deconv, tfs.universe, clinical, trait, trait.positive, partition,
          metric, stack, feature.selection, seed = random.seed,
          file_name = file.name, return = return
        )

        # Save result as RDS file with unique identifier based on iteration (random seed)
        saveRDS(list(result = result, seed = random.seed),
                file = file.path(paste0(folder, "/ML_result_", random.seed, ".rds")))

        gc()

      }, error = function(e) {

        # Save error information as RDS file with random seed identifier
        saveRDS(list(result = NULL, error = e$message, seed = random.seed),
                file = file.path(paste0(folder, "/ML_result_", random.seed, ".rds")))

        gc()
      })

    }

    #Stop cluster after all runs
    parallel::stopCluster(cl)
    unregister_dopar() #Stop Dopar from running in the background

  }

  message("Analysis is done!")

  message("ML models are saved in Results/ folder")

}

compute.bootstrap.TME.features.ML = function(raw.counts, normalized = F, TME, clinical, trait, trait.positive, partition = 0.8, metric = "Accuracy", iterations, feature.selection = F, stack, workers = NULL, file.name = NULL, return = F){

  folder = paste0("Results/", file.name)
  dir.create(file.path(getwd(), folder))

  if(is.null(iterations) == T){
    stop("No iterations specified, please set a number")
  }else{
    if(is.null(workers)==T){
      num_cores <- parallel::detectCores() - 1
    }else{
      num_cores <- workers
    }

    cl = parallel::makeCluster(num_cores) #Forking just copy the R session in its current state. - makeCluster() all must be exported (copied) to the clusters, which can add some overhead
    doParallel::registerDoParallel(cl)

    message("\nRunning ", iterations, " splits for training and test using ", num_cores, " cores")

    # Run foreach loop using each random seed directly
    foreach::foreach(iteration = seq_len(iterations), random.seed = sample.int(100000, iterations)) %dopar% {

      # Use absolute path for the source file to avoid path issues
      source("src/environment_set.R")

      # Run foreach loop using each random seed directly
      tryCatch({
        # Compute the result with the current random seed
        result <- compute.TME.features.ML(
          raw.counts, normalized, TME, clinical, trait, trait.positive, partition,
          metric, stack, feature.selection, seed = random.seed,
          file_name = file.name, return = return
        )

        # Save result as RDS file with unique identifier based on iteration (random seed)
        saveRDS(list(result = result, seed = random.seed),
                file = file.path(paste0(folder, "/ML_result_", random.seed, ".rds")))

        gc()

      }, error = function(e) {

        # Save error information as RDS file with random seed identifier
        saveRDS(list(result = NULL, error = e$message, seed = random.seed),
                file = file.path(paste0(folder, "/ML_result_", random.seed, ".rds")))

        gc()
      })

    }

    #Stop cluster after all runs
    parallel::stopCluster(cl)
    unregister_dopar() #Stop Dopar from running in the background

  }

  message("Analysis is done!")

  message("ML models are saved in Results/ folder")

}

compute.bootstrap.LODO.ML = function(raw.counts, normalized = F, deconv, clinical, trait, trait.positive, trait.out, out, metric = "Accuracy", iterations, feature.selection = F, stack, workers = NULL, file.name = NULL, return = F){

  folder = paste0("Results/", file.name)
  dir.create(file.path(getwd(), folder))

  if(is.null(iterations) == T){
    stop("No iterations specified, please set a number")
  }else{
    if(is.null(workers)==T){
      num_cores <- parallel::detectCores() - 1
    }else{
      num_cores <- workers
    }

    cl = parallel::makeCluster(num_cores) #Forking just copy the R session in its current state. - makeCluster() all must be exported (copied) to the clusters, which can add some overhead
    doParallel::registerDoParallel(cl)

    message("\nRunning ", iterations, " splits for training and test using ", num_cores, " cores")

    tfs.universe = decoupleR::get_collectri(organism = 'human', split_complexes = F)

    # Run foreach loop using each random seed directly
    foreach::foreach(iteration = seq_len(iterations), random.seed = sample.int(100000, iterations)) %dopar% {

      # Use absolute path for the source file to avoid path issues
      source("src/environment_set.R")

      # Run foreach loop using each random seed directly
      tryCatch({
        # Compute the result with the current random seed
        result <- compute.LODO.ML(
          raw.counts, normalized, deconv, tfs.universe, clinical, trait, trait.positive,
          trait.out, out, metric, stack, feature.selection,
          file_name = file.name, seed = random.seed, return = return
        )

        # Save result as RDS file with unique identifier based on iteration (random seed)
        saveRDS(list(result = result, seed = random.seed),
                file = file.path(paste0(folder, "/ML_result_", random.seed, ".rds")))

        gc()

      }, error = function(e) {

        # Save error information as RDS file with random seed identifier
        saveRDS(list(result = NULL, error = e$message, seed = random.seed),
                file = file.path(paste0(folder, "/ML_result_", random.seed, ".rds")))

        gc()
      })

    }

    #Stop cluster after all runs
    parallel::stopCluster(cl)
    unregister_dopar() #Stop Dopar from running in the background

  }

  message("Analysis is done!")

  message("ML models are saved in Results/ folder")

}

compute.bootstrap.features.ML = function(train, test, clinical, trait, trait.positive, metric = "Accuracy", iterations, feature.selection = F, stack, workers = NULL, file.name = NULL, return = F){

  folder = paste0("Results/ML_models_",file.name)

  dir.create(file.path(getwd(), folder))

  if(is.null(iterations) == T){
    stop("No iterations specified, please set a number")
  }else{
    if(is.null(workers)==T){
      num_cores <- parallel::detectCores() - 1
    }else{
      num_cores <- workers
    }

    cl = parallel::makeCluster(num_cores) #Forking just copy the R session in its current state. - makeCluster() all must be exported (copied) to the clusters, which can add some overhead
    doParallel::registerDoParallel(cl)

    message("\nRunning ", iterations, " splits for training and test using ", num_cores, " cores")

    # Run foreach loop using each random seed directly
    foreach::foreach(iteration = seq_len(iterations), random.seed = sample.int(100000, iterations)) %dopar% {

      # Use absolute path for the source file to avoid path issues
      source("src/environment_set.R")

      # Run foreach loop using each random seed directly
      tryCatch({
        # Compute the result with the current random seed
        result <- compute.features.ML(
          train, test, clinical, trait, trait.positive,
          metric, stack, feature.selection, seed = random.seed,
          file_name = file.name, return = return
        )

        # Save result as RDS file with unique identifier based on iteration (random seed)
        saveRDS(list(result = result, seed = random.seed),
                file = file.path(paste0(folder, "/ML_result_", random.seed, ".rds")))

      }, error = function(e) {

        # Save error information as RDS file with random seed identifier
        saveRDS(list(result = NULL, error = e$message, seed = random.seed),
                file = file.path(paste0(folder, "/ML_result_", random.seed, ".rds")))
      })

    }

    #Stop cluster after all runs
    parallel::stopCluster(cl)
    unregister_dopar() #Stop Dopar from running in the background

  }

  message("Analysis is done!")

  message("ML models are saved in Results/ML_models folder")

}

compute.bootstrap.deconvolution.ML = function(deconv, clinical, trait, trait.positive, partition = 0.8, metric = "Accuracy", iterations, feature.selection = F, stack, workers = NULL, file.name = NULL, return = F){

  folder = paste0("Results/", file.name)

  dir.create(file.path(getwd(), folder))

  if(is.null(iterations) == T){
    stop("No iterations specified, please set a number")
  }else{
    if(is.null(workers)==T){
      num_cores <- parallel::detectCores() - 1
    }else{
      num_cores <- workers
    }

    cl = parallel::makeCluster(num_cores) #Forking just copy the R session in its current state. - makeCluster() all must be exported (copied) to the clusters, which can add some overhead
    doParallel::registerDoParallel(cl)

    message("\nRunning ", iterations, " splits for training and test using ", num_cores, " cores")

    # Run foreach loop using each random seed directly
    foreach::foreach(iteration = seq_len(iterations), random.seed = sample.int(100000, iterations)) %dopar% {

      # Use absolute path for the source file to avoid path issues
      source("src/environment_set.R")

      # Run foreach loop using each random seed directly
      tryCatch({
        # Compute the result with the current random seed
        result <- compute.deconvolution.ML(
          deconv, clinical, trait, trait.positive, partition,
          metric, stack, feature.selection, seed = random.seed,
          file_name = file.name, return = return
        )

        # Save result as RDS file with unique identifier based on iteration (random seed)
        saveRDS(list(result = result, seed = random.seed),
                file = file.path(paste0(folder, "/ML_result_", random.seed, ".rds")))

      }, error = function(e) {

        # Save error information as RDS file with random seed identifier
        saveRDS(list(result = NULL, error = e$message, seed = random.seed),
                file = file.path(paste0(folder, "/ML_result_", random.seed, ".rds")))
      })

    }

    #Stop cluster after all runs
    parallel::stopCluster(cl)
    unregister_dopar() #Stop Dopar from running in the background

  }

  message("Analysis is done!")

  message("ML models are saved in Results/ML_models folder")

}

compute.bootstrap.raw.deconvolution.ML = function(deconv, clinical, trait, trait.positive, partition = 0.8, metric = "Accuracy", iterations, feature.selection = F, stack, workers = NULL, file.name = NULL, return = F){

  folder = paste0("Results/ML_models_", file.name)
  dir.create(file.path(getwd(), folder))

  if(is.null(iterations) == T){
    stop("No iterations specified, please set a number")
  }else{
    if(is.null(workers)==T){
      num_cores <- parallel::detectCores() - 1
    }else{
      num_cores <- workers
    }

    cl = parallel::makeCluster(num_cores) #Forking just copy the R session in its current state. - makeCluster() all must be exported (copied) to the clusters, which can add some overhead
    doParallel::registerDoParallel(cl)

    message("\nRunning ", iterations, " splits for training and test using ", num_cores, " cores")

    # Run foreach loop using each random seed directly
    foreach::foreach(iteration = seq_len(iterations), random.seed = sample.int(100000, iterations)) %dopar% {

      # Use absolute path for the source file to avoid path issues
      source("src/environment_set.R")

      # Run foreach loop using each random seed directly
      tryCatch({
        # Compute the result with the current random seed
        result <- compute.raw.deconvolution.ML(
          deconv, clinical, trait, trait.positive, partition,
          metric, stack, feature.selection, seed = random.seed,
          file_name = file.name, return = return
        )

        # Save result as RDS file with unique identifier based on iteration (random seed)
        saveRDS(list(result = result, seed = random.seed),
                file = file.path(paste0(folder, "/ML_result_", random.seed, ".rds")))

      }, error = function(e) {

        # Save error information as RDS file with random seed identifier
        saveRDS(list(result = NULL, error = e$message, seed = random.seed),
                file = file.path(paste0(folder, "/ML_result_", random.seed, ".rds")))
      })

    }

    #Stop cluster after all runs
    parallel::stopCluster(cl)
    unregister_dopar() #Stop Dopar from running in the background

  }

  message("Analysis is done!")

  message("ML models are saved in Results/ML_models folder")

}

#Get boxplot from ML models from a single folder
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

#Get a sequence of boxplots from ML models from different folders
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

compute_cv_accuracy = function(models, file_name = NULL, base_models = F, return = T){

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


plot_cv_metrics <- function(models, file_name = NULL) {

  # Helper function to process metrics: Extract resample matrix
  process_metric <- function(metric_name) {
    metric_data <- map_dfr(names(models), function(model_name) {
      models[[model_name]]$resample %>%
        dplyr::mutate(model = model_name)
    })

    # Compute median and SD for the metric
    res_metric <- metric_data %>%
      dplyr::group_by(model) %>%
      dplyr::summarise(
        Mean = median(.data[[metric_name]], na.rm = TRUE),
        SD = sd(.data[[metric_name]], na.rm = TRUE)
      ) %>%
      dplyr::arrange(desc(Mean))

    # Plot results
    grDevices::pdf(paste0("Results/", metric_name, "_CV_methods_", file_name, ".pdf"), width = 10)
    print(
      ggplot2::ggplot(res_metric, ggplot2::aes(x = model, y = Mean, fill = model)) +
        ggplot2::geom_bar(stat = "identity", position = ggplot2::position_dodge(), width = 0.6) +
        ggplot2::geom_errorbar(ggplot2::aes(ymin = Mean - SD, ymax = Mean + SD),
                      width = 0.2, position = ggplot2::position_dodge(0.6)) +
        ggplot2::labs(title = paste("Performance of Models -", metric_name),
             x = "Model",
             y = paste("Median", metric_name)) +
        ggplot2::theme_minimal() +
        ggplot2::theme(legend.position = "none") +
        ggplot2::scale_y_continuous(breaks = seq(0, 1, by = 0.05))
    )
    grDevices::dev.off()
  }

  # Process both F1 and MCC metrics
  process_metric("F1")
  process_metric("MCC")
}

compute_cv_AUC = function(models, file_name = NULL, base_models = F, AUC_type = "AUROC", return = T){

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

  if(base_models == T){
    cat("Choosing base models for stacking.......................................\n\n")
    base_models = choose_base_models(models, metric = AUC_type)
    cat("Models chosen are:", paste0(base_models, collapse = ", "), "\n\n")
    return(list("AUROC" = res_auroc, "AUPRC" = res_auprc, "Top_model" = top_model, "Base_models" = base_models))
  }else{
    return(list("AUROC" = res_auroc, "AUPRC" = res_auprc, "Top_model" = top_model))
  }

}

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

#Matthews Correlation Coefficient
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

#Take sensitivities values based on values of specificities
get_sensitivity = function(x, data){
  data %>%
    dplyr::filter(specificity - x >= 0)%>% #Take specificity values above threshold x
    dplyr::top_n(sensitivity, n=1) %>% #Take highest sensitivity from that threshold
    dplyr::mutate(specificity = x, fpr = 1-x) %>% #Define sensitivity based on the specified threshold
    dplyr::distinct() #If multiple thresholds have same sensitivity values take only one
}

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
  meta_importance = varImp(meta_learner, scale = F)$importance %>%
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

compute.prediction.stacked = function(super.learner, test_data, target, ml.models, base.models){

  #Learning from simple meta-learner
  base_predictions = list()
  for (i in 1:length(base.models)) {
    base_predictions[[i]] = stats::predict(ml.models[[base.models[i]]], test_data, type = "prob")$yes
    names(base_predictions)[i] = base.models[i]
  }

  base_predictions = do.call(cbind, base_predictions)

  prediction_simple = data.frame(stats::predict(super.learner, base_predictions, type = "prob"))

  # #Learning from simple meta-learner
  # all_predictions = list()
  # for (i in 1:length(ml.models)) {
  #   all_predictions[[i]] = predict(ml.models[[i]], test_data, type = "prob")$yes
  #   names(all_predictions)[i] = names(ml.models)[i]
  # }
  #
  # all_predictions = do.call(cbind, all_predictions)
  #
  # prediction_all = data.frame(predict(super.learner[["all"]], all_predictions, type = "prob"))
  #
  #Metrics

  #Meta-learner simple
  sens_spec_simple = get_sensitivity_specificity(prediction_simple, target, "Meta-learner_simple")
  auroc_simple = calculate_auroc(sens_spec_simple$fpr, sens_spec_simple$Sensitivity)
  auprc_simple = calculate_auprc(sens_spec_simple$Recall, sens_spec_simple$Precision)
  #Meta-learner all
  # sens_spec_all = get_sensitivity_specificity(prediction_all, target, "Meta-learner_all")
  # auc_all = calculate_auc(sens_spec_all$fpr, sens_spec_all$sensitivity)

  #Not returning all (discarded)

  return(list(Metrics = sens_spec_simple, AUC = list("AUROC" = auroc_simple, "AUPRC" = auprc_simple), Predictions = prediction_simple))

}

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

calculate_precision <- function(metrics, target) {
  confusion_values <- calculate_confusion_values(metrics, target)
  TP <- confusion_values$TP
  FP <- confusion_values$FP

  precision <- TP / (TP + FP)

  return(precision)
}

calculate_recall <- function(metrics, target) {
  confusion_values <- calculate_confusion_values(metrics, target)
  TP <- confusion_values$TP
  FN <- confusion_values$FN

  # Calculate recall (sensitivity)
  recall <- TP / (TP + FN)

  return(recall)
}

#Get ROC and Recall curve
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

identify.cell.signatures = function(model, deconvolution, var_importance, n_top = 20, sign, file.name){

  # Extract deconvolution subgroups composition per ML model
  deconv_subgroups = list()
  contador = 1
  subgroups = model[["Deconvolution_subgroups"]][["Deconvolution groups - Linear-based correlation"]]
  for (i in 1:length(subgroups)) {
    if(length(subgroups[[i]])!=0){ #Whether a specific cell type does not contains subgroups
      for (j in 1:length(subgroups[[i]])) {
        deconv_subgroups[[contador]] = subgroups[[i]][[j]]
        names(deconv_subgroups)[contador] = names(subgroups[[i]])[j]
        contador = contador + 1
      }
    }
  }

  # Take top features based on variable importance
  top <- var_importance %>%
    tidyr::pivot_longer(cols = dplyr::everything(), names_to = "feature", values_to = "shap_value") %>%
    dplyr::group_by(feature) %>%
    dplyr::summarise(mean_shap = mean(shap_value, na.rm = TRUE)) %>% #Average SHAP values across samples
    dplyr::mutate(direction = ifelse(mean_shap > 0, "Increase", "Decrease")) %>% #Give direction
    dplyr::filter(direction == sign) %>%
    {
      if (nrow(.) < n_top) {
        warning("Not enough features for selecting n_top = ", n_top, " features in model. Selecting all available features.\n")
        .  # Use all available rows
      } else {
        dplyr::top_n(., n_top, wt = mean_shap)  # Select the top n_top features
      }
    } %>%
    dplyr::arrange(desc(mean_shap)) %>% #Order in decreasing order
    dplyr::top_n(n_top, wt = mean_shap) %>% #Select n_top features
    dplyr::pull(feature)

  # Initialize presence matrix for cells
  cells_types = extract_cells(colnames(deconvolution))
  presence_matrix <- data.frame(matrix(data = 0, nrow = length(top), ncol = length(cells_types)))
  colnames(presence_matrix) <- cells_types
  rownames(presence_matrix) = top

  contador = 1 #Iterator for features inside each top_features
  cell.groups_top = list() #Vector to save top cell groups

  # Extract cell composition from top features
  for (feature in top) {
    composition = model[["Cell_groups"]][[2]][[feature]]
    cell.groups_top[[contador]] = composition
    contador = contador + 1
  }

  # Update presence_matrix for this model
  row_index_cell <- 1 #Initialize row number
  for (j in seq_along(cell.groups_top)) {
    features <- cell.groups_top[[j]] #Extract cell group j from ML model i
    for (cell_feature in features) { #Iterate over features in cell group j ML model i
      cells <- get_all_cells(cell_feature, deconv_subgroups) #Use recursive function to extract all nested cell features from different subgroup levels
      cells_types = extract_cells(cells)
      presence_matrix[row_index_cell, cells_types] <- 1 #Set 1 if cell feature is present in subgroup
    }
    row_index_cell <- row_index_cell + 1
  }

  # Remove cells with no presence in any group
  remove = which(colSums(presence_matrix) == 0) #Identify cell features with no presence in any group
  if(length(remove)>0){
    presence_matrix = presence_matrix[,-remove] #Remove features
  }

  # Identify clusters of cell combinations
  jaccard_dist = stats::dist(presence_matrix, method="binary") #Jaccard distance
  hc <- stats::hclust(jaccard_dist, method = "ward.D2") #Hierarchical clustering of distance
  matrica <- as.matrix(jaccard_dist) #Convert distance object to matrix
  silhouette <- factoextra::fviz_nbclust(matrica, FUNcluster = hcut, method = "silhouette", k.max = 6) #Identify k clusters from matrix
  k_cluster = as.numeric(silhouette$data$clusters[which.max(silhouette$data$y)]) #Extract k value
  sub_grp <- stats::cutree(hc, k = k_cluster) #Obtain k cluster composition

  p = pheatmap::pheatmap(matrica,
                         cluster_rows = hc,
                         cluster_cols = hc,
                         show_rownames = TRUE,
                         show_colnames = TRUE,
                         fontsize = 8,
                         border_color = NA,
                         color = grDevices::hcl.colors(20, palette = "PRGn"),
                         main = "Jaccard Distance Matrix Heatmap")

  grDevices::pdf(paste0("Results/Cell_combinations_jaccard_", file.name, ".pdf"), width = 8)
  print(p)
  grDevices::dev.off()

  #### Determine which features are important for the clustering
  baseline_quality <- compute_silhouette(sub_grp, matrica) #baseline quality of default clustering
  feature_impacts <- numeric(ncol(presence_matrix)) #Initialize vector to store feature impacts

  ### Features permutation to verify clustering quality

  n_bootstrap = 100 #Number of bootstraps
  for (feature_idx in seq_len(ncol(presence_matrix))) {

    # Perform bootstrapping with different seeds
    for (i in 1:n_bootstrap) {
      set.seed(sample.int(100000, 1)) # Set a different random seed for each bootstrap iteration

      bootstrap_impacts <- numeric(n_bootstrap) # Initialize vector to hold impacts across bootstrap iterations

      permuted_matrix <- presence_matrix # Obtain default presence matrix
      permuted_matrix[, feature_idx] <- sample(permuted_matrix[, feature_idx]) # Permute the feature values (only the current feature is shuffled)

      # Recompute Jaccard distance and clustering for the permuted matrix
      permuted_dist = stats::dist(permuted_matrix, method="binary")
      permuted_hc <- stats::hclust(permuted_dist, method = "ward.D2")
      permuted_sub_grp <- stats::cutree(permuted_hc, k = k_cluster)

      # Recompute clustering quality for the permuted matrix
      permuted_quality <- compute_silhouette(permuted_sub_grp, as.matrix(permuted_dist))

      # Measure impact as the change in clustering quality
      bootstrap_impacts[i] <- baseline_quality - permuted_quality
    }

    # Store the mean feature impact across all bootstrap iterations
    feature_impacts[feature_idx] <- mean(bootstrap_impacts)

  }

  # Create a dataframe of feature impacts
  feature_importance <- data.frame(
    Feature = colnames(presence_matrix),
    Impact = feature_impacts
  )

  # Sort features by impact
  feature_importance <- feature_importance[order(-feature_importance$Impact), ]

  # Plot feature importance
  p = ggplot2::ggplot(feature_importance, ggplot2::aes(x = stats::reorder(Feature, -Impact), y = Impact)) +
    ggplot2::geom_bar(stat = "identity", fill = "steelblue") +
    ggplot2::theme_minimal() +
    ggplot2::labs(title = "Feature Importance Based on Clustering Impact",
                  x = "Feature",
                  y = "Impact on Clustering Quality") +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))

  grDevices::pdf(paste0("Results/Clustering_feature_importance_", file.name, ".pdf"), width = 8)
  print(p)
  grDevices::dev.off()

  # Keep only features with a positive impact in clustering quality
  features = feature_importance %>%
    dplyr::filter(Impact > 0) %>%
    dplyr::pull(Feature)

  # Keep only features with positive impact in clustering
  presence_matrix_important = presence_matrix[,features]

  # Remove cell groups with 0 or 1 cells only (because all the features were discard/no important for clustering) and
  idx <- which(rowSums(presence_matrix_important) %in% c(0, 1))

  if(length(idx)>0){
    presence_matrix_important = presence_matrix_important[-idx,]
  }

  ############################################################################# Plot presence scores per feature

  #Calculate scores of presence for each feature
  scores = colSums(presence_matrix_important)/nrow(presence_matrix_important)
  top_scores_df <- utils::stack(scores)

  #Plot scores of the cell features across ML models
  p = ggplot2::ggplot(top_scores_df, ggplot2::aes(y = values, x = stats::reorder(ind, -values, decreasing = F))) +
    ggplot2::geom_bar(stat = "identity", fill = "skyblue", width = 0.6) +
    ggplot2::theme_minimal() +
    ggplot2::labs(title = paste0("Cell features presence scores across ", nrow(presence_matrix_important), " predictive cell groups"),
         x = "Cell features",
         y = "Presence scores") +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
                   plot.margin = ggplot2::margin(t = 10, r = 10, b = 30, l = 60)) +
    ggplot2::scale_y_continuous(labels = scales::comma)

  grDevices::pdf(paste0("Results/Cell_features_scores_", file.name, ".pdf"), width = 8)
  print(p)
  grDevices::dev.off()

  #Extract TF modules
  colors_groups = c()
  for(i in top){
    colors_groups = c(colors_groups, extract_colors(names(model$TF_network$`TFs per module`),i))
  }
  top_colors_df = data.frame(prop.table(table(colors_groups)))

  #Plot scores of the colors across ML models
  p = ggplot2::ggplot(top_colors_df, ggplot2::aes(y = Freq, x = stats::reorder(colors_groups, -Freq, decreasing = F))) +
    ggplot2::geom_bar(stat = "identity", fill = "skyblue", width = 0.6) +
    ggplot2::theme_minimal() +
    ggplot2::labs(title = paste0("TF modules presence scores across ", nrow(presence_matrix_important), " predictive cell groups"),
         x = "TF modules",
         y = "Presence scores") +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
                   plot.margin = ggplot2::margin(t = 10, r = 10, b = 30, l = 60)) +
    ggplot2::scale_y_continuous(labels = scales::comma)

  grDevices::pdf(paste0("Results/TF_modules_scores_", file.name, ".pdf"), width = 8)
  print(p)
  grDevices::dev.off()

}

identify.predictive.cell.group = function(folder_path, deconvolution, AUC = 0.7, metric = "ROC", n_top = 10, scores_threshold = 0.5, n_cells = 5, file.name){

  #Read ML models
  res <- list.files(folder_path, pattern = "\\.rds$", full.names = TRUE)
  cell_subgroups <- list() #List to fill cell subgroups
  res.cell.groups.top = list()
  res.tfs.top = list()
  ml_models <- 0 #Iterator to fill the number of valid ML models

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
      valid_indices[i] <- auc_value >= AUC
    }
  }

  # Filter files based on valid indices
  res <- res[valid_indices]
  if (length(res) == 0) {
    stop("No ML models with AUC above ", AUC, " were found. Try with another value.")
  }else{
    cat(length(res), "ML models found with AUC >", AUC, "!\n")
  }

  # Initialize presence matrix
  cells_types = extract_cells(colnames(deconvolution))
  presence_matrix <- data.frame(matrix(data = 0, nrow = n_top * length(res), ncol = length(cells_types)))
  colnames(presence_matrix) <- cells_types

  # Initialize TF matrix
  tfs_matrix = data.frame(matrix(data = 0, nrow = n_top*length(res), ncol = ncol(tfs)))
  colnames(tfs_matrix) <- colnames(tfs)

  row_index_cell <- 1
  row_index_tfs <- 1

  for (file in res) {

    model <- readRDS(file)

    #Extract deconvolution subgroups composition for each ML model
    deconv_subgroups = list()
    contador = 1
    subgroups = model[["result"]][["Deconvolution_subgroups"]]
    for (i in 1:length(subgroups)) {
      if(length(subgroups[[i]])!=0){ #Whether a specific cell type does not contains subgroups
        for (j in 1:length(subgroups[[i]])) {
          deconv_subgroups[[contador]] = subgroups[[i]][[j]]
          names(deconv_subgroups)[contador] = names(subgroups[[i]])[j]
          contador = contador + 1
        }
      }
    }
    cell_subgroups[[ml_models + 1]] <- deconv_subgroups #Save subgroups from each ML model

    #Align feature importance
    model[["result"]][["Variable_importance"]] = feature.importance.alignment(model)

    #Extract network
    tf_network = model[["result"]][["TF_network"]]

    #Extract top cell groups for each model
    top <- model[["result"]][["Variable_importance"]] %>%
      dplyr::filter(final_importance > 0) %>% #Filter features with no-importance/importance to classify positive class
      {
        if (nrow(.) < n_top) {
          warning("Not enough features for selecting n_top = ", n_top, " features in model ", file, ". Selecting all available features.\n")
          .  # Use all available rows
        } else {
          dplyr::top_n(., n_top, wt = final_importance)  # Select the top n_top features
        }
      } %>%
      dplyr::pull(features)

    n_features = length(top)
    contador = 1 #Iterator for features inside each top_features
    cell.groups_top = list()
    tfs_top = list()

    #Extract cell composition from top features
    for (feature in top) {

      #Cells composition
      composition = model[["result"]][["Cell_groups"]][[2]][[feature]]
      cell.groups_top[[contador]] = composition

      #TFs composition
      colors = extract_colors(unique(tf_network$`TFs colors`), feature)
      tfs_top[[contador]] = do.call(c, tf_network$`TFs per module`[colors])

      contador = contador + 1
    }

    res.cell.groups.top[[ml_models + 1]] <- cell.groups_top
    res.tfs.top[[ml_models + 1]] <- tfs_top

    #Check if available features were indeed n_top (sometimes not enough features for selecting n_top), if not adjust the matrix dynamically
    if(n_features < n_top){
      rows_to_delete = n_top - n_features
      rows = seq(nrow(presence_matrix) - rows_to_delete + 1, nrow(presence_matrix))
      presence_matrix = presence_matrix[-rows,]
      tfs_matrix = tfs_matrix[-rows,]
    }else if(n_features > n_top){
      rows_to_add = n_features - n_top
      extra_rows <- data.frame(matrix(data = 0, nrow = rows_to_add, ncol = ncol(presence_matrix)))
      colnames(extra_rows) = colnames(presence_matrix)
      presence_matrix <- rbind(presence_matrix, extra_rows)
      tfs_matrix = rbind(tfs_matrix, extra_rows)
    }

    # Update presence_matrix for this model
    for (j in seq_along(cell.groups_top)) {
      features <- cell.groups_top[[j]] #Extract cell group j from ML model i
      for (cell_feature in features) { #Iterate over features in cell group j ML model i
        cells <- get_all_cells(cell_feature, deconv_subgroups) #Use recursive function to extract all nested cell features from different subgroup levels
        cells_types = extract_cells(cells)
        presence_matrix[row_index_cell, cells_types] <- 1 #Set 1 if cell feature is present in subgroup
      }
      row_index_cell <- row_index_cell + 1
    }

    # Update tfs_matrix for this model
    for (j in seq_along(tfs_top)) {
      markers <- tfs_top[[j]] #Extract cell group j from ML model i
      tfs_matrix[row_index_tfs, markers] <- 1 #Set 1 if TF is present in cell group top
      row_index_tfs <- row_index_tfs + 1
    }

    ml_models <- ml_models + 1

  }

  #### Identify clusters of cell combinations
  remove = which(colSums(presence_matrix) == 0) #Identify cell features with no presence in any group
  if(length(remove)>0){
    presence_matrix = presence_matrix[,-remove] #Remove features
  }

  remove = which(colSums(tfs_matrix) == 0) #Identify TFs with no presence in any group
  if(length(remove)>0){
    tfs_matrix = tfs_matrix[,-remove] #Remove features
  }


  # Cell features
  jaccard_dist = stats::dist(presence_matrix, method="binary") #Jaccard distance
  hc <- stats::hclust(jaccard_dist, method = "ward.D2") #Hierarchical clustering of distance
  matrica <- as.matrix(jaccard_dist) #Convert distance object to matrix
  silhouette <- factoextra::fviz_nbclust(matrica, FUNcluster = hcut, method = "silhouette", k.max = 6) #Identify k clusters from matrix
  k_cluster = as.numeric(silhouette$data$clusters[which.max(silhouette$data$y)]) #Extract k value
  sub_grp <- stats::cutree(hc, k = k_cluster) #Obtain k cluster composition

  p = pheatmap::pheatmap(matrica,
                         cluster_rows = hc,
                         cluster_cols = hc,
                         show_rownames = TRUE,
                         show_colnames = TRUE,
                         fontsize = 8,
                         border_color = NA,
                         color = grDevices::hcl.colors(20, palette = "PRGn"),
                         main = "Jaccard Distance Matrix Heatmap")

  # pdf(paste0("Results/Cell_combinations_jaccard_", file.name, ".pdf"), width = 8)
  # print(p)
  # dev.off()

  #### Determine which features are important for the clustering
  baseline_quality <- compute_silhouette(sub_grp, matrica) #baseline quality of default clustering
  feature_impacts <- numeric(ncol(presence_matrix)) #Initialize vector to store feature impacts

  ### Features permutation to verify clustering quality

  n_bootstrap = 100 #Number of bootstraps
  for (feature_idx in seq_len(ncol(presence_matrix))) {

    # Perform bootstrapping with different seeds
    for (i in 1:n_bootstrap) {
      set.seed(sample.int(100000, 1)) # Set a different random seed for each bootstrap iteration

      bootstrap_impacts <- numeric(n_bootstrap) # Initialize vector to hold impacts across bootstrap iterations

      permuted_matrix <- presence_matrix # Obtain default presence matrix
      permuted_matrix[, feature_idx] <- sample(permuted_matrix[, feature_idx]) # Permute the feature values (only the current feature is shuffled)

      # Recompute Jaccard distance and clustering for the permuted matrix
      permuted_dist = stats::dist(permuted_matrix, method="binary")
      permuted_hc <- stats::hclust(permuted_dist, method = "ward.D2")
      permuted_sub_grp <- stats::cutree(permuted_hc, k = k_cluster)

      # Recompute clustering quality for the permuted matrix
      permuted_quality <- compute_silhouette(permuted_sub_grp, as.matrix(permuted_dist))

      # Measure impact as the change in clustering quality
      bootstrap_impacts[i] <- baseline_quality - permuted_quality
    }

    # Store the mean feature impact across all bootstrap iterations
    feature_impacts[feature_idx] <- mean(bootstrap_impacts)

  }

  # Create a dataframe of feature impacts
  feature_importance <- data.frame(
    Feature = colnames(presence_matrix),
    Impact = feature_impacts
  )

  # Sort features by impact
  #feature_importance <- feature_importance[order(-feature_importance$Impact), ]

  # Plot feature importance
  # ggplot(feature_importance, aes(x = reorder(Feature, -Impact), y = Impact)) +
  #   geom_bar(stat = "identity", fill = "steelblue") +
  #   theme_minimal() +
  #   labs(title = "Feature Importance Based on Clustering Impact",
  #        x = "Feature",
  #        y = "Impact on Clustering Quality") +
  #   theme(axis.text.x = element_text(angle = 45, hjust = 1))

  # Keep only features with a positive impact in clustering quality
  features = feature_importance %>%
    dplyr::filter(Impact > 0) %>%
    dplyr::pull(Feature)

  # Keep only features with positive impact in clustering
  presence_matrix_important = presence_matrix[,features]

  # Remove cell groups with 0 or 1 cells only (because all the features were discard/no important for clustering) and
  idx <- which(rowSums(presence_matrix_important) %in% c(0, 1))

  if(length(idx)>0){
    presence_matrix_important = presence_matrix_important[-idx,]
  }

  ############################################################################# Plot presence scores per feature

  #Calculate scores of presence for each feature
  scores = colSums(presence_matrix_important)/nrow(presence_matrix_important)

  #Select top scores above threshold
  top_scores_df <- utils::stack(scores)

  #Plot scores of the cell features across ML models
  p = ggplot2::ggplot(top_scores_df, ggplot2::aes(y = values, x = reorder(ind, -values, decreasing = F))) +
    ggplot2::geom_bar(stat = "identity", fill = "skyblue", width = 0.6) +
    ggplot2::theme_minimal() +
    ggplot2::labs(title = paste0("Cell features presence scores across ", nrow(presence_matrix_important), " predictive cell groups"),
                  subtitle = paste0(ml_models, " ML models (AUC > ", AUC, ")"),
                  x = "Cell features",
                  y = "Presence scores") +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
                   plot.margin = ggplot2::margin(t = 10, r = 10, b = 30, l = 60)) +
    ggplot2::scale_y_continuous(labels = scales::comma)

  grDevices::pdf(paste0("Results/Cell_features_scores_", file.name, ".pdf"), width = 8)
  print(p)
  grDevices::dev.off()

  cat("Picking", n_cells, "cells types.........................................................\n")

  # Pick top_n cells
  cells = top_scores_df %>%
    dplyr::arrange(desc(values)) %>%
    dplyr::slice(1:n_cells) %>%
    dplyr::pull(ind) %>%
    as.character()

  print(cells)

  #Calculate scores of presence for each TF
  scores_tfs = colSums(tfs_matrix)/nrow(tfs_matrix)

  #Select top scores above threshold
  top_scores_df <- utils::stack(scores_tfs[scores_tfs > 0.8])

  entrz <- AnnotationDbi::select(org.Hs.eg.db, keys = as.character(top_scores_df$ind), columns = "ENTREZID", keytype = "SYMBOL") #Change to EntrezID

  paths_df = paths[,1:2]
  colnames(paths_df) = c("TERM", "GENE")

  # Perform ORA using enricher
  ora_results <- clusterProfiler::enricher(
    gene = as.character(top_scores_df$ind),
    TERM2GENE = paths_df,
    pvalueCutoff = 0.05
  )

  grDevices::pdf(paste0("Results/Enrichment_TFs_CGs_", file.name, ".pdf"), width = 5)
  graphics::dotplot(ora_results) + ggtitle("ORA - Progeny pathways")
  grDevices::dev.off()

  ############################################################################# Plot frequency scores per cell group (to be checked)

  # cell.composition = list()
  # for (i in 1:length(res.cell.groups.top)) {
  #   cell.composition = append(cell.composition, res.cell.groups.top[[i]])
  # }

  # #Calculate frequency of equal cell groups composition
  # row_strings = apply(presence_matrix_important, 1, paste, collapse = "|")
  # composition_counts = table(row_strings)
  #
  # #Extract composition of sorted cell groups and fill frequency matrix
  # composition = list()
  # frequency = data.frame(matrix(0, nrow = length(composition_counts), ncol = 2))
  # colnames(frequency) = c("Cell.group", "Frequency")
  #
  # for (i in 1:length(composition_counts)) {
  #   vec = as.numeric(unlist(strsplit(names(composition_counts[i]), "|", fixed = T)))
  #   composition[[i]] = colnames(presence_matrix_important)[vec == 1]
  #   names(composition)[i] = paste0("Cell.group.composition_", i)
  #   frequency[i,2] = composition_counts[[i]]/nrow(presence_matrix_important)
  #   frequency[i,1] = paste0("Cell.group.composition_", i)
  # }
  #
  # # Discard groups below threshold
  # frequent_groups = frequency %>%
  #   filter(Frequency > frequency_threshold)
  #
  # if(nrow(frequent_groups)>0){
  #
  #   #Plot frequency of the cell groups across ML models
  #   p = ggplot(frequent_groups, aes(x = Frequency, y = reorder(Cell.group, -Frequency, decreasing = T))) +
  #     geom_bar(stat = "identity", fill = "skyblue") +
  #     theme_minimal() +
  #     labs(title = paste0("Cell groups frequency across ", nrow(presence_matrix_important), " predictive cell groups"),
  #          x = "Frequency",
  #          y = "Cell groups") +
  #     theme(axis.text.y = element_text(angle = 0, hjust = 1)) +  # Ensure y-axis labels are horizontal
  #     scale_x_continuous(labels = scales::comma)  # Add comma formatting for x-axis labels if needed
  #
  #   pdf(paste0("Results/Cell_groups_frequency_", file.name, ".pdf"))
  #   print(p)
  #   dev.off()
  #
  #   cat(nrow(frequent_groups), "cell groups were found appearing more than", frequency_threshold, "times.\n")
  #
  # }else{
  #   cat("There is no cell groups appearing more than", frequency_threshold, "times. Try with another value.\n")
  # }
  #
  return(cells)
}

identify.predictive.deconvolution = function(folder_path, deconvolution, AUC = 0.7, metric = "ROC", n_top = 10, scores_threshold = 0.5, n_cells = 5, file.name){

  #Read ML models
  res <- list.files(folder_path, pattern = "\\.rds$", full.names = TRUE)
  cell_subgroups <- list() #List to fill cell subgroups
  res.cell.groups.top = list()
  ml_models <- 0 #Iterator to fill the number of valid ML models

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
      valid_indices[i] <- auc_value >= AUC
    }
  }

  # Filter files based on valid indices
  res <- res[valid_indices]
  if (length(res) == 0) {
    stop("No ML models with AUC above ", AUC, " were found. Try with another value.")
  }else{
    cat(length(res), "ML models found with AUC >", AUC, "!\n")
  }

  # Initialize presence matrix
  cells_types = extract_cells(colnames(deconvolution))
  presence_matrix <- data.frame(matrix(data = 0, nrow = n_top * length(res), ncol = length(cells_types)))
  colnames(presence_matrix) <- cells_types

  row_index <- 1

  for (file in res) {

    model <- readRDS(file)

    #Extract deconvolution subgroups composition for each ML model
    deconv_subgroups = list()
    contador = 1
    subgroups = model[["result"]][["Deconvolution_subgroups"]]
    for (i in 1:length(subgroups)) {
      if(length(subgroups[[i]])!=0){ #Whether a specific cell type does not contains subgroups
        for (j in 1:length(subgroups[[i]])) {
          deconv_subgroups[[contador]] = subgroups[[i]][[j]]
          names(deconv_subgroups)[contador] = names(subgroups[[i]])[j]
          contador = contador + 1
        }
      }
    }
    cell_subgroups[[ml_models + 1]] <- deconv_subgroups #Save subgroups from each ML model

    #Align feature importance
    model[["result"]][["Variable_importance"]] = feature.importance.alignment.deconvolution(model, dt[[1]])

    #Extract top cell groups for each model
    top <- model[["result"]][["Variable_importance"]] %>%
      dplyr::filter(final_importance > 0) %>% #Filter features with no-importance/importance to classify positive class
      {
        if (nrow(.) < n_top) {
          warning("Not enough features for selecting n_top = ", n_top, " features in model ", file, ". Selecting all available features.\n")
          .  # Use all available rows
        } else {
          dplyr::top_n(., n_top, wt = final_importance)  # Select the top n_top features
        }
      } %>%
      dplyr::pull(features)

    n_features = length(top)

    #Check if available features were indeed n_top, if not adjust the matrix dynamically
    if(n_features < n_top){
      rows_to_delete = n_top - n_features
      rows = seq(nrow(presence_matrix) - rows_to_delete + 1, nrow(presence_matrix))
      presence_matrix = presence_matrix[-rows,]
    }else if(n_features > n_top){
      rows_to_add = n_features - n_top
      extra_rows <- data.frame(matrix(data = 0, nrow = rows_to_add, ncol = ncol(presence_matrix)))
      colnames(extra_rows) = colnames(presence_matrix)
      presence_matrix <- rbind(presence_matrix, extra_rows)
    }

    # Update presence_matrix for this model
    for (cell_feature in top) { #Iterate over features in cell group j ML model i
        cells <- get_all_cells(cell_feature, deconv_subgroups) #Use recursive function to extract all nested cell features from different subgroup levels
        cells_types = extract_cells(cells)
        presence_matrix[row_index, cells_types] <- 1 #Set 1 if cell feature is present in subgroup
    }

    row_index <- row_index + 1

    ml_models <- ml_models + 1

  }

  #### Identify clusters of cell combinations
  remove = which(colSums(presence_matrix) == 0) #Identify cell features with no presence in any group
  if(length(remove)>0){
    presence_matrix = presence_matrix[,-remove] #Remove features
  }

  jaccard_dist = stats::dist(presence_matrix, method="binary") #Jaccard distance
  hc <- stats::hclust(jaccard_dist, method = "ward.D2") #Hierarchical clustering of distance
  matrica <- as.matrix(jaccard_dist) #Convert distance object to matrix
  silhouette <- factoextra::fviz_nbclust(matrica, FUNcluster = hcut, method = "silhouette", k.max = 6) #Identify k clusters from matrix
  k_cluster = as.numeric(silhouette$data$clusters[which.max(silhouette$data$y)]) #Extract k value
  sub_grp <- stats::cutree(hc, k = k_cluster) #Obtain k cluster composition

  p = pheatmap::pheatmap(matrica,
               cluster_rows = hc,
               cluster_cols = hc,
               show_rownames = TRUE,
               show_colnames = TRUE,
               fontsize = 8,
               border_color = NA,
               color = grDevices::hcl.colors(20, palette = "PRGn"),
               main = "Jaccard Distance Matrix Heatmap")

  # pdf(paste0("Results/Cell_combinations_jaccard_", file.name, ".pdf"), width = 8)
  # print(p)
  # dev.off()

  #### Determine which features are important for the clustering
  baseline_quality <- compute_silhouette(sub_grp, matrica) #baseline quality of default clustering
  feature_impacts <- numeric(ncol(presence_matrix)) #Initialize vector to store feature impacts

  ### Features permutation to verify clustering quality

  n_bootstrap = 100 #Number of bootstraps
  for (feature_idx in seq_len(ncol(presence_matrix))) {

    # Perform bootstrapping with different seeds
    for (i in 1:n_bootstrap) {
      seed = sample.int(100000, 1) # Set a different random seed for each bootstrap iteration

      bootstrap_impacts <- numeric(n_bootstrap) # Initialize vector to hold impacts across bootstrap iterations

      permuted_matrix <- presence_matrix # Obtain default presence matrix
      permuted_matrix[, feature_idx] <- sample(permuted_matrix[, feature_idx]) # Permute the feature values (only the current feature is shuffled)

      # Recompute Jaccard distance and clustering for the permuted matrix
      permuted_dist = stats::dist(permuted_matrix, method="binary")
      permuted_hc <- stats::hclust(permuted_dist, method = "ward.D2")
      permuted_sub_grp <- stats::cutree(permuted_hc, k = k_cluster)

      # Recompute clustering quality for the permuted matrix
      permuted_quality <- compute_silhouette(permuted_sub_grp, as.matrix(permuted_dist))

      # Measure impact as the change in clustering quality
      bootstrap_impacts[i] <- baseline_quality - permuted_quality
    }

    # Store the mean feature impact across all bootstrap iterations
    feature_impacts[feature_idx] <- mean(bootstrap_impacts)

  }

  # Create a dataframe of feature impacts
  feature_importance <- data.frame(
    Feature = colnames(presence_matrix),
    Impact = feature_impacts
  )

  # Keep only features with a positive impact in clustering quality
  features = feature_importance %>%
    dplyr::filter(Impact > 0) %>%
    dplyr::pull(Feature)

  # Keep only features with positive impact in clustering
  presence_matrix_important = presence_matrix[,features]

  # Remove cell groups with 0 and 1 features (because all the features were discard/no important for clustering) and
  idx <- which(rowSums(presence_matrix_important) %in% c(0, 1))

  if(length(idx)>0){
    presence_matrix_important = presence_matrix_important[-idx,]
  }

  ############################################################################# Plot presence scores per feature

  #Calculate scores of presence for each feature
  scores = colSums(presence_matrix_important)/nrow(presence_matrix_important)

  #Select top scores above threshold
  top_scores_df <- utils::stack(scores)

  #Plot scores of the cell features across ML models
  p = ggplot2::ggplot(top_scores_df, ggplot2::aes(y = values, x = stats::reorder(ind, -values, decreasing = F))) +
    ggplot2::geom_bar(stat = "identity", fill = "skyblue", width = 0.6) +
    ggplot2::theme_minimal() +
    ggplot2::labs(title = paste0("Cell features presence scores across ", nrow(presence_matrix_important), " predictive cell groups"),
         subtitle = paste0(ml_models, " ML models (AUC > ", AUC, ")"),
         x = "Cell features",
         y = "Presence scores") +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
          plot.margin = ggplot2::margin(t = 10, r = 10, b = 30, l = 60)) +
    ggplot2::scale_y_continuous(labels = scales::comma)

  grDevices::pdf(paste0("Results/Cell_features_scores_", file.name, ".pdf"), width = 8)
  print(p)
  grDevices::dev.off()

  cat("Picking", n_cells, "cells types.........................................................\n")

  # Pick top_n cells
  cells = top_scores_df %>%
    dplyr::arrange(desc(values)) %>%
    dplyr::slice(1:n_cells) %>%
    dplyr::pull(ind) %>%
    as.character()

  print(cells)

  return(cells)
}


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

feature.importance.alignment.deconvolution = function(model, deconv){

  #positive and negative class are defined based on the factor() from trait (this have been defined in compute.ML() function already)
  importance = model[["result"]][["Variable_importance"]]
  trait = model[["result"]][["Model"]][["trainingData"]]

  features_values = deconv[rownames(deconv) %in% rownames(trait), ]


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

compute.independent.test = function(ml_trained, deconvolution_test, coldata_test, trait, trait.positive, maximize_threshold = "Accuracy", file.name = NULL){

  #### Testing set
  testing_set = compute.test.set(ml_trained$Deconvolution_subgroups, ml_trained$Cell_groups, names(ml_trained[["Cell_groups"]][[2]]), deconvolution_test)
  target = coldata_test %>%
    dplyr::mutate(target = ifelse(coldata_test[,trait] == trait.positive, "yes", "no")) %>%
    dplyr::pull(target) %>%
    factor(levels = c("no", "yes"))

  ##### ML prediction
  model = ml_trained$Model$Model
  #prediction = compute.prediction.stacked(model, testing_set, target, ml_trained$Model$ML_models, ml_trained$Model$Base_models)
  prediction = compute.prediction(model, testing_set, target, file.name, maximize = maximize_threshold)
  auc_roc_score = prediction[["AUC"]][["AUROC"]]
  auc_prc_score = prediction[["AUC"]][["AUPRC"]]
  metrics = prediction[["Metrics"]]
  predictions = prediction[["Predictions"]]

  get_curves(metrics, "Specificity", "Sensitivity", "Recall", "Precision", "model", auc_roc_score, auc_prc_score, file.name)

  return(prediction)
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
      source("~/Documents/CellTFusion_paper/user_projects/src/machine_learning.R")
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

compute.variable.importance = function(model, stacking = F, n_cores = 2){

  if(stacking == T){
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

###### Check folds stratification
# # Count class distribution in each fold
# check_stratification <- function(multifolds, target) {
#   class_ratios <- lapply(multifolds, function(fold_idx) {
#     table(target[fold_idx]) / length(fold_idx)  # Compute class proportions
#   })
#
#   return(class_ratios)
# }
#
# # Apply the function to check stratification
# class_ratios_per_fold <- check_stratification(multifolds, train_data$target)
#
# library(ggplot2)
#
# # Convert class_ratios_per_fold to a dataframe
# df <- do.call(rbind, lapply(1:length(class_ratios_per_fold), function(i) {
#   data.frame(Fold = names(class_ratios_per_fold)[i],
#              Class = names(class_ratios_per_fold[[i]]),
#              Proportion = as.numeric(class_ratios_per_fold[[i]]))
# }))
#
# # Plot class proportions
# ggplot(df, aes(x = Fold, y = Proportion, fill = Class)) +
#   geom_bar(stat = "identity", position = "dodge") +
#   theme_minimal() +
#   labs(title = "Class Distribution Across Folds", x = "Fold", y = "Proportion")
