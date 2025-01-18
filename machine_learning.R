
##Machine learning pipeline functions adapted for CellTFusion (https://github.com/VeraPancaldiLab/CellTFusion)

##Basic function to compute boruta algorithm
compute.boruta <- function(data, seed, fix = TRUE) {
  
  set.seed(seed)
  boruta_output <- Boruta(target ~ ., data = data, doTrace = 0)
  
  if (fix) {
    roughFixMod <- TentativeRoughFix(boruta_output)
    boruta_output <- roughFixMod
  }
  
  imps <- attStats(boruta_output)
  decision <- as.character(imps$decision)
  
  res <- imps %>%
    data.frame() %>%
    rownames_to_column("Variable") %>%
    dplyr::select(-decision)

  
  return(list(res, decision))
}

##Merge results from boruta iterations 
merge_boruta_results = function(importance_values, decisions, file_name, iterations, threshold, return = T){
  
  ### Construct matrix of importance
  combined_importance <- do.call(rbind, importance_values)
  combined_results_long <- combined_importance %>% #Matrix for plotting
    pivot_longer(cols = meanImp, names_to = "Measure", values_to = "Value")
  
  median_df <- combined_importance %>% #Calculate the median for each column, grouped by the variable name
    group_by(Variable) %>%
    dplyr::summarize(across(everything(), \(x) median(x, na.rm = TRUE)))
  
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
    arrange(meanImp) %>%
    pull(Variable)
  
  # For result 
  median_df$Decision = "Rejected"
  median_df$Decision[which(median_df$Variable %in% confirmed_vars)] = "Confirmed"
  median_df$Decision[which(median_df$Variable %in% tentative_vars)] = "Tentative"
  
  # Plot variable importance boxplots
  if(return){
    pdf(paste0("Results/Boruta_variable_importance_", file_name, ".pdf"), width = 8, height = 12)
    print(ggplot(combined_results_long, aes(x = factor(Variable, levels = mean_order), y = Value, fill = Decision)) +
            geom_bar(stat = "identity", position = "dodge") +
            theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
            coord_flip() +
            labs(x = "Features", y = "Importance", title = paste0("Variable Importance by Boruta after ", iterations, " bootstraps\n", file_name)) +
            scale_fill_manual(values = c("Confirmed" = "green", "Tentative" = "yellow", "Rejected" = "red")) +
            facet_wrap(~ Measure, scales = "free_y"))
    dev.off()
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
        num_cores <- detectCores() - 1  
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
      
      res <- foreach(seed = sample.int(100000, iterations)) %dopar% {
        
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
    folds <- createFolds(data$target, k = k_folds, returnTrain = TRUE, list = TRUE)
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
compute.k_fold_CV = function(model, k_folds, n_rep, stacking = F, metric = "Accuracy", boruta, boruta_iterations = NULL, fix_boruta = NULL, tentative = F, boruta_threshold = NULL, file_name = NULL, return){
  
  if(!(metric %in% c("AUC","Accuracy"))){
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
          mutate(target = model$target)
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
          mutate(target = model$target) 
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

  ######### Stratify K fold cross-validation 
  #folds <- createFolds(train_data[,'target'], k = k_folds, returnTrain = T, list = T) #this for single folds
  multifolds <- createMultiFolds(train_data[,'target'], k = k_folds, times = n_rep) #repeated folds
  trainControl <- trainControl(index = multifolds, method="repeatedcv", number=k_folds, repeats=n_rep, verboseIter = F, allowParallel = F, classProbs = TRUE, savePredictions=T)
  
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
  fit.treebag <- train(target~., data = train_data, method = "treebag", metric = "Accuracy",trControl = trainControl) 
  
  ################## RF
  require(randomForest)
  fit.rf <- train(target~., data = train_data, method = "rf", metric = "Accuracy",trControl = trainControl)
  
  ################## C5.0
  require(C50)
  fit.c50 <- train(target~., data = train_data, method = "C5.0", metric = "Accuracy",trControl = trainControl)
  
  ################## LG - Logistic Regression
  fit.glm <- train(target~., data = train_data, method="glm", metric="Accuracy",trControl=trainControl)
  
  ################## LDA - Linear Discriminate Analysis
  fit.lda <- train(target~., data = train_data, method="lda", metric="Accuracy",trControl=trainControl)
  
  ################## GLMNET - Regularized Logistic Regression (Elastic net)
  fit.glmnet <- train(target~., data = train_data, method="glmnet", metric="Accuracy",trControl=trainControl)
  
  ################## KNN - k-Nearest Neighbors 
  fit.knn <- train(target~., data = train_data, method="knn", metric="Accuracy",trControl=trainControl)
  
  ################## CART - Classification and Regression Trees (CART), 
  fit.cart <- train(target~., data = train_data, method="rpart", metric="Accuracy",trControl=trainControl)
  
  # NB - Naive Bayes (NB) 
  #Grid = expand.grid(usekernel=TRUE,adjust=1,fL=c(0,0.2,0.5,0.8,1))
  #fit.nb <- train(target~., data = train_data, method="nb", metric="Accuracy",trControl=trainControl, tuneGrid=Grid)
  
  ################## Regularized Lasso
  fit.lasso <- train(target~., data = train_data, method="glmnet", metric="Accuracy",trControl=trainControl, tuneGrid = expand.grid(alpha = 1, lambda = seq(0.001, 1, length = 20)))
  
  ################## Ridge regression
  fit.ridge <- train(target~., data = train_data, method="glmnet", metric="Accuracy",trControl=trainControl, tuneGrid = expand.grid(alpha = 0, lambda = seq(0.001, 1, length = 20)))
  
  ################## Support Vector Machine with Radial Kernel
  fit.svm_radial <- train(target ~ ., data = train_data, method = "svmRadial", metric = "Accuracy", trControl = trainControl)
  
  ################## Support Vector Machine with Linear Kernel
  fit.svm_linear <- train(target ~ ., data = train_data, method = "svmLinear", metric = "Accuracy", trControl = trainControl)
  
  ####### Optimized based on metric (only AUC or Accuracy available)
  if(metric == "AUC"){
    
    ################################################Bagged CART
    
    ## Integrate AUCs into prediction matrix
    fit.treebag$pred = fit.treebag$pred %>% 
      group_by(Resample) %>%
      mutate(AUC = calculate_auc_resample(obs, yes)) %>% #Calculate resamples AUC scores 
      ungroup() 
    
    ## Integrate AUCs into resamples matrix
    auc = c()
    for (i in 1:nrow(fit.treebag$resample)) {
      auc_val = fit.treebag$pred %>%
        filter(Resample == fit.treebag$resample$Resample[i]) %>%
        pull(AUC) %>% #AUC per resample is the same
        unique()
      
      auc = c(auc, auc_val)
    }
    
    fit.treebag$resample = fit.treebag$resample %>%
      mutate(AUC = auc) %>%
      select(AUC, everything())
    
    ## Integrate average CV AUCs into results
    fit.treebag$results = fit.treebag$results %>%
      mutate(AUC = mean(fit.treebag$resample$AUC))
    
    ################################################Random Forest
    
    ## Integrate AUCs into prediction matrix
    fit.rf$pred = fit.rf$pred %>%
      group_by(Resample, mtry) %>% #Parameters for tunning
      mutate(AUC = calculate_auc_resample(obs, yes)) %>%
      ungroup()
    
    ## Integrate AUCs into results per parameter 
    auc_values = fit.rf$pred %>%
      group_by(mtry) %>%
      summarise(AUC = mean(AUC, na.rm = TRUE), .groups = 'drop')
    
    fit.rf[["results"]] <- fit.rf[["results"]] %>%
      left_join(auc_values, by = "mtry")
    
    #Tuning parameter (select combination with top AUC)
    tune = which.max(fit.rf$results$AUC)
    fit.rf$bestTune = fit.rf$bestTune %>%
      mutate(mtry = fit.rf$results$mtry[tune])
    
    #Configure resamples to have the AUCs only using tuned parameter
    fit.rf$resample = fit.rf$resample[order(fit.rf$resample$Resample),] #Order resamples (just in case) to match with correct AUCs from prediction object 
    auc = c()
    for (i in 1:nrow(fit.rf$resample)){
      auc_val = fit.rf$pred %>%
        filter(mtry == as.numeric(fit.rf$bestTune),
               Resample == fit.rf$resample$Resample[i]) %>%
        pull(AUC) %>% #AUC per resample is the same
        unique()
      auc = c(auc, auc_val)
    }
    
    fit.rf$resample = fit.rf$resample %>%
      mutate(AUC = auc) 
    
    ################################################C5.0
    
    ## Integrate AUCs into prediction matrix
    fit.c50$pred = fit.c50$pred %>%
      group_by(trials, model, winnow) %>% #Parameters for tunning
      mutate(AUC = calculate_auc_resample(obs, yes)) %>%
      ungroup()
    
    ## Integrate AUCs into results per parameter 
    auc_values = fit.c50$pred %>%
      group_by(trials, model, winnow) %>%
      summarise(AUC = mean(AUC, na.rm = TRUE), .groups = 'drop')
    
    fit.c50[["results"]] <- fit.c50[["results"]] %>%
      left_join(auc_values, by = c("trials", "model", "winnow"))
    
    #Tuning parameter (select combination with top AUC)
    tune = which.max(fit.c50$results$AUC)
    fit.c50$bestTune = fit.c50$bestTune %>%
      mutate(trials = fit.c50$results$trials[tune],
             model = fit.c50$results$model[tune],
             winnow = fit.c50$results$winnow[tune])
    
    #Configure resamples to have the AUCs only using tuned parameter
    fit.c50$resample = fit.c50$resample[order(fit.c50$resample$Resample),] #Order resamples (just in case) to match with correct AUCs from prediction object 
    auc = c()
    for (i in 1:nrow(fit.c50$resample)){
      auc_val = fit.c50$pred %>%
        filter(trials == as.numeric(fit.c50$bestTune$trials),
               model == as.character(fit.c50$bestTune$model),
               winnow == as.character(fit.c50$bestTune$winnow),
               Resample == fit.c50$resample$Resample[i]) %>%
        pull(AUC) %>% #AUC per resample is the same
        unique()
      auc = c(auc, auc_val)
    }
    
    fit.c50$resample = fit.c50$resample %>%
      mutate(AUC = auc) 
    
    ################################################LG
    
    ## Integrate AUCs into prediction matrix
    fit.glm$pred = fit.glm$pred %>% #Calculate resamples AUC scores 
      group_by(Resample) %>%
      mutate(AUC = calculate_auc_resample(obs, yes)) %>% 
      ungroup() 
    
    ## Integrate AUCs into resamples matrix
    auc = c()
    for (i in 1:nrow(fit.glm$resample)) {
      auc_val = fit.glm$pred %>%
        filter(Resample == fit.glm$resample$Resample[i]) %>%
        pull(AUC) %>% #AUC per resample is the same
        unique()
      
      auc = c(auc, auc_val)
    }
    
    fit.glm$resample = fit.glm$resample %>%
      mutate(AUC = auc) 
    
    ## Integrate average CV AUCs into results
    fit.glm$results = fit.glm$results %>%
      mutate(AUC = mean(fit.glm$resample$AUC))
    
    ################################################LDA
    
    ## Integrate AUCs into prediction matrix
    fit.lda$pred = fit.lda$pred %>% #Calculate resamples AUC scores 
      group_by(Resample) %>%
      mutate(AUC = calculate_auc_resample(obs, yes)) %>% 
      ungroup() 
    
    ## Integrate AUCs into resamples matrix
    auc = c()
    for (i in 1:nrow(fit.lda$resample)) {
      auc_val = fit.lda$pred %>%
        filter(Resample == fit.lda$resample$Resample[i]) %>%
        pull(AUC) %>% #AUC per resample is the same
        unique()
      
      auc = c(auc, auc_val)
    }
    
    fit.lda$resample = fit.lda$resample %>%
      mutate(AUC = auc) 
    
    ## Integrate average CV AUCs into results
    fit.lda$results = fit.lda$results %>%
      mutate(AUC = mean(fit.lda$resample$AUC))
    
    ################################################GLMNET
    
    ## Integrate AUCs into prediction matrix
    fit.glmnet$pred = fit.glmnet$pred %>% #Calculate resamples AUC scores 
      group_by(Resample, alpha, lambda) %>%
      mutate(AUC = calculate_auc_resample(obs, yes)) %>% 
      ungroup() 
    
    ## Integrate AUCs into results per parameter 
    auc_values = fit.glmnet$pred %>%
      group_by(alpha, lambda) %>%
      summarise(AUC = mean(AUC, na.rm = TRUE), .groups = 'drop')
    
    fit.glmnet[["results"]] <- fit.glmnet[["results"]] %>%
      left_join(auc_values, by = c("alpha", "lambda"))
    
    #Tuning parameter (select combination with top AUC)
    tune = which.max(fit.glmnet$results$AUC)
    fit.glmnet$bestTune = fit.glmnet$bestTune %>%
      mutate(alpha = fit.glmnet$results$alpha[tune],
             lambda = fit.glmnet$results$lambda[tune])
    
    #Configure resamples to have the AUCs only using tuned parameter
    fit.glmnet$resample = fit.glmnet$resample[order(fit.glmnet$resample$Resample),] #Order resamples (just in case) to match with correct AUCs from prediction object 
    auc = c()
    for (i in 1:nrow(fit.glmnet$resample)){
      auc_val = fit.glmnet$pred %>%
        filter(alpha == as.numeric(fit.glmnet$bestTune$alpha),
               lambda == as.numeric(fit.glmnet$bestTune$lambda),
               Resample == fit.glmnet$resample$Resample[i]) %>%
        pull(AUC) %>% #AUC per resample is the same
        unique()
      auc = c(auc, auc_val)
    }
    
    fit.glmnet$resample = fit.glmnet$resample %>%
      mutate(AUC = auc) 
    
    ################################################KNN
    
    ## Integrate AUCs into prediction matrix
    fit.knn$pred = fit.knn$pred %>% #Calculate resamples AUC scores 
      group_by(Resample, k) %>%
      mutate(AUC = calculate_auc_resample(obs, yes)) %>% 
      ungroup() 
    
    ## Integrate AUCs into results per parameter 
    auc_values = fit.knn$pred %>%
      group_by(k) %>%
      summarise(AUC = mean(AUC, na.rm = TRUE), .groups = 'drop')
    
    fit.knn[["results"]] <- fit.knn[["results"]] %>%
      left_join(auc_values, by = "k")
    
    #Tuning parameter (select combination with top AUC)
    tune = which.max(fit.knn$results$AUC)
    fit.knn$bestTune = fit.knn$bestTune %>%
      mutate(k = fit.knn$results$k[tune])
    
    #Configure resamples to have the AUCs only using tuned parameter
    fit.knn$resample = fit.knn$resample[order(fit.knn$resample$Resample),] #Order resamples (just in case) to match with correct AUCs from prediction object 
    auc = c()
    for (i in 1:nrow(fit.knn$resample)){
      auc_val = fit.knn$pred %>%
        filter(k == as.numeric(fit.knn$bestTune$k),
               Resample == fit.knn$resample$Resample[i]) %>%
        pull(AUC) %>% #AUC per resample is the same
        unique()
      auc = c(auc, auc_val)
    }
    
    fit.knn$resample = fit.knn$resample %>%
      mutate(AUC = auc)
    
    ################################################CART
    
    ## Integrate AUCs into prediction matrix
    fit.cart$pred = fit.cart$pred %>% #Calculate resamples AUC scores 
      group_by(Resample, cp) %>%
      mutate(AUC = calculate_auc_resample(obs, yes)) %>% 
      ungroup() 
    
    ## Integrate AUCs into results per parameter 
    auc_values = fit.cart$pred %>%
      group_by(cp) %>%
      summarise(AUC = mean(AUC, na.rm = TRUE), .groups = 'drop')
    
    fit.cart[["results"]] <- fit.cart[["results"]] %>%
      left_join(auc_values, by = "cp")
    
    #Tuning parameter (select combination with top AUC)
    tune = which.max(fit.cart$results$AUC)
    fit.cart$bestTune = fit.cart$bestTune %>%
      mutate(cp = fit.cart$results$cp[tune])
    
    #Configure resamples to have the AUCs only using tuned parameter
    fit.cart$resample = fit.cart$resample[order(fit.cart$resample$Resample),] #Order resamples (just in case) to match with correct AUCs from prediction object 
    auc = c()
    for (i in 1:nrow(fit.cart$resample)){
      auc_val = fit.cart$pred %>%
        filter(cp == as.numeric(fit.cart$bestTune$cp),
               Resample == fit.cart$resample$Resample[i]) %>%
        pull(AUC) %>% #AUC per resample is the same
        unique()
      auc = c(auc, auc_val)
    }
    
    fit.cart$resample = fit.cart$resample %>%
      mutate(AUC = auc) 
    
    ################################################Regularized Lasso
    
    ## Integrate AUCs into prediction matrix
    fit.lasso$pred = fit.lasso$pred %>% #Calculate resamples AUC scores 
      group_by(Resample, alpha, lambda) %>%
      mutate(AUC = calculate_auc_resample(obs, yes)) %>% 
      ungroup() 
    
    ## Integrate AUCs into results per parameter 
    auc_values = fit.lasso$pred %>%
      group_by(alpha, lambda) %>%
      summarise(AUC = mean(AUC, na.rm = TRUE), .groups = 'drop')
    
    fit.lasso[["results"]] <- fit.lasso[["results"]] %>%
      left_join(auc_values, by = c("alpha", "lambda"))
    
    #Tuning parameter (select combination with top AUC)
    tune = which.max(fit.lasso$results$AUC)
    fit.lasso$bestTune = fit.lasso$bestTune %>%
      mutate(alpha = fit.lasso$results$alpha[tune],
             lambda = fit.lasso$results$lambda[tune])
    
    #Configure resamples to have the AUCs only using tuned parameter
    fit.lasso$resample = fit.lasso$resample[order(fit.lasso$resample$Resample),] #Order resamples (just in case) to match with correct AUCs from prediction object 
    auc = c()
    for (i in 1:nrow(fit.lasso$resample)){
      auc_val = fit.lasso$pred %>%
        filter(alpha == as.numeric(fit.lasso$bestTune$alpha),
               lambda == as.numeric(fit.lasso$bestTune$lambda),
               Resample == fit.lasso$resample$Resample[i]) %>%
        pull(AUC) %>% #AUC per resample is the same
        unique()
      auc = c(auc, auc_val)
    }
    
    fit.lasso$resample = fit.lasso$resample %>%
      mutate(AUC = auc) 
    
    ################################################Ridge regression
    
    ## Integrate AUCs into prediction matrix
    fit.ridge$pred = fit.ridge$pred %>% #Calculate resamples AUC scores 
      group_by(Resample, alpha, lambda) %>%
      mutate(AUC = calculate_auc_resample(obs, yes)) %>% 
      ungroup() 
    
    ## Integrate AUCs into results per parameter 
    auc_values = fit.ridge$pred %>%
      group_by(alpha, lambda) %>%
      summarise(AUC = mean(AUC, na.rm = TRUE), .groups = 'drop')
    
    fit.ridge[["results"]] <- fit.ridge[["results"]] %>%
      left_join(auc_values, by = c("alpha", "lambda"))
    
    #Tuning parameter (select combination with top AUC)
    tune = which.max(fit.ridge$results$AUC)
    fit.ridge$bestTune = fit.ridge$bestTune %>%
      mutate(alpha = fit.ridge$results$alpha[tune],
             lambda = fit.ridge$results$lambda[tune])
    
    #Configure resamples to have the AUCs only using tuned parameter
    fit.ridge$resample = fit.ridge$resample[order(fit.ridge$resample$Resample),] #Order resamples (just in case) to match with correct AUCs from prediction object 
    auc = c()
    for (i in 1:nrow(fit.ridge$resample)){
      auc_val = fit.ridge$pred %>%
        filter(alpha == as.numeric(fit.ridge$bestTune$alpha),
               lambda == as.numeric(fit.ridge$bestTune$lambda),
               Resample == fit.ridge$resample$Resample[i]) %>%
        pull(AUC) %>% #AUC per resample is the same
        unique()
      auc = c(auc, auc_val)
    }
    
    fit.ridge$resample = fit.ridge$resample %>%
      mutate(AUC = auc) 
    
    ################################################SVM radial
    
    ## Integrate AUCs into prediction matrix
    fit.svm_radial$pred = fit.svm_radial$pred %>% #Calculate resamples AUC scores 
      group_by(Resample, sigma, C) %>%
      mutate(AUC = calculate_auc_resample(obs, yes)) %>% 
      ungroup() 
    
    ## Integrate AUCs into results per parameter 
    auc_values = fit.svm_radial$pred %>%
      group_by(sigma, C) %>%
      summarise(AUC = mean(AUC, na.rm = TRUE), .groups = 'drop')
    
    fit.svm_radial[["results"]] <- fit.svm_radial[["results"]] %>%
      left_join(auc_values, by = c("sigma", "C"))
    
    #Tuning parameter (select combination with top AUC)
    tune = which.max(fit.svm_radial$results$AUC)
    fit.svm_radial$bestTune = fit.svm_radial$bestTune %>%
      mutate(sigma = fit.svm_radial$results$sigma[tune],
             C = fit.svm_radial$results$C[tune])
    
    #Configure resamples to have the AUCs only using tuned parameter
    fit.svm_radial$resample = fit.svm_radial$resample[order(fit.svm_radial$resample$Resample),] #Order resamples (just in case) to match with correct AUCs from prediction object 
    auc = c()
    for (i in 1:nrow(fit.svm_radial$resample)){
      auc_val = fit.svm_radial$pred %>%
        filter(sigma == as.numeric(fit.svm_radial$bestTune$sigma),
               C == as.numeric(fit.svm_radial$bestTune$C),
               Resample == fit.svm_radial$resample$Resample[i]) %>%
        pull(AUC) %>% #AUC per resample is the same
        unique()
      auc = c(auc, auc_val)
    }
    
    fit.svm_radial$resample = fit.svm_radial$resample %>%
      mutate(AUC = auc)
    
    ################################################SVM linear
    
    ## Integrate AUCs into prediction matrix
    fit.svm_linear$pred = fit.svm_linear$pred %>% #Calculate resamples AUC scores 
      group_by(Resample, C) %>%
      mutate(AUC = calculate_auc_resample(obs, yes)) %>% 
      ungroup() 
    
    ## Integrate AUCs into results per parameter 
    auc_values = fit.svm_linear$pred %>%
      group_by(C) %>%
      summarise(AUC = mean(AUC, na.rm = TRUE), .groups = 'drop')
    
    fit.svm_linear[["results"]] <- fit.svm_linear[["results"]] %>%
      left_join(auc_values, by = "C")
    
    #Tuning parameter (select combination with top AUC)
    tune = which.max(fit.svm_linear$results$AUC)
    fit.svm_linear$bestTune = fit.svm_linear$bestTune %>%
      mutate(C = fit.svm_linear$results$C[tune])
    
    #Configure resamples to have the AUCs only using tuned parameter
    fit.svm_linear$resample = fit.svm_linear$resample[order(fit.svm_linear$resample$Resample),] #Order resamples (just in case) to match with correct AUCs from prediction object 
    auc = c()
    for (i in 1:nrow(fit.svm_linear$resample)){
      auc_val = fit.svm_linear$pred %>%
        filter(C == as.numeric(fit.svm_linear$bestTune$C),
               Resample == fit.svm_linear$resample$Resample[i]) %>%
        pull(AUC) %>% #AUC per resample is the same
        unique()
      auc = c(auc, auc_val)
    }
    
    fit.svm_linear$resample = fit.svm_linear$resample %>%
      mutate(AUC = auc) 
    
  }
  
  ###Prediction with best tuned hyper-parameters
  
  ###Bagged CART
  
  predictions.bag <- data.frame(predict(fit.treebag, newdata = train_data, type = "prob")) %>% #Predictions using tuned model
    dplyr::select(yes) %>%
    dplyr::rename(BAG = yes) 
  
  ###Random Forest
  
  predictions.rf = data.frame(predict(fit.rf, newdata = train_data, type = "prob"))[,"yes", drop=F]  %>%
    dplyr::select(yes) %>%
    dplyr::rename(RF = yes) #Predictions of model (already ordered)
  
  ###C5.0
  
  predictions.c50 = data.frame(predict(fit.c50$finalModel, newdata = train_data, type = "prob"))[,"yes", drop=F]  %>% 
    dplyr::select(yes) %>%
    dplyr::rename(C50 = yes)  #Predictions of model (already ordered)
  
  ### LG
  
  predictions.glm = predict(fit.glm, newdata = train_data, type = "prob")[,"yes", drop=F]  %>%
    dplyr::select(yes) %>%
    dplyr::rename(GLM = yes)  #Predictions of model (already ordered)
  
  ### LDA
  
  predictions.lda = predict(fit.lda, newdata = train_data, type = "prob")[,"yes", drop=F]  %>%
    dplyr::select(yes) %>%
    dplyr::rename(LDA = yes)  #Predictions of model (already ordered)
  
  ### GLMNET
  
  predictions.glmnet = predict(fit.glmnet, newdata = train_data, type = "prob")[,"yes", drop=F]  %>%
    dplyr::select(yes) %>%
    dplyr::rename(GLMNET = yes)  #Predictions of model (already ordered)
  
  ### KNN
  
  predictions.knn = predict(fit.knn, newdata = train_data, type = "prob")[,"yes", drop=F]  %>%
    dplyr::select(yes) %>%
    dplyr::rename(KNN = yes) #Predictions of model (already ordered)
  
  ## CART
  
  predictions.cart = predict(fit.cart, newdata = train_data, type = "prob")[,"yes", drop=F]  %>%
    dplyr::select(yes) %>%
    dplyr::rename(CART = yes)  #Predictions of model (already ordered)
  
  ## Regularized Lasso
  
  predictions.lasso = predict(fit.lasso, newdata = train_data, type = "prob")[,"yes", drop=F]  %>%
    dplyr::select(yes) %>%
    dplyr::rename(LASSO = yes)  #Predictions of model (already ordered)
  
  ## Ridge regression
  
  predictions.ridge = predict(fit.ridge, newdata = train_data, type = "prob")[,"yes", drop=F]  %>%
    dplyr::select(yes) %>%
    dplyr::rename(RIDGE = yes)  #Predictions of model (already ordered)
  
  ## SVM radial
  
  predictions.svm_radial = predict(fit.svm_radial, newdata = train_data, type = "prob")[,"yes", drop=F]  %>%
    dplyr::select(yes) %>%
    dplyr::rename(SVM_radial = yes)  #Predictions of model (already ordered)
  
  ## SVM linear
  
  predictions.svm_linear = predict(fit.svm_linear, newdata = train_data, type = "prob")[,"yes", drop=F]  %>%
    dplyr::select(yes) %>%
    dplyr::rename(SVM_linear = yes)  #Predictions of model (already ordered)
  
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
                          SVM_linear = fit.svm_linear)
  
  
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
                           SVM_linear = predictions.svm_linear)
  
  #Remove models with same predictions across samples (not able to make distinction)
  model_predictions <- lapply(model_predictions, function(df) {
    df = df %>%
      select(where(~ n_distinct(.) > 1))
    
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
    }else if(metric == "AUC"){
      base_models = compute_cv_AUC(ensembleResults, base_models = T, file_name = file_name, return = return)
    }
    
    #Save variable importance of each base model
    importance = list()
    for (i in 1:length(base_models$Base_models)) {
      importance[[i]] = varImp(ensembleResults[[base_models$Base_models[i]]], scale = F)
      names(importance)[i] = base_models$Base_models[i]
    }
    
    features_predictions = model_predictions %>%
      t() %>%
      data.frame() %>%
      rownames_to_column("Models") %>%
      filter(grepl(paste0("\\b(", paste(base_models$Base_models, collapse = "|"), ")\\b"), Models)) %>%
      column_to_rownames("Models") %>%
      t() %>%
      data.frame()
    
    meta_features = cbind(features_predictions, "true_label" = train_data$target) 
    
    meta_learner <- train(true_label ~ ., data = meta_features, method = "glmnet", trControl = trainControl) #Staking based on simple logistic regression
    
    #Base models using ALL ML models 
    meta_features_all = cbind(model_predictions, "true_label" = train_data$target) 
    
    meta_learner_all <- train(true_label ~ ., data = meta_features_all, method = "glmnet", trControl = trainControl) #Staking based on simple logistic regression
    
    cat("Meta-learners ML model based on GLM\n")
    output = list("Features" = features, "Meta_learner" = meta_learner, "Base_models" = base_models$Base_models, "ML_models" = ensembleResults, "Variable_importance" = importance)
    
  }else{
    features = colnames(train_data)[colnames(train_data) != "target"] #Extract features used for model training
    
    #Top model with best accuracy or AUC
    if(metric == "Accuracy"){
      metrics = compute_cv_accuracy(ensembleResults, file_name, file_name = file_name, return = return)
    }else if(metric == "AUC"){
      metrics = compute_cv_AUC(ensembleResults, file_name, file_name = file_name, return = return)
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
compute.ML = function(raw.counts, normalized = F, deconv, tf.universe, paths.universe, clinical, trait, trait.positive, partition, metric = "Accuracy", stack, feature.selection = F, deconv_methods = c("Quantiseq", "CBSX", "Epidish", "DeconRNASeq", "DWLS"), doParallel = F, workers = NULL, seed, file_name = NULL, return = F){
  
  set.seed(seed)   
  
  # Do stratified partition 
  index = createDataPartition(clinical[,trait], times = 1, p = partition, list = FALSE) 
  
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
  #deconv_train = compute.deconvolution(raw.counts_train, normalized = normalized, methods = deconv_methods, doParallel = doParallel, workers = workers, credentials.mail = "marcelo.hurtado@inserm.fr", credentials.token = "734212f6ad77fc4eea2bdb502792f294", return = return)
  
  # Test cohort
  traitData_test = clinical[-index, ]
  raw.counts_test = raw.counts[,-index]
  counts.normalized_test = norm.counts[,-index]
  deconv_test = deconv[-index,]
  tfs_test = compute.TFs.activity(counts.normalized_test, universe = tf.universe)
  #deconv_test = compute.deconvolution(raw.counts_test, normalized = normalized, methods = deconv_methods, doParallel = doParallel, workers = workers, credentials.mail = "marcelo.hurtado@inserm.fr", credentials.token = "734212f6ad77fc4eea2bdb502792f294", return = return)
  
  ###############################################################################################################################################################################
  
  # CellTFusion
  set.seed(seed)
  network = compute.WTCNA(tfs_train, corr_mod = 0.7, clustering.method = "ward.D2", return = return) 
  pathways = compute.pathway.activity(counts.normalized_train, paths = paths.universe)
  dt = compute.deconvolution.analysis(deconv_train, corr = 0.7, seed = 123, return = return, file_name = file_name) #Here we try to put the same seed just in case we have the same high_corr features, to ensure the same are choose for easier comparison between ML models (this choosing is random anyway so it doesn't affect)
  corr_modules = compute.modules.relationship(network[[1]], dt[[1]], return = T, plot = return)
  cell_dendrograms = identify.cell.groups(corr_modules, height = 20, return = return)
  cell.groups = cell.groups.computation(dt[[1]], cell.dendrograms = cell_dendrograms, corr_modules[[1]], network, return = return) #Identify cell groups with specific cut 
  
  ###############################################################################################################################################################################
  
  ####################################################Training

  #Set training set
  train_data = cell.groups[[1]] %>%
    data.frame() %>%
    mutate(Trait = traitData_train[,trait],
           target = as.factor(ifelse(Trait == trait.positive, 'yes', 'no'))) %>%
    dplyr::select(-Trait)
  
  train_data$target <- factor(train_data$target, levels = c("no", "yes"))  # Order, just in case to ensure positive class is not well defined

  #Cross-validation training (5 k-folds and 100 repetitions)
  training = compute.k_fold_CV(train_data, k_folds = 5, n_rep = 100, metric = metric, stacking = stack, boruta = feature.selection, boruta_iterations = 100, fix_boruta = F, boruta_threshold = 0.8, file_name = file_name, return= return)
  
  ####################################################Predicting
  if(length(training)!=0){
    cell_groups = cell.groups #Save cell groups scores per partition
    network = network #Save TF network per partition
    features = training[["Features"]] #Save selected features per partition (only useful if Boruta = T - needs to be improve it)
    deconvolution_subgroups = dt[["Deconvolution groups - Linear-based correlation"]] #Save cell subgroups  
    ####################### Testing set
    testing_set = compute_cell_groups_signatures(dt, network, cell.groups, features, deconv_test, tfs_test) #Cell groups projection
    #Extract target variable
    target = traitData_test %>%
      mutate(target = ifelse(traitData_test[,trait] == trait.positive, "yes", "no")) %>%
      pull(target)
    
    target = factor(target, levels = c("no", "yes"))
    
    if(stack){
      model = training[["Meta_learner"]]
      var_importance = calculate_feature_importance_stacking(training[["Variable_importance"]], training[["Base_models"]], model)
      prediction = compute.prediction.stacked(model, testing_set, target, training[["ML_models"]], training[["Base_models"]])
      
    }else{
      model = training[["Model"]] #Save best ML model based on the Accuracy/AUC from CV per partition
      var_importance = varImp(model, scale = F) #Retrieve variable importance
      prediction = compute.prediction(model, testing_set, target)
    }
    
    auc_roc_score = prediction[["AUC"]][["AUROC"]]
    auc_prc_score = prediction[["AUC"]][["AUPRC"]]
    
    metrics = prediction[["Metrics"]]
    predictions = prediction[["Predictions"]]
    
    if(return == T){
      get_curves(metrics, "specificity", "sensitivity", "recall", "precision", "model", auc_roc_score, auc_prc_score, file_name)
    }
    
    rm(pathways, tfs.modules.clusters, dt, corr_modules, cell_dendrograms, cell.groups,
       traitData_train, counts.normalized_train, tfs_train, deconv_train,
       traitData_test, counts.normalized_test, deconv_test, tfs_test, 
       clinical, norm.counts) #Remove variables
    
    gc() #Clean garbage
    
    return(list(Model = model, Features = features, Variable_importance = var_importance, Cell_groups = cell_groups, Deconvolution_subgroups = deconvolution_subgroups, TF_network = network, Prediction_metrics = metrics, AUC = list(AUROC = auc_roc_score, AUPRC = auc_prc_score), Prediction = predictions))
  }else{  #No features are selected as predictive
    
    rm(pathways, tfs.modules.clusters, dt, corr_modules, cell_dendrograms, cell.groups,
       traitData_train, counts.normalized_train, tfs_train, deconv_train,
       traitData_test, counts.normalized_test, deconv_test, tfs_test, 
       clinical, norm.counts) #Remove variables
    
    gc() #Clean garbage
    
    message("No features selected as predictive after Boruta runs. No model returned.")
    
    return(NULL)
  }
  
}

compute.deconvolution.ML = function(deconv, clinical, trait, trait.positive, partition, metric = "Accuracy", stack, feature.selection = F, seed, file_name = NULL, return = F){
  
  set.seed(seed)   
  
  # Do stratified partition 
  index = createDataPartition(clinical[,trait], times = 1, p = partition, list = FALSE) 
  
  # Train cohort
  traitData_train = clinical[index, ]
  deconv_train = deconv[index,]
  dt = compute.deconvolution.analysis(deconv_train, corr = 0.7, seed = 123, return = return, file_name = file_name) 
  
  # Test cohort
  traitData_test = clinical[-index, ]
  deconv_test = deconv[-index,]
  
  
  ###############################################################################################################################################################################
  
  ####################################################Training
  
  #Set training set
  train_data = dt[[1]] %>%
    data.frame() %>%
    mutate(Trait = traitData_train[,trait],
           target = as.factor(ifelse(Trait == trait.positive, 'yes', 'no'))) %>%
    dplyr::select(-Trait)
  
  train_data$target <- factor(train_data$target, levels = c("no", "yes"))  # Order, just in case to ensure positive class is not well defined
  
  #Cross-validation training (5 k-folds and 100 repetitions)
  training = compute.k_fold_CV(train_data, k_folds = 5, n_rep = 100, metric = metric, stacking = stack, boruta = feature.selection, boruta_iterations = 100, fix_boruta = F, boruta_threshold = 0.8, file_name = file_name, return= return)
  
  ####################################################Predicting
  if(length(training)!=0){
    features = training[["Features"]] #Save selected features per partition (only useful if Boruta = T - needs to be improve it)
    deconvolution_subgroups = dt[["Deconvolution groups - Linear-based correlation"]] #Save cell subgroups  
    
    ####################### Testing set
    testing_set = replicate_deconvolution_subgroups(dt, features, deconv_test) #Replicate subgroups
    
    #Extract target variable
    target = traitData_test %>%
      mutate(target = ifelse(traitData_test[,trait] == trait.positive, "yes", "no")) %>%
      pull(target)
    
    target = factor(target, levels = c("no", "yes"))
    
    if(stack){
      model = training[["Meta_learner"]]
      var_importance = calculate_feature_importance_stacking(training[["Variable_importance"]], training[["Base_models"]], model)
      prediction = compute.prediction.stacked(model, testing_set, target, training[["ML_models"]], training[["Base_models"]])
      
    }else{
      model = training[["Model"]] #Save best ML model based on the Accuracy/AUC from CV per partition
      var_importance = varImp(model, scale = F) #Retrieve variable importance
      prediction = compute.prediction(model, testing_set, target)
    }
    
    auc_roc_score = prediction[["AUC"]][["AUROC"]]
    auc_prc_score = prediction[["AUC"]][["AUPRC"]]
    
    metrics = prediction[["Metrics"]]
    predictions = prediction[["Predictions"]]
    
    if(return == T){
      get_curves(metrics, "specificity", "sensitivity", "recall", "precision", "model", auc_roc_score, auc_prc_score, file_name)
    }
    
    
    return(list(Model = model, Features = features, Variable_importance = var_importance, Deconvolution_subgroups = deconvolution_subgroups, Prediction_metrics = metrics, AUC = list(AUROC = auc_roc_score, AUPRC = auc_prc_score), Prediction = predictions))
  }else{  #No features are selected as predictive
    
    message("No features selected as predictive after Boruta runs. No model returned.")
    
    return(NULL)
  }
  
}

compute.raw.deconvolution.ML = function(deconv, clinical, trait, trait.positive, partition, metric = "Accuracy", stack, feature.selection = F, seed, file_name = NULL, return = F){
  
  set.seed(seed)   
  
  # Do stratified partition 
  index = createDataPartition(clinical[,trait], times = 1, p = partition, list = FALSE) 
  
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
    mutate(Trait = traitData_train[,trait],
           target = as.factor(ifelse(Trait == trait.positive, 'yes', 'no'))) %>%
    dplyr::select(-Trait)
  
  train_data$target <- factor(train_data$target, levels = c("no", "yes"))  # Order, just in case to ensure positive class is not well defined
  
  #Cross-validation training (5 k-folds and 100 repetitions)
  training = compute.k_fold_CV(train_data, k_folds = 5, n_rep = 100, metric = metric, stacking = stack, boruta = feature.selection, boruta_iterations = 100, fix_boruta = F, boruta_threshold = 0.8, file_name = file_name, return= return)
  
  ####################################################Predicting
  if(length(training)!=0){
    features = training[["Features"]] #Save selected features per partition (only useful if Boruta = T - needs to be improve it)

    ####################### Testing set
    testing_set = deconv_test[,features] #Replicate subgroups
    
    #Extract target variable
    target = traitData_test %>%
      mutate(target = ifelse(traitData_test[,trait] == trait.positive, "yes", "no")) %>%
      pull(target)
    
    target = factor(target, levels = c("no", "yes"))
    
    if(stack){
      model = training[["Meta_learner"]]
      var_importance = calculate_feature_importance_stacking(training[["Variable_importance"]], training[["Base_models"]], model)
      prediction = compute.prediction.stacked(model, testing_set, target, training[["ML_models"]], training[["Base_models"]])
      
    }else{
      model = training[["Model"]] #Save best ML model based on the Accuracy/AUC from CV per partition
      var_importance = varImp(model, scale = F) #Retrieve variable importance
      prediction = compute.prediction(model, testing_set, target)
    }
    
    auc_roc_score = prediction[["AUC"]][["AUROC"]]
    auc_prc_score = prediction[["AUC"]][["AUPRC"]]
    
    metrics = prediction[["Metrics"]]
    predictions = prediction[["Predictions"]]
    
    if(return == T){
      get_curves(metrics, "specificity", "sensitivity", "recall", "precision", "model", auc_roc_score, auc_prc_score, file_name)
    }
    
    
    return(list(Model = model, Features = features, Variable_importance = var_importance, Prediction_metrics = metrics, AUC = list(AUROC = auc_roc_score, AUPRC = auc_prc_score), Prediction = predictions))
  }else{  #No features are selected as predictive
    
    message("No features selected as predictive after Boruta runs. No model returned.")
    
    return(NULL)
  }
  
}

#Main ML Pipeline for doing Leaving-one-dataset-out (LODO)

compute.LODO.ML = function(raw.counts, normalized = F, deconv, tf.universe, paths.universe, clinical, trait, trait.positive, trait.out, out, metric = "Accuracy", stack, feature.selection = T, doParallel = F, workers = NULL, deconv_methods = c("Quantiseq", "CBSX", "Epidish", "DeconRNASeq", "DWLS"), file_name = NULL, return = T){
  
  clinical = clinical %>%
    mutate(trait.out = clinical[,trait.out]) ##Just a way to make work "filter" cause it does not allow "" variables (might change after)
  
  # Normalize counts
  if(normalized == T){
    norm.counts = data.frame(ADImpute::NormalizeTPM(raw.counts)) 
  }else{ #If they are already normalized
    norm.counts = raw.counts
  }
  
  # Test cohort
  traitData_test = clinical %>%
    filter(trait.out == out)
  raw.counts_test = raw.counts[,colnames(raw.counts)%in%rownames(traitData_test)]
  counts.normalized_test = norm.counts[,colnames(norm.counts)%in%rownames(traitData_test)]
  
  tfs_test = compute.TFs.activity(counts.normalized_test, universe = tf.universe)
  #deconv_test = compute.deconvolution(raw.counts_test, normalized = normalized, methods = deconv_methods, doParallel = doParallel, workers = workers, credentials.mail = "marcelo.hurtado@inserm.fr", credentials.token = "734212f6ad77fc4eea2bdb502792f294")
  deconv_test = deconv[rownames(deconv)%in%rownames(traitData_test),]
  
  # Train cohort
  traitData_train = clinical %>%
    filter(trait.out != out)
  raw.counts_train = raw.counts[,colnames(raw.counts)%in%rownames(traitData_train)]
  counts.normalized_train = norm.counts[,colnames(norm.counts)%in%rownames(traitData_train)]
  
  tfs_train = compute.TFs.activity(counts.normalized_train, universe = tf.universe)
  #deconv_train = compute.deconvolution(raw.counts_train, normalized = normalized, methods = deconv_methods, doParallel = doParallel, workers = workers, credentials.mail = "marcelo.hurtado@inserm.fr", credentials.token = "734212f6ad77fc4eea2bdb502792f294")
  deconv_train = deconv[rownames(deconv)%in%rownames(traitData_train),]
  
  ###############################################################################################################################################################################
  
  # CellTFusion
  network = compute.WTCNA(tfs_train, corr_mod = 0.8, clustering.method = "ward.D2", return = return) 
  pathways = compute.pathway.activity(counts.normalized_train, paths = paths.universe)
  tfs.modules.clusters = compute.TF.network.classification(network, pathways, return = return)
  dt = compute.deconvolution.analysis(deconv_train, corr = 0.7, seed = 123)
  corr_modules = compute.modules.relationship(network[[1]], dt[[1]], return = T, plot = return)
  cell_dendrograms = identify.cell.groups(corr_modules, tfs.modules.clusters, height = 20, return = return)
  cell.groups = cell.groups.computation(dt[[1]], tfs.module.network = network, cell.dendrograms = cell_dendrograms) #Identify cell groups with specific cut 
  
  ###############################################################################################################################################################################
  
  ####################################################Training
  
  #Set training set
  train_data = cell.groups[[1]] %>%
    data.frame() %>%
    mutate(Trait = traitData_train[,trait],
           target = as.factor(ifelse(Trait == trait.positive, 'yes', 'no'))) %>%
    dplyr::select(-Trait)
  
  train_data$target <- factor(train_data$target, levels = c("no", "yes"))  # Just in case positive class is not well defined
  
  #Cross-validation training (5 k-folds and 100 repetitions)
  training = compute.k_fold_CV(train_data, k_folds = 5, n_rep = 100,  metric = metric, stacking = stack, boruta = feature.selection, boruta_iterations = 100, fix_boruta = F, boruta_threshold = 0.8, file_name = file_name, return= return)
  
  ####################################################Predicting
  if(length(training)!=0){
    cell_groups = cell.groups #Save cell groups scores per partition
    features = training[["Features"]] #Save selected features per partition
    deconvolution_subgroups = dt[["Deconvolution groups - Linear-based correlation"]] #Save cell subgroups  
    ####################### Testing set
    testing_set = compute_cell_groups_signatures(dt, network, cell.groups, features, deconv_test, tfs_test) #Compute features in testing set
    #Extract target variable
    target = traitData_test %>%
      mutate(target = ifelse(traitData_test[,trait] == trait.positive, "yes", "no")) %>%
      pull(target) 
    
    target = factor(target, levels = c("no", "yes"))
    
    if(stack){
      model = training[["Meta_learner"]]
      prediction = compute.prediction.stacked(model, testing_set, target, training[["ML_models"]], training[["Base_models"]])
    }else{
      model = training[["Model"]] #Save best ML model based on the Accuracy from CV per partition
      prediction = compute.prediction(model, testing_set, target)
    }
    
    auc_roc_score = prediction[["AUC"]][["AUROC"]]
    auc_prc_score = prediction[["AUC"]][["AUPRC"]]
    
    metrics = prediction[["Metrics"]]
    predictions = prediction[["Predictions"]]
    
    if(return == T){
      get_curves(metrics, "specificity", "sensitivity", "recall", "precision", "model", auc_roc_score, auc_prc_score, file_name)
    }
    
    rm(network, pathways, tfs.modules.clusters, dt, corr_modules, cell_dendrograms, cell.groups,
       traitData_train, counts.normalized_train, tfs_train, deconv_train,
       traitData_test, counts.normalized_test, deconv_test, tfs_test, 
       clinical, norm.counts) #Remove variables
    
    gc() #Clean garbage
    
    return(list(Model = model, Features = features, Cell_groups = cell_groups, Deconvolution_subgroups = deconvolution_subgroups, Prediction_metrics = metrics, AUC = list(AUROC = auc_roc_score, AUPRC = auc_prc_score), Prediction = predictions))
  }else{  #No features are selected as predictive
    
    rm(network, pathways, tfs.modules.clusters, dt, corr_modules, cell_dendrograms, cell.groups,
       traitData_train, counts.normalized_train, tfs_train, deconv_train,
       traitData_test, counts.normalized_test, deconv_test, tfs_test, 
       clinical, norm.counts) #Remove variables
    
    gc() #Clean garbage
    
    return(NULL)
  }
  
}

unregister_dopar <- function() {
  env <- foreach:::.foreachGlobals
  rm(list=ls(name=env), pos=env)
  gc()
}

#Main ML pipeline for several partitions
compute.bootstrap.ML = function(raw.counts, normalized = F, deconv, clinical, trait, trait.positive, partition = 0.8, metric = "Accuracy", iterations, feature.selection = F, stack, deconv_methods = c("Quantiseq", "CBSX", "Epidish", "DeconRNASeq", "DWLS"), workers = NULL, file.name = NULL, return = F){
  
  folder = paste0("Results/ML_models_", file.name)
  dir.create(file.path(getwd(), folder))
  
  if(is.null(iterations) == T){
    stop("No iterations specified, please set a number")
  }else{
    if(is.null(workers)==T){
      num_cores <- detectCores() - 1
    }else{
      num_cores <- workers
    }
    
    cl = parallel::makeCluster(num_cores) #Forking just copy the R session in its current state. - makeCluster() all must be exported (copied) to the clusters, which can add some overhead
    doParallel::registerDoParallel(cl)

    message("\nRunning ", iterations, " splits for training and test using ", num_cores, " cores")

    tfs.universe = decoupleR::get_collectri(organism = 'human', split_complexes = F)
    pathways.universe <- get_progeny(organism = 'human', top = 500)
    
    # Run foreach loop using each random seed directly
    foreach(iteration = seq_len(iterations), random.seed = sample.int(100000, iterations)) %dopar% {
      
      # Use absolute path for the source file to avoid path issues
      source("src/environment_set.R")
      
      # Run foreach loop using each random seed directly
      tryCatch({
        # Compute the result with the current random seed
        result <- compute.ML(
          raw.counts, normalized, deconv, tfs.universe, pathways.universe, clinical, trait, trait.positive, partition, 
          metric, stack, feature.selection, doParallel = F, workers = NULL, 
          seed = random.seed, deconv_methods = deconv_methods, 
          file_name = file.name, return = return
        )
        
        # Save result as RDS file with unique identifier based on iteration (random seed)
        saveRDS(list(result = result, seed = random.seed), 
                file = file.path(paste0(folder, "/ML_result_", iteration, ".rds")))
        
      }, error = function(e) {

        # Save error information as RDS file with random seed identifier
        saveRDS(list(result = NULL, error = e$message, seed = random.seed), 
                file = file.path(paste0(folder, "/ML_result_", iteration, ".rds")))
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
  
  folder = paste0("Results/ML_models_",file.name)
  
  dir.create(file.path(getwd(), folder))
  
  if(is.null(iterations) == T){
    stop("No iterations specified, please set a number")
  }else{
    if(is.null(workers)==T){
      num_cores <- detectCores() - 1
    }else{
      num_cores <- workers
    }
    
    cl = parallel::makeCluster(num_cores) #Forking just copy the R session in its current state. - makeCluster() all must be exported (copied) to the clusters, which can add some overhead
    doParallel::registerDoParallel(cl)
    
    message("\nRunning ", iterations, " splits for training and test using ", num_cores, " cores")
    
    # Run foreach loop using each random seed directly
    foreach(iteration = seq_len(iterations), random.seed = sample.int(100000, iterations)) %dopar% {
      
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
                file = file.path(paste0(folder, "/ML_result_", iteration, ".rds")))
        
      }, error = function(e) {
        
        # Save error information as RDS file with random seed identifier
        saveRDS(list(result = NULL, error = e$message, seed = random.seed), 
                file = file.path(paste0(folder, "/ML_result_", iteration, ".rds")))
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
      num_cores <- detectCores() - 1
    }else{
      num_cores <- workers
    }
    
    cl = parallel::makeCluster(num_cores) #Forking just copy the R session in its current state. - makeCluster() all must be exported (copied) to the clusters, which can add some overhead
    doParallel::registerDoParallel(cl)
    
    message("\nRunning ", iterations, " splits for training and test using ", num_cores, " cores")
    
    # Run foreach loop using each random seed directly
    foreach(iteration = seq_len(iterations), random.seed = sample.int(100000, iterations)) %dopar% {
      
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
                file = file.path(paste0(folder, "/ML_result_", iteration, ".rds")))
        
      }, error = function(e) {
        
        # Save error information as RDS file with random seed identifier
        saveRDS(list(result = NULL, error = e$message, seed = random.seed), 
                file = file.path(paste0(folder, "/ML_result_", iteration, ".rds")))
      })
      
    }
    
    #Stop cluster after all runs
    parallel::stopCluster(cl)
    unregister_dopar() #Stop Dopar from running in the background
    
  }
  
  message("Analysis is done!")
  
  message("ML models are saved in Results/ML_models folder")
  
}


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
    group_by(Cohort) %>%
    dplyr::summarize(medianAUROC = median(AUC_roc))
  
  median_auc_prc = cumulative_data %>%
    group_by(Cohort) %>%
    dplyr::summarize(medianAUPRC = median(AUC_prc))
  
  # Plot boxplot with mean AUC annotations
  plot_roc = ggplot(cumulative_data, aes(x = Cohort, y = AUC_roc, fill = Cohort)) +
    geom_boxplot() +
    labs(title = paste0("Distribution of AUROC values across ", iterations, " splits"),
         x = "Model",
         y = "AUROC") +
    theme_minimal() +
    theme(legend.position = "right") +
    geom_text(data = median_auc_roc, aes(x = Cohort, y = max(cumulative_data$AUC_roc), 
                                   label = paste("Median AUROC:", round(medianAUROC, 3))),
              size = 4, color = "black", vjust = -0.5)
  
  pdf(paste0("Results/Boxplot_AUROC_performance_", file.name, ".pdf"))
  print(plot_roc)
  dev.off()
  
  plot_prc = ggplot(cumulative_data, aes(x = Cohort, y = AUC_prc, fill = Cohort)) +
    geom_boxplot() +
    labs(title = paste0("Distribution of AUPRC values across ", iterations, " splits"),
         x = "Model",
         y = "AUPRC") +
    theme_minimal() +
    theme(legend.position = "right") +
    geom_text(data = median_auc_prc, aes(x = Cohort, y = max(cumulative_data$AUC_prc), 
                                       label = paste("Median AUPRC:", round(medianAUPRC, 3))),
              size = 4, color = "black", vjust = -0.5)
  
  pdf(paste0("Results/Boxplot_AUPRC_performance_", file.name, ".pdf"))
  print(plot_prc)
  dev.off()
  
}

compute_cv_accuracy = function(models, file_name = NULL, base_models = F, return = T){
  
  #Bind accuracy values from each model
  accuracy = list()
  for (i in 1:length(models)){
    accuracy[[i]] = models[[i]]$resample %>% 
      mutate(model = names(models)[i])
    names(accuracy)[i] = names(models)[i]
  }
  accuracy_data = do.call(rbind, accuracy)
  
  #Retrieve top model based on accuracy
  res_accuracy <- accuracy_data %>%
    group_by(model) %>%
    summarise(Accuracy = mean(Accuracy)) %>%
    arrange(desc(Accuracy)) 
  
  top_model = res_accuracy %>%
    slice(1) %>%
    pull(model)
  
  if(return){
    pdf(paste0("Results/Accuracy_CV_methods_", file_name, ".pdf"), width = 10)
    plot(ggplot(accuracy_data, aes(x = model, y = Accuracy, fill = model)) +
           geom_boxplot() +
           labs(title = "Distribution of Accuracy Values by Model",
                x = "Model",
                y = "Accuracy") +
           theme_minimal() +
           theme(legend.position = "none"))
    dev.off()
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

compute_cv_AUC = function(models, file_name = NULL, base_models = F, return = T){
  
  #Bind AUC values from each model
  auc = list()
  for (i in 1:length(models)){
    auc[[i]] = models[[i]]$resample %>% #we use the resample matrix and not directly the results matrix as some have hyperparameters so we will need to define best on the tuned parameter (=more code) - resample matrix is made based on the best tuning
      mutate(model = names(models)[i])
    names(auc)[i] = names(models)[i]
  }
  auc_data = do.call(rbind, auc)
  
  #Retrieve top model based on accuracy
  res_auc <- auc_data %>%
    group_by(model) %>%
    summarise(AUC = mean(AUC))  %>%
    arrange(desc(AUC)) 
  
  top_model = res_auc %>%
    slice(1) %>%
    pull(model)
  
  if(return){
    pdf(paste0("Results/AUC_CV_methods_", file_name, ".pdf"), width = 10)
    plot(ggplot(auc_data, aes(x = model, y = AUC, fill = model)) +
           geom_boxplot() +
           labs(title = "Distribution of AUC scores by Model",
                x = "Model",
                y = "AUC") +
           theme_minimal() +
           theme(legend.position = "none"))
    dev.off()
  }
  
  if(base_models == T){
    cat("Choosing base models for stacking.......................................\n\n")
    base_models = choose_base_models(models, metric = "AUC")
    cat("Models chosen are:", paste0(base_models, collapse = ", "), "\n\n")
    return(list("AUC" = res_auc, "Top_model" = top_model, "Base_models" = base_models))
  }else{
    return(list("AUC" = res_auc, "Top_model" = top_model))
  }
  
}

choose_base_models = function(models, metric = "Accuracy"){
  
  #Bind accuracy values from each model
  resample_df = list()
  for (i in 1:length(models)){
    resample_df[[i]] = models[[i]]$resample %>% 
      mutate(model = names(models)[i])
    names(resample_df)[i] = names(models)[i]
  }
  resample_df = do.call(rbind, resample_df)
  
  if(metric == "Accuracy"){
    #Prepare data frame for ploting
    resample_df <- resample_df %>%
      group_by(model) %>%
      summarise(Accuracy = mean(Accuracy)) 
  }else if(metric == "AUC"){
    #Prepare data frame for ploting
    resample_df <- resample_df %>%
      group_by(model) %>%
      summarise(AUC = mean(AUC)) 
  }
  
  resample_df <- resample_df %>%
    mutate(Category = case_when(
      model %in% c("BAG", "C50", "CART", "RF") ~ "Tree-based Methods",
      model %in% c("GLM", "LDA", "GLMNET", "LASSO", "RIDGE") ~ "Linear Models",
      model %in% c("KNN", "SVM_linear", "SVM_radial") ~ "Instance-based Methods",
      TRUE ~ "Other"  # In case there are models not in the above lists
    ))
  
  if(metric == "Accuracy"){
    groupped_df <- resample_df %>%
      group_by(Category) %>%
      filter(Accuracy == max(Accuracy)) %>%
      ungroup()     
  }else if(metric == "AUC"){
    groupped_df <- resample_df %>%
      group_by(Category) %>%
      filter(AUC == max(AUC)) %>%
      ungroup()     
  }else{
    stop("No metric defined")
  }
  
  #Retrieve top model based on accuracy/auc
  base_models <- groupped_df %>%
    pull(model)
  
  return(base_models)
}

calculate_auc_resample = function(obs, pred){
  
  prob_obs = data.frame("yes" = pred, "obs" = obs)
  
  prob_obs = prob_obs %>%
    arrange(desc(pred)) %>% #need to be arrange for apply cumulative sum
    mutate(is_yes = (obs == "yes"),
           tp = cumsum(is_yes), #true positive above the threshold - cumulative sum to refer to the threshold 
           fp = cumsum(!is_yes), #false positive above the threshold - cumulative sum to refer to the threshold
           fpr = fp/sum(obs == 'no'),
           tpr = tp/sum(obs == 'yes'))
  
  auc_value = calculate_auc(prob_obs$fpr, prob_obs$tpr)
  
  return(auc_value)
}

get_sensitivity_specificity = function(predictions, observed, ml.model){
  prob_obs = bind_cols(predictions, observed = observed) 
  
  prob_obs = prob_obs %>%
    arrange(desc(yes)) %>% #need to be arrange for apply cumulative sum
    mutate(is_yes = (observed == "yes"),
           tp = cumsum(is_yes), #true positive above the threshold - cumulative sum to refer to the threshold 
           fp = cumsum(!is_yes), #false positive above the threshold - cumulative sum to refer to the threshold
           sensitivity = tp/sum(observed == 'yes'),
           fpr = fp/sum(observed == 'no'),
           specificity = 1 - fpr) %>%
    select(sensitivity, specificity, fpr) %>%
    mutate(model = ml.model)
  
  # starts_at_zero <- any(prob_obs$sensitivity == 0 & prob_obs$fpr == 0)
  
  # ##Add dummy row if it doesnt start at 0
  # if(!starts_at_zero){
  #   dummy_row <- data.frame(
  #     sensitivity = 0,
  #     specificity = 1,
  #     fpr = 0,
  #     model = ml.model
  #   )
  # 
  #   prob_obs = rbind(dummy_row, prob_obs)
  # }
  
  prob_obs = prob_obs %>%
    mutate(Accuracy = calculate_accuracy(., observed),
           precision = calculate_precision(., observed),
           recall = calculate_recall(., observed)) 
  
  
  return(prob_obs)
  
}

#Take sensitivities values based on values of specificities
get_sensitivity = function(x, data){
  data %>%
    filter(specificity - x >= 0)%>% #Take specificity values above threshold x
    top_n(sensitivity, n=1) %>% #Take highest sensitivity from that threshold
    mutate(specificity = x, fpr = 1-x) %>% #Define sensitivity based on the specified threshold
    distinct() #If multiple thresholds have same sensitivity values take only one
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
        rownames_to_column("features") %>%
        dplyr::select(features, yes) %>% #Take only importance for positive class
        dplyr::rename(importance = yes)
    }else{
      base_importance_list[[i]] = base_importance[[base_models[i]]][["importance"]] %>%
        rownames_to_column("features") %>%
        dplyr::rename(importance = Overall)
    }
    names(base_importance_list)[i] = base_models[i]
  }
  
  #Combine all base model importances in one data frame and add the model name
  combined_importance <- bind_rows(
    lapply(names(base_importance_list), function(model) {
      base_importance_list[[model]] %>%
        data.frame() %>%
        mutate(model = model)
    })
  )
  
  #Calculate base-models importance for the meta-learner 
  meta_importance = varImp(meta_learner, scale = F)$importance %>%
    rownames_to_column("model")
  
  #Normalize the meta-learner's importance scores so they sum to 1 
  meta_importance$Overall <- meta_importance$Overall / sum(meta_importance$Overall)
  
  #Combine features importance within base models with the overall importance for meta-learner
  combined_importance <- combined_importance %>%
    left_join(meta_importance, by = "model") %>%
    mutate(weighted_importance = importance * Overall) # importance is from base, Overall is from meta
  
  #Sum the weighted importance by feature across all models
  final_importance <- combined_importance %>%
    group_by(features) %>%
    summarise(final_importance = sum(weighted_importance, na.rm = TRUE)) %>%
    arrange(desc(final_importance))
  
  return(final_importance)
  
}

compute.prediction = function(model, test_data, target){
  
  cat("Predicting target variable using provided ML model")
  
  features <- colnames(test_data)
  are_equal = setequal(model[["coefnames"]], features)
  if(are_equal == T){
    #Predict target variable
    predict <- data.frame(predict(model, test_data, type = "prob"))
    #Get metrics
    sens_spec = get_sensitivity_specificity(predict, target, model$method) 
    auroc = calculate_auroc(sens_spec$fpr, sens_spec$sensitivity)
    auprc = calculate_auprc(sens_spec$recall, sens_spec$precision)

    return(list(Metrics = sens_spec, AUC = list("AUROC" = auroc, "AUPRC" = auprc), Predictions = predict))
  }else{
    message("Testing set does not count with the same features as model")
  }
}

compute.prediction.stacked = function(super.learner, test_data, target, ml.models, base.models){
  
  #Learning from simple meta-learner
  base_predictions = list()
  for (i in 1:length(base.models)) {
    base_predictions[[i]] = predict(ml.models[[base.models[i]]], test_data, type = "prob")$yes
    names(base_predictions)[i] = base.models[i]
  }
  
  base_predictions = do.call(cbind, base_predictions)
  
  prediction_simple = data.frame(predict(super.learner, base_predictions, type = "prob")) 
  
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
  auroc_simple = calculate_auroc(sens_spec_simple$fpr, sens_spec_simple$sensitivity)
  auprc_simple = calculate_auprc(sens_spec_simple$recall, sens_spec_simple$precision)
  #Meta-learner all
  # sens_spec_all = get_sensitivity_specificity(prediction_all, target, "Meta-learner_all") 
  # auc_all = calculate_auc(sens_spec_all$fpr, sens_spec_all$sensitivity)
  
  #Not returning all (discarded)
  
  return(list(Metrics = sens_spec_simple, AUC = list("AUROC" = auroc_simple, "AUPRC" = auprc_simple), Predictions = prediction_simple))    
  
}

calculate_accuracy <- function(metrics, target) {
  sensitivity = metrics[,"sensitivity"]
  specificity = metrics[,"specificity"]
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
  sensitivity <- metrics[,"sensitivity"]
  specificity <- metrics[,"specificity"]
  
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

get_curves = function(data, spec, sens, reca, prec, color, auc_roc, auc_prc, file.name){
  
  data = data %>%
    mutate(specificity = data[,spec],
           sensitivity = data[,sens],
           recall = data[,reca],
           precision = data[,prec],
           color = data[,color])
  
  #Add AUC scores to data frame
  data$color.roc <- paste(data$color, "\n(AUC-ROC =", round(auc_roc, 2), ")\n")
  
  # Plot the ROC curves
  roc = ggplot(data = data, aes(x = 1- specificity, y = sensitivity, color = color.roc)) +
    geom_line() +  
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey") +
    labs(title = "ROC Curve", x = "1 - Specificity", y = "Sensitivity") +
    theme_minimal() +
    theme(legend.title = element_blank())
  
  #Add AUC scores to data frame
  data$color.prc <- paste(data$color, "\n(AUC-PRC =", round(auc_prc, 2), ")\n")
  
  # Plot recall curves
  recall = ggplot(data = data, aes(x = recall, y = precision, color = color.prc)) +
      geom_line() +
      labs(title = "Precision-Recall Curve", x = "Recall", y = "Precision") +
      ylim(0, 1) +  
      theme_minimal() +
      theme(legend.title = element_blank())

  pdf(paste0("Results/ROC_curve_", file.name, ".pdf"))
  print(roc)
  dev.off()
  
  pdf(paste0("Results/Recall_curve_", file.name, ".pdf"))
  print(recall)
  dev.off()

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
      filter(final_importance > 0) %>% #Filter features with no-importance/importance to classify positive class
      {
        if (nrow(.) < n_top) {
          warning("Not enough features for selecting n_top = ", n_top, " features in model ", file, ". Selecting all available features.\n")
          .  # Use all available rows
        } else {
          top_n(., n_top, wt = final_importance)  # Select the top n_top features
        }
      } %>%
      pull(features)
    
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
    
    #Check if available features were indeed n_top, if not adjust the matrix dynamically
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
  jaccard_dist = dist(presence_matrix, method="binary") #Jaccard distance
  hc <- hclust(jaccard_dist, method = "ward.D2") #Hierarchical clustering of distance
  matrica <- as.matrix(jaccard_dist) #Convert distance object to matrix 
  silhouette <- fviz_nbclust(matrica, FUNcluster = hcut, method = "silhouette", k.max = 6) #Identify k clusters from matrix
  k_cluster = as.numeric(silhouette$data$clusters[which.max(silhouette$data$y)]) #Extract k value
  sub_grp <- cutree(hc, k = k_cluster) #Obtain k cluster composition
  
  p = pheatmap(matrica,
               cluster_rows = hc,
               cluster_cols = hc,
               show_rownames = TRUE,
               show_colnames = TRUE,
               fontsize = 8,
               border_color = NA,
               color = hcl.colors(20, palette = "PRGn"),
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
      permuted_dist = dist(permuted_matrix, method="binary") 
      permuted_hc <- hclust(permuted_dist, method = "ward.D2")
      permuted_sub_grp <- cutree(permuted_hc, k = k_cluster)
      
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
    filter(Impact > 0) %>%
    pull(Feature)
  
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
  top_scores_df <- stack(scores)
    
  #Plot scores of the cell features across ML models
  p = ggplot(top_scores_df, aes(y = values, x = reorder(ind, -values, decreasing = F))) +
    geom_bar(stat = "identity", fill = "skyblue", width = 0.6) +
    theme_minimal() + 
    labs(title = paste0("Cell features presence scores across ", nrow(presence_matrix_important), " predictive cell groups"), 
         subtitle = paste0(ml_models, " ML models (AUC > ", AUC, ")"),
         x = "Cell features", 
         y = "Presence scores") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          plot.margin = ggplot2::margin(t = 10, r = 10, b = 30, l = 60)) +  
    scale_y_continuous(labels = scales::comma)  
  
  pdf(paste0("Results/Cell_features_scores_", file.name, ".pdf"), width = 8)
  print(p)
  dev.off()
  
  cat("Picking", n_cells, "cells types.........................................................\n")
  
  # Pick top_n cells 
  cells = top_scores_df %>% 
    arrange(desc(values)) %>%
    slice(1:n_cells) %>%
    pull(ind) %>%
    as.character()
    
  print(cells)
  
  #Calculate scores of presence for each TF
  scores_tfs = colSums(tfs_matrix)/nrow(tfs_matrix)

  #Select top scores above threshold
  top_scores_df <- stack(scores_tfs[scores_tfs > 0.8])
  
  entrz <- AnnotationDbi::select(org.Hs.eg.db, keys = as.character(top_scores_df$ind), columns = "ENTREZID", keytype = "SYMBOL") #Change to EntrezID

  paths_df = paths[,1:2]
  colnames(paths_df) = c("TERM", "GENE")
  
  # Perform ORA using enricher
  ora_results <- clusterProfiler::enricher(
    gene = as.character(top_scores_df$ind),
    TERM2GENE = paths_df,
    pvalueCutoff = 0.05
  )
  
  pdf(paste0("Results/Enrichment_TFs_CGs_", file.name, ".pdf"), width = 5)
  dotplot(ora_results) + ggtitle("ORA - Progeny pathways")
  dev.off()
  
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
      filter(final_importance > 0) %>% #Filter features with no-importance/importance to classify positive class
      {
        if (nrow(.) < n_top) {
          warning("Not enough features for selecting n_top = ", n_top, " features in model ", file, ". Selecting all available features.\n")
          .  # Use all available rows
        } else {
          top_n(., n_top, wt = final_importance)  # Select the top n_top features
        }
      } %>%
      pull(features)
    
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
  
  jaccard_dist = dist(presence_matrix, method="binary") #Jaccard distance
  hc <- hclust(jaccard_dist, method = "ward.D2") #Hierarchical clustering of distance
  matrica <- as.matrix(jaccard_dist) #Convert distance object to matrix 
  silhouette <- fviz_nbclust(matrica, FUNcluster = hcut, method = "silhouette", k.max = 6) #Identify k clusters from matrix
  k_cluster = as.numeric(silhouette$data$clusters[which.max(silhouette$data$y)]) #Extract k value
  sub_grp <- cutree(hc, k = k_cluster) #Obtain k cluster composition
  
  p = pheatmap(matrica,
               cluster_rows = hc,
               cluster_cols = hc,
               show_rownames = TRUE,
               show_colnames = TRUE,
               fontsize = 8,
               border_color = NA,
               color = hcl.colors(20, palette = "PRGn"),
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
      permuted_dist = dist(permuted_matrix, method="binary") 
      permuted_hc <- hclust(permuted_dist, method = "ward.D2")
      permuted_sub_grp <- cutree(permuted_hc, k = k_cluster)
      
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
    filter(Impact > 0) %>%
    pull(Feature)
  
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
  top_scores_df <- stack(scores)
  
  #Plot scores of the cell features across ML models
  p = ggplot(top_scores_df, aes(y = values, x = reorder(ind, -values, decreasing = F))) +
    geom_bar(stat = "identity", fill = "skyblue", width = 0.6) +
    theme_minimal() + 
    labs(title = paste0("Cell features presence scores across ", nrow(presence_matrix_important), " predictive cell groups"), 
         subtitle = paste0(ml_models, " ML models (AUC > ", AUC, ")"),
         x = "Cell features", 
         y = "Presence scores") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          plot.margin = ggplot2::margin(t = 10, r = 10, b = 30, l = 60)) +  
    scale_y_continuous(labels = scales::comma)  
  
  pdf(paste0("Results/Cell_features_scores_", file.name, ".pdf"), width = 8)
  print(p)
  dev.off()
  
  cat("Picking", n_cells, "cells types.........................................................\n")
  
  # Pick top_n cells 
  cells = top_scores_df %>% 
    arrange(desc(values)) %>%
    slice(1:n_cells) %>%
    pull(ind) %>%
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
    logreg = glm(trait$.outcome ~ features_values[,i], family = binomial) #Calculate logistic regression using cell groups values and trait outcome from training
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
    logreg = glm(trait$.outcome ~ features_values[,i], family = binomial) #Calculate logistic regression using cell groups values and trait outcome from training
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
