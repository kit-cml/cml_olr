# Function to calculate sigmoid
sigmoid <- function(x){
  return(1.0/(1.0 + exp(-x)))
}

# Function to calculate variable B
calculate_B1 <- function(alpha1, alpha2, beta1, beta2, feature1, feature2) {
  z1 <- alpha1 - beta1 * feature1 - beta2 * feature2
  z2 <- alpha2 - beta1 * feature1 - beta2 * feature2
  return(sigmoid(z1) - 0.5 * sigmoid(z2))
}

calculate_B2 <- function(alpha1, alpha2, beta1, beta2, feature1, feature2) {
  z1 <- alpha1 - beta1 * feature1 - beta2 * feature2
  z2 <- alpha2 - beta1 * feature1 - beta2 * feature2
  return(1.0 + sigmoid(z1) - 2.0 * sigmoid(z2))
}

# Function to calculate sigmoid
sigmoid <- function(x){
  return(1.0/(1.0 + exp(-x)))
}

# Function to calculate variable B
calculate_B1 <- function(alpha1, alpha2, beta1, beta2, feature1, feature2) {
  z1 <- alpha1 - beta1 * feature1 - beta2 * feature2
  z2 <- alpha2 - beta1 * feature1 - beta2 * feature2
  return(sigmoid(z1) - 0.5 * sigmoid(z2))
}

calculate_B2 <- function(alpha1, alpha2, beta1, beta2, feature1, feature2) {
  z1 <- alpha1 - beta1 * feature1 - beta2 * feature2
  z2 <- alpha2 - beta1 * feature1 - beta2 * feature2
  return(1.0 + sigmoid(z1) - 2.0 * sigmoid(z2))
}

run_all <- function(results_folder = 'results',
                    training = 'trainingdf',
                    testing = 'testingdf',
                    features_vector = 'features',
                    units_vector = 'units',
                    num_tests = 10000,
                    is_normalized = FALSE){
  # Determine if the analysis involves a single feature or multiple features
  is_single <- length(features) == 1
  
  # Construct a log file name based on the number of features
  log_filename <- if (!is_single) {
    paste(features_vector, collapse="-")
  } else {
    features_vector
  }
  logfile <- file(paste(results_folder,"/",paste(log_filename, "log.txt", sep="_"), sep = ""), open="wt")
  
  # Dynamically select columns for the analysis based on features
  # Include common columns like "label", "drug_name", "Sample_ID"
  common_cols <- c("label", "drug_name", "Sample_ID")
  training_cols <- as.character(c(features_vector, common_cols))
  testing_cols <- as.character(c(features_vector, common_cols))
  
  # Subset training and testing dataframes based on selected columns
  ordinaldf <- training[, training_cols]
  ordinaldf_test <- testing[, testing_cols]

  # Construct dataset for training and testing
  levels <- c("low", "intermediate", "high")
  values <- c(1, 2, 3)
  ordinaldf$label <- as.integer(factor(ordinaldf$label, levels = levels, labels = values))
  ordinaldf_test$label <- as.integer(factor(ordinaldf_test$label, levels = levels, labels = values))
  ordinaldf$label <- as.factor(ordinaldf$label)
  ordinaldf_test$label <- as.factor(ordinaldf_test$label)

  # Remove errors
  ordinaldf <- na.omit(ordinaldf)
  ordinaldf_test <- na.omit(ordinaldf_test)
  
  # Normalize the data
  dimension <- length(features_vector)
  if (is_normalized) {
    # Calculate mean and SD for the first m columns of ordinaldf
    means <- sapply(ordinaldf[1:dimension], mean)
    sds <- sapply(ordinaldf[1:dimension], sd)
    tempdf <- ordinaldf
    tempdf[1:dimension] <- mapply(function(x, mean, sd) (x - mean) / sd, ordinaldf[1:dimension], means, sds)
    ordinaldf <- tempdf
    tempdf <- ordinaldf_test
    tempdf[1:dimension] <- mapply(function(x, mean, sd) (x - mean) / sd, ordinaldf_test[1:dimension], means, sds)
    ordinaldf_test <- tempdf
  }
  
  # Fit the ordinal logistic regression model
  # Prepare the formula for logistic regression based on the number of features
  formula_string <- paste("label ~", paste(features_vector[1:dimension], collapse = " + "))
  formula <- as.formula(formula_string)

  # Specify the number of attempts
  max_attempts <- 10000

  # Flag to check if the model has converged
  converged <- FALSE

  # First trial without start option
  tryCatch({
    # Attempt to fit the model without specifying start
    mod <- polr(formula, data = ordinaldf, Hess = TRUE)
    converged <- TRUE
  }, error = function(e) {
    write(sprintf("First trial without 'start' failed. Retrying..."),logfile)
  })

  # If the first trial fails, attempt to fit with random starting values
  attempts <- as.integer(0)
  if (!converged) {
    for (attemp_idx in 1:max_attempts) {
      # Generate random starting values
      start_values <- rnorm(dimension + 2, mean = 0.0, sd = 100.0)

      # Attempt to fit the model with random starting values
      tryCatch({
        mod <- polr(formula, data = ordinaldf, start = start_values, Hess = TRUE)
        converged <- TRUE
        break
      }, error = function(e) {
        attempts <- attempts + 1
      })
    }
  }

  # If all attempts fail, print a message
  if (!converged) {
    write(sprintf("Could not fit the model after %d attempts",max_attempts),logfile)

  } else {
    # Make predictions for training dataset
    predicted_labels <- predict(mod, newdata = ordinaldf, type = "class")

    # Variables from Ordinal Logistic Regression model
    alphas <- as.numeric(mod$zeta)  # Alpha values
    betas <- as.numeric(mod$coefficients)  # Beta coefficients for all features
    
    if (!is_single) {
      # Apply TMS function row-wise efficiently
      ordinaldf$TMS <- apply(ordinaldf[, 1:dimension], 1, function(row) TMS(alphas, betas, row))
      ordinaldf_test$TMS <- apply(ordinaldf_test[, 1:dimension], 1, function(row) TMS(alphas, betas, row))
    }
  }

  if (dimension == 2) {
    scatterplotfun(data = ordinaldf,
                   mod = mod,
                   feature_1 = as.character(features_vector[1]),
                   feature_2 = as.character(features_vector[2]),
                   unit_1 = as.character(units_vector[1]),
                   unit_2 = as.character(units_vector[2]),
                   is_converged = converged,
                   is_training = TRUE,
                   is_legend = FALSE,
                   results_folder = results_folder,
                   is_normalized = is_normalized)
    scatterplotfun(data = ordinaldf_test,
                   mod = mod,
                   feature_1 = as.character(features_vector[1]),
                   feature_2 = as.character(features_vector[2]),
                   unit_1 = as.character(units_vector[1]),
                   unit_2 = as.character(units_vector[2]),
                   is_converged = converged,
                   is_training = FALSE,
                   is_legend = FALSE,
                   results_folder = results_folder,
                   is_normalized = is_normalized)
  }
  # Constants of the TMS
  if (converged) {
    # if (!is_single) {
      th1 <- (alphas[2] - alphas[1] ) / 2.0 + log(1.0 - 2.0 * exp(alphas[1] - alphas[2]))
      th2 <- -th1

    # Plot TMS for training dataset
    label_colors <- c("low" = "green", "intermediate" = "blue", "high" = "red")
    ordinaldf$risk <- ifelse(ordinaldf$label == 1, "low",
                             ifelse(ordinaldf$label == 2, "intermediate", "high"))

    ordinaldf_test$risk <- ifelse(ordinaldf_test$label == 1, "low",
                                  ifelse(ordinaldf_test$label == 2, "intermediate", "high"))
    # if (!is_single) {
      tms_name <- "TMS"
      filename_training <- paste(results_folder,"/",paste(log_filename, "training_dataset_tms.jpg", sep="_"), sep = "")
      filename_testing <- paste(results_folder,"/",paste(log_filename, "testing_dataset_tms.jpg", sep="_"), sep = "")
      # Plot TMS for training dataset
      tmsplotfun(data = ordinaldf,
                 th1 = th1,
                 th2 = th2,
                 label_colors = label_colors,
                 title = "Training dataset",
                 file_name = filename_training,
                 tms_name = "TMS")
      tmsplotfun(data = ordinaldf_test,
                 th1 = th1,
                 th2 = th2,
                 label_colors = label_colors,
                 title = "Testing dataset",
                 file_name = filename_testing,
                 tms_name = "TMS")
  } else {
    return(NA)
  }

  # Preallocate the pmeasures dataframe
  pmeasures <- data.frame()

  # Preallocate the drug test predictions
  drugs_predict <- data.frame(
    azimilide = numeric(num_tests),
    disopyramide = numeric(num_tests),
    ibutilide = numeric(num_tests),
    vandetanib = numeric(num_tests),
    astemizole = numeric(num_tests),
    clarithromycin = numeric(num_tests),
    clozapine = numeric(num_tests),
    domperidone = numeric(num_tests),
    droperidol = numeric(num_tests),
    pimozide = numeric(num_tests),
    risperidone = numeric(num_tests),
    loratadine = numeric(num_tests),
    metoprolol = numeric(num_tests),
    nifedipine = numeric(num_tests),
    nitrendipine = numeric(num_tests),
    tamoxifen = numeric(num_tests)
  )

  # Preallocate the drug TSM calculation
  drugs_tms <- data.frame(
    azimilide = numeric(num_tests),
    disopyramide = numeric(num_tests),
    ibutilide = numeric(num_tests),
    vandetanib = numeric(num_tests),
    astemizole = numeric(num_tests),
    clarithromycin = numeric(num_tests),
    clozapine = numeric(num_tests),
    domperidone = numeric(num_tests),
    droperidol = numeric(num_tests),
    pimozide = numeric(num_tests),
    risperidone = numeric(num_tests),
    loratadine = numeric(num_tests),
    metoprolol = numeric(num_tests),
    nifedipine = numeric(num_tests),
    nitrendipine = numeric(num_tests),
    tamoxifen = numeric(num_tests),
    bepridil = numeric(num_tests),
    dofetilide = numeric(num_tests),
    quinidine = numeric(num_tests),
    sotalol = numeric(num_tests),
    chlorpromazine = numeric(num_tests),
    cisapride = numeric(num_tests),
    ondansetron = numeric(num_tests),
    terfenadine = numeric(num_tests),
    diltiazem = numeric(num_tests),
    mexiletine = numeric(num_tests),
    ranolazine = numeric(num_tests),
    verapamil = numeric(num_tests)
  )

  # Iterate through the num_tests tests
  for (i in 1:num_tests) {
    # Create a subset with 16 testing drug names
    selected_drugs <- unique(ordinaldf_test$drug_name)
    sampled_rows <- data.frame()
    for (drug in selected_drugs) {
      drug_data <- ordinaldf_test[ordinaldf_test$drug_name == drug, ]
      sampled_row <- drug_data[sample(nrow(drug_data), 1), ]
      sampled_rows <- rbind(sampled_rows, sampled_row)
    }
    test_data <- sampled_rows
    test_data["is_training"] <- 0

    # Create a subset with 12 training drug names
    selected_drugs <- unique(ordinaldf$drug_name)
    sampled_rows <- data.frame()
    for (drug in selected_drugs) {
      drug_data <- ordinaldf[ordinaldf$drug_name == drug, ]
      sampled_row <- drug_data[sample(nrow(drug_data), 1), ]
      sampled_rows <- rbind(sampled_rows, sampled_row)
    }
    training_data <- sampled_rows
    training_data["is_training"] <- 1

    # Calculate all pmeasures
    temppmeasures <- pmeasuresfun(training_data = training_data,
                              test_data = test_data,
                              mod = mod,
                              label_values = values,
                              tms_name = tms_name,
                              is_whole = FALSE)

    # Store the calculated pmeasures
    pmeasures <- rbind(pmeasures,temppmeasures)

    # Fill the drug prediction row by row
    predicted_labels_test <- predict(mod, newdata = test_data, type = "class")
    new_test_data <- cbind(test_data,predicted_labels_test)
    for (drug in test_data$drug_name) {
      temp_data <- subset(new_test_data,new_test_data$drug_name == drug)
      drugs_predict[i, drug] <- as.numeric(temp_data$label == temp_data$predicted_labels_test)
    }

    if (!is_single) {
      # Fill the drug TMS row by row
      for (drug in test_data$drug_name) {
        drugs_tms[i, drug] <- test_data$TMS[test_data$drug_name == drug]
      }
      for (drug in training_data$drug_name) {
        drugs_tms[i, drug] <- training_data$TMS[training_data$drug_name == drug]
      }
    }
  } # i-th test

  # Calculate the classification error from the whole training dataset
  predicted_labels_training <- predict(mod, newdata = ordinaldf, type = "class")
  pred_err_training <- (abs(as.integer(predicted_labels_training) - as.integer(ordinaldf$label)))

  # Calculate the classification error from the whole testing dataset
  predicted_labels_testing <- predict(mod, newdata = ordinaldf_test, type = "class")
  pred_err_testing <- (abs(as.integer(predicted_labels_testing) - as.integer(ordinaldf_test$label)))

  # Calculate the pmeasures for the whole dataset
  ordinaldf["is_training"] = as.integer(1)
  ordinaldf_test["is_training"] = as.integer(0)
  pmeasures_training_data <- pmeasuresfun(training_data = ordinaldf,
                                   test_data = ordinaldf,
                                   mod = mod,
                                   label_values = values,
                                   tms_name = tms_name,
                                   is_whole = TRUE)
  pmeasures_testing_data <- pmeasuresfun(training_data = ordinaldf,
                                     test_data = ordinaldf_test,
                                     mod = mod,
                                     label_values = values,
                                     tms_name = tms_name,
                                     is_whole = TRUE)

  if (!is_single) {
    # Save the pmeasures dataframe
    write.csv(pmeasures, paste(results_folder,"/",paste(log_filename,"pmeasures.csv",sep = "_"), sep = ""), row.names = FALSE)

    # Save the drug_predict dataframe
    write.csv(drugs_predict, paste(results_folder,"/",paste(log_filename,"drugs_predict.csv",sep = "_"), sep = ""), row.names = FALSE)

    # Save the drug_tms dataframe
    write.csv(drugs_tms, paste(results_folder,"/",paste(log_filename,"drugs_tms.csv",sep = "_"), sep = ""), row.names = FALSE)
  } else {
    # Save the pmeasures dataframe
    write.csv(pmeasures, paste(results_folder,"/",paste(log_filename,"pmeasures.csv",sep = "_"), sep = ""), row.names = FALSE)
    
    # Save the drug_predict dataframe
    write.csv(drugs_predict, paste(results_folder,"/",paste(log_filename,"drugs_predict.csv",sep = "_"), sep = ""), row.names = FALSE)
  }

  writepmeasuresfun(pmeasures = pmeasures,
                  pmeasures_training_data = pmeasures_training_data,
                  pmeasures_testing_data = pmeasures_testing_data,
                  pred_error_training = pred_err_training,
                  pred_error_testing = pred_err_testing,
                  logfile = logfile,
                  is_single = is_single,
                  th1 = th1,
                  th2 = th2)

  # Close the connection to the text file
  close(logfile)

  # Rank score for the OLR model
  rank_score <- rankscorefun(pmeasures = pmeasures,
                             pred_error = pred_err_testing,
                             is_normalized = TRUE)

  # Store summary of pmeasures to summarydf
  feature_pair_name <- log_filename

  summarydf <- data.frame(
    Feature_Pair = feature_pair_name,
    Accuracy_Class_1 = quantile(pmeasures$Accuracy_Class_1, 0.025),
    Accuracy_Class_2 = quantile(pmeasures$Accuracy_Class_2, 0.025),
    Accuracy_Class_3 = quantile(pmeasures$Accuracy_Class_3, 0.025),
    AUC_Class_1 = quantile(pmeasures$AUC_Class_1, 0.025),
    AUC_Class_2 = quantile(pmeasures$AUC_Class_2, 0.025),
    AUC_Class_3 = quantile(pmeasures$AUC_Class_3, 0.025),
    Sensitivity_class_1 = quantile(pmeasures$Sensitivity_class_1, 0.025),
    Sensitivity_class_2 = quantile(pmeasures$Sensitivity_class_2, 0.025),
    Sensitivity_class_3 = quantile(pmeasures$Sensitivity_class_3, 0.025),
    Specificity_class_1 = quantile(pmeasures$Specificity_class_1, 0.025),
    Specificity_class_2 = quantile(pmeasures$Specificity_class_2, 0.025),
    Specificity_class_3 = quantile(pmeasures$Specificity_class_3, 0.025),
    LR_positive_class_1 = quantile(pmeasures$LR_positive_class_1, 0.025),
    LR_positive_class_2 = quantile(pmeasures$LR_positive_class_2, 0.025),
    LR_positive_class_3 = quantile(pmeasures$LR_positive_class_3, 0.025),
    LR_negative_class_1 = quantile(pmeasures$LR_negative_class_1, 0.975),
    LR_negative_class_2 = quantile(pmeasures$LR_negative_class_2, 0.975),
    LR_negative_class_3 = quantile(pmeasures$LR_negative_class_3, 0.975),
    F1score_class_1 = quantile(pmeasures$F1score_class_1, 0.025),
    F1score_class_2 = quantile(pmeasures$F1score_class_2, 0.025),
    F1score_class_3 = quantile(pmeasures$F1score_class_3, 0.025),
    Classification_error = mean(pred_err_testing) + 1.96 * sd(pred_err_testing) / sqrt(length(pred_err_testing)),
    Pairwise_classification_accuracy = quantile(pmeasures$Pairwise, 0.025),
    Rank_score = rank_score
  )
  
  # Add Alphas
  for (i in 1:length(alphas)) {
    summarydf[[paste0("Alpha_", i)]] <- alphas[i]
  }
  
  # Add Betas
  for (i in 1:length(betas)) {
    summarydf[[paste0("Beta_", i)]] <- betas[i]
  }
  
  return(summarydf)
}

TMS <- function(alphas, betas, features_row){
  z1 <- alphas[1] - sum(betas * features_row)
  z2 <- alphas[2] - sum(betas * features_row)

  return((z1 + z2) / 2.0)
}

pairwisefun<- function(fulltable){
  cmb <- combn(seq_len(nrow(fulltable)), 2)
  mergedtable<-cbind(fulltable[cmb[1,],], fulltable[cmb[2,],])
  validpairidx <- (mergedtable[,4]!=mergedtable[,10])&(!mergedtable[,6]|!mergedtable[,12])
  correctidx1 <- ((mergedtable[,4]>mergedtable[,10])&(mergedtable[,5]>mergedtable[,11]))|((mergedtable[,4]<mergedtable[,10])&(mergedtable[,5]<mergedtable[,11])) #when predicted class are different
  correctidx2 <- (mergedtable[,5]==3)&(mergedtable[,11]==3)&(((mergedtable[,4]>mergedtable[,10])&(mergedtable[,3]<mergedtable[,9]))|((mergedtable[,4]<mergedtable[,10])&(mergedtable[,3]>mergedtable[,9]))) #when predicted class are both high
  correctidx3 <- (mergedtable[,5]==1)&(mergedtable[,11]==1)&(((mergedtable[,4]>mergedtable[,10])&(mergedtable[,3]<mergedtable[,9]))|((mergedtable[,4]<mergedtable[,10])&(mergedtable[,3]>mergedtable[,9]))) #when predicted class are both low
  correctidx4 <- (mergedtable[,5]==2)&(mergedtable[,11]==2)&(((mergedtable[,4]>mergedtable[,10])&(mergedtable[,3]<mergedtable[,9]))|((mergedtable[,4]<mergedtable[,10])&(mergedtable[,3]>mergedtable[,9]))) #when predicted class are both intermediate
  correctidx <- correctidx1|correctidx2|correctidx3|correctidx4
  sum(validpairidx&correctidx)/sum(validpairidx)
}

aucrocfun <- function(data, mod, label_values){
  auc_scores <- c()
  for (class_label in label_values) {
    actual <- as.integer(data$label == class_label)
    predicted_prob <- predict(mod, newdata = data, type = "probs")
    predicted_prob <- predicted_prob[, class_label]
    roc_obj <- roc(actual,
                   predicted_prob,
                   direction = "<",
                   quiet = TRUE)
    auc_score <- auc(roc_obj)
    auc_scores <- c(auc_scores, auc_score)
  }
  return(auc_scores)
}

pmeasuresfun <- function(training_data, test_data, mod, label_values, tms_name = "TMS", is_whole){
  # Calculate the accuracy and AUC for each class
  predicted_labels_test <- predict(mod, newdata = test_data, type = "class")
  
  # Create confusion matrix
  confusion_matrix <- table(predicted_labels_test, test_data$label)
  
  # Calculate the accuracy for each class
  tp_class_1 = confusion_matrix[1,1]
  tn_class_1 = confusion_matrix[2,2] + confusion_matrix[2,3] + confusion_matrix[3,2] + confusion_matrix[3,3]
  fp_class_1 = confusion_matrix[1,2] + confusion_matrix[1,3]
  fn_class_1 = confusion_matrix[2,1] + confusion_matrix[3,1]
  
  tp_class_2 = confusion_matrix[2,2]
  tn_class_2 = confusion_matrix[1,1] + confusion_matrix[1,3] + confusion_matrix[3,1] + confusion_matrix[3,3]
  fp_class_2 = confusion_matrix[2,1] + confusion_matrix[2,3]
  fn_class_2 = confusion_matrix[1,2] + confusion_matrix[3,2]
  
  tp_class_3 = confusion_matrix[3,3]
  tn_class_3 = confusion_matrix[1,1] + confusion_matrix[1,2] + confusion_matrix[2,1] + confusion_matrix[2,2]
  fp_class_3 = confusion_matrix[3,1] + confusion_matrix[3,2]
  fn_class_3 = confusion_matrix[1,3] + confusion_matrix[2,3]
  
  f1score_class_1 <- 2.0 * tp_class_1 / (2.0 * tp_class_1 + fp_class_1 + fn_class_1)
  f1score_class_2 <- 2.0 * tp_class_2 / (2.0 * tp_class_2 + fp_class_2 + fn_class_2)
  f1score_class_3 <- 2.0 * tp_class_3 / (2.0 * tp_class_3 + fp_class_3 + fn_class_3)
  
  accuracy_class_1 <- (tp_class_1 + tn_class_1) / sum(confusion_matrix)
  accuracy_class_2 <- (tp_class_2 + tn_class_2) / sum(confusion_matrix)
  accuracy_class_3 <- (tp_class_3 + tn_class_3) / sum(confusion_matrix)
  
  sensitivity_class_1 <- tp_class_1 / (tp_class_1 + fn_class_1)
  sensitivity_class_2 <- tp_class_2 / (tp_class_2 + fn_class_2)
  sensitivity_class_3 <- tp_class_3 / (tp_class_3 + fn_class_3)
  
  specificity_class_1 <- tn_class_1 / (tn_class_1 + fp_class_1)
  specificity_class_2 <- tn_class_2 / (tn_class_2 + fp_class_2)
  specificity_class_3 <- tn_class_3 / (tn_class_3 + fp_class_3)
  
  # Add random number to LR+ and LR- to prevent zero devision
  u <- 1e-6
  sd <- 1e-12
  
  lr_positive_class_1 <- (sensitivity_class_1 + rnorm(1, mean = u, sd = sd)) / (1 - specificity_class_1 + rnorm(1, mean = u, sd = sd))
  lr_positive_class_2 <- (sensitivity_class_2 + rnorm(1, mean = u, sd = sd)) / (1 - specificity_class_2 + rnorm(1, mean = u, sd = sd))
  lr_positive_class_3 <- (sensitivity_class_3 + rnorm(1, mean = u, sd = sd)) / (1 - specificity_class_3 + rnorm(1, mean = u, sd = sd))
  
  lr_negative_class_1 <- (1 - sensitivity_class_1 + rnorm(1, mean = u, sd = sd)) / (specificity_class_1 + rnorm(1, mean = u, sd = sd))
  lr_negative_class_2 <- (1 - sensitivity_class_2 + rnorm(1, mean = u, sd = sd)) / (specificity_class_2 + rnorm(1, mean = u, sd = sd))
  lr_negative_class_3 <- (1 - sensitivity_class_3 + rnorm(1, mean = u, sd = sd)) / (specificity_class_3 + rnorm(1, mean = u, sd = sd))
  
  auc_scores <- aucrocfun(test_data, mod, label_values)
  
  # Calculate the pairwise classification error
  complete_data <- rbind(test_data,training_data)
  complete_data["pred"] <- predict(mod, newdata = complete_data, type = "class")
  complete_data["risk_label"] <- as.integer(complete_data$label)
  complete_data["risk_pred"] <- as.integer(complete_data$pred)
  complete_data <- complete_data[,c("Sample_ID","drug_name",tms_name,"risk_label","risk_pred","is_training")]
  if (is_whole) {
    pairwise <- NA
  } else{
    pairwise <- pairwisefun(complete_data)
  }
  
  # Fill the pmeasures row by row
  pmeasures <- data.frame(
    TP_class_1 = tp_class_1,
    TP_class_2 = tp_class_2,
    TP_class_3 = tp_class_3,
    TN_class_1 = tn_class_1,
    TN_class_2 = tn_class_2,
    TN_class_3 = tn_class_3,
    FP_class_1 = fp_class_1,
    FP_class_2 = fp_class_2,
    FP_class_3 = fp_class_3,
    FN_class_1 = fn_class_1,
    FN_class_2 = fn_class_2,
    FN_class_3 = fn_class_3,
    F1score_class_1 = f1score_class_1,
    F1score_class_2 = f1score_class_2,
    F1score_class_3 = f1score_class_3,
    Accuracy_Class_1 = accuracy_class_1,
    Accuracy_Class_2 = accuracy_class_2,
    Accuracy_Class_3 = accuracy_class_3,
    AUC_Class_1 = auc_scores[1],
    AUC_Class_2 = auc_scores[2],
    AUC_Class_3 = auc_scores[3],
    Sensitivity_class_1 = sensitivity_class_1,
    Sensitivity_class_2 = sensitivity_class_2,
    Sensitivity_class_3 = sensitivity_class_3,
    Specificity_class_1 = specificity_class_1,
    Specificity_class_2 = specificity_class_2,
    Specificity_class_3 = specificity_class_3,
    LR_positive_class_1 = lr_positive_class_1,
    LR_positive_class_2 = lr_positive_class_2,
    LR_positive_class_3 = lr_positive_class_3,
    LR_negative_class_1 = lr_negative_class_1,
    LR_negative_class_2 = lr_negative_class_2,
    LR_negative_class_3 = lr_negative_class_3,
    Pairwise = pairwise
  )
  
  return(pmeasures)
}

tmsplotfun <- function(data, th1, th2, label_colors, title, file_name, tms_name, tms_unit){
  data$drug_name <- factor(data$drug_name, levels = unique(data$drug_name[order(data$label)]))
  tms <- tms_name
  plot <- ggplot(data, aes_string(x = tms_name, y = "drug_name", fill = "risk")) +
    geom_boxplot(color = "black", width = 0.5, size = 0.2, outlier.size = 0.5, outlier.shape = NA) +
    labs(title = title, x = tms, y = "") +
    geom_vline(xintercept = th1, linetype = "dashed", color = "blue", size = 1)  +
    geom_vline(xintercept = th2, linetype = "dashed", color = "red", size = 1)  +
    scale_fill_manual(values = label_colors) + # Set the fill colors
    theme(plot.title = element_text(size = 20), # Title font size
          # Change axis title font sizes
          axis.title.x = element_text(size = 14), # X axis title font size
          axis.title.y = element_text(size = 14), # Y axis title font size
          # Change axis text font sizes
          axis.text.x = element_text(size = 12), # X axis text font size
          axis.text.y = element_text(size = 12), # Y axis text font size
          # Change legend title and text font sizes
          legend.title = element_text(size = 10), # Legend title font size
          legend.text = element_text(size = 8) # Legend text font size
    )
  ggsave(file_name, plot, width = 8, height = 6, dpi = 900)
}

scatterplotfun <- function(data, 
                           mod = NA, 
                           feature_1, 
                           feature_2, 
                           unit_1, 
                           unit_2, 
                           is_converged, 
                           is_training, 
                           is_legend, 
                           results_folder,
                           is_normalized){
  # Check the column index
  data <- data.frame(data)
  idx_model_1 <- as.integer(which(colnames(data) == feature_1))
  idx_model_2 <- as.integer(which(colnames(data) == feature_2))
  idx_label <- as.integer(which(colnames(data) == "label"))
  
  if (is_converged) {
    # Variables from Ordinal Logistic Regression model
    alpha1 <- as.numeric(mod$zeta[1])
    alpha2 <- as.numeric(mod$zeta[2])
    beta1 <- as.numeric(mod$coefficients[1])
    beta2 <- as.numeric(mod$coefficients[2])
    
    # Some descriptions of decision boundaries
    m <- - beta1 / beta2
    c1 <- - 1.0 / beta2 * (- alpha1 + log(1.0 - 2.0 * exp(alpha1 - alpha2)))
    c2 <- 1.0 / beta2 * (alpha1 + log(exp(alpha2 - alpha1) - 2.0))
  }
  
  # Create a meshgrid for contour plotting testing dataset
  x <- seq(min(data[,idx_model_1]), max(data[,idx_model_1]), length.out = 100)
  y <- seq(min(data[,idx_model_2]), max(data[,idx_model_2]), length.out = 100)
  if (is_converged) {
    z1 <- outer(x, y, Vectorize(function(x, y) calculate_B1(alpha1, alpha2, beta1, beta2, x, y)))
    z2 <- outer(x, y, Vectorize(function(x, y) calculate_B2(alpha1, alpha2, beta1, beta2, x, y)))
  }
  if (is_training) {
    title <- "Training dataset"
    file_name <- "training"
  } else {
    title <- "Testing dataset"
    file_name <- "testing"
  }
  if (is_normalized) {
    unit_1 <- ""
    unit_2 <- ""
  }
  jpeg(paste(results_folder,"/",paste(feature_1,feature_2,file_name,"dataset.jpg",sep = "_"), sep = ""),quality = 100, units = "in", width = 5, height = 5, res = 900)
  plot(data[,idx_model_1], 
       data[,idx_model_2], 
       xlab = paste(feature_1, unit_1, sep = " "), 
       ylab = paste(feature_2, unit_2, sep = " "),
       main = title,
       cex.axis = 1.5, 
       cex.lab = 1.5, 
       cex.main = 1.5, 
       cex = 0.5)
  if (is_legend) {
    legend("bottomright", legend = c("Low", "Intermediate", "High"), fill = c("green", "blue", "red"))
  }
  points(data[,idx_model_1][data$label == "1"], data[,idx_model_2][data$label == "1"], col = "green", cex = 0.5)
  points(data[,idx_model_1][data$label == "2"], data[,idx_model_2][data$label == "2"], col = "blue", cex = 0.5)
  points(data[,idx_model_1][data$label == "3"], data[,idx_model_2][data$label == "3"], col = "red", cex = 0.5)
  if (is_converged) {
    abline(a = c1, b = m, col = "blue", lty = 2, lwd = 2)  # Replace 'a' with the intercept if needed
    abline(a = c2, b = m, col = "red", lty = 2, lwd = 2)  # Replace 'a' with the intercept if needed
  }
  dev.off()
}

pairsdfinitfun <- function(features, units, dimension) {
  if (dimension > length(features)) {
    stop("Dimension cannot be greater than the number of features.")
  }

  # Calculate all possible combinations of features and units
  feature_combinations <- combn(features, dimension, simplify = FALSE)
  unit_combinations <- combn(units, dimension, simplify = FALSE)

  # Initialize an empty list to store data frames for each combination
  pairsdf_list <- list()

  # Loop through each combination and create a data frame
  for (i in seq_along(feature_combinations)) {
    feature_combination <- feature_combinations[[i]]
    unit_combination <- unit_combinations[[i]]

    # Create a data frame with the feature names and units for this combination
    feature_names <- paste0("feature_", seq_len(dimension))
    unit_names <- paste0("unit_", seq_len(dimension))
    tempdf <- setNames(data.frame(t(feature_combination)), feature_names)
    unitdf <- setNames(data.frame(t(unit_combination)), unit_names)

    # Combine features and units into one data frame
    combinedf <- cbind(tempdf, unitdf)

    # Add the combined data frame to the list
    pairsdf_list[[i]] <- combinedf
  }

  # Combine all data frames in the list into a single data frame
  pairsdf <- do.call("rbind", pairsdf_list)

  return(pairsdf)
}

solvepair <- function(results_folder,
                      pair_id,
                      filepath_training,
                      filepath_testing,
                      features_vector,
                      units_vector,
                      num_tests, 
                      is_normalized){
  # Print the current pair_id and features being processed
  print(paste(c(pair_id, features_vector), collapse = "_"))

  # Read in the training and testing datasets
  training <- read_csv(filepath_training, show_col_types = FALSE)
  testing <- read_csv(filepath_testing, show_col_types = FALSE)

  # Perform all
  resultdf <- run_all(results_folder = results_folder,
                      training = training,
                      testing = testing,
                      features_vector = features_vector,
                      units_vector = units_vector,
                      num_tests = num_tests,
                      is_normalized = is_normalized)
  return(resultdf)
}

writepmeasuresfun <- function(pmeasures, 
                            pmeasures_training_data, 
                            pmeasures_testing_data, 
                            pred_error_training,
                            pred_error_testing,
                            logfile,
                            is_single = FALSE,
                            th1 = NA,
                            th2 = NA){
  
  # Print the Thresholds into logfile
  write("=======================",logfile)
  write(sprintf("TMS thresholds"),logfile)
  write("=======================",logfile)
  write(paste(sprintf('Threshold_1: %.4f ', th1)), logfile)
  write(paste(sprintf('Threshold_2: %.4f ', th2)), logfile)

  # Print the pmeasures dataframe into logfile
  write("=======================",logfile)
  write(sprintf("Performance measures for %d test",num_tests),logfile)
  write("=======================",logfile)
  write(sprintf('AUC of low risk: %.4f (%.4f, %.4f)',
                median(pmeasures$AUC_Class_1),
                quantile(pmeasures$AUC_Class_1, 0.025),
                quantile(pmeasures$AUC_Class_1, 0.975)),
        logfile)
  write(sprintf('AUC of intermediate risk: %.4f (%.4f, %.4f)',
                median(pmeasures$AUC_Class_2),
                quantile(pmeasures$AUC_Class_2, 0.025),
                quantile(pmeasures$AUC_Class_2, 0.975)),
        logfile)
  write(sprintf('AUC of high risk: %.4f (%.4f, %.4f)',
                median(pmeasures$AUC_Class_3),
                quantile(pmeasures$AUC_Class_3, 0.025),
                quantile(pmeasures$AUC_Class_3, 0.975)),
        logfile)
  write(sprintf('Accuracy of low risk: %.4f (%.4f, %.4f)',
                median(pmeasures$Accuracy_Class_1),
                quantile(pmeasures$Accuracy_Class_1, 0.025),
                quantile(pmeasures$Accuracy_Class_1, 0.975)),
        logfile)
  write(sprintf('Accuracy of intermediate risk: %.4f (%.4f, %.4f)',
                median(pmeasures$Accuracy_Class_2),
                quantile(pmeasures$Accuracy_Class_2, 0.025),
                quantile(pmeasures$Accuracy_Class_2, 0.975)),
        logfile)
  write(sprintf('Accuracy of high risk: %.4f (%.4f, %.4f)',
                median(pmeasures$Accuracy_Class_3),
                quantile(pmeasures$Accuracy_Class_3, 0.025),
                quantile(pmeasures$Accuracy_Class_3, 0.975)),
        logfile)
  write(sprintf('Sensitivity of low risk: %.4f (%.4f, %.4f)',
                median(pmeasures$Sensitivity_class_1),
                quantile(pmeasures$Sensitivity_class_1, 0.025),
                quantile(pmeasures$Sensitivity_class_1, 0.975)),
        logfile)
  write(sprintf('Sensitivity of intermediate risk: %.4f (%.4f, %.4f)',
                median(pmeasures$Sensitivity_class_2),
                quantile(pmeasures$Sensitivity_class_2, 0.025),
                quantile(pmeasures$Sensitivity_class_2, 0.975)),
        logfile)
  write(sprintf('Sensitivity of high risk: %.4f (%.4f, %.4f)',
                median(pmeasures$Sensitivity_class_3),
                quantile(pmeasures$Sensitivity_class_3, 0.025),
                quantile(pmeasures$Sensitivity_class_3, 0.975)),
        logfile)
  write(sprintf('Specificity of low risk: %.4f (%.4f, %.4f)',
                median(pmeasures$Specificity_class_1),
                quantile(pmeasures$Specificity_class_1, 0.025),
                quantile(pmeasures$Specificity_class_1, 0.975)),
        logfile)
  write(sprintf('Specificity of intermediate risk: %.4f (%.4f, %.4f)',
                median(pmeasures$Specificity_class_2),
                quantile(pmeasures$Specificity_class_2, 0.025),
                quantile(pmeasures$Specificity_class_2, 0.975)),
        logfile)
  write(sprintf('Specificity of high risk: %.4f (%.4f, %.4f)',
                median(pmeasures$Specificity_class_3),
                quantile(pmeasures$Specificity_class_3, 0.025),
                quantile(pmeasures$Specificity_class_3, 0.975)),
        logfile)
  write(sprintf('LR+ of low risk: %.4f (%.4f, %.4f)',
                median(pmeasures$LR_positive_class_1),
                quantile(pmeasures$LR_positive_class_1, 0.025),
                quantile(pmeasures$LR_positive_class_1, 0.975)),
        logfile)
  write(sprintf('LR+ of intermediate risk: %.4f (%.4f, %.4f)',
                median(pmeasures$LR_positive_class_2),
                quantile(pmeasures$LR_positive_class_2, 0.025),
                quantile(pmeasures$LR_positive_class_2, 0.975)),
        logfile)
  write(sprintf('LR+ of high risk: %.4f (%.4f, %.4f)',
                median(pmeasures$LR_positive_class_3),
                quantile(pmeasures$LR_positive_class_3, 0.025),
                quantile(pmeasures$LR_positive_class_3, 0.975)),
        logfile)
  write(sprintf('LR- of low risk: %.4f (%.4f, %.4f)',
                median(pmeasures$LR_negative_class_1),
                quantile(pmeasures$LR_negative_class_1, 0.025),
                quantile(pmeasures$LR_negative_class_1, 0.975)),
        logfile)
  write(sprintf('LR- of intermediate risk: %.4f (%.4f, %.4f)',
                median(pmeasures$LR_negative_class_2),
                quantile(pmeasures$LR_negative_class_2, 0.025),
                quantile(pmeasures$LR_negative_class_2, 0.975)),
        logfile)
  write(sprintf('LR- of high risk: %.4f (%.4f, %.4f)',
                median(pmeasures$LR_negative_class_3),
                quantile(pmeasures$LR_negative_class_3, 0.025),
                quantile(pmeasures$LR_negative_class_3, 0.975)),
        logfile)
  write(sprintf('F1score of low risk: %.4f (%.4f, %.4f)',
                median(pmeasures$F1score_class_1),
                quantile(pmeasures$F1score_class_1, 0.025),
                quantile(pmeasures$F1score_class_1, 0.975)),
        logfile)
  write(sprintf('F1score of intermediate risk: %.4f (%.4f, %.4f)',
                median(pmeasures$F1score_class_2),
                quantile(pmeasures$F1score_class_2, 0.025),
                quantile(pmeasures$F1score_class_2, 0.975)),
        logfile)
  write(sprintf('F1score of high risk: %.4f (%.4f, %.4f)',
                median(pmeasures$F1score_class_3),
                quantile(pmeasures$F1score_class_3, 0.025),
                quantile(pmeasures$F1score_class_3, 0.975)),
        logfile)
  write(sprintf('Classification error: %.4f (%.4f, %.4f)',
                mean(pred_error_testing),
                mean(pred_error_testing) - 1.96 * sd(pred_error_testing) / sqrt(length(pred_error_testing)),
                mean(pred_error_testing) + 1.96 * sd(pred_error_testing) / sqrt(length(pred_error_testing))),
        logfile)
  write(sprintf('Pairwise classification accuracy: %.4f (%.4f, %.4f)',
                median(pmeasures$Pairwise),
                quantile(pmeasures$Pairwise, 0.025),
                quantile(pmeasures$Pairwise, 0.975)),
        logfile)
  
  write("==============================",logfile)
  write("Performance measures for whole training data",logfile)
  write("==============================",logfile)
  write(sprintf('AUC of low risk: %.4f',pmeasures_training_data$AUC_Class_1),logfile)
  write(sprintf('AUC of intermediate risk: %.4f',pmeasures_training_data$AUC_Class_2),logfile)
  write(sprintf('AUC of high risk: %.4f',pmeasures_training_data$AUC_Class_3),logfile)
  write(sprintf('Accuracy of low risk: %.4f',pmeasures_training_data$Accuracy_Class_1),logfile)
  write(sprintf('Accuracy of intermediate risk: %.4f',pmeasures_training_data$Accuracy_Class_2),logfile)
  write(sprintf('Accuracy of high risk: %.4f',pmeasures_training_data$Accuracy_Class_3),logfile)
  write(sprintf('Sensitivity of low risk: %.4f',pmeasures_training_data$Sensitivity_class_1),logfile)
  write(sprintf('Sensitivity of intermediate risk: %.4f',pmeasures_training_data$Sensitivity_class_2),logfile)
  write(sprintf('Sensitivity of high risk: %.4f',pmeasures_training_data$Sensitivity_class_3),logfile)
  write(sprintf('Specificity of low risk: %.4f',pmeasures_training_data$Specificity_class_1),logfile)
  write(sprintf('Specificity of intermediate risk: %.4f',pmeasures_training_data$Specificity_class_2),logfile)
  write(sprintf('Specificity of high risk: %.4f',pmeasures_training_data$Specificity_class_3),logfile)
  write(sprintf('LR+ of low risk: %.4f',pmeasures_training_data$LR_positive_class_1),logfile)
  write(sprintf('LR+ of intermediate risk: %.4f',pmeasures_training_data$LR_positive_class_2),logfile)
  write(sprintf('LR+ of high risk: %.4f',pmeasures_training_data$LR_positive_class_3),logfile)
  write(sprintf('LR- of low risk: %.4f',pmeasures_training_data$LR_negative_class_1),logfile)
  write(sprintf('LR- of intermediate risk: %.4f',pmeasures_training_data$LR_negative_class_2),logfile)
  write(sprintf('LR- of high risk: %.4f',pmeasures_training_data$LR_negative_class_3),logfile)
  write(sprintf('F1score of low risk: %.4f',pmeasures_training_data$F1score_class_1),logfile)
  write(sprintf('F1score of intermediate risk: %.4f',pmeasures_training_data$F1score_class_2),logfile)
  write(sprintf('F1score of high risk: %.4f',pmeasures_training_data$F1score_class_3),logfile)
  write(sprintf('Classification error: %.4f (%.4f, %.4f)',
                mean(pred_error_training),
                mean(pred_error_training) - 1.96 * sd(pred_error_training) / sqrt(length(pred_error_training)),
                mean(pred_error_training) + 1.96 * sd(pred_error_training) / sqrt(length(pred_error_training))),
        logfile)
  write(sprintf('Pairwise classification accuracy: %.4f',pmeasures_training_data$Pairwise),logfile)
  
  write("==============================",logfile)
  write("Performance measures for whole testing data",logfile)
  write("==============================",logfile)
  write(sprintf('AUC of low risk: %.4f',pmeasures_testing_data$AUC_Class_1),logfile)
  write(sprintf('AUC of intermediate risk: %.4f',pmeasures_testing_data$AUC_Class_2),logfile)
  write(sprintf('AUC of high risk: %.4f',pmeasures_testing_data$AUC_Class_3),logfile)
  write(sprintf('Accuracy of low risk: %.4f',pmeasures_testing_data$Accuracy_Class_1),logfile)
  write(sprintf('Accuracy of intermediate risk: %.4f',pmeasures_testing_data$Accuracy_Class_2),logfile)
  write(sprintf('Accuracy of high risk: %.4f',pmeasures_testing_data$Accuracy_Class_3),logfile)
  write(sprintf('Sensitivity of low risk: %.4f',pmeasures_testing_data$Sensitivity_class_1),logfile)
  write(sprintf('Sensitivity of intermediate risk: %.4f',pmeasures_testing_data$Sensitivity_class_2),logfile)
  write(sprintf('Sensitivity of high risk: %.4f',pmeasures_testing_data$Sensitivity_class_3),logfile)
  write(sprintf('Specificity of low risk: %.4f',pmeasures_testing_data$Specificity_class_1),logfile)
  write(sprintf('Specificity of intermediate risk: %.4f',pmeasures_testing_data$Specificity_class_2),logfile)
  write(sprintf('Specificity of high risk: %.4f',pmeasures_testing_data$Specificity_class_3),logfile)
  write(sprintf('LR+ of low risk: %.4f',pmeasures_testing_data$LR_positive_class_1),logfile)
  write(sprintf('LR+ of intermediate risk: %.4f',pmeasures_testing_data$LR_positive_class_2),logfile)
  write(sprintf('LR+ of high risk: %.4f',pmeasures_testing_data$LR_positive_class_3),logfile)
  write(sprintf('LR- of low risk: %.4f',pmeasures_testing_data$LR_negative_class_1),logfile)
  write(sprintf('LR- of intermediate risk: %.4f',pmeasures_testing_data$LR_negative_class_2),logfile)
  write(sprintf('LR- of high risk: %.4f',pmeasures_testing_data$LR_negative_class_3),logfile)
  write(sprintf('F1score of low risk: %.4f',pmeasures_testing_data$F1score_class_1),logfile)
  write(sprintf('F1score of intermediate risk: %.4f',pmeasures_testing_data$F1score_class_2),logfile)
  write(sprintf('F1score of high risk: %.4f',pmeasures_testing_data$F1score_class_3),logfile)
  write(sprintf('Pairwise classification accuracy: %.4f',pmeasures_testing_data$Pairwise),logfile)
}

rankscorefun <- function(pmeasures,
                         pred_error,
                         is_normalized = TRUE){
  Performance_measures <- c("AUC_Class_1", "AUC_Class_3",
                            "LR_positive_class_1", "LR_positive_class_3",
                            "LR_negative_class_1", "LR_negative_class_3",
                            "Pairwise","Classification_error")
  Performance_levels <- c("Excellent_performance", 
                          "Good_performance",
                          "Minimally_acceptable_performance",
                          "Not_acceptabel")
  
  # A dataframe for the weights of performance measures
  pm_df <- data.frame(
    Performance_measure = Performance_measures,
    Weight = c(1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 2.0, 2.0)
  )

  # A dataframe for the weights of performance levels
  pl_df <- data.frame(
    Performance_level = Performance_levels,
    Weight = c(3.0, 2.0, 1.0, 0.0)
  )
  if (is_normalized) {
    pm_df$Weight <- pm_df$Weight / sum(pm_df$Weight)
    pl_df$Weight <- pl_df$Weight / max(pl_df$Weight)
  }
  
  # Initialize the performance dataframe for the model
  model_df <- data.frame(
    Performance_measure = Performance_measures,
    Performance_level_weight = c(NA, NA, NA, NA, NA, NA, NA, NA)
  )
  
  # Check the performance level for each performance measure 
  # by looking at the 95% confidence interval that match the CiPA's criteria
  
  # AUC_Class_1
  score_to_check <- quantile(pmeasures$AUC_Class_1, 0.025)
  if (score_to_check < 0.7) {
    model_df[model_df$Performance_measure == "AUC_Class_1",]$Performance_level_weight <- pl_df[pl_df$Performance_level == "Not_acceptabel",]$Weight
  } else if (score_to_check >= 0.7 & score_to_check < 0.8) {
    model_df[model_df$Performance_measure == "AUC_Class_1",]$Performance_level_weight <- pl_df[pl_df$Performance_level == "Minimally_acceptable_performance",]$Weight
  } else if (score_to_check >= 0.8 & score_to_check < 0.9) {
    model_df[model_df$Performance_measure == "AUC_Class_1",]$Performance_level_weight <- pl_df[pl_df$Performance_level == "Good_performance",]$Weight
  } else {
    model_df[model_df$Performance_measure == "AUC_Class_1",]$Performance_level_weight <- pl_df[pl_df$Performance_level == "Excellent_performance",]$Weight
  }
  
  # AUC_Class_3
  score_to_check <- quantile(pmeasures$AUC_Class_3, 0.025)
  if (score_to_check < 0.7) {
    model_df[model_df$Performance_measure == "AUC_Class_3",]$Performance_level_weight <- pl_df[pl_df$Performance_level == "Not_acceptabel",]$Weight
  } else if (score_to_check >= 0.7 & score_to_check < 0.8) {
    model_df[model_df$Performance_measure == "AUC_Class_3",]$Performance_level_weight <- pl_df[pl_df$Performance_level == "Minimally_acceptable_performance",]$Weight
  } else if (score_to_check >= 0.8 & score_to_check < 0.9) {
    model_df[model_df$Performance_measure == "AUC_Class_3",]$Performance_level_weight <- pl_df[pl_df$Performance_level == "Good_performance",]$Weight
  } else {
    model_df[model_df$Performance_measure == "AUC_Class_3",]$Performance_level_weight <- pl_df[pl_df$Performance_level == "Excellent_performance",]$Weight
  }
  
  # LR_positive_class_1
  score_to_check <- quantile(pmeasures$LR_positive_class_1, 0.025)
  if (score_to_check < 2.0) {
    model_df[model_df$Performance_measure == "LR_positive_class_1",]$Performance_level_weight <- pl_df[pl_df$Performance_level == "Not_acceptabel",]$Weight
  } else if (score_to_check >= 2.0 & score_to_check < 5.0) {
    model_df[model_df$Performance_measure == "LR_positive_class_1",]$Performance_level_weight <- pl_df[pl_df$Performance_level == "Minimally_acceptable_performance",]$Weight
  } else if (score_to_check >= 5.0 & score_to_check < 10.0) {
    model_df[model_df$Performance_measure == "LR_positive_class_1",]$Performance_level_weight <- pl_df[pl_df$Performance_level == "Good_performance",]$Weight
  } else {
    model_df[model_df$Performance_measure == "LR_positive_class_1",]$Performance_level_weight <- pl_df[pl_df$Performance_level == "Excellent_performance",]$Weight
  }
  
  # LR_positive_class_3
  score_to_check <- quantile(pmeasures$LR_positive_class_3, 0.025)
  if (score_to_check < 2.0) {
    model_df[model_df$Performance_measure == "LR_positive_class_3",]$Performance_level_weight <- pl_df[pl_df$Performance_level == "Not_acceptabel",]$Weight
  } else if (score_to_check >= 2.0 & score_to_check < 5.0) {
    model_df[model_df$Performance_measure == "LR_positive_class_3",]$Performance_level_weight <- pl_df[pl_df$Performance_level == "Minimally_acceptable_performance",]$Weight
  } else if (score_to_check >= 5.0 & score_to_check < 10.0) {
    model_df[model_df$Performance_measure == "LR_positive_class_3",]$Performance_level_weight <- pl_df[pl_df$Performance_level == "Good_performance",]$Weight
  } else {
    model_df[model_df$Performance_measure == "LR_positive_class_3",]$Performance_level_weight <- pl_df[pl_df$Performance_level == "Excellent_performance",]$Weight
  }
  
  # LR_negative_class_1
  score_to_check <- quantile(pmeasures$LR_negative_class_1, 0.975)
  if (score_to_check > 0.5) {
    model_df[model_df$Performance_measure == "LR_negative_class_1",]$Performance_level_weight <- pl_df[pl_df$Performance_level == "Not_acceptabel",]$Weight
  } else if (score_to_check <= 0.5 & score_to_check > 0.2) {
    model_df[model_df$Performance_measure == "LR_negative_class_1",]$Performance_level_weight <- pl_df[pl_df$Performance_level == "Minimally_acceptable_performance",]$Weight
  } else if (score_to_check <= 0.2 & score_to_check > 0.1) {
    model_df[model_df$Performance_measure == "LR_negative_class_1",]$Performance_level_weight <- pl_df[pl_df$Performance_level == "Good_performance",]$Weight
  } else {
    model_df[model_df$Performance_measure == "LR_negative_class_1",]$Performance_level_weight <- pl_df[pl_df$Performance_level == "Excellent_performance",]$Weight
  }
  
  # LR_negative_class_3
  score_to_check <- quantile(pmeasures$LR_negative_class_3, 0.975)
  if (score_to_check > 0.5) {
    model_df[model_df$Performance_measure == "LR_negative_class_3",]$Performance_level_weight <- pl_df[pl_df$Performance_level == "Not_acceptabel",]$Weight
  } else if (score_to_check <= 0.5 & score_to_check > 0.2) {
    model_df[model_df$Performance_measure == "LR_negative_class_3",]$Performance_level_weight <- pl_df[pl_df$Performance_level == "Minimally_acceptable_performance",]$Weight
  } else if (score_to_check <= 0.2 & score_to_check > 0.1) {
    model_df[model_df$Performance_measure == "LR_negative_class_3",]$Performance_level_weight <- pl_df[pl_df$Performance_level == "Good_performance",]$Weight
  } else {
    model_df[model_df$Performance_measure == "LR_negative_class_3",]$Performance_level_weight <- pl_df[pl_df$Performance_level == "Excellent_performance",]$Weight
  }
  
  # Pairwise
  score_to_check <- quantile(pmeasures$Pairwise, 0.025)
  if (score_to_check < 0.7) {
    model_df[model_df$Performance_measure == "Pairwise",]$Performance_level_weight <- pl_df[pl_df$Performance_level == "Not_acceptabel",]$Weight
  } else if (score_to_check >= 0.7 & score_to_check < 0.8) {
    model_df[model_df$Performance_measure == "Pairwise",]$Performance_level_weight <- pl_df[pl_df$Performance_level == "Minimally_acceptable_performance",]$Weight
  } else if (score_to_check >= 0.8 & score_to_check < 0.9) {
    model_df[model_df$Performance_measure == "Pairwise",]$Performance_level_weight <- pl_df[pl_df$Performance_level == "Good_performance",]$Weight
  } else {
    model_df[model_df$Performance_measure == "Pairwise",]$Performance_level_weight <- pl_df[pl_df$Performance_level == "Excellent_performance",]$Weight
  }
  
  # Classification_error
  score_to_check <- mean(pred_error) + 1.96 * sd(pred_error) / sqrt(length(pred_error))
  if (score_to_check > 1.0) {
    model_df[model_df$Performance_measure == "Classification_error",]$Performance_level_weight <- pl_df[pl_df$Performance_level == "Not_acceptabel",]$Weight
  } else if (score_to_check <= 1.0 & score_to_check > 0.5) {
    model_df[model_df$Performance_measure == "Classification_error",]$Performance_level_weight <- pl_df[pl_df$Performance_level == "Minimally_acceptable_performance",]$Weight
  } else if (score_to_check <= 0.5 & score_to_check > 0.3) {
    model_df[model_df$Performance_measure == "Classification_error",]$Performance_level_weight <- pl_df[pl_df$Performance_level == "Good_performance",]$Weight
  } else {
    model_df[model_df$Performance_measure == "Classification_error",]$Performance_level_weight <- pl_df[pl_df$Performance_level == "Excellent_performance",]$Weight
  }
  
  # Calculate rank score
  rank_score <- model_df$Performance_level_weight %*% pm_df$Weight
  
  return(as.numeric(rank_score))
  }