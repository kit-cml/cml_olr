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

# Ordinal logistic regression
run_all <- function(training_1 = 'filepath_training_1',testing_1 = 'filepath_testing_1',
                    training_2 = 'filepath_training_2',testing_2 = 'filepath_testing_2',
                    feature_1 = 'feature_1',feature_2 = 'feature_2',
                    cell_model_1 = 'cell_model_1',cell_model_2 = 'cell_model_2',
                    num_tests = 10000){
  # Rename the feature according to the cell model
  feature_model_1 <- paste(feature_1,cell_model_1,sep = "_")
  feature_model_2 <- paste(feature_2,cell_model_2,sep = "_")
  
  # Log file
  logfile <- file(paste(feature_model_1,feature_model_2,"log.txt",sep = "_"),open = "wt")
  
  # Load training data and change colnames
  trainingdf_1 <- training_1[,c(feature_1,"label","drug_name","Sample_ID")]
  colnames(trainingdf_1)[1] <- feature_model_1
  trainingdf_2 <- training_2[,c(feature_2,"label","drug_name","Sample_ID")]
  colnames(trainingdf_2)[1] <- feature_model_2
  ordinaldf <- merge(trainingdf_1,trainingdf_2, by = c("Sample_ID","drug_name","label"))
  
  testingdf_1 <- testing_1[,c(feature_1,"label","drug_name","Sample_ID")]
  colnames(testingdf_1)[1] <- feature_model_1
  testingdf_2 <- testing_2[,c(feature_2,"label","drug_name","Sample_ID")]
  colnames(testingdf_2)[1] <- feature_model_2
  ordinaldf_test <- merge(testingdf_1,testingdf_2, by = c("Sample_ID","drug_name","label"))
  
  # Construct dataset for training and testing
  levels <- c("low", "intermediate", "high")
  values <- c(1, 2, 3)
  ordinaldf$label <- as.integer(factor(ordinaldf$label, levels = levels, labels = values))
  ordinaldf_test$label <- as.integer(factor(ordinaldf_test$label, levels = levels, labels = values))
  ordinaldf$label <- as.factor(ordinaldf$label)
  ordinaldf_test$label <- as.factor(ordinaldf_test$label)
  
  # Check the column index
  idx_model_1 <- as.integer(which(colnames(ordinaldf) == feature_model_1))
  idx_model_2 <- as.integer(which(colnames(ordinaldf) == feature_model_2))
  idx_label <- as.integer(which(colnames(ordinaldf) == "label"))
  
  # Remove errors
  ordinaldf <- na.omit(ordinaldf)
  ordinaldf_test <- na.omit(ordinaldf_test)
  
  # Fit the ordinal logistic regression model
  formula_string <- paste("label~", feature_model_1, "+", feature_model_2, sep = "")
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
      start_values <- rnorm(4, mean = 0.0, sd = 100.0)
      
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
    alpha1 <- as.numeric(mod$zeta[1])
    alpha2 <- as.numeric(mod$zeta[2])
    beta1 <- as.numeric(mod$coefficients[1])
    beta2 <- as.numeric(mod$coefficients[2])
    
    # Some descriptions of decision boundaries
    m <- - beta1 / beta2
    c1 <- - 1.0 / beta2 * (- alpha1 + log(1.0 - 2.0 * exp(alpha1 - alpha2)))
    c2 <- 1.0 / beta2 * (alpha1 + log(exp(alpha2 - alpha1) - 2.0))
    
    # Calculate the TMS for training and testing dataset
    ordinaldf["TMS"] <- TMS(alpha1,alpha2,beta1,beta2,ordinaldf[,feature_model_1],ordinaldf[,feature_model_2])
    ordinaldf_test["TMS"] <- TMS(alpha1,alpha2,beta1,beta2,ordinaldf_test[,feature_model_1],ordinaldf_test[,feature_model_2])
  }
 
  scatterplotfun(data = ordinaldf,
                 mod = mod,
                 feature_model_1 = feature_model_1, 
                 feature_model_2 = feature_model_2, 
                 is_converged = converged, 
                 is_training = TRUE, 
                 is_legend = FALSE)
  scatterplotfun(data = ordinaldf_test,
                 mod = mod,
                 feature_model_1 = feature_model_1, 
                 feature_model_2 = feature_model_2, 
                 is_converged = converged, 
                 is_training = FALSE, 
                 is_legend = FALSE)

  # Constants of the TMS
  if (converged) {
    c1_tms <- c1 * beta2
    c2_tms <- c2 * beta2
    th1 <- (c2_tms - c1_tms) / 2.0
    th2 <- (c1_tms - c2_tms) / 2.0
    
    # Plot TMS for training dataset
    label_colors <- c("green", "blue", "red")
    tmsplotfun(data = ordinaldf, 
               th1 = th1, 
               th2 = th2, 
               label_colors = label_colors, 
               title = "Training dataset",
               file_name = paste(feature_model_1,feature_model_2,"training_dataset_tms.jpg",sep = "_"),
               tms_name = "TMS")

    # Plot TMS for testing dataset
    tmsplotfun(data = ordinaldf_test,
               th1 = th1, 
               th2 = th2, 
               label_colors = label_colors, 
               title = "Testing dataset",
               file_name = paste(feature_model_1,feature_model_2,"testing_dataset_tms.jpg",sep = "_"), 
               tms_name = "TMS")
  } else {
    return(NA)
  }
  
  # Preallocate the metrics dataframe
  metrics <- data.frame()
  
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
    
    # Calculate all metrics 
    tempmetrics <- metricsfun(training_data = training_data, 
                              test_data = test_data, 
                              mod = mod, 
                              label_values = values, 
                              tms_name = "TMS",
                              is_whole = FALSE)
    
    # Store the calculated metrics 
    metrics <- rbind(metrics,tempmetrics)

    # Fill the drug prediction row by row
    predicted_labels_test <- predict(mod, newdata = test_data, type = "class")
    new_test_data <- cbind(test_data,predicted_labels_test)
    for (drug in test_data$drug_name) {
      temp_data <- subset(new_test_data,new_test_data$drug_name == drug)
      drugs_predict[i, drug] <- as.numeric(temp_data$label == temp_data$predicted_labels_test)      
    }
    
    # Fill the drug TMS row by row
    for (drug in test_data$drug_name) {
      f1 <- test_data[,feature_model_1][test_data$drug_name == drug]
      f2 <- test_data[,feature_model_2][test_data$drug_name == drug]
      drugs_tms[i, drug] <- TMS(alpha1,alpha2,beta1,beta2,f1,f2)
    }
    for (drug in training_data$drug_name) {
      f1 <- training_data[,feature_model_1][training_data$drug_name == drug]
      f2 <- training_data[,feature_model_2][training_data$drug_name == drug]
      drugs_tms[i, drug] <- TMS(alpha1,alpha2,beta1,beta2,f1,f2)
    }
  } # i-th test
  
  # Calculate the classification error from the whole testing dataset
  predicted_labels <- predict(mod, newdata = ordinaldf_test, type = "class")
  pred_err <- (abs(as.integer(predicted_labels) - as.integer(ordinaldf_test$label)))
  
  # Calculate the metrics for the whole dataset
  ordinaldf["is_training"] = as.integer(1)
  ordinaldf_test["is_training"] = as.integer(0)
  metrics_whole_data <- metricsfun(training_data = ordinaldf, 
                                   test_data = ordinaldf_test,
                                   mod = mod,
                                   label_values = values,
                                   tms_name = "TMS",
                                   is_whole = TRUE)
  
  # Save the metrics dataframe
  write.csv(metrics, paste(feature_model_1,feature_model_2,"metrics.csv",sep = "_"), row.names = FALSE)
  
  # Save the drug_predict dataframe
  write.csv(drugs_predict, paste(feature_model_1,feature_model_2,"drugs_predict.csv",sep = "_"), row.names = FALSE)
  
  # Save the drug_tms dataframe
  write.csv(drugs_tms, paste(feature_model_1,feature_model_2,"drugs_tms.csv",sep = "_"), row.names = FALSE)
  
  # Print the metrics dataframe into logfile
  write("=======================",logfile)
  write(sprintf("Metrics for %d test",num_tests),logfile)
  write("=======================",logfile)
  write(sprintf('AUC of low risk: %.2f (%.2f, %.2f)',
                      median(metrics$AUC_Class_1),
                      quantile(metrics$AUC_Class_1, 0.025),
                      quantile(metrics$AUC_Class_1, 0.975)),
        logfile)
  write(sprintf('AUC of intermediate risk: %.2f (%.2f, %.2f)',
                median(metrics$AUC_Class_2),
                quantile(metrics$AUC_Class_2, 0.025),
                quantile(metrics$AUC_Class_2, 0.975)),
        logfile)
  write(sprintf('AUC of high risk: %.2f (%.2f, %.2f)',
                median(metrics$AUC_Class_3),
                quantile(metrics$AUC_Class_3, 0.025),
                quantile(metrics$AUC_Class_3, 0.975)),
        logfile)
  write(sprintf('Accuracy of low risk: %.2f (%.2f, %.2f)',
                median(metrics$Accuracy_Class_1),
                quantile(metrics$Accuracy_Class_1, 0.025),
                quantile(metrics$Accuracy_Class_1, 0.975)),
        logfile)
  write(sprintf('Accuracy of intermediate risk: %.2f (%.2f, %.2f)',
                median(metrics$Accuracy_Class_2),
                quantile(metrics$Accuracy_Class_2, 0.025),
                quantile(metrics$Accuracy_Class_2, 0.975)),
        logfile)
  write(sprintf('Accuracy of high risk: %.2f (%.2f, %.2f)',
                median(metrics$Accuracy_Class_3),
                quantile(metrics$Accuracy_Class_3, 0.025),
                quantile(metrics$Accuracy_Class_3, 0.975)),
        logfile)
  write(sprintf('Sensitivity of low risk: %.2f (%.2f, %.2f)',
                median(metrics$Sensitivity_class_1),
                quantile(metrics$Sensitivity_class_1, 0.025),
                quantile(metrics$Sensitivity_class_1, 0.975)),
        logfile)
  write(sprintf('Sensitivity of intermediate risk: %.2f (%.2f, %.2f)',
                median(metrics$Sensitivity_class_2),
                quantile(metrics$Sensitivity_class_2, 0.025),
                quantile(metrics$Sensitivity_class_2, 0.975)),
        logfile)
  write(sprintf('Sensitivity of high risk: %.2f (%.2f, %.2f)',
                median(metrics$Sensitivity_class_3),
                quantile(metrics$Sensitivity_class_3, 0.025),
                quantile(metrics$Sensitivity_class_3, 0.975)),
        logfile)
  write(sprintf('Specificity of low risk: %.2f (%.2f, %.2f)',
                median(metrics$Specificity_class_1),
                quantile(metrics$Specificity_class_1, 0.025),
                quantile(metrics$Specificity_class_1, 0.975)),
        logfile)
  write(sprintf('Specificity of intermediate risk: %.2f (%.2f, %.2f)',
                median(metrics$Specificity_class_2),
                quantile(metrics$Specificity_class_2, 0.025),
                quantile(metrics$Specificity_class_2, 0.975)),
        logfile)
  write(sprintf('Specificity of high risk: %.2f (%.2f, %.2f)',
                median(metrics$Specificity_class_3),
                quantile(metrics$Specificity_class_3, 0.025),
                quantile(metrics$Specificity_class_3, 0.975)),
        logfile)
  write(sprintf('LR+ of low risk: %.2f (%.2f, %.2f)',
                median(metrics$LR_positive_class_1),
                quantile(metrics$LR_positive_class_1, 0.025),
                quantile(metrics$LR_positive_class_1, 0.975)),
        logfile)
  write(sprintf('LR+ of intermediate risk: %.2f (%.2f, %.2f)',
                median(metrics$LR_positive_class_2),
                quantile(metrics$LR_positive_class_2, 0.025),
                quantile(metrics$LR_positive_class_2, 0.975)),
        logfile)
  write(sprintf('LR+ of high risk: %.2f (%.2f, %.2f)',
                median(metrics$LR_positive_class_3),
                quantile(metrics$LR_positive_class_3, 0.025),
                quantile(metrics$LR_positive_class_3, 0.975)),
        logfile)
  write(sprintf('LR- of low risk: %.2f (%.2f, %.2f)',
                median(metrics$LR_negative_class_1),
                quantile(metrics$LR_negative_class_1, 0.025),
                quantile(metrics$LR_negative_class_1, 0.975)),
        logfile)
  write(sprintf('LR- of intermediate risk: %.2f (%.2f, %.2f)',
                median(metrics$LR_negative_class_2),
                quantile(metrics$LR_negative_class_2, 0.025),
                quantile(metrics$LR_negative_class_2, 0.975)),
        logfile)
  write(sprintf('LR- of high risk: %.2f (%.2f, %.2f)',
                median(metrics$LR_negative_class_3),
                quantile(metrics$LR_negative_class_3, 0.025),
                quantile(metrics$LR_negative_class_3, 0.975)),
        logfile)
  write(sprintf('F1score of low risk: %.2f (%.2f, %.2f)',
                median(metrics$F1score_class_1),
                quantile(metrics$F1score_class_1, 0.025),
                quantile(metrics$F1score_class_1, 0.975)),
        logfile)
  write(sprintf('F1score of intermediate risk: %.2f (%.2f, %.2f)',
                median(metrics$F1score_class_2),
                quantile(metrics$F1score_class_2, 0.025),
                quantile(metrics$F1score_class_2, 0.975)),
        logfile)
  write(sprintf('F1score of high risk: %.2f (%.2f, %.2f)',
                median(metrics$F1score_class_3),
                quantile(metrics$F1score_class_3, 0.025),
                quantile(metrics$F1score_class_3, 0.975)),
        logfile)
  write(sprintf('Classification error: %.2f (%.2f, %.2f)',
                mean(pred_err),
                mean(pred_err) - 1.96 * sd(pred_err) / sqrt(nrow(ordinaldf_test)),
                mean(pred_err) + 1.96 * sd(pred_err) / sqrt(nrow(ordinaldf_test))),
        logfile)
  write(sprintf('Pairwise classification accuracy: %.2f (%.2f, %.2f)',
                median(metrics$Pairwise),
                quantile(metrics$Pairwise, 0.025),
                quantile(metrics$Pairwise, 0.975)),
        logfile)

  write("==============================",logfile)
  write("Metrics for whole testing data",logfile)
  write("==============================",logfile)
  write(sprintf('AUC of low risk: %.2f',metrics_whole_data$AUC_Class_1),logfile)
  write(sprintf('AUC of intermediate risk: %.2f',metrics_whole_data$AUC_Class_2),logfile)
  write(sprintf('AUC of high risk: %.2f',metrics_whole_data$AUC_Class_3),logfile)
  write(sprintf('Accuracy of low risk: %.2f',metrics_whole_data$Accuracy_Class_1),logfile)
  write(sprintf('Accuracy of intermediate risk: %.2f',metrics_whole_data$Accuracy_Class_2),logfile)
  write(sprintf('Accuracy of high risk: %.2f',metrics_whole_data$Accuracy_Class_3),logfile)
  write(sprintf('Sensitivity of low risk: %.2f',metrics_whole_data$Sensitivity_class_1),logfile)
  write(sprintf('Sensitivity of intermediate risk: %.2f',metrics_whole_data$Sensitivity_class_2),logfile)
  write(sprintf('Sensitivity of high risk: %.2f',metrics_whole_data$Sensitivity_class_3),logfile)
  write(sprintf('Specificity of low risk: %.2f',metrics_whole_data$Specificity_class_1),logfile)
  write(sprintf('Specificity of intermediate risk: %.2f',metrics_whole_data$Specificity_class_2),logfile)
  write(sprintf('Specificity of high risk: %.2f',metrics_whole_data$Specificity_class_3),logfile)
  write(sprintf('LR+ of low risk: %.2f',metrics_whole_data$LR_positive_class_1),logfile)
  write(sprintf('LR+ of intermediate risk: %.2f',metrics_whole_data$LR_positive_class_2),logfile)
  write(sprintf('LR+ of high risk: %.2f',metrics_whole_data$LR_positive_class_3),logfile)
  write(sprintf('LR- of low risk: %.2f',metrics_whole_data$LR_negative_class_1),logfile)
  write(sprintf('LR- of intermediate risk: %.2f',metrics_whole_data$LR_negative_class_2),logfile)
  write(sprintf('LR- of high risk: %.2f',metrics_whole_data$LR_negative_class_3),logfile)
  write(sprintf('F1score of low risk: %.2f',metrics_whole_data$F1score_class_1),logfile)
  write(sprintf('F1score of intermediate risk: %.2f',metrics_whole_data$F1score_class_2),logfile)
  write(sprintf('F1score of high risk: %.2f',metrics_whole_data$F1score_class_3),logfile)
  write(sprintf('Pairwise classification accuracy: %.2f',metrics_whole_data$Pairwise),logfile)

  # Close the connection to the text file
  close(logfile) 
  
  # Store summary of metrics to summarydf
  summarydf <- data.frame(
    Feature_Pair = paste(feature_model_1,feature_model_2,sep = "_"),
    Alpha_1 = alpha1,
    Alpha_2 = alpha2,
    Beta_1 = beta1,
    Beta_2 = beta2,
    Accuracy_Class_1 = median(metrics$Accuracy_Class_1),
    Accuracy_Class_2 = median(metrics$Accuracy_Class_2),
    Accuracy_Class_3 = median(metrics$Accuracy_Class_3),
    AUC_Class_1 = median(metrics$AUC_Class_1),
    AUC_Class_2 = median(metrics$AUC_Class_2),
    AUC_Class_3 = median(metrics$AUC_Class_3),
    Sensitivity_class_1 = median(metrics$Sensitivity_class_1),
    Sensitivity_class_2 = median(metrics$Sensitivity_class_2),
    Sensitivity_class_3 = median(metrics$Sensitivity_class_3),
    Specificity_class_1 = median(metrics$Specificity_class_1),
    Specificity_class_2 = median(metrics$Specificity_class_2),
    Specificity_class_3 = median(metrics$Specificity_class_3),
    LR_positive_class_1 = median(metrics$LR_positive_class_1),
    LR_positive_class_2 = median(metrics$LR_positive_class_2),
    LR_positive_class_3 = median(metrics$LR_positive_class_3),
    LR_negative_class_1 = median(metrics$LR_negative_class_1),
    LR_negative_class_2 = median(metrics$LR_negative_class_2),
    LR_negative_class_3 = median(metrics$LR_negative_class_3),
    F1score_class_1 = median(metrics$F1score_class_1),
    F1score_class_2 = median(metrics$F1score_class_2),
    F1score_class_3 = median(metrics$F1score_class_3),
    Classification_error = mean(pred_err),
    Pairwise_classification_accuracy = median(metrics$Pairwise)
  )
  
  return(summarydf)
}

TMS <- function(alpha1,alpha2,beta1,beta2,feature1,feature2){
  z1 <- alpha1 - beta1 * feature1 - beta2 * feature2
  z2 <- alpha2 - beta1 * feature1 - beta2 * feature2
  
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

metricsfun <- function(training_data, test_data, mod, label_values, tms_name = "TMS", is_whole){
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
  tn_class_2 = confusion_matrix[1,2] + confusion_matrix[1,3] + confusion_matrix[3,1] + confusion_matrix[3,3]
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
  
  # Fill the metrics row by row
  metrics <- data.frame(
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
  
  return(metrics)
}

tmsplotfun <- function(data, th1, th2, label_colors, title, file_name, tms_name){
  data$drug_name <- factor(data$drug_name, levels = unique(data$drug_name[order(data$label)]))
  plot <- ggplot(data, aes_string(x = tms_name, y = "drug_name", fill = "label")) +
    geom_boxplot(color = "black", width = 0.5) +
    labs(title = title, x = tms_name, y = "") +
    geom_vline(xintercept = th1, linetype = "dashed", color = "blue", size = 1)  +
    geom_vline(xintercept = th2, linetype = "dashed", color = "red", size = 1)  +
    scale_fill_manual(values = label_colors)  # Set the fill colors
  ggsave(file_name, plot, width = 8, height = 6, dpi = 300)
}

scatterplotfun <- function(data, mod = NA, feature_model_1, feature_model_2, is_converged, is_training, is_legend){
  # Check the column index
  idx_model_1 <- as.integer(which(colnames(data) == feature_model_1))
  idx_model_2 <- as.integer(which(colnames(data) == feature_model_2))
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
  jpeg(paste(feature_model_1,feature_model_2,file_name,"dataset.jpg",sep = "_"),quality = 100, units = "in", width = 5, height = 5, res = 300)
  plot(data[,idx_model_1], 
       data[,idx_model_2], 
       xlab = feature_model_1, 
       ylab = feature_model_2,
       main = title,
       cex.axis = 1.5, 
       cex.lab = 1.5, 
       cex.main = 1.5)
  if (is_legend) {
    legend("bottomright", legend = c("Low", "Intermediate", "High"), fill = c("green", "blue", "red"))
  }
  points(data[,idx_model_1][data$label == "1"], data[,idx_model_2][data$label == "1"], col = "green")
  points(data[,idx_model_1][data$label == "2"], data[,idx_model_2][data$label == "2"], col = "blue")
  points(data[,idx_model_1][data$label == "3"], data[,idx_model_2][data$label == "3"], col = "red")
  if (is_converged) {
    abline(a = c1, b = m, col = "blue", lty = 2, lwd = 2)  # Replace 'a' with the intercept if needed
    abline(a = c2, b = m, col = "red", lty = 2, lwd = 2)  # Replace 'a' with the intercept if needed
  }
  dev.off()
}

pairsdfinitfun <- function(features, cell_models_1, cell_models_2){
  pairsdf <- data.frame()
  for (cell_model_1 in cell_models_1){
    for (cell_model_2 in cell_models_2){
      if (cell_model_1 == cell_model_2){
        for (i in 1:(length(features) - 1)) {
          for (j in (i + 1):length(features)) {
            feature_1 <- features[i]
            feature_2 <- features[j]
            tempdf <- data.frame(cell_model_1 = cell_model_1,
                                 cell_model_2 = cell_model_2,
                                 feature_1 = feature_1,
                                 feature_2 = feature_2)
            pairsdf <- rbind(pairsdf,tempdf)
          }# feature_1
        } # feature_2
      } else {
        for (feature_1 in features){
          for (feature_2 in features){
            tempdf <- data.frame(cell_model_1 = cell_model_1,
                                 cell_model_2 = cell_model_2,
                                 feature_1 = feature_1,
                                 feature_2 = feature_2)
            pairsdf <- rbind(pairsdf,tempdf)
          } # feature_2
        } # feature_1
      }
    } # cell_model_2
  } # cell_model_1
  
  return(pairsdf)
}

solvecellpair <- function(pair_id,
                          filepath_training_ORD,
                          filepath_testing_ORD,
                          filepath_training_Tomek,
                          filepath_testing_Tomek,
                          feature_1,
                          feature_2,
                          cell_model_1,
                          cell_model_2,
                          num_tests){
  feature_model_1 <- paste(feature_1,cell_model_1,sep = "_")
  feature_model_2 <- paste(feature_2,cell_model_2,sep = "_")
  print(paste(pair_id,feature_model_1,feature_model_2,sep = "_"))
  
  # Prepare the raw training and testing dataset
  if (cell_model_1 == "ORD"){
    training_1 <- read_csv(filepath_training_ORD, show_col_types = FALSE)
    testing_1 <- read_csv(filepath_testing_ORD, show_col_types = FALSE)
  } else if(cell_model_1 == "Tomek"){
    training_1 <- read_csv(filepath_training_Tomek, show_col_types = FALSE)
    testing_1 <- read_csv(filepath_testing_Tomek, show_col_types = FALSE)
  }
  if (cell_model_2 == "ORD"){
    training_2 <- read_csv(filepath_training_ORD, show_col_types = FALSE)
    testing_2 <- read_csv(filepath_testing_ORD, show_col_types = FALSE)
  } else if(cell_model_2 == "Tomek"){
    training_2 <- read_csv(filepath_training_Tomek, show_col_types = FALSE)
    testing_2 <- read_csv(filepath_testing_Tomek, show_col_types = FALSE)
  }
  # Perform all
  resultdf <- run_all(training_1,testing_1,
                    training_2,testing_2,
                    feature_1,feature_2,
                    cell_model_1,cell_model_2,
                    num_tests)
  return(resultdf)
}