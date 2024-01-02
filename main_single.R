library(readr)
library(MASS)
library(pROC)
library(ggplot2)
source("functions.r")

# Input
features <- c("dVm_dt_Repol",
              "dVm_dt_Max",
              "Vm_Peak",	
              "APD90",
              "APD50",
              "APDtri",
              "Vm_Resting",
              "Ca_Peak",
              "Ca_Diastole",
              "CaD90",
              "CaD50",
              "Catri",
              "qNet_Vm_repol_90",
              "qNet_CL",
              "qInward")
filepath_training <- "data_single/AVG_ORD_4ch_training dataset.csv"
filepath_testing <- "data_single/AVG_ORD_4ch_testing dataset.csv"
filepath_drug_list <- "data_single/drug_data_28.csv"
num_tests <- 10000

# Loop over the features
for (feature in features) {
  # Search for drug list 
  drug_list <- read.csv(filepath_drug_list)
  training_drugs <- drug_list["drug"][drug_list["test.training"]=="training"]
  testing_drugs <- drug_list["drug"][drug_list["test.training"]=="test"]
  raw_training <- read.csv(filepath_training)
  raw_testing <- read.csv(filepath_testing)
  trainingdf <- raw_training[raw_training$drug_name %in% training_drugs,]
  testingdf <- raw_testing[raw_testing$drug_name %in% testing_drugs,]
  
  # Log file
  logfile <- file(paste(feature,"log.txt",sep = "_"),open = "wt")
  
  # Add label for each data
  trainingdf["label"] = ""
  for (drug in training_drugs){
    trainingdf["label"][trainingdf["drug_name"] == drug] = drug_list["risk"][drug_list["drug"] == drug]
  }
  testingdf["label"] = ""
  for (drug in testing_drugs){
    testingdf["label"][testingdf["drug_name"] == drug] = drug_list["risk"][drug_list["drug"] == drug]
  }
  
  # Modify the label
  levels <- c("low", "intermediate", "high")
  label_values <- c(1, 2, 3)
  trainingdf$label <- as.integer(factor(trainingdf$label, levels = levels, labels = label_values))
  trainingdf$label <- as.factor(trainingdf$label)
  testingdf$label <- as.integer(factor(testingdf$label, levels = levels, labels = label_values))
  testingdf$label <- as.factor(testingdf$label)
  
  # Reorder the "Drug" factor levels based on the order of "label" levels
  trainingdf$drug_name <- factor(trainingdf$drug_name, levels = unique(trainingdf$drug_name[order(trainingdf$label)]))
  testingdf$drug_name <- factor(testingdf$drug_name, levels = unique(testingdf$drug_name[order(testingdf$label)]))
  
  # Fit the ordinal logistic regression model
  formula_string <- paste("label~", feature, sep = "")
  formula <- as.formula(formula_string)
  mod <- polr(formula, data = trainingdf, Hess = TRUE)
  
  # Make predictions
  predicted_labels <- predict(mod, newdata = trainingdf, type = "class")
  
  # Variables from Ordinal Logistic Regression model
  alpha_1 <- as.numeric(mod$zeta[1])
  alpha_2 <- as.numeric(mod$zeta[2])
  beta <- as.numeric(mod$coefficients[1])
  
  # Thresholds of the feature
  threshold_1 <- -log(exp(-alpha_1) - 2.0 * exp(-alpha_2))/beta
  threshold_2 <- -log(exp(-alpha_1) * exp(-alpha_2) / (exp(-alpha_1) - 2.0 * exp(-alpha_2))) / beta
  
  # Plot the training dataset
  label_colors <- c("green", "blue", "red")
  tmsplotfun(data = trainingdf, 
             th1 = threshold_1, 
             th2 = threshold_2, 
             label_colors = label_colors, 
             title = "Training dataset",
             file_name = paste(feature,"training_dataset.jpg",sep = "_"),
             tms_name = feature)
  
  # Plot the testing dataset
  tmsplotfun(data = testingdf, 
             th1 = threshold_1, 
             th2 = threshold_2, 
             label_colors = label_colors, 
             title = "Testing dataset",
             file_name = paste(feature,"testing_dataset.jpg",sep = "_"),
             tms_name = feature)
  
  # Preallocate the metrics dataframe
  metrics <- data.frame()
  
  # Iterate through the 10,000 tests
  for (i in 1:num_tests) {
    # Create a subset with 16 random drug names
    selected_drugs <- testing_drugs
    sampled_rows <- data.frame()
    for (drug in selected_drugs) {
      drug_data <- testingdf[testingdf$drug_name == drug, ]
      sampled_row <- drug_data[sample(nrow(drug_data), 1), ]
      sampled_rows <- rbind(sampled_rows, sampled_row)
    }
    test_data <- sampled_rows
    test_data["is_training"] <- 0
    
    # Create a subset with 12 training drug names
    selected_drugs <- training_drugs
    sampled_rows <- data.frame()
    for (drug in selected_drugs) {
      drug_data <- trainingdf[trainingdf$drug_name == drug, ]
      sampled_row <- drug_data[sample(nrow(drug_data), 1), ]
      sampled_rows <- rbind(sampled_rows, sampled_row)
    }
    training_data <- sampled_rows
    training_data["is_training"] <- 1 
    
    # Calculate all metrics 
    tempmetrics <- metricsfun(training_data = training_data, 
                              test_data = test_data, 
                              mod = mod, 
                              label_values = label_values, 
                              tms_name = feature,
                              is_whole = FALSE)
    
    # Store the calculated metrics 
    metrics <- rbind(metrics,tempmetrics)
  } # i-th test
  
  # Calculate the classification error from the whole testing dataset
  predicted_labels <- predict(mod, newdata = testingdf, type = "class")
  pred_err <- (abs(as.integer(predicted_labels) - as.integer(testingdf$label)))
  
  # Calculate the metrics for the whole dataset
  trainingdf["is_training"] = as.integer(1)
  testingdf["is_training"] = as.integer(0)
  metrics_whole_data <- metricsfun(training_data = trainingdf, 
                                   test_data = testingdf,
                                   mod = mod,
                                   label_values = label_values,
                                   tms_name = feature,
                                   is_whole = TRUE)
  
  # Save the metrics dataframe
  write.csv(metrics, paste(feature,"metrics.csv",sep = "_"), row.names = FALSE)
  
  # Print the Thresholds into logfile
  write("=======================",logfile)
  write(sprintf("TMS thresholds"),logfile)
  write("=======================",logfile)
  write(paste(sprintf('Threshold_1: %.4f ', threshold_1)), logfile)
  write(paste(sprintf('Threshold_2: %.4f ', threshold_2)), logfile)
  
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
                mean(pred_err) - 1.96 * sd(pred_err) / sqrt(nrow(testingdf)),
                mean(pred_err) + 1.96 * sd(pred_err) / sqrt(nrow(testingdf))),
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
}