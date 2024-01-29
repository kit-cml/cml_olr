library(MASS)

features <- c("FPD",
              "FPDc",
              "Beat_Period",
              "Amplitude")
# features <- c("qNet", "dvdtmax", "vmax", "vrest", "APD50", "APD90", "max_dv", "camax", "carest", "CaTD50", "CaTD90")

# File paths
# filepath_training <- "C:/Users/USER-PC/OneDrive - Universitas Airlangga/5. CML/4. Data_Ali_E/5. MLR_ORD_TOMEK/data/5. dynamic_hERG_li/metrics_li_training_avg.csv"
# filepath_testing <- "C:/Users/USER-PC/OneDrive - Universitas Airlangga/5. CML/4. Data_Ali_E/5. MLR_ORD_TOMEK/data/5. dynamic_hERG_li/metrics_li_testing_avg.csv"
filepath_training <- "C:/Users/USER-PC/OneDrive - Universitas Airlangga/5. CML/2. RESEARCH/26. MEA/data/14_cells_28_CiPA.csv"
filepath_testing <- "C:/Users/USER-PC/OneDrive - Universitas Airlangga/5. CML/2. RESEARCH/26. MEA/data/14_cells_28_CiPA.csv"

# Open the connection to logfile
sink("LDA_and_PCA_log.txt")

# Load and preprocess data
raw_training <- na.omit(read.csv(filepath_training))
raw_testing <- na.omit(read.csv(filepath_testing))
raw_training$label <- as.factor(raw_training$label)
raw_testing$label <- as.factor(raw_testing$label)

print("===========================================================")
print("The correlation matrix for raw training data")
print("===========================================================")
print(cor(raw_training[,features]))

print("===========================================================")
print("The correlation matrix for raw testing data")
print("===========================================================")
print(cor(raw_testing[,features]))

# Normalize training and testing data
cols_to_normalize <- features
mean_vals <- sapply(raw_training[cols_to_normalize], mean)
training <- raw_training
testing <- raw_testing
training[cols_to_normalize] <- sweep(training[cols_to_normalize], 2, mean_vals, "/")
testing[cols_to_normalize] <- sweep(testing[cols_to_normalize], 2, mean_vals, "/")

# Perform LDA on training data
formula_string <- paste("label ~", paste(cols_to_normalize, collapse = " + "))
mod <- lda(as.formula(formula_string), data = training)
print("===========================================================")
print("The LDA model")
print("===========================================================")
print(mod)

# Print the LD1 and LD2
print("===========================================================")
print("The correlation matrix of LD1 and LD2 for training dataset")
print("===========================================================")
pred <- predict(mod)$x[,1:2]
print(cor(pred))
write.csv(data.frame(pred, 
                     label = training$label, 
                     drug_name = training$drug_name, 
                     Sample_ID = training$Sample_ID),
          "resultLDA_training.csv", 
          row.names = FALSE)
print("===========================================================")
print("The correlation matrix of LD1 and LD2 for testing dataset")
print("===========================================================")
pred <- predict(mod, newdata = testing)$x[,1:2]
print(cor(pred))
write.csv(data.frame(pred,
                     label = testing$label, 
                     drug_name = testing$drug_name, 
                     Sample_ID = testing$Sample_ID),
          "resultLDA_testing.csv", 
          row.names = FALSE)

# Perform PCA on normalized training data using princomp
pca_result <- princomp(training[cols_to_normalize], cor = FALSE)
print("===========================================================")
print("The PCA model")
print("===========================================================")
print(pca_result)
print(summary(pca_result))

# Print the components
print("===========================================================")
print("The correlation matrix of PCA components for training dataset")
print("===========================================================")
pred <- predict(pca_result)
print(cor(pred))
write.csv(data.frame(pred, 
                     label = training$label, 
                     drug_name = training$drug_name, 
                     Sample_ID = training$Sample_ID),
          "resultPCA_training.csv",
          row.names = FALSE)
print("===========================================================")
print("The correlation matrix of PCA components for testing dataset")
print("===========================================================")
pred <- predict(pca_result, newdata = testing)
print(cor(pred))
write.csv(data.frame(pred,
                     label = testing$label, 
                     drug_name = testing$drug_name, 
                     Sample_ID = testing$Sample_ID),
          "resultPCA_testing.csv", 
          row.names = FALSE)

# Close the connection to the logfile
sink()