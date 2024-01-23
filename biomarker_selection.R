library(MASS)

# File paths
filepath_training <- "C:/Users/USER-PC/OneDrive - Universitas Airlangga/5. CML/4. Data_Ali_E/5. MLR_ORD_TOMEK/data/5. dynamic_hERG_li/metrics_li_training_avg.csv"
filepath_testing <- "C:/Users/USER-PC/OneDrive - Universitas Airlangga/5. CML/4. Data_Ali_E/5. MLR_ORD_TOMEK/data/5. dynamic_hERG_li/metrics_li_testing_avg.csv"

# Load and preprocess data
raw_training <- na.omit(read.csv(filepath_training))
raw_testing <- na.omit(read.csv(filepath_testing))
raw_training$label <- as.factor(raw_training$label)
raw_testing$label <- as.factor(raw_testing$label)

# Normalize training and testing data
cols_to_normalize <- c("qNet", "dvdtmax", "vmax", "vrest", "APD50", "APD90", "max_dv", "camax", "carest", "CaTD50", "CaTD90")
mean_vals <- sapply(raw_training[cols_to_normalize], mean)
training <- raw_training
testing <- raw_testing
training[cols_to_normalize] <- sweep(training[cols_to_normalize], 2, mean_vals, "/")
testing[cols_to_normalize] <- sweep(testing[cols_to_normalize], 2, mean_vals, "/")

# Perform LDA on training data
formula_string <- paste("label ~", paste(cols_to_normalize, collapse = " + "))
mod <- lda(as.formula(formula_string), data = training)
print(mod)

# Collecting LD1 and LD2
resultLDA <- rbind(
  data.frame(predict(mod)$x[,1:2], label = training$label, is_training = 1, drug_name = training$drug_name, Sample_ID = training$Sample_ID),
  data.frame(predict(mod, newdata = testing)$x[,1:2], label = testing$label, is_training = 0, drug_name = testing$drug_name, Sample_ID = testing$Sample_ID)
)
write.csv(resultLDA, "resultLDA.csv", row.names = FALSE)

# Perform PCA on normalized training data using princomp
pca_result <- princomp(training[cols_to_normalize], cor = FALSE)
print(pca_result)

# Collect the components
resultPCA <- rbind(
  data.frame(predict(pca_result), label = training$label, is_training = 1, drug_name = training$drug_name, Sample_ID = training$Sample_ID),
  data.frame(predict(pca_result, newdata = testing), label = testing$label, is_training = 0, drug_name = testing$drug_name, Sample_ID = testing$Sample_ID)
)
write.csv(resultPCA, "resultPCA.csv", row.names = FALSE)