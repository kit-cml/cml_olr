library(readr)
library(MASS)
library(pROC)
library(ggplot2)
library(foreach)
library(doParallel)
source("functions.r")

# Declare features
# features <- c("Comp.1",
#              "Comp.7")
# features <- c("LD1","LD2")
# features <- c("carest",
features <- c("qNet",
              "dvdtmax",
              "vmax",
              "vrest",
              "APD50",
              "APD90",
              "max_dv",
              "camax",
              "carest",
              "CaTD50",
              "CaTD90")
# units <- c("","")
# units <- c("","")
# units <- c("(mM)",
units <- c("(nC/uF)",
           "(mV/ms)",
           "(mV)",
          "(mV)",
          "(ms)",
          "(ms)",
          "(mV/ms)",
          "(mM)",
          "(mM)",
          "(ms)",
           "(ms)")

# Declare the file paths
# filepath_training <- "C:/Users/USER-PC/OneDrive - Universitas Airlangga/5. CML/4. Data_Ali_E/5. MLR_ORD_TOMEK/0. Rscript_Multi_OLR/Results_LDA_PCA/resultPCA_training.csv"
# filepath_testing <- "C:/Users/USER-PC/OneDrive - Universitas Airlangga/5. CML/4. Data_Ali_E/5. MLR_ORD_TOMEK/0. Rscript_Multi_OLR/Results_LDA_PCA/resultPCA_testing.csv"
# filepath_training <- "C:/Users/USER-PC/OneDrive - Universitas Airlangga/5. CML/4. Data_Ali_E/5. MLR_ORD_TOMEK/0. Rscript_Multi_OLR/Results_LDA_PCA/resultLDA_training.csv"
# filepath_testing <- "C:/Users/USER-PC/OneDrive - Universitas Airlangga/5. CML/4. Data_Ali_E/5. MLR_ORD_TOMEK/0. Rscript_Multi_OLR/Results_LDA_PCA/resultLDA_testing.csv"
# filepath_training <- "C:/Users/USER-PC/OneDrive - Universitas Airlangga/5. CML/4. Data_Ali_E/5. MLR_ORD_TOMEK/data/4. dynamic_hERG_chantest/metrics_chantest_training_avg.csv"
# filepath_testing <- "C:/Users/USER-PC/OneDrive - Universitas Airlangga/5. CML/4. Data_Ali_E/5. MLR_ORD_TOMEK/data/4. dynamic_hERG_chantest/metrics_chantest_testing_avg.csv"
filepath_training <- "C:/Users/USER-PC/OneDrive - Universitas Airlangga/5. CML/4. Data_Ali_E/5. MLR_ORD_TOMEK/data/5. dynamic_hERG_li/metrics_li_training_avg.csv"
filepath_testing <- "C:/Users/USER-PC/OneDrive - Universitas Airlangga/5. CML/4. Data_Ali_E/5. MLR_ORD_TOMEK/data/5. dynamic_hERG_li/metrics_li_testing_avg.csv"

# Set the number of tests
num_tests <- 10000

# Set feature dimension
dimension <- 11

# Create pairsdf with all unique combinations
pairsdf <- pairsdfinitfun(features = features, units = units, dimension = dimension)

# Create the results folder
results_folder <- "results"

# Choose whether data needs to be normalized
is_normalized <- TRUE

# Check if the folder exists
if (!dir.exists(results_folder)) {
  # The folder does not exist, so create it
  dir.create(results_folder)
  cat("Folder created:", results_folder, "\n")
} else {
  # The folder already exists
  # List all files in the folder
  files <- list.files(path = results_folder, full.names = TRUE)
  # Remove all files in the folder
  if (length(files) > 0) {
    file.remove(files)
    cat("All files in the folder have been removed.\n")
  } else {
    cat("The folder is already empty.\n")
  }
}

# Register parallel backend
numCores <- 1
cl <- makeCluster(numCores)
registerDoParallel(cl)

# Execute the tasks in parallel
summarydf <- foreach(pair_id = 1:nrow(pairsdf),
                     .combine = 'rbind',
                     .packages = c("readr","MASS","pROC","ggplot2")) %dopar% {
  # Select the row corresponding to pair_id
  pair_row <- pairsdf[pair_id, ]
 
  # Extract vectors of features and units
  features_vector <- pair_row[grep("feature_", names(pair_row))]
  units_vector <- pair_row[grep("unit_", names(pair_row))]
  
  result <- solvepair(results_folder = results_folder,
                      pair_id = pair_id,
                      filepath_training = filepath_training,
                      filepath_testing = filepath_testing,
                      features_vector = features_vector,
                      units_vector = units_vector,
                      num_tests = num_tests,
                      is_normalized = is_normalized)

  # Return the result
  result
}

# Stop the parallel cluster
stopCluster(cl)

# Save the summarydf dataframe
write.csv(summarydf, paste(results_folder,"/","summary.csv", sep = ""), row.names = FALSE)
