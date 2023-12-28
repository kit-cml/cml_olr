library(readr)
library(MASS)
library(pROC)
library(ggplot2)
library(foreach)
library(doParallel)
source("functions.r")

# Declare features
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

# Declare cell models 
cell_models_1 <- c("ORD")
cell_models_2 <- c("ORD")

# Declare the file paths
filepath_training_ORD <- "../data/4. dynamic_hERG_chantest/metrics_chantest_training_4_cmax.csv"
filepath_testing_ORD <- "../data/4. dynamic_hERG_chantest/metrics_chantest_testing_4_cmax.csv"
filepath_training_Tomek <- "../data/0. Tomek/Tomek_train_rev.csv"
filepath_testing_Tomek <- "../data/0. Tomek/Tomek_test_rev.csv"

# Set the number of tests
num_tests <- 10

# Create pairsdf with all unique combinations
pairsdf <- pairsdfinitfun(features = features, 
                          cell_models_1 = cell_models_1, 
                          cell_models_2 = cell_models_2)

# Register parallel backend
numCores <- 3
cl <- makeCluster(numCores)
registerDoParallel(cl)

# Execute the tasks in parallel
summarydf <- foreach(pair_id = 1:nrow(pairsdf),
                     .combine = 'rbind',
                     .packages = c("readr","MASS","pROC","ggplot2")) %dopar% {
  pair <- pairsdf[pair_id,]
  result <- solvecellpair(pair_id = pair_id,
                          filepath_training_ORD = filepath_training_ORD,
                          filepath_testing_ORD = filepath_testing_ORD,
                          filepath_training_Tomek = filepath_training_Tomek,
                          filepath_testing_Tomek = filepath_testing_Tomek,
                          feature_1 = pair$feature_1,
                          feature_2 = pair$feature_2,
                          cell_model_1 = pair$cell_model_1,
                          cell_model_2 = pair$cell_model_2,
                          num_tests = num_tests)
  # Return the result
  result
}

# Stop the parallel cluster
stopCluster(cl)

# Save the summarydf dataframe
write.csv(summarydf, "summary.csv", row.names = FALSE)
