library(readr)
library(MASS)
library(pROC)
library(ggplot2)
library(foreach)
library(doParallel)
source("functions.r")

# Input
features <- c("qNet")
              # "CaTD90")
# features <- c("qNet",
#               "dvdtmax",
#               "vmax",
#               "vrest",
#               "APD50",
#               "APD90",
#               "max_dv",
#               "camax",
#               "carest",
#               "CaTD50",
#               "CaTD90")
units <- c("(nC/uF)")
# units <- c("(nC/uF)",
#            "(mV/ms)",
#            "(mV)",
#            "(mV)",
#            "(ms)",
#            "(ms)",
#            "(mV/ms)",
#            "(mM)",
#            "(mM)",
#            "(ms)",
#            "(ms)")

filepath_training <- "C:/Users/USER-PC/OneDrive - Universitas Airlangga/5. CML/4. Data_Ali_E/5. MLR_ORD_TOMEK/data/4. dynamic_hERG_chantest/metrics_chantest_training_avg.csv"
filepath_testing <- "C:/Users/USER-PC/OneDrive - Universitas Airlangga/5. CML/4. Data_Ali_E/5. MLR_ORD_TOMEK/data/4. dynamic_hERG_chantest/metrics_chantest_testing_avg.csv"

num_tests <- 1

# Load training and testing data
training <- data.frame(read.csv(filepath_training))
testing <- data.frame(read.csv(filepath_testing))

# Register parallel backend
numCores <- 1
cl <- makeCluster(numCores)
registerDoParallel(cl)

# Loop over the features
summarydf <- foreach(feature_id = 1:length(features),
                     .combine = 'rbind',
                     .packages = c("readr","MASS","pROC","ggplot2")) %dopar% {
                      result <- run_all(training = training,
                                        testing = testing,
                                        feature_1 = features[feature_id],
                                        unit_1 = units[feature_id],
                                        num_tests = num_tests,
                                        is_single = TRUE)

                      # Return the result
                      result
                     }

# Stop the parallel cluster
stopCluster(cl)

# Save the summarydf dataframe
write.csv(summarydf, "summary.csv", row.names = FALSE)