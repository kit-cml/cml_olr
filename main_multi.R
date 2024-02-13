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
features <- c("carest",
# features <- c("qNet",
#               "dvdtmax")
              # "vmax",
              # "vrest",
              # "APD50",
              # "APD90",
              # "max_dv",
              # "camax",
              # "carest")
              # "CaTD50",
              "CaTD90")
# units <- c("","")
# units <- c("","")
units <- c("(mM)",
# units <- c("(nC/uF)",
#            "(mV/ms)")
           # "(mV)",
           # "(mV)",
           # "(ms)",
           # "(ms)",
           # "(mV/ms)",
           # "(mM)",
           # "(mM)")
           # "(ms)",
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
num_tests <- 1

# Create pairsdf with all unique combinations
pairsdf <- pairsdfinitfun(features = features, units = units)

# Register parallel backend
numCores <- 1
cl <- makeCluster(numCores)
registerDoParallel(cl)

# Execute the tasks in parallel
summarydf <- foreach(pair_id = 1:nrow(pairsdf),
                     .combine = 'rbind',
                     .packages = c("readr","MASS","pROC","ggplot2")) %dopar% {
  pair <- pairsdf[pair_id,]
  result <- solvepair(pair_id = pair_id,
                          filepath_training = filepath_training,
                          filepath_testing = filepath_testing,
                          feature_1 = pair$feature_1,
                          feature_2 = pair$feature_2,
                          unit_1 = pair$unit_1,
                          unit_2 = pair$unit_2,
                          num_tests = num_tests)
  # Return the result
  result
}

# Stop the parallel cluster
stopCluster(cl)

# Save the summarydf dataframe
write.csv(summarydf, "summary.csv", row.names = FALSE)
