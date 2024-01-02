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
num_tests <- 10
summarydf <- data.frame()

# Loop over the features
for (feature in features) {
  # Load training and testing data
  raw_training <- read.csv(filepath_training)
  raw_testing <- read.csv(filepath_testing)
  tempdf <- run_all(training_1 = raw_training,
                    testing_1 = raw_testing,
                    feature_1 = feature,
                    num_tests = num_tests,
                    is_single = TRUE)
  summarydf <- rbind(summarydf,tempdf)
}

# Save the summarydf dataframe
write.csv(summarydf, "summary.csv", row.names = FALSE)