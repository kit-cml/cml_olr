library(MASS)

filepath_training <- "C:/Users/USER-PC/OneDrive - Universitas Airlangga/5. CML/4. Data_Ali_E/5. MLR_ORD_TOMEK/data/5. dynamic_hERG_li/metrics_li_training_avg.csv"
filepath_testing <- "C:/Users/USER-PC/OneDrive - Universitas Airlangga/5. CML/4. Data_Ali_E/5. MLR_ORD_TOMEK/data/5. dynamic_hERG_li/metrics_li_testing_avg.csv"

# Load training and testing data
raw_training <- read.csv(filepath_training)
raw_testing <- read.csv(filepath_testing)
raw_training <- na.omit(raw_training)
raw_testing <- na.omit(raw_testing)
raw_training$label <- as.factor(raw_training$label)
raw_testing$label <- as.factor(raw_testing$label)

# Normalize training and testing data
cols_to_normalize <- c("qNet",
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
training <- raw_training
testing <- raw_testing
for (col in cols_to_normalize) {
  training[[col]] <- training[[col]] / mean(raw_training[[col]])
  testing[[col]] <- testing[[col]] / mean(raw_training[[col]])
}

# Perform LDA on training data
formula_string <- "label ~"
for (col in cols_to_normalize){
  formula_string <- paste(formula_string,"+",col,sep = " ")
}
formula <- as.formula(formula_string)
# mod <- lda(label ~ qNet + dvdtmax + vmax + vrest + APD50 + APD90 + max_dv + camax + carest + CaTD50 + CaTD90,
mod <- lda(formula,
          data = training)

# Calculate training accuracy
pred <- predict(mod)$class
print(mean(pred == training$label))

# Confusion matrix for training data
table(pred,training$label)

# Plot for training data
pred <- predict(mod)
plot(pred$x[,1:2], 
     col = training$label,
     main = "Training data")
legend("bottomright", legend = c("Low", "Intermediate", "High"), fill = c("green", "blue", "red"))
points(pred$x[,1][training$label == "low"], pred$x[,2][training$label == "low"], col = "green")
points(pred$x[,1][training$label == "intermediate"], pred$x[,2][training$label == "intermediate"], col = "blue")
points(pred$x[,1][training$label == "high"], pred$x[,2][training$label == "high"], col = "red")

# Calculate testing accuracy
pred <- predict(mod, newdata = testing)$class
print(mean(pred == testing$label))

# Confusion matrix for testing data
table(pred,testing$label)

# Plot for testing data
pred <- predict(mod, newdata = testing)
plot(pred$x[,1:2], col = testing$label,
     main = "Testing data")
legend("bottomright", legend = c("Low", "Intermediate", "High"), fill = c("green", "blue", "red"))
points(pred$x[,1][testing$label == "low"], pred$x[,2][testing$label == "low"], col = "green")
points(pred$x[,1][testing$label == "intermediate"], pred$x[,2][testing$label == "intermediate"], col = "blue")
points(pred$x[,1][testing$label == "high"], pred$x[,2][testing$label == "high"], col = "red")
