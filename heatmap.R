library(pheatmap)
library(RColorBrewer)

# File paths
filepath_training <- "C:/Users/USER-PC/OneDrive - Universitas Airlangga/5. CML/4. Data_Ali_E/5. MLR_ORD_TOMEK/data/5. dynamic_hERG_li/metrics_li_training_avg.csv"
filepath_testing <- "C:/Users/USER-PC/OneDrive - Universitas Airlangga/5. CML/4. Data_Ali_E/5. MLR_ORD_TOMEK/data/5. dynamic_hERG_li/metrics_li_testing_avg.csv"

# Function to load and preprocess data
load_and_process_data <- function(filepath) {
  data <- na.omit(read.csv(filepath))
  data$label <- as.factor(data$label)
  data
}

# Function to normalize data
normalize_data <- function(data, cols_to_normalize) {
  data[cols_to_normalize] <- scale(data[cols_to_normalize])
  data
}

# Load and preprocess data
raw_training <- load_and_process_data(filepath_training)
raw_testing <- load_and_process_data(filepath_testing)

# Features to normalize
features <- c("qNet", "dvdtmax", "vmax", "vrest", "APD50", "APD90", "max_dv", "camax", "carest", "CaTD50", "CaTD90")

# Normalize training and testing data
training <- normalize_data(raw_training, features)
testing <- normalize_data(raw_testing, features)

# Make annotation
row_annotation <- training[, c("drug_name", "label")]
drug_names <- unique(row_annotation$drug_name)
drug_colors <- brewer.pal(min(12, length(drug_names)), "Paired")
drug_annotation <- setNames(drug_colors, drug_names)
annotation_color <- list(
  drug_name = drug_annotation,
  label = c("low" = "green", "intermediate" = "blue", "high" = "red")
)

# Create heatmap
create_heatmap <- function(data, annotation, filename) {
  jpeg(filename, quality = 100, units = "in", width = 5, height = 5, res = 300)
  pheatmap(data[, features],
           main = "Heatmap",
           show_rownames = FALSE, show_colnames = TRUE,
           cluster_rows = TRUE, cluster_cols = TRUE,
           clustering_distance_rows = "euclidean",
           clustering_distance_cols = "euclidean",
           clustering_method = "ward.D2",
           color = brewer.pal(5, "Blues"),
           annotation_row = annotation,
           annotation_colors = annotation_color,
           annotation_names_row = FALSE,
           annotation_names_col = FALSE,
           angle_col = 45,
           border_color = "grey")
  dev.off()
}

# Generate heatmap for training data
create_heatmap(training, row_annotation, "heatmap_training.jpg")


# library(pheatmap)
# library(RColorBrewer)
# 
# features <- c("qNet", "dvdtmax", "vmax", "vrest", "APD50", "APD90", "max_dv", "camax", "carest", "CaTD50", "CaTD90")
# 
# # File paths
# filepath_training <- "C:/Users/USER-PC/OneDrive - Universitas Airlangga/5. CML/4. Data_Ali_E/5. MLR_ORD_TOMEK/data/5. dynamic_hERG_li/metrics_li_training_avg.csv"
# filepath_testing <- "C:/Users/USER-PC/OneDrive - Universitas Airlangga/5. CML/4. Data_Ali_E/5. MLR_ORD_TOMEK/data/5. dynamic_hERG_li/metrics_li_testing_avg.csv"
# 
# # Load and preprocess data
# raw_training <- na.omit(read.csv(filepath_training))
# raw_testing <- na.omit(read.csv(filepath_testing))
# raw_training$label <- as.factor(raw_training$label)
# raw_testing$label <- as.factor(raw_testing$label)
# 
# # Normalize training and testing data
# cols_to_normalize <- features
# mean_vals <- sapply(raw_training[cols_to_normalize], mean)
# training <- raw_training
# testing <- raw_testing
# for (col in cols_to_normalize) {
#   training[[col]] <- (training[[col]] - mean(training[[col]])) / sd(training[[col]])
#   testing[[col]] <- (testing[[col]] - mean(testing[[col]])) / sd(testing[[col]])
# }
# # training[cols_to_normalize] <- sweep(training[cols_to_normalize], 2, mean_vals, "/")
# # testing[cols_to_normalize] <- sweep(testing[cols_to_normalize], 2, mean_vals, "/")
# 
# # Make annotation
# row_annotation <- training[,c("drug_name","label")]
# drug_names <- unique(row_annotation$drug_name)
# drug_colors <- brewer.pal(12,"Paired")
# drug_annotation <- setNames(drug_colors,drug_names)
# annotation_color <- list(
#   drug_name = drug_annotation,
#   label = c("low" = "green",
#             "intermediate" = "blue",
#             "high" = "red")
# )
# 
# jpeg(paste("heatmap_training.jpg"),quality = 100, units = "in", width = 5, height = 5, res = 300)
# pheatmap(training[,features],
#          main = "Training dataset",
#          show_rownames = FALSE, show_colnames = TRUE,
#          cluster_rows = TRUE, cluster_cols = TRUE,
#          clustering_distance_rows = "euclidean",
#          clustering_distance_cols = "euclidean",
#          clustering_method = "ward.D",
#          color = brewer.pal(5,"Blues"),
#          annotation_row = row_annotation,
#          annotation_colors = annotation_color,
#          annotation_names_row = FALSE,
#          annotation_names_col = FALSE,
#          angle_col = 45,
#          border_color = "grey")
# dev.off()
