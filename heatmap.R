library(pheatmap)
library(RColorBrewer)
library(dplyr)

# File paths
filepath_training <- "C:/Users/USER-PC/OneDrive - Universitas Airlangga/5. CML/4. Data_Ali_E/5. MLR_ORD_TOMEK/data/5. dynamic_hERG_li/metrics_li_training_avg.csv"
filepath_testing <- "C:/Users/USER-PC/OneDrive - Universitas Airlangga/5. CML/4. Data_Ali_E/5. MLR_ORD_TOMEK/data/5. dynamic_hERG_li/metrics_li_testing_avg.csv"

# Load and preprocess data function
load_and_process_data <- function(filepath) {
  data <- na.omit(read.csv(filepath))
  data$label <- as.factor(data$label)
  data
}

# Normalize data function
normalize_data <- function(data, cols_to_normalize) {
  data[cols_to_normalize] <- scale(data[cols_to_normalize])
  data
}

# Load and preprocess data
raw_training <- load_and_process_data(filepath_training)
raw_testing <- load_and_process_data(filepath_testing)
combined <- rbind(raw_training, raw_testing)

# Features to normalize
# features <- c("qNet", "dvdtmax", "vmax", "vrest", "APD50", "APD90", "max_dv", "camax", "carest", "CaTD50", "CaTD90")
features <- c("carest","CaTD90")

# Normalize data
training <- normalize_data(raw_training, features)
testing <- normalize_data(raw_testing, features)
combined_normalized <- normalize_data(combined, features)

# Sort the data according to drug's label
training <- training %>% arrange(label)
testing <- testing %>% arrange(label)
combined_normalized < combined_normalized %>% arrange(label)

# Sample 100 samples from combined dataframes
sampled_combined <- subset(combined_normalized,combined_normalized$Sample_ID <= 100)
sampled_training <- subset(training,training$Sample_ID <= 500)
sampled_testing <- subset(testing,testing$Sample_ID <= 500)

# Function to create a custom color palette for drugs
create_custom_color_palette <- function(num_colors) {
  if (num_colors <= 12) {
    return(brewer.pal(num_colors, "Paired"))
  } else {
    # Combine colors from different Brewer palettes and/or custom colors
    additional_colors <- c("#e6194b", "#3cb44b", "#ffe119", "#4363d8", "#f58231", "#911eb4", "#46f0f0", "#f032e6", "#bcf60c", "#fabebe", "#008080", "#e6beff", "#9a6324", "#fffac8", "#800000", "#aaffc3", "#808000", "#ffd8b1", "#000075", "#808080", "#ffffff", "#000000")
    return(c(brewer.pal(12, "Paired"), additional_colors[1:(num_colors - 12)]))
  }
}

# Create heatmap function
create_heatmap <- function(data, annotation, annotation_colors, filename, map_title) {
  jpeg(filename, quality = 100, units = "in", width = 5, height = 7, res = 300)
  blues_gradient <- colorRampPalette(brewer.pal(9, "Blues"))(100) # Adjust the number of colors as needed
  pheatmap(data[, features],
           main = map_title,
           show_rownames = FALSE, show_colnames = TRUE,
           cluster_rows = FALSE, cluster_cols = FALSE,
           # clustering_distance_rows = "euclidean",
           # clustering_distance_cols = "euclidean",
           # clustering_method = "complete",
           color = blues_gradient,
           annotation_row = annotation,
           annotation_colors = annotation_colors,
           annotation_names_row = FALSE,
           annotation_names_col = FALSE,
           annotation_legend_spacing = 0.5,
           fontsize = 8,
           angle_col = 45,
           border_color = "grey")
  dev.off()
}

# Make annotation for training data
row_annotation <- sampled_training[, c("drug_name", "label")]
drug_names <- unique(row_annotation$drug_name)
# drug_colors <- rainbow(length(drug_names))
drug_colors <- create_custom_color_palette(12)
drug_annotation <- setNames(drug_colors, drug_names)
annotation_color <- list(
  drug_name = drug_annotation,
  label = c("low" = "green", "intermediate" = "blue", "high" = "red")
)

# Generate heatmap for training data
create_heatmap(sampled_training,
               row_annotation,
               annotation_color,
               "heatmap_training.jpg",
               "Training dataset")

# Make annotation for testing data
row_annotation <- sampled_testing[, c("drug_name", "label")]
drug_names <- unique(row_annotation$drug_name)
# drug_colors <- rainbow(length(drug_names))
drug_colors <- create_custom_color_palette(16)
drug_annotation <- setNames(drug_colors, drug_names)
annotation_color <- list(
  drug_name = drug_annotation,
  label = c("low" = "green", "intermediate" = "blue", "high" = "red")
)

# Generate heatmap for testing data
create_heatmap(sampled_testing,
               row_annotation,
               annotation_color,
               "heatmap_testing.jpg",
               "Testing dataset")

# # Make annotation for combined data
# row_annotation <- sampled_combined[, c("drug_name", "label")]
# drug_names <- unique(row_annotation$drug_name)
# drug_colors <- rainbow(length(drug_names))
# drug_annotation <- setNames(drug_colors, drug_names)
# annotation_color <- list(
#   drug_name = drug_annotation,
#   label = c("low" = "green", "intermediate" = "blue", "high" = "red")
# )

# # Generate heatmap for combined data
# create_heatmap(sampled_combined,
#                row_annotation,
#                annotation_color,
#                "heatmap_combined.jpg",
#                "Combined dataset")