library(readr)

# Initialize the result folder and features
result_folder <- "Accepted_models_manual"
features <- c("qNet", "dvdtmax", "vmax", "vrest", "APD50", "APD90", "max_dv", "camax", "carest", "CaTD50", "CaTD90")

# Initialize a dataframe to store the betas
feature_counts <- data.frame(Feature = features)
feature_betas_avg <- data.frame(Feature = features)
feature_betas_min <- data.frame(Feature = features)
feature_betas_max <- data.frame(Feature = features)
feature_rank_scores_avg <- data.frame(Feature = features)
feature_rank_scores_min <- data.frame(Feature = features)
feature_rank_scores_max <- data.frame(Feature = features)

max_scenarios <- 11

# Calculation for feature occurrences
for (i in 1:max_scenarios) {
  summary_file <- paste(result_folder,"/summary_",i,".csv",sep="")
  summary_df <- read_csv(summary_file, show_col_types = FALSE)
  summary_df <- na.omit(summary_df)
  feature_pairs <- summary_df$Feature_Pair
  col_name <- as.character(i)
  feature_counts[col_name] <- 0
  for (long_string in feature_pairs) {
    for (j in seq_along(features)) {
      feature_counts[[col_name]][j] <- feature_counts[[col_name]][j] + as.numeric(grepl(features[j], long_string))
    }
  }
  feature_counts[col_name] <- feature_counts[col_name] / nrow(summary_df)
}
write.csv(feature_counts,"Biomarker_occurrences.csv",row.names = FALSE)

# Calculation for betas
for (i in 1:max_scenarios) {
  summary_file <- paste(result_folder, "/summary_", i, ".csv", sep="")
  summary_df <- read_csv(summary_file, show_col_types = FALSE)
  summary_df <- na.omit(summary_df)
  feature_betas_avg[paste0(i)] <- 0
  feature_betas_min[paste0(i)] <- 1000
  feature_betas_max[paste0(i)] <- 0
  for (row in 1:nrow(summary_df)) {
    long_string <- summary_df$Feature_Pair[row]
    betas <- abs(unlist(summary_df[row, grep("^Beta", names(summary_df))]))  # Extract beta columns
    betas <- betas / max(betas)
    feature_names_in_pair <- strsplit(long_string, "-")[[1]]
    for (feature_idx in seq_along(feature_names_in_pair)) {
      feature_name <- feature_names_in_pair[feature_idx]
      feature_row <- which(feature_betas_avg$Feature == feature_name)
      if (length(feature_row) > 0 && !is.na(betas[feature_idx])) {
        feature_betas_avg[feature_row, paste0(i)] <- feature_betas_avg[feature_row, paste0(i)] + as.numeric(betas[feature_idx])
        feature_betas_min[feature_row, paste0(i)] <- min(feature_betas_min[feature_row, paste0(i)],as.numeric(betas[feature_idx]))
        feature_betas_max[feature_row, paste0(i)] <- max(feature_betas_max[feature_row, paste0(i)],as.numeric(betas[feature_idx]))
      }
    }
  }
  for (row in 1:nrow(feature_betas_avg)){
    if (feature_counts[row,paste0(i)] > 0.0) {
      feature_betas_avg[row,paste0(i)] <- feature_betas_avg[row,paste0(i)]/(nrow(summary_df) * feature_counts[row,paste0(i)])
    }
    if (feature_betas_min[row,paste0(i)] == 1000) {
      feature_betas_min[row,paste0(i)] <- 0
    }
  }
}
write.csv(feature_betas_avg, "Feature_Betas_avg.csv", row.names = FALSE)
write.csv(feature_betas_min, "Feature_Betas_min.csv", row.names = FALSE)
write.csv(feature_betas_max, "Feature_Betas_max.csv", row.names = FALSE)

# Calculation of rank scores
for (i in 1:max_scenarios) {
  summary_file <- paste(result_folder,"/summary_",i,".csv",sep="")
  summary_df <- read_csv(summary_file, show_col_types = FALSE)
  summary_df <- na.omit(summary_df)
  feature_rank_scores_avg[paste0(i)] <- 0
  feature_rank_scores_min[paste0(i)] <- 1000
  feature_rank_scores_max[paste0(i)] <- 0
  for (row in 1:nrow(summary_df)) {
    long_string <- summary_df$Feature_Pair[row]
    rank_score <- abs(unlist(summary_df[row, grep("^Rank_score", names(summary_df))]))  # Extract rank score column
    feature_names_in_pair <- strsplit(long_string, "-")[[1]]
    for (feature_idx in seq_along(feature_names_in_pair)) {
      feature_name <- feature_names_in_pair[feature_idx]
      feature_row <- which(feature_rank_scores_avg$Feature == feature_name)
      if (length(feature_row) > 0 && !is.na(rank_score)) {
        feature_rank_scores_avg[feature_row, paste0(i)] <- feature_rank_scores_avg[feature_row, paste0(i)] + as.numeric(rank_score)
        feature_rank_scores_min[feature_row, paste0(i)] <- min(feature_rank_scores_min[feature_row, paste0(i)],as.numeric(rank_score))
        feature_rank_scores_max[feature_row, paste0(i)] <- max(feature_rank_scores_max[feature_row, paste0(i)],as.numeric(rank_score))
      }
    }
  }
  for (row in 1:nrow(feature_rank_scores_avg)){
    if (feature_counts[row,paste0(i)] > 0.0) {
      feature_rank_scores_avg[row,paste0(i)] <- feature_rank_scores_avg[row,paste0(i)]/(nrow(summary_df) * feature_counts[row,paste0(i)])
    }
    if (feature_rank_scores_min[row,paste0(i)] == 1000) {
      feature_rank_scores_min[row,paste0(i)] <- 0
    }
  }
}
write.csv(feature_rank_scores_avg,"Biomarker_rank_scores_avg.csv",row.names = FALSE)
write.csv(feature_rank_scores_min,"Biomarker_rank_scores_min.csv",row.names = FALSE)
write.csv(feature_rank_scores_max,"Biomarker_rank_scores_max.csv",row.names = FALSE)