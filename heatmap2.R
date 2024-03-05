library(readr)
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
result_folder <- "C:/Users/USER-PC/OneDrive - Universitas Airlangga/5. CML/4. Data_Ali_E/5. MLR_ORD_TOMEK/0. Rscript_Multi_OLR/Accepted_models"
feature_counts <- data.frame(Feature = features)
for (i in 1:9) {
  summary_file <- paste(result_folder,"/summary_",i,".csv",sep="")
  summary_df <- read_csv(summary_file, show_col_types = FALSE)
  summary_df <- na.omit(summary_df)
  feature_pairs <- summary_df$Feature_Pair
  col_name <- paste0("Count_", i)
  feature_counts[col_name] <- 0
  for (long_string in feature_pairs) {
    for (j in seq_along(features)) {
      feature_counts[[col_name]][j] <- feature_counts[[col_name]][j] + as.numeric(grepl(features[j], long_string))
    }
  }
  feature_counts[col_name] <- feature_counts[col_name] / max(feature_counts[col_name])
}
