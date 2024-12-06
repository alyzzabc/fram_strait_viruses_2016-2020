library(vegan)
library(tidyr)
library(zoo)

# Load dataframes
caudo_df <- read.table("/path/to/normalized/abundances/caudoviricetes",sep='\t', header=T)
mega_df <- read.table("/path/to/normalized/abundances/megaviricetes",sep='\t', header=T)
others_df <- read.table("/path/to/normalized/abundances/others",sep='\t', header=T)

dates <- read.table("/Users/acalayag/Downloads/metadata.date.txt",check.names=F,header=T)

# Filter dates for the "WSC" location
dates_wsc <- dates[dates$location == "WSC", ]

# Ensure dates are in Date format
dates_wsc$date <- as.Date(dates_wsc$date)

# Function to preprocess and get sorted and normalized matrices
get_sorted_and_normalized_matrix <- function(df, dates) {
  # Ensure 'CVPG' is in the expected column
  if (!"CVPG" %in% colnames(df)) {
    stop("CVPG column not found in dataframe")
  }
  
  b <- df[df$Breadth_coverage > 0.25,]
  merged_df <- merge(b, dates, by.x="Library_name", by.y="short")
  
  # Debug: Check the structure of merged_df
  print(head(merged_df))
  
  # Assuming CVPG is indeed in the correct column
  tryCatch({
    merged_df_matrix <- spread(merged_df[,c(1,3,8)], Library_name, CV_per_genome) #edit this depending on CVPG or CV_per_genome needs
  }, error = function(e) {
    stop("Error in spreading merged_df: ", e$message)
  })
  
  merged_df_matrix[is.na(merged_df_matrix)] <- 0
  rownames(merged_df_matrix) <- merged_df_matrix[,1]
  merged_df_matrix_2 <- merged_df_matrix[,-1]
  short_ids <- dates$short[dates$location == "WSC"]
  wsc <- merged_df_matrix_2[,colnames(merged_df_matrix_2) %in% short_ids]
  wsc_dates <- dates[dates$short %in% colnames(wsc),]
  wsc_dates_ordered <- wsc_dates[order(wsc_dates$date),]
  wsc_sorted_by_date <- wsc[,order(match(colnames(wsc), wsc_dates_ordered$short))]
  
  # Normalize the sorted matrix
  wsc_normalized <- t(t(wsc_sorted_by_date) / rowSums(t(wsc_sorted_by_date)))
  
  # Return both sorted (non-normalized) and normalized matrices
  list(sorted_matrix=wsc_sorted_by_date, normalized_matrix=wsc_normalized, dates_ordered=wsc_dates_ordered)
}


# Processing for caudoviricetes
caudo_result <- get_sorted_and_normalized_matrix(caudo_df, dates)
# Processing for megaviricetes
mega_result <- get_sorted_and_normalized_matrix(mega_df, dates)
# Processing for others
others_result <- get_sorted_and_normalized_matrix(others_df, dates)


#Plot using ggplot2
library(ggplot2)
library(gridExtra)

plot_with_averages_ggplot <- function(result, title) {
  # Calculate distances and matrix dissimilarities
  dates_dist <- as.vector(vegdist(as.Date(result$dates_ordered$date), method='euclidean'))
  matrix_dist <- 1 - as.vector(vegdist(t(result$normalized_matrix), method='bray'))
  
  # Combine and order data by dates_dist
  data_comb <- data.frame(dates_dist, matrix_dist)
  data_comb <- data_comb[order(data_comb$dates_dist), ]
  
  # Define intervals for running averages
  breaks <- seq(min(data_comb$dates_dist, na.rm = TRUE), max(data_comb$dates_dist, na.rm = TRUE), by = 30)
  data_comb$intervals <- cut(data_comb$dates_dist, breaks = breaks, include.lowest = TRUE, right = FALSE)
  
  # Aggregate data by intervals
  means_data <- aggregate(cbind(dates_dist, matrix_dist) ~ intervals, data = data_comb, FUN = mean)
  
  # Define x-axis ticks every 6 months (180 days)
  x_breaks <- seq(0, max(data_comb$dates_dist, na.rm = TRUE), by = 180)
  x_labels <- seq(0, length(x_breaks) - 1) * 6  # Generate labels for each 6-month interval starting from 0 months
  
  # Create the base ggplot
  p <- ggplot(data_comb, aes(x = dates_dist, y = matrix_dist)) +
    geom_point(color = "grey", size = 2, alpha = 0.6) +
    geom_point(data = means_data, aes(x = dates_dist, y = matrix_dist), color = "red", size = 2, shape = 19) +
    ggtitle(title) +
    xlab("Months Separating Samples") +
    ylab("Bray-Curtis Dissimilarity") +
    theme_bw() +
    scale_x_continuous(breaks = x_breaks, labels = x_labels)  # Custom x-ticks and labels every 6 months starting from 0
  
  # Print the plot
  print(p)
}


plot1 <- plot_with_averages_ggplot(caudo_result, "Caudoviricetes")
plot2 <- plot_with_averages_ggplot(mega_result, "Megaviricetes")
plot3 <- plot_with_averages_ggplot(others_result, "Others")

# Arrange plots
grid.arrange(plot1, plot2, plot3, nrow = 3)

# Plot with average, cutoff of >0.1
plot_with_averages_ggplot_2 <- function(result, title) {
  # Load required library
  library(vegan)
  library(ggplot2)
  
  # Step 1: Filter the matrix and the dates based on column sums
  column_sums <- colSums(result$sorted_matrix)
  filtered_matrix <- result$sorted_matrix[, column_sums >= 0.01]
  print("Filtered matrix:")
  print(filtered_matrix) # Add this line to inspect the filtered matrix
  filtered_normalized_matrix <- result$normalized_matrix[, colnames(filtered_matrix)]
  filtered_dates <- result$dates_ordered[result$dates_ordered$short %in% colnames(filtered_normalized_matrix),]
  
  # Ensure vectors match in length before computing distances
  if(ncol(filtered_normalized_matrix) != length(filtered_dates$date)) {
    stop("Mismatch in the number of columns between normalized matrix and dates.")
  }
  
  # Step 2: Calculate distances using Bray-Curtis for the matrix and Euclidean for the dates
  dates_dist <- as.vector(vegdist(as.Date(filtered_dates$date), method='euclidean'))
  matrix_dist <- 1 - as.vector(vegdist(t(filtered_normalized_matrix), method="bray"))
  
  # Step 3: Combine and order data by dates_dist
  data_comb <- data.frame(dates_dist, matrix_dist)
  data_comb <- data_comb[order(data_comb$dates_dist),]
  
  # Define intervals using the original method of aggregation
  breaks_original <- seq(min(data_comb$dates_dist, na.rm = TRUE), max(data_comb$dates_dist, na.rm = TRUE), by = 30)
  data_comb$intervals_original <- cut(data_comb$dates_dist, breaks = breaks_original, include.lowest = TRUE, right = FALSE)
  means_data_original <- aggregate(cbind(dates_dist, matrix_dist) ~ intervals_original, data = data_comb, FUN = mean)
  
  # Define new intervals for x-axis ticks every 6 months (180 days)
  x_breaks <- seq(0, max(data_comb$dates_dist, na.rm = TRUE), by = 180)
  x_labels <- seq(0, length(x_breaks) - 1) * 6  # Generate labels for each 6-month interval starting from 0 months
  
  # Create the base ggplot
  p <- ggplot(data_comb, aes(x = dates_dist, y = matrix_dist)) +
    geom_point(color = "grey", size = 2, alpha = 0.6) +
    geom_point(data = means_data_original, aes(x = dates_dist, y = matrix_dist), color = "red", size = 3, shape = 19) +
    ggtitle(title) +
    xlab("Months Separating Samples") +
    ylab("Bray-Curtis Similarity") +
    theme_bw() +
    scale_x_continuous(breaks = x_breaks, labels = x_labels) + # Custom x-ticks and labels every 6 months starting from 0
    scale_y_continuous(limits = c(0, 0.6)) # Set y-axis range
  # Print the plot
  print(p)
}


# Call the function
plt1 <- plot_with_averages_ggplot_2(caudo_result, "Caudoviricetes")
plt2 <- plot_with_averages_ggplot_2(mega_result, "Megaviricetes")
plt3 <- plot_with_averages_ggplot_2(others_result, "Others")

grid.arrange(plt1, plt2, plt3 nrow = 3)