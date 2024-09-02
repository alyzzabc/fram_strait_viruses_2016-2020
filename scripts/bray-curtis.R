library(tidyr)
library(vegan)

#cVPR
par(mfrow=c(4,1), oma=c(2,2,2,1))
par(mar=c(2, 2, 2, 2))

calculate_running_averages <- function(breadth_values, data_file, metadata_file) {
  # Read the data
  dat <- read.table(data_file, header = TRUE)[,-2]
  dates <- read.table(metadata_file, check.names = FALSE, header = TRUE)
  
  for (breadth in breadth_values) {
    print(paste("Processing breadth coverage:", breadth))
    
    # Filter and merge the data
    b <- dat[dat$Breadth_coverage > breadth,]
    merged_df <- merge(b, dates, by.x = "Library_name", by.y = "short")
    merged_df.matrix <- spread(merged_df[,c(1,2,7)], Library_name, CV_per_genome)
    merged_df.matrix[is.na(merged_df.matrix)] <- 0
    rownames(merged_df.matrix) <- merged_df.matrix[,1]
    merged_df.matrix.2 <- merged_df.matrix[,-1]
    
    # Subset the data for WSC location
    short.wsc.ids <- dates$short[dates$location == "WSC"]
    wsc <- merged_df.matrix.2[,colnames(merged_df.matrix.2) %in% short.wsc.ids]
    wsc.dates <- dates[dates$short %in% colnames(wsc),]
    wsc.dates.ordered <- wsc.dates[order(wsc.dates$date),]
    wsc.sorted.by.date <- wsc[,order(match(colnames(wsc), wsc.dates.ordered$short))]
    wsc.sorted.by.date.2 <- t(t(wsc.sorted.by.date) / rowSums(t(wsc.sorted.by.date), na.rm = TRUE))
    
    # Calculate distances
    dates_dist <- as.vector(vegdist(as.Date(wsc.dates.ordered$date), method = 'euclidean'))
    matrix_dist <- 1 - as.vector(vegdist(t(wsc.sorted.by.date.2), method = "bray"))
    
    # Combine and order data by dates_dist
    data_comb <- data.frame(dates_dist, matrix_dist)
    data_comb <- data_comb[order(data_comb$dates_dist), ]
    
    # Plot the original data
    plot(data_comb$dates_dist, data_comb$matrix_dist, main = paste("Breadth", breadth), 
         xlab = "Euclidean Distance of Dates", ylab = "Bray-Curtis Similarity",
         xaxt = 'n', 
         ylim = c(0, 0.8))
    axis(1, at = seq(0, 1400, by = 180), labels = seq(0, 1400 / 180 * 6, by = 6))  # Custom x-axis for months
    
    # Calculate running averages
    breaks <- seq(min(data_comb$dates_dist, na.rm = TRUE), 
                  max(data_comb$dates_dist, na.rm = TRUE), by = 30)
    
    data_comb$intervals <- cut(data_comb$dates_dist, breaks = breaks, 
                               include.lowest = TRUE, right = FALSE)
    
    # Aggregate data by intervals
    means_data <- aggregate(matrix_dist ~ intervals, data = data_comb, FUN = mean)
    
    # Ensure no NA intervals are present
    means_data <- means_data[!is.na(means_data$intervals), ]
    
    # Simplified conversion of intervals to numeric
    interval_centers <- sapply(levels(data_comb$intervals), function(x) {
      bounds <- unlist(strsplit(gsub("\\(|\\[|\\]|\\)", "", x), ","))
      bounds <- as.numeric(bounds)
      if (length(bounds) == 2) {
        return(mean(bounds))
      } else {
        return(NA)
      }
    })
    
    # Remove NAs from interval_centers and means_data
    valid_indices <- !is.na(interval_centers)
    interval_centers <- interval_centers[valid_indices]
    means_data <- means_data[valid_indices, ]
    
    print("Valid interval centers:")
    print(interval_centers)
    print("Valid means_data:")
    print(means_data)
    
    # Add running average points
    if (nrow(means_data) > 0) {
      points(interval_centers, means_data$matrix_dist, col = 'red', pch = 19)
    }
  }
}

# Call the function with the given breadth coverage values
breadth_values <- c(0, 0.1, 0.25, 0.5)
calculate_running_averages(breadth_values,  "/path/to/normalized/abundances", "/path/to/metadata")

mtext("Read-based, cVPR", outer = TRUE, cex = 1.5, line = 0)

#cVGB
par(mfrow=c(4,1), oma=c(2,2,2,1))
par(mar=c(2, 2, 2, 2))

calculate_running_averages <- function(breadth_values, data_file, metadata_file) {
  # Read the data
  dat <- read.table(data_file, header = TRUE)[,-2]
  dates <- read.table(metadata_file, check.names = FALSE, header = TRUE)
  
  for (breadth in breadth_values) {
    print(paste("Processing breadth coverage:", breadth))
    
    # Filter and merge the data
    b <- dat[dat$Breadth_coverage > breadth,]
    merged_df <- merge(b, dates, by.x = "Library_name", by.y = "short")
    merged_df.matrix <- spread(merged_df[,c(1,2,6)], Library_name, CVPG)
    merged_df.matrix[is.na(merged_df.matrix)] <- 0
    rownames(merged_df.matrix) <- merged_df.matrix[,1]
    merged_df.matrix.2 <- merged_df.matrix[,-1]
    
    # Subset the data for WSC location
    short.wsc.ids <- dates$short[dates$location == "WSC"]
    wsc <- merged_df.matrix.2[,colnames(merged_df.matrix.2) %in% short.wsc.ids]
    wsc.dates <- dates[dates$short %in% colnames(wsc),]
    wsc.dates.ordered <- wsc.dates[order(wsc.dates$date),]
    wsc.sorted.by.date <- wsc[,order(match(colnames(wsc), wsc.dates.ordered$short))]
    wsc.sorted.by.date.2 <- t(t(wsc.sorted.by.date) / rowSums(t(wsc.sorted.by.date), na.rm = TRUE))
    
    # Calculate distances
    dates_dist <- as.vector(vegdist(as.Date(wsc.dates.ordered$date), method = 'euclidean'))
    matrix_dist <- 1 - as.vector(vegdist(t(wsc.sorted.by.date.2), method = "bray"))
    
    # Combine and order data by dates_dist
    data_comb <- data.frame(dates_dist, matrix_dist)
    data_comb <- data_comb[order(data_comb$dates_dist), ]
    
    # Plot the original data
    plot(data_comb$dates_dist, data_comb$matrix_dist, main = paste("Breadth", breadth), 
         xlab = "Euclidean Distance of Dates", ylab = "Bray-Curtis Similarity",
         xaxt = 'n', 
         ylim = c(0, 0.8))
    axis(1, at = seq(0, 1400, by = 180), labels = seq(0, 1400 / 180 * 6, by = 6))  # Custom x-axis for months
    
    # Calculate running averages
    breaks <- seq(min(data_comb$dates_dist, na.rm = TRUE), 
                  max(data_comb$dates_dist, na.rm = TRUE), by = 30)
    
    data_comb$intervals <- cut(data_comb$dates_dist, breaks = breaks, 
                               include.lowest = TRUE, right = FALSE)
    
    # Aggregate data by intervals
    means_data <- aggregate(matrix_dist ~ intervals, data = data_comb, FUN = mean)
    
    # Ensure no NA intervals are present
    means_data <- means_data[!is.na(means_data$intervals), ]
    
    # Simplified conversion of intervals to numeric
    interval_centers <- sapply(levels(data_comb$intervals), function(x) {
      bounds <- unlist(strsplit(gsub("\\(|\\[|\\]|\\)", "", x), ","))
      bounds <- as.numeric(bounds)
      if (length(bounds) == 2) {
        return(mean(bounds))
      } else {
        return(NA)
      }
    })
    
    # Remove NAs from interval_centers and means_data
    valid_indices <- !is.na(interval_centers)
    interval_centers <- interval_centers[valid_indices]
    means_data <- means_data[valid_indices, ]
    
    print("Valid interval centers:")
    print(interval_centers)
    print("Valid means_data:")
    print(means_data)
    
    # Add running average points
    if (nrow(means_data) > 0) {
      points(interval_centers, means_data$matrix_dist, col = 'red', pch = 19)
    }
  }
}

# Call the function with the given breadth coverage values
breadth_values <- c(0, 0.1, 0.25, 0.5)
calculate_running_averages(breadth_values, "/path/to/normalized/abundances", "/path/to/metadata")

mtext("Read-based, cVGB", outer = TRUE, cex = 1.5, line = 0)
