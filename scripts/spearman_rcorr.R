library(Hmisc)
library(reshape2)
library(ggplot2)
library(tidyr)
library(vegan)
library(dplyr)

dat <- read.table("/path/to/normalized/abundance")
b <- dat[dat$Breadth_coverage > 0.25,]
merged_df <- merge(b, metadata, by.x="Library_name", by.y="short")
merged_df.matrix <- spread(merged_df[,c(1,2,7)], Library_name, CV_per_genome)
merged_df.matrix[is.na(merged_df.matrix)] <- 0
rownames(merged_df.matrix) <- merged_df.matrix[,1]
merged_df.matrix.2 <- merged_df.matrix[,-1]
short.wsc.ids <- dates$short[dates$location == "WSC"]
wsc <- merged_df.matrix.2[,colnames(merged_df.matrix.2) %in% short.wsc.ids]
wsc.dates <- dates[dates$short %in% colnames(wsc),]
wsc.dates.ordered <- wsc.dates[order(wsc.dates$date),]
wsc.sorted.by.date <- wsc[,order(match(colnames(wsc), wsc.dates.ordered$short))]
wsc.sorted.by.date.2 <- t(t(wsc.sorted.by.date)/rowSums(t(wsc.sorted.by.date)))
wsc.sorted.by.date.3 <- wsc.sorted.by.date.2[rowSums(wsc.sorted.by.date.2)>0.005,]

network.meta <- read.csv("/path/to/network/file")
merged_data <- merge(network.meta, wsc.sorted.by.date.3, by.x="original_id", by.y="row.names")
merged_data_sor <- merged_data[order(merged_data$LouvainLabelD),]
result <- aggregate(merged_data_sor$temp, by=list(LouvainLabelD=testing_sor$LouvainLabelD), FUN=mean)
major.module.sums <- data.frame(result[result$LouvainLabelD %in% c(6, 8, 0, 2, 10, 11),])
rownames(major.module.sums) <- major.module.sums$LouvainLabelD
major.module.sums <- major.module.sums[,-1]
desired_order <- c(6,8,2,0,10,11) # Explicitly reorder rows in the desired order
major.module.sums <- major.module.sums[as.character(desired_order),]


# Select relevant environmental variables from metadata
metadata <- read.csv("/path/to/metadata")
env_data <- metadata[, c("temp", "sal", "daylight", "AW_frac", "PW_frac", "MLD", "iceDist",
                          "icePast", "O2_conc", "O2_sat", "PAR_satellite", "Richness",
                          "Shannon_diversity", "Evenness", "CV_per_genome_sum")]  # replace 'other_vars' with actual names of other environmental columns

major.module.sums_transposed <- t(major.module.sums)

# Combine the data frames
comb_df <- cbind(env_data, major.module.sums_transposed)
comb_df[] <- lapply(comb_df, function(x) as.numeric(as.character(x))) 
comb_df[] <- lapply(comb_df, function(x) { # Replace NAs with column-wise median (or mean, or any strategy that suits your data)
  if(any(is.na(x))) {
    x[is.na(x)] <- median(x, na.rm = TRUE)  # Ensure no NAs remain
  }
  x
})

sapply(comb_df, function(x) sum(is.na(x))) # Check for any remaining NAs
comb_df <- data.frame(lapply(comb_df, function(x) as.numeric(x)), stringsAsFactors = FALSE) # Ensure the dataframe is entirely numeric

# Correlation analysis
rcorr_results <- rcorr(as.matrix(comb_df), type = "spearman")

# Desired order of rows/columns specified
desired_order <- c("X6", "X8", "X2", "X0", "X11", "X10", "MLD", "sal", "O2_sat", "O2_conc", "PAR_satellite",
                  "iceDist", "CV_per_genome_sum", "temp", "icePast","PW_frac")

# Ensure that all elements in desired_order exist in your data frame
if (!all(desired_order %in% colnames(rcorr_results$r))) {
  stop("One or more elements in the desired order are not present in the correlation matrix.")
}

# Reorder correlation and P-value matrices
reordered_corr <- rcorr_results$r[desired_order, desired_order]
reordered_p <- rcorr_results$P[desired_order, desired_order]

# Melt correlation matrix
melted_cor <- melt(rcorr_results$r)
colnames(melted_cor) <- c("Variable1", "Variable2", "Correlation")

# Melt p-value matrix
melted_p <- melt(rcorr_results$P)
colnames(melted_p) <- c("Variable1", "Variable2", "PValue")

# Merge correlation and p-value data frames
melted_data <- merge(melted_cor, melted_p, by = c("Variable1", "Variable2"))

# Convert factors to characters to avoid comparison issues
melted_data$Variable1 <- as.character(melted_data$Variable1)
melted_data$Variable2 <- as.character(melted_data$Variable2)

# Function to add significance labels based on p-values
significance_label <- function(p) {
  # Check for NA values first
  if (is.na(p)) {
    return("")
  }
  
  if (p < 0.0005) {
    return("***")
  } else if (p < 0.005) {
    return("**")
  } else if (p < 0.05) {
    return("*")
  } else {
    return("")
  }
}

# Apply significance label function to p-values
melted_data$Significance <- sapply(melted_data$PValue, significance_label)

# Plotting
ggplot(melted_data, aes(x = Variable2, y = Variable1, fill = Correlation)) +
  geom_tile() +
  geom_text(aes(label = Significance), color = "black") +
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", midpoint = 0, 
                       limit = c(-1, 1), space = "Lab", name="Spearman\nCorrelation") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title = element_blank()) +
  labs(fill = "Correlation") +
  coord_fixed()
