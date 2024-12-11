library(Hmisc)
library(reshape2)
library(ggplot2)
library(tidyr)
library(vegan)
library(dplyr)

dat <- read.table("/path/to/normalized/abundances/file",header=T)
metadata <- read.csv("/path/to/metadata/file",check.names=F,header=T)
b <- dat[dat$Breadth_coverage > 0.25,]
merged_df <- merge(b, metadata, by.x="RAS_id", by.y="RAS_id_2")
merged_df.matrix <- spread(merged_df[,c(1,3,8)], RAS_id, CV_per_genome)
merged_df.matrix[is.na(merged_df.matrix)] <- 0
rownames(merged_df.matrix) <- merged_df.matrix[,1]
merged_df.matrix.2 <- merged_df.matrix[,-1]
short.wsc.ids <- metadata$RAS_id_2[metadata$location == "WSC"]
wsc <- merged_df.matrix.2[,colnames(merged_df.matrix.2) %in% short.wsc.ids]
wsc.dates <- metadata[metadata$RAS_id_2 %in% colnames(wsc),]
wsc.dates.ordered <- wsc.dates[order(wsc.dates$date),]
wsc.sorted.by.date <- wsc[,order(match(colnames(wsc), wsc.dates.ordered$RAS_id))]
wsc.sorted.by.date.2 <- t(t(wsc.sorted.by.date)/rowSums(t(wsc.sorted.by.date)))
abundance_matrix <- wsc.sorted.by.date.2

#Prepare asv matrix -- need to remove other F4 samples
asv <- read.table("/path/to/ASV/table/",header=T)
names(asv)[-1] <- gsub("^X0", "", names(asv)[-1])
names(asv)[-1] <- gsub("^X", "", names(asv)[-1])
rownames(asv) <- asv$ASV_name

asv.b <- asv[colnames(asv) %in% colnames(abundance_matrix)]
# Get the column names of asv.b

asv_b_colnames <- colnames(asv.b)
# Subset and order metadata
metadata.ordered <- metadata[order(match(metadata$RAS_id_2, colnames(asv.b))), ]
metadata.2 <- metadata.ordered[metadata.ordered$RAS_id_2 %in% colnames(asv.b),]
rownames(metadata.2) <- metadata.2$RAS_id_2 # Set sample identifiers as row names in metadata
metadata.2 <- metadata.2[, -which(names(metadata.2) == "RAS_id_2")] # Remove the identifier column as it's now row names

#Extract environmental parameters
temp <- metadata2$temp[metadata2$location == "WSC"]
sal <- metadata2$sal[metadata2$location == "WSC"]
daylight <- metadata2$daylight[metadata2$location == "WSC"]
AW_frac <- metadata2$AW_frac[metadata2$location == "WSC"]
PW_frac <- metadata2$PW_frac[metadata2$location == "WSC"]
MLD <- metadata2$MLD[metadata2$location == "WSC"]
iceDist <- metadata2$iceDist[metadata2$location == "WSC"]
icePast <- metadata2$icePast[metadata2$location == "WSC"]
O2_conc_frac <- metadata2$O2_conc[metadata2$location == "WSC"]
O2_sat <- metadata2$O2_sat[metadata2$location == "WSC"]
PAR <- metadata2$PAR_satellite[metadata2$location == "WSC"]
Richness <- metadata2$Richness[metadata2$location == "WSC"]
Shannon_div <- metadata2$Shannon_diversity[metadata2$location == "WSC"]
Evenness <- metadata2$Evenness[metadata2$location == "WSC"]
CV_per_genome_sum <- metadata2$CV_per_genome_sum[metadata2$location == "WSC"]
# Subset and order virus data
abundance_matrix.matched.and.ordered <- abundance_matrix[, colnames(asv.b), drop = FALSE]

# Transposing abundance matrices to match the dimensions for correlation
virus_transposed <- t(abundance_matrix.matched.and.ordered)
asv_transposed <- t(asv.b)

# Select relevant environmental variables from metadata.2 (adjust column names as needed)
env_data <- metadata.2[, c("temp", "sal", "daylight", "AW_frac", "PW_frac", "MLD",
                           "O2_conc", "O2_sat", "PAR_satellite", "Richness",
                          "Shannon_diversity", "Evenness", "CV_per_genome_sum")]  # replace 'other_vars' with actual names of other environmental columns

major.module.sums.2_transposed <- t(major.module.sums.2)
major.module.sums.2_transposed <- major.module.sums.2_transposed[row.names(major.module.sums.2_transposed) != "07_2017_F4_4", ]

# Combine the data frames
#comb_df <- cbind(env_data, virus_transposed[,1:10], asv_transposed[,1:10])
comb_df <- cbind(env_data, major.module.sums.2_transposed)
comb_df[] <- lapply(comb_df, function(x) as.numeric(as.character(x))) # Convert all elements of comb_df to numeric to avoid type mismatch errors
str(comb_df) # Check structure to ensure conversion was successful
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

# Ensure desired order is set
desired_order <- c("X10", "X11", "X0", "X2", "X8", "X6", "MLD", "sal", "O2_sat", "O2_conc", "PAR_satellite",
                   "iceDist", "CV_per_genome_sum", "temp", "icePast","PW_frac")

# Ensure that all elements in desired_order exist in your data frame
if (!all(desired_order %in% colnames(rcorr_results$r))) {
  stop("One or more elements in the desired order are not present in the correlation matrix.")
}

# Reorder correlation and p-value matrices
reordered_corr <- rcorr_results$r[desired_order, desired_order]
reordered_p <- rcorr_results$P[desired_order, desired_order]

# Melt the reordered matrices
melted_cor <- melt(reordered_corr)
colnames(melted_cor) <- c("Variable1", "Variable2", "Correlation")

melted_p <- melt(reordered_p)
colnames(melted_p) <- c("Variable1", "Variable2", "PValue")
melted_p$AdjustedPValue <- p.adjust(melted_p$PValue, method = "BH")

# Merge correlation and adjusted p-value data frames
melted_data <- merge(melted_cor, melted_p[, c("Variable1", "Variable2", "PValue", "AdjustedPValue")], by = c("Variable1", "Variable2"))

# Define significance labels based on adjusted p-values
significance_label <- function(p) {
  if (is.na(p)) {
    return("")
  } else if (p < 0.0005) {
    return("***")
  } else if (p < 0.005) {
    return("**")
  } else if (p < 0.05) {
    return("*")
  } else {
    return("")
  }
}
melted_data$Significance <- sapply(melted_data$AdjustedPValue, significance_label)

# Plot the reordered full correlation matrix
ggplot(melted_data, aes(x = Variable2, y = Variable1, fill = Correlation)) +
  geom_tile() +
  geom_text(aes(label = Significance), color = "black") +
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", midpoint = 0, limit = c(-1, 1), space = "Lab", name="Spearman\nCorrelation") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title = element_blank()) +
  labs(fill = "Correlation") +
  coord_fixed()

