library(vegan)

dates <- read.table("/path/to dates")
metadata <- read.csv("/path/to/metadata")

#Load major.module.sums from spearman_rcorr.R

# Create a named vector for mapping
name_mapping <- setNames(dates$RAS_id_GOVID, dates$short)

# Check if all current column names in major.module.sums.2 exist in the name_mapping
missing_names <- !colnames(major.module.sums) %in% names(name_mapping)
if(any(missing_names)) {
  warning("The following column names are missing in 'dates':", paste(colnames(major.module.sums)[missing_names], collapse=", "))
}

# Replace column names using the mapping
# Filter out columns not found in name_mapping to avoid NA
filtered_mapping <- name_mapping[colnames(major.module.sums)]
colnames(major.module.sums) <- filtered_mapping

# Filter metadata to only include matching RAS_ids
filtered_metadata <- metadata[metadata$RAS_id %in% colnames(major.module.sums), ]

# Create a factor with levels in the order of the gene expression data columns
filtered_metadata$RAS_id <- factor(filtered_metadata$RAS_id, levels = colnames(major.module.sums))

# Reorder metadata based on this factor
filtered_metadata <- filtered_metadata[order(filtered_metadata$RAS_id), ]

#Extract environmental parameters
temperature_vector <- as.matrix(dist(filtered_metadata$temp[filtered_metadata$location == "WSC"]))
salinity_vector <- as.matrix(dist(filtered_metadata$sal[filtered_metadata$location == "WSC"]))
daylight_vector <- as.matrix(dist(filtered_metadata$daylight[filtered_metadata$location == "WSC"]))
AW_frac_vector <- as.matrix(dist(filtered_metadata$AW_frac[filtered_metadata$location == "WSC"]))
PW_frac_vector <- as.matrix(dist(filtered_metadata$PW_frac[filtered_metadata$location == "WSC"]))
MLD_vector <- as.matrix(dist(filtered_metadata$MLD[filtered_metadata$location == "WSC"]))
iceDist_vector <- as.matrix(dist(filtered_metadata$iceDist[filtered_metadata$location == "WSC"]))
icePast_vector <- as.matrix(dist(filtered_metadata$icePast[filtered_metadata$location == "WSC"]))
O2_conc_vector <- as.matrix(dist(filtered_metadata$O2_conc[filtered_metadata$location == "WSC"]))
O2_sat_vector <- as.matrix(dist(filtered_metadata$O2_sat[filtered_metadata$location == "WSC"]))
PAR_satellite_vector <- as.matrix(dist(filtered_metadata$PAR_satellite[filtered_metadata$location == "WSC"]))
Richness_vector <- as.matrix(dist(filtered_metadata$Richness[filtered_metadata$location == "WSC"]))
Shannon_diversity_vector <- as.matrix(dist(filtered_metadata$Shannon_diversity[filtered_metadata$location == "WSC"]))
Evenness_vector <- as.matrix(dist(filtered_metadata$Evenness[filtered_metadata$location == "WSC"]))
CV_per_genome_sum_vector <- as.matrix(dist(metadata2$CV_per_genome_sum[metadata2$location == "WSC"]))

split_and_transpose_to_plain_vector <- function(df) {
  # Create a list to store each individual transposed vector
  list_of_transposed_vectors <- list()
  
  # Loop over the rows of the dataframe
  for (i in 1:nrow(df)) {
    # Extract the row into a new dataframe and transpose it
    single_df <- df[i, , drop = FALSE]
    transposed_df <- t(single_df)
    
    # Convert the transposed dataframe to a numeric vector
    transposed_vector <- as.numeric(transposed_df[, 1])
    
    # Create a dynamic name for the transposed vector
    name <- paste0("module_", rownames(df)[i])
    
    # Assign the plain transposed vector to the list with a name
    list_of_transposed_vectors[[name]] <- transposed_vector
  }
  
  # Return the list of plain transposed vectors
  return(list_of_transposed_vectors)
}

# Call the function
list_of_transposed_vectors <- split_and_transpose_to_plain_vector(major.module.sums)

# Access a specific transposed module vector like this
list_of_transposed_vectors$module_11


module_6_dist <- as.matrix(vegdist(list_of_transposed_vectors$module_6, method = "euclidean"))
module_8_dist <- as.matrix(vegdist(list_of_transposed_vectors$module_8, method = "euclidean"))
module_2_dist <- as.matrix(vegdist(list_of_transposed_vectors$module_2, method = "euclidean"))
module_0_dist <- as.matrix(vegdist(list_of_transposed_vectors$module_0, method = "euclidean"))
module_11_dist <- as.matrix(vegdist(list_of_transposed_vectors$module_11, method = "euclidean"))
module_10_dist <- as.matrix(vegdist(list_of_transposed_vectors$module_10, method = "euclidean"))

#Prepare ASV matrix
asv <- read.table("/path/to/asv/table")
date_lookup <- setNames(metadata$date, metadata$RAS_id_GOVID)
asv_dist <- vegdist(t(asv), method = "bray")
asv_dist_matrix <- as.matrix(asv_dist)

datasets <- c("temperature_vector", "salinity_vector", "daylight_vector", "MLD_vector",
              "iceDist_vector","icePast_vector", "PAR_satellite_vector", "asv_dist_matrix", 
              "viral_dist_adjusted", "PW_frac_vector", "major_modules_dist",
              "Richness_vector","Shannon_diversity_vector","Evenness_vector", 
               "CV_per_genome_sum_vector", "module_6_dist","module_8_dist", 
              "module_0_dist", "module_2_dist", "module_10_dist", "module_11_dist")

# Generate combinations of datasets
dataset_combinations <- combn(datasets, 2, simplify = TRUE)

# Perform Mantel test for each combination
mantel_results <- lapply(seq(ncol(dataset_combinations)), function(i) {
  dataset1 <- dataset_combinations[1, i]
  dataset2 <- dataset_combinations[2, i]
  result <- mantel(get(dataset1), get(dataset2), method = "pearson")
  return(result)
})

# Combine results into a named list
names <- apply(dataset_combinations, 2, paste, collapse = ":")
mantel_results <- setNames(mantel_results, names)

# Extract Mantel statistics
mantel_statistics <- lapply(mantel_results, function(result) result$statistic)

# Convert to data frame
mantel_statistics_df <- as.data.frame(unlist(mantel_statistics))
colnames(mantel_statistics_df) <- "Mantel_Statistic"

# Display the data frame
mantel_statistics_df <- mantel_statistics_df[order(mantel_statistics_df$Mantel_Statistic,decreasing=T), , drop = FALSE]
mantel_statistics_df