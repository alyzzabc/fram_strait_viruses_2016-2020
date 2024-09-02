library(dplyr)
library(tidyr)
library(vegan)
library(ggplot2)
library(RColorBrewer)

# Prepare abundance_matrix
dat <- read.table("/path/to/normalized/abundance")
metadata <- read.csv("/path/to/metadata")
b <- dat[dat$Breadth_coverage > 0.25,]
merged_df <- merge(b, metadata, by.x="Library_name", by.y="short")
merged_df.matrix <- spread(merged_df[,c(1,3,8)], Library_name, CV_per_genome)
merged_df.matrix[is.na(merged_df.matrix)] <- 0
rownames(merged_df.matrix) <- merged_df.matrix[,1]
merged_df.matrix.2 <- merged_df.matrix[,-1]
short.wsc.ids <- metadata$short[metadata$location_x == "WSC"]
wsc <- merged_df.matrix.2[,colnames(merged_df.matrix.2) %in% short.wsc.ids]
wsc.dates <- metadata[metadata$short %in% colnames(wsc),]

wsc.dates.ordered <- wsc.dates[order(wsc.dates$date_x),]
wsc.sorted.by.date <- wsc[,order(match(colnames(wsc), wsc.dates.ordered$short))]
wsc.sorted.by.date.2 <- t(t(wsc.sorted.by.date)/rowSums(t(wsc.sorted.by.date)))
abundance_matrix_1 <- wsc.sorted.by.date.2

#Extract environmental parameters
temperature_vector <- metadata$temp
PW_frac_vector <- metadata$PW_frac
MLD_vector <- metadata$MLD
iceDist_vector <- metadata$iceDist
icePast_vector <- metadata$icePast
PAR_satellite_vector <- metadata$PAR_satellite

env_data <- data.frame(
  Temperature = temperature_vector,
  PW_fraction = PW_frac_vector,
  MLD = MLD_vector,
  iceDist = iceDist_vector,
  icePast = icePast_vector,
  PAR_satellite = PAR_satellite_vector
)

# Filter out vOTUs with zero total abundance across all time points
abundance_matrix_filtered <- abundance_matrix_1[rowSums(abundance_matrix_1) > 0, ]
abundance_matrix_transposed <- t(abundance_matrix_filtered)

cca_result <- cca(abundance_matrix_transposed ~ Temperature +
                    PW_fraction + MLD + iceDist + icePast + PAR_satellite, data = env_data)

# Summarize the results
summary(cca_result)

# Plotting
# Create a dataframe from the CCA species scores
virus_scores <- as.data.frame(scores(cca_result, display = "species"))

# Create a dataframe from the CCA site scores
site_scores <- as.data.frame(scores(cca_result, display = "sites"))
site_scores$dates_x <- wsc.dates$date_x

# Convert wsc.dates$date_x to date format
wsc.dates$date_x <- as.Date(wsc.dates$date_x)
months <- format(wsc.dates$date_x, "%m")
paired_colors <- c("#981243","#B94250","#D96B57","#ECA067","#ECD088","#F7EEB5","#ECF5BD","#D6E7AC",
                   "#A2D1B2","#70ACA2","#3B81B2", "#5D538B")
colors_by_month <- setNames(paired_colors, sprintf("%02d", 1:12))
site_scores$month <- months

colors_by_month_2 <- c("01" = "#000000", "02" = "#004949", "03" = "#009292", "04" = "#B6DBFF", "05" = "#006DDB", 
                  "06" = "#490092", "07" = "#B66DFF", "08" = "#FF6DB6", "09" = "#920000", "10" = "#924900", 
                  "11" = "#DB6D00", "12" = "#4CBB17") 
# Plot using ggplot2
# Extract inertia values for CCA1 and CCA2
inertia_values <- c(0.5710, 0.4630)

# Calculate the proportion of total variance for each axis
total_constrained_inertia <- 1.7317  # Total inertia for constrained axes
inertia_proportion <- inertia_values / total_constrained_inertia

# Format the axis labels with inertia values
x_label <- paste0("CCA1 (", round(inertia_proportion[1] * 100, 2), "%)")
y_label <- paste0("CCA2 (", round(inertia_proportion[2] * 100, 2), "%)")

# Plot using ggplot2
p <- ggplot() +
  geom_point(data = virus_scores, aes(x = CCA1, y = CCA2), size = 0.5, color = "#cccccc") +
  geom_point(data = site_scores, aes(x = CCA1, y = CCA2, color = factor(month)), shape = 17, size = 5) +
  scale_color_manual(values = colors_by_month_2) +
  theme_minimal() +
  labs(color = 'Month') +
  xlab(x_label) +
  ylab(y_label)

p + 
  geom_segment(data = vectors, aes(x = 0, y = 0, xend = CCA1, yend = CCA2), 
               arrow = arrow(length = unit(0.2, "inches")), color = "black") +
  geom_text(data = vectors, aes(x = CCA1, y = CCA2, label = variable), 
            vjust = -1, size = 3, hjust = 1.1)