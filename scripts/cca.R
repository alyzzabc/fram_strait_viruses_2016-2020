library(dplyr)
library(tidyr)
library(vegan)
library(ggplot2)
library(RColorBrewer)

# Prepare abundance_matrix
dat <- read.table("/path/to/normalized/abundance")
metadata <- read.csv("/path/to/metadata")
b <- dat[dat$Breadth_coverage > 0.25,]
merged_df <- merge(b, metadata, by.x = "Library_name", by.y = "short")
merged_df.matrix <- spread(merged_df[, c(1, 3, 8)], Library_name, CV_per_genome)
merged_df.matrix[is.na(merged_df.matrix)] <- 0
rownames(merged_df.matrix) <- merged_df.matrix[, 1]
merged_df.matrix.2 <- merged_df.matrix[, -1]
short.wsc.ids <- metadata$short[metadata$location_x == "WSC"]
wsc <- merged_df.matrix.2[, colnames(merged_df.matrix.2) %in% short.wsc.ids]
wsc.dates <- metadata[metadata$short %in% colnames(wsc),]

# Order by date
wsc.dates$date_x <- as.Date(wsc.dates$date_x)
wsc.dates.ordered <- wsc.dates[order(wsc.dates$date_x),]
wsc.sorted.by.date <- wsc[, order(match(colnames(wsc), wsc.dates.ordered$short))]
abundance_matrix_1 <- t(t(wsc.sorted.by.date) / rowSums(t(wsc.sorted.by.date)))
abundance_matrix_filtered <- abundance_matrix_1[rowSums(abundance_matrix_1) > 0, ]
abundance_matrix_transposed <- t(abundance_matrix_filtered)

# Prepare environmental data
temperature_vector <- metadata$temp
PW_frac_vector <- metadata$PW_frac
MLD_vector <- metadata$MLD
PAR_satellite_vector <- metadata$PAR_satellite

env_data <- data.frame(
  Temperature = temperature_vector,
  PW_fraction = PW_frac_vector,
  MLD = MLD_vector,
  PAR_satellite = PAR_satellite_vector
)

# Run CCA
cca_result <- cca(abundance_matrix_transposed ~ Temperature + PW_fraction + MLD + PAR_satellite, data = env_data)

# Prepare data for ggplot
vectors <- as.data.frame(scores(cca_result, display = "bp"))
virus_scores <- as.data.frame(scores(cca_result, display = "species"))
site_scores <- as.data.frame(scores(cca_result, display = "sites"))
site_scores$date_x <- wsc.dates$date_x

# Extract year and month
site_scores$year <- format(site_scores$date_x, "%Y")
site_scores$month <- format(site_scores$date_x, "%m")

# Create color palette for months and shapes for years
colors_by_month <- c("01" = "#000000", "02" = "#004949", "03" = "#009292", "04" = "#B6DBFF", "05" = "#006DDB", 
                     "06" = "#490092", "07" = "#B66DFF", "08" = "#FF6DB6", "09" = "#920000", "10" = "#924900", 
                     "11" = "#DB6D00", "12" = "#4CBB17")
year_shapes <- c("2016" = 21, "2017" = 22, "2018" = 23, "2019" = 24, "2020" = 25)

# Calculate axis labels with inertia values
inertia_values <- c(0.5699, 0.4627)
total_constrained_inertia <- 1.452
inertia_proportion <- inertia_values / total_constrained_inertia
x_label <- paste0("CCA1 (", round(inertia_proportion[1] * 100, 2), "%)")
y_label <- paste0("CCA2 (", round(inertia_proportion[2] * 100, 2), "%)")

# Plot using ggplot2 with shape and fill aesthetics
p <- ggplot() +
  geom_point(data = virus_scores, aes(x = CCA1, y = CCA2), size = 0.5, color = "#cccccc") +  # Virus points
  geom_point(data = site_scores, aes(x = CCA1, y = CCA2, color = factor(month), fill = factor(month), shape = factor(year)), size = 5) +  # Site points by month and year
  scale_color_manual(values = colors_by_month, name = "Month") +  
  scale_fill_manual(values = colors_by_month, name = "Month") +   
  scale_shape_manual(values = year_shapes, name = "Year") +       
  theme_minimal() +
  labs(x = x_label, y = y_label)

# Add CCA vectors
p <- p + 
  geom_segment(data = vectors, aes(x = 0, y = 0, xend = CCA1, yend = CCA2), 
               arrow = arrow(length = unit(0.2, "inches")), color = "black") +
  geom_text(data = vectors, aes(x = CCA1, y = CCA2, label = variable), 
            vjust = -1, size = 3, hjust = 1.1)

p