#####

### Script for finding structure in data - sPLS analysis

#####

### Load libraries

suppressPackageStartupMessages(library("vegan"))
suppressPackageStartupMessages(library("reshape2"))
suppressPackageStartupMessages(library("dplyr"))
suppressPackageStartupMessages(library("tidyr"))
suppressPackageStartupMessages(library("ggplot2"))
suppressPackageStartupMessages(library("RColorBrewer"))
suppressPackageStartupMessages(library("patchwork"))

### Define input and output directories

input_dir=('XXXXX')
output_dir=('XXXXXX')
dir.create(output_dir)

#####

### Import data

#####
# Import mapped read data
virus_raw_count=read.table(file=paste0(input_dir,"FRAM_RAS_read_CVPG_mapped_counts.txt"), sep="\t",
                            check.names=F, header=T)  %>%
                            rename_with(~ "Mapped_count", contains("count", ignore.case = TRUE))

# Import metadata
virus_meta=read.table(file=paste0(input_dir,"FRAM_ALL_months_and_moorings.txt"), sep="\t",
                            check.names=F, header=T)

### Define an output prefix based on dataset
output_prefix=('FRAM_RAS_read')
options(scipen = 999999)

### Filter those with breadth of coverage < 0.25 and reformat dataframe to wide using Mapped_contig_count as the values
virus_mapped_count_filt_wide = virus_raw_count %>%
    filter(., Breadth_coverage > 0.25) %>%
    select(Read_name,RAS_id,Mapped_count) %>%
    dcast(Read_name~RAS_id, value.var="Mapped_count", data=.) %>%
    tibble::column_to_rownames(., var="Read_name") %>%
    replace(is.na(.), 0)

### Assess how the richness of viruses detected changes with the number of reads mapped to viruses

# FUNCTION to rarefy mapped read counts, starting with the smallest mapped read count in the data and going to the largest and calculating richness at increments of 50 reads

### Assess how the eveness of viruses detected changes with the number of reads mapped to viruses

# FUNCTION to rarefy mapped read counts, starting with the smallest mapped read count in the data and going to the largest and calculating richness 
# at increments of 50 reads

calculate_incremental_alpha <- function(counts_table){
  
  counts_table <- as.data.frame(counts_table)
  
  results_list <- list()
    
  count_values <- seq(25, 16000, by = 50)
  
  for (count in count_values) {
    # Check which samples have counts less than the current rarefaction level
    sample_counts <- colSums(counts_table)
    too_low <- sample_counts < count
    
    # Rarefy the counts for samples that can be rarefied
    rarefied_abund <- t(rrarefy(t(as.matrix(counts_table[, !too_low])), count))
    
    # Initialize the metrics vectors with NA for samples that can't be rarefied
    R <- rep(NA, ncol(counts_table))
    H <- rep(NA, ncol(counts_table))
    E <- rep(NA, ncol(counts_table))
    
    # Calculate richness, Shannon diversity, and evenness for samples that were rarefied
    R[!too_low] <- specnumber(t(rarefied_abund))
    H[!too_low] <- diversity(t(rarefied_abund), "shannon")
    E[!too_low] <- H[!too_low] / log(R[!too_low])
    
    alpha_div_df <- data.frame(Sample = colnames(counts_table), 
                               Richness = R, 
                               Shannon_diversity = H, 
                               Evenness = E,
                               Sequencing_depth = count)
      
    results_list[[length(results_list) + 1]] <- alpha_div_df
  }
  
  combined_results <- bind_rows(results_list) %>%
    mutate(across(c(Richness, Shannon_diversity, Evenness), ~ format(., scientific = FALSE, trim = TRUE)))
  
  return(combined_results)
}

# Apply the function
virus_incremental_alpha <- calculate_incremental_alpha(virus_mapped_count_filt_wide)

# Export tables
write.table(virus_incremental_richness_df, file=paste0(output_dir,output_prefix,"_richness_rarefaction_values.txt"), sep="\t", row.names=F, quote=F)

# Combine incremental richness with metadata
virus_incremental_alpha_w_meta = virus_incremental_alpha %>%
    left_join(virus_meta, by=c("RAS_id" = "RAS_id_alicia")) %>%
    select(-RAS_id.y)

### Visualize rarefaction of richness

### Colour palette
month_cols = c("Jan" = "#000000", "Feb" = "#004949", "Mar" = "#009292", "Apr" = "#b6dbff", "May" = "#006ddb", "Jun" = "#490092", "Jul" = "#b66dff", "Aug" = "#ff6db6",
    "Sep" = "#920000", "Oct" = "#924900", "Nov" = "#db6d00", "Dec" = "#4CBB17")

### Visualise virus richness across read mappings in a rarefaction-style plot
virus_incremental_richness_rarefaction = virus_incremental_alpha_w_meta %>%
    filter(., Location == "WSC") %>%
    ggplot(., aes(x = as.numeric(Sequencing_depth), y = as.numeric(Richness))) + 
    geom_line(aes(group = RAS_id, colour = Month)) + 
    theme_bw() + 
    labs(x = "Number of viral reads", y = "Richness of vOTUs") + 
    scale_colour_manual(values = month_cols) + 
    theme(legend.position = "bottom",
        legend.title = element_blank(),
        axis.title.y = element_text(size = 11),
        axis.text.y = element_text(size = 11),
        axis.title.x = element_text(size = 11),
        axis.text.x = element_text(size = 11),
        title = element_text(size = 12),
        strip.background.x = element_rect(fill = "white", colour = "black"),
        strip.text.x = element_text(size = 12)) +  
    guides(colour=guide_legend(nrow=2,byrow=TRUE, override.aes = list(size= 8)))

# Export figure
pdf(file = paste0(output_dir,output_prefix,"_viruses_richness_rarefaction.pdf"), height=6, width = 7)
virus_incremental_richness_rarefaction
dev.off()

### Calculate the slope of the richness rarefaction for each sample to determine capturing of diversity

calculate_slope <- function(df) {
  # Remove rows with NA values in any of the metrics or Sequencing_depth
  df <- df %>% 
    filter(., !grepl("NA",Richness)) %>%
    filter(., !grepl("NA",Shannon_diversity)) %>%
    filter(., !grepl("NA",Evenness))
  
  # If there are less than 2 data points after removing NAs, return NAs for all slopes
  if (nrow(df) < 2) {
    return(data.frame(Richness_Slope = NA, Evenness_Slope = NA, Shannon_Slope = NA))
  }
  
  # Fit linear models for each metric
  richness_model <- lm(Richness ~ Sequencing_depth, data = df)
  evenness_model <- lm(Evenness ~ Sequencing_depth, data = df)
  shannon_model <- lm(Shannon_diversity ~ Sequencing_depth, data = df)
  
  # Extract slopes (coefficients of Sequencing_depth)
  richness_slope <- coef(richness_model)["Sequencing_depth"]
  evenness_slope <- coef(evenness_model)["Sequencing_depth"]
  shannon_slope <- coef(shannon_model)["Sequencing_depth"]
  
  # Return a data frame with the slopes
  return(data.frame(Richness_Slope = richness_slope, 
                    Evenness_Slope = evenness_slope, 
                    Shannon_Slope = shannon_slope))
}

# Apply the function to each sample
slopes_df <- virus_incremental_alpha %>%
  group_by(Sample) %>%
  group_modify(~ calculate_slope(.x)) %>%
  rename(., RAS_id = Sample)

# Export slope information
write.table(slopes_df, file=paste0(output_dir,output_prefix,"_richness_rarefaction_slope_values.txt"), sep="\t", row.names=F, quote=F)

### Combine the slope results with information on the month of sampling, mooring and location of mooring
virus_richness_slopes_w_meta = slopes_df %>%
    left_join(virus_meta, by=c("RAS_id" = "RAS_id_alicia")) %>%
    select(-RAS_id.y) %>%
    mutate(Month = factor(Month, levels = c("Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec")))

### Visualise the slope of the richness vs mapped read counts grouped by month of sampling
virus_richness_vs_depth_slope_boxplot = virus_richness_slopes_w_meta %>%
    filter(., Location == "WSC") %>%
    ggplot(., aes(x = Month, y = Richness_Slope)) + 
    geom_boxplot(aes(group = Month, fill = Month), colour = "black") + 
    geom_point(colour = "black", size = 1) + 
    theme_bw() + 
    labs(y = "Slope of Richness") + 
    scale_fill_manual(values = month_cols) + 
    theme(legend.position = "none",
        legend.title = element_blank(),
        axis.title.y = element_text(size = 11),
        axis.text.y = element_text(size = 11),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        title = element_text(size = 12),
        strip.background.x = element_rect(fill = "white", colour = "black"),
        strip.text.x = element_text(size = 12)) +
        guides(colour=guide_legend(nrow=2,byrow=TRUE, override.aes = list(size= 8)))

# Export the figure
pdf(file = paste0(output_dir,output_prefix,"_viruses_richness_slope_boxplots.pdf"), height=6, width = 8)
virus_richness_vs_depth_slope_boxplot
dev.off()

### Combine richness rarefaction and slopes into a single figure and export
virus_incremental_richness_rarefaction_noleg = virus_incremental_richness_rarefaction + theme(legend.position="none")
virus_richness_vs_depth_slope_boxplot_noleg = virus_richness_vs_depth_slope_boxplot + theme(legend.position="none")

pdf(file = paste0(output_dir,output_prefix,"_viruses_richness_rarefaction_and_slope_boxplots.pdf"), height=6, width = 9)
virus_incremental_richness_rarefaction_noleg+virus_richness_vs_depth_slope_boxplot_noleg
dev.off()

######

### The slopes of the richness rarefaction curves calculated above indidcated how much of the diversity is being captured
### at the achieved sequencing depth

### To compare the alpha diversity of vOTUs across months, we will iteratively subsample the dataset to a defined number of viral reads

sort(colSums(virus_mapped_count_filt_wide))
# 500 viral reads determined to be a good cut-off to capture sufficient diversity while excluding only a few, low sequencing depth samples

calculate_alpha_div <- function(counts_table, subsample_depth){
  
  counts_table <- as.data.frame(counts_table)
  
  results_list <- list()
  
  for (i in 1:100) {
    rarefied_abund <- t(rrarefy(t(as.matrix(counts_table)), subsample_depth))
    
    R <- specnumber(t(rarefied_abund))
    H <- diversity(t(rarefied_abund), "shannon")
    E <- H/log(R)
    
    alpha_div_df <- data.frame(Sample = colnames(counts_table), 
                               Richness = R, 
                               Shannon_diversity = H, 
                               Evenness = E)
    
    results_list[[i]] <- alpha_div_df
  }
  
  combined_results <- bind_rows(results_list)

  
  summary_stats <- combined_results %>%
    group_by(Sample) %>%
    summarise(
      Mean_Richness = mean(Richness),
      SD_Richness = sd(Richness),
      Mean_Shannon_diversity = mean(Shannon_diversity),
      SD_Shannon_diversity = sd(Shannon_diversity),
      Mean_Evenness = mean(Evenness),
      SD_Evenness = sd(Evenness)
    )
  
  return(summary_stats)
}

virus_alpha_df <- calculate_incremental_alpha(virus_mapped_count_filt_wide, 500)

# Export table
write.table(virus_alpha_df, file=paste0(output_dir,output_prefix,"_rarefied_alpha_diversity_values.txt"), sep="\t", row.names=F, quote=F)


### Combine the subsampled alpha diversity dataframe with information on the month of sampling, mooring and location of mooring
virus_alpha_df_w_meta = virus_alpha_df %>%
    left_join(virus_meta, by=c("RAS_id" = "RAS_id_alicia")) %>%
    select(-RAS_id.y) %>%
    mutate(Month = factor(Month, levels = c("Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec")))


#### Visualise mean richness per month

# Calculate mean richness and standard deviation per month and visualise as a barplot
virus_rarefied_richness_barplot = virus_alpha_df_w_meta %>%
    select(Month, Mean_Richness, SD_Richness) %>%
    group_by(Month) %>%
    mutate(Richness = mean(Mean_Richness),
    SD = sd(Mean_Richness+SD_Richness)) %>%
    select(-Mean_Richness,-SD_Richness) %>%
    unique() %>%
    ggplot(., aes(x = Month, y = Richness)) + 
    geom_bar(aes(fill = Month), position = "identity", stat = "identity", colour = "black", show.legend=F) + 
    geom_errorbar( aes(ymin=Richness-SD, ymax=Richness+SD), width=0.4, colour="black", alpha=1) +
    theme_bw() + 
    labs(y = "mean vOTU richness at rarefied depth") + 
    scale_fill_manual(values = month_cols) + 
    theme(legend.position = "bottom",
        legend.title = element_blank(),
        axis.title.y = element_text(size = 11),
        axis.text.y = element_text(size = 11),
        axis.title.x = element_blank(),
        axis.text.x = element_text(size = 11, angle = 45, hjust=1),
        title = element_text(size = 12)) +
    guides(fill=guide_legend(nrow=2,byrow=TRUE, override.aes = list(size= 8)))

### Combine with rarefaction curve and export
pdf(file = paste0(output_dir,output_prefix,"_viruses_rarefied_richness_and_rarefaction.pdf"), height=6, width = 9)
virus_incremental_richness_rarefaction_noleg+virus_rarefied_richness_barplot
dev.off()