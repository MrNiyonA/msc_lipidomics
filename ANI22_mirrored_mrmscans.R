file_paths <- c(
  `181_16` = "C:/Users/alecn/Desktop/msc project_lipodomics/MS data/ANI21/ANI21_18116_redo-ANI20_IDA_ 18_1 _ 16.mzXML",
  `16_181` = "C:/Users/alecn/Desktop/msc project_lipodomics/MS data/ANI21/ANI21_16181_redo-ANI20_IDA_ 16 _ 18_1 .mzXML"  )

library(data.table)
library(dplyr)
library(ggplot2)
library(plotly)
library(mzR)

select_precursor <- function(available_precursors) {
  cat("Available precursors:\n")
  print(sort(available_precursors))
  precursor_input <- as.numeric(readline("Enter one available precursor m/z: "))
  precursor <- available_precursors[which.min(abs(available_precursors - precursor_input))]
  return(precursor)  # Explicit return
}

# Function to process and average spectra within a time window
process_averaged_spectrum <- function(msfile_handle, header_table, center_scan, time_window = 5) {
  # Get retention time of center scan
  center_time <- header_table$retentionTime[header_table$acquisitionNum == center_scan]
  
  # Find scans within time window
  scan_nums <- header_table %>%
    filter(between(retentionTime, center_time - time_window, center_time + time_window)) %>%
    pull(acquisitionNum)
  
  # Extract all spectra
  spectra_list <- lapply(scan_nums, function(scan) {
    spec <- mzR::peaks(msfile_handle, scan)
    if (nrow(spec) == 0) return(NULL)
    spec <- as.data.frame(spec)
    colnames(spec) <- c("mz", "intensity")
    spec$scan <- scan
    spec
  })
  
  # Remove NULL spectra
  spectra_list <- Filter(Negate(is.null), spectra_list)
  
  # Combine all spectra
  all_spectra <- bind_rows(spectra_list)
  
  # Create averaged spectrum (mean intensity per m/z)
  averaged_spectrum <- all_spectra %>%
    group_by(mz) %>%
    summarise(intensity = mean(intensity, na.rm = TRUE), .groups = 'drop')
  
  return(averaged_spectrum)
}

# Function to average spectra across replicates
process_replicates <- function(file_group, time_window = 5) {
  all_spectra <- list()
  
  for (file in file_group) {
    msfile <- mzR::openMSfile(file)
    header <- header(msfile)
    center_scan <- header$acquisitionNum[which.max(header$totIonCurrent)]
    
    # Get averaged spectrum for this replicate
    avg_spectrum <- process_averaged_spectrum(msfile, header, center_scan, time_window)
    all_spectra[[file]] <- avg_spectrum
    
    close(msfile)
  }
  
  # Combine all replicates and average
  combined <- bind_rows(all_spectra, .id = "file") %>%
    group_by(mz) %>%
    summarise(intensity = mean(intensity, na.rm = TRUE), .groups = "drop")
  
  return(combined)
}

# Group your files by configuration
config_181_16 <- file_paths[grep("181_16", names(file_paths))]
config_16_181 <- file_paths[grep("16_181", names(file_paths))]

# Process all replicates for each configuration
avg_16_181 <- process_replicates(config_16_181)  # Combined 16:0/18:1 replicates
avg_181_16 <- process_replicates(config_181_16)  # Combined 18:1/16:0 replicates


# First you need to open the files and get header tables
centroid_msfile_handle <- mzR::openMSfile(file_paths["16_181"])
profile_header_table <- header(centroid_msfile_handle)

centroid_msfile_handle2 <- mzR::openMSfile(file_paths["181_16"])
profile_header_table2 <- header(centroid_msfile_handle2)

# Select precursor (you need to have precursor values available)
#select_precursor(unique(c(profile_header_table$precursorMZ, profile_header_table2$precursorMZ)))

precursor <- 760.578

# Find scans with highest TIC (you'll need to define EAD_num and CID_num)
PC16181 <- profile_header_table$acquisitionNum[which.max(profile_header_table$totIonCurrent)]
PC18116 <- profile_header_table2$acquisitionNum[which.max(profile_header_table2$totIonCurrent)]

# Process both 16_181 and 18:1/16 spectra
avg_16_181 <- process_averaged_spectrum(centroid_msfile_handle, profile_header_table, PC16181)
avg_181_16 <- process_averaged_spectrum(centroid_msfile_handle2, profile_header_table2, PC18116)

# Create a line plot comparison
create_line_plot <- function(df1, df2, precursor) {
  # Mirror 181_16 spectrum
  df2$intensity <- -df2$intensity
  
  # Combine data with source indicator
  plot_data <- bind_rows(
    mutate(df1, technique = "16_181"),
    mutate(df2, technique = "181_16")
  )
  
  ggplot(plot_data, aes(x = mz, y = intensity, color = technique)) +
    geom_line(linewidth = 0.5) +
    scale_color_manual(values = c("16_181" = "red4", "181_16" = "blue4")) +
    ggtitle(paste("Averaged Spectra Comparison for Precursor", round(precursor, 4))) +
    labs(x = "m/z", y = "Relative Intensity (%)", 
         caption = "Averaged over ±5 seconds around highest TIC scan") +
    geom_hline(yintercept = 0, color = "black") +
    theme_minimal() +
    theme(legend.position = "top",
          plot.caption = element_text(hjust = 0.5, face = "italic"))+
    geom_hline(yintercept = 0, color = "black") +
    scale_y_continuous(
      breaks = seq(-100, 100, by = 1),
      minor_breaks = seq(-100, 100, by = 1),  
      limits = c(-100, 100)
    )+
    scale_x_continuous(
      breaks = seq(0, 1000, by = 50), 
      minor_breaks = seq(0, 1000, by = 5),  # Now matches limits
      limits = c(0, 1000))+
  geom_text(
    data = label_data,
    aes(label = sprintf("%.2f", mz)),
    size = 3.5,
    position = position_dodge(width = 1.0),
    vjust = ifelse(label_data$intensity > 0, -50, 50),
    color = "black")
}

# Generate and display plot
p <- create_line_plot(avg_16_181, avg_181_16, precursor)
ggplotly(p, tooltip = c("x", "y")) %>% 
  layout(hovermode = "x unified",
         legend = list(orientation = "h", x = 0.4, y = 1.1))




# Create mirrored difference spectra
create_difference_spectra <- function(df1, df2) {
  # Join by mz
  combined <- full_join(df1, df2, by = "mz", suffix = c("_16_181", "_181_16")) %>%
    mutate(
      intensity_16_181 = ifelse(is.na(intensity_16_181), 0, intensity_16_181),
      intensity_181_16 = ifelse(is.na(intensity_181_16), 0, intensity_181_16),
      diff = intensity_16_181 - intensity_181_16
    )
  
  # Split based on sign of difference
  diff_16_181vs181_16 <- combined %>%
    filter(diff > 0) %>%
    select(mz, intensity = diff)
  
  diff_181_16vs16_181 <- combined %>%
    filter(diff < 0) %>%
    mutate(intensity = -diff) %>%  # Flip for mirroring
    select(mz, intensity)
  
  return(list(
    diff_16_181vs181_16 = diff_16_181vs181_16,
    diff_181_16vs16_181 = diff_181_16vs16_181
  ))
}

# Generate the difference spectra
diff_spectra <- create_difference_spectra(avg_16_181, avg_181_16)

# Plot mirrored difference spectra
create_mirrored_difference_plot <- function(diff1, diff2, precursor) {
  # Flip second spectrum for mirroring
  diff2$intensity <- -diff2$intensity
  
  plot_data <- bind_rows(
    mutate(diff1, source = "16_181 > 181_16"),
    mutate(diff2, source = "181_16 > 16_181")
  )
  
  ggplot(plot_data, aes(x = mz, y = intensity, color = source)) +
    geom_line(linewidth = 0.5) +
    scale_color_manual(values = c("16_181 > 181_16" = "red", "181_16 > 16_181" = "blue")) +
    geom_hline(yintercept = 0, color = "black") +
    ggtitle(paste("Mirrored Intensity Differences for Precursor", round(precursor, 4))) +
    labs(x = "m/z", y = "Intensity Difference",
         caption = "Red: Higher in 16_181 | Blue: Higher in 181_16") +
    theme_minimal() +
    theme(legend.position = "top",
          plot.caption = element_text(hjust = 0.5, face = "italic"))
}

# Display plot
p_mirrored_diff <- create_mirrored_difference_plot(
  diff_spectra$diff_16_181vs181_16,
  diff_spectra$diff_181_16vs16_181,
  precursor
)

ggplotly(p_mirrored_diff, tooltip = c("x", "y")) %>%
  layout(hovermode = "x unified",
         legend = list(orientation = "h", x = 0.3, y = 1.1))
