file_paths <- c(
  '16_181ANI23' = "C:/Users/alecn/Desktop/msc project_lipodomics/MS data/ANI23_sn1_sn2/ANI23_16_181-ANI23_16181.mzXML",
  '182cis' = "C:/Users/alecn/Desktop/msc project_lipodomics/MS data/ANI22_db/ANI22_182cis.mzXML",
  'bl_empty' = "C:/Users/alecn/Desktop/msc project_lipodomics/MS data/ANI27/ANI27_rerun_BL.mzXML"
)

library(data.table)
library(dplyr)
library(ggplot2)
library(plotly)
library(mzR)

# Function to select precursor
select_precursor <- function(available_precursors) {
  cat("Available precursors:\n")
  print(sort(available_precursors))
  precursor_input <- as.numeric(readline("Enter one available precursor m/z: "))
  precursor <- available_precursors[which.min(abs(available_precursors - precursor_input))]
  return(precursor)
}

# Function to select retention time
select_rt <- function(header_table) {
  cat("Retention time range (minutes):\n")
  cat(sprintf("Min: %.2f, Max: %.2f\n", min(header_table$retentionTime)/60, max(header_table$retentionTime)/60))
  
  rt_input <- as.numeric(readline("Enter retention time (minutes): ")) * 60  # Convert to seconds
  closest_scan <- header_table$acquisitionNum[which.min(abs(header_table$retentionTime - rt_input))]
  
  cat(sprintf("Selected scan: %d at %.2f min\n", closest_scan, header_table$retentionTime[header_table$acquisitionNum == closest_scan]/60))
  return(closest_scan)
}

# Function to process and average spectra around a given RT
process_averaged_spectrum <- function(msfile_handle, header_table, target_rt, precursor, time_window = 5, mz_tol = 1.00) {
  # Find scans within time window AND matching precursor m/z
  scan_nums <- header_table %>%
    filter(
      between(retentionTime, target_rt - time_window, target_rt + time_window),
      abs(precursorMZ - precursor) <= mz_tol
    ) %>%
    pull(acquisitionNum)
  
  if (length(scan_nums) == 0) {
    warning("No scans found with precursor m/z ", precursor, " ± ", mz_tol, " in RT window ", 
            target_rt - time_window, " to ", target_rt + time_window, " seconds.")
    return(data.frame(mz = numeric(), intensity = numeric()))
  }
  
  # Extract and average spectra
  spectra_list <- lapply(scan_nums, function(scan) {
    spec <- mzR::peaks(msfile_handle, scan)
    if (nrow(spec) == 0) return(NULL)
    data.frame(mz = spec[, 1], intensity = spec[, 2], scan = scan)
  })
  
  # Combine and average
  averaged_spectrum <- bind_rows(Filter(Negate(is.null), spectra_list)) %>%
    group_by(mz) %>%
    summarise(intensity = mean(intensity, na.rm = TRUE), .groups = 'drop')
  
  return(averaged_spectrum)
}

# Open file and get header for cis PC analysis
msfile <- mzR::openMSfile(file_paths["bl_empty"])  # cis PC file
header <- header(msfile)

# Select precursor for PC (PC 18:2)
cat("\nSelect precursor for PC (PC 18:2):\n")
precursor <- select_precursor(unique(header$precursorMZ))

# Select retention time
cat("\nSelect retention time for cis PC (PC 18:2):\n")
scan <- select_rt(header)
target_rt <- header$retentionTime[header$acquisitionNum == scan]

# Process spectrum
avg_spectrum <- process_averaged_spectrum(msfile, header, target_rt, precursor)

# Close file handle
close(msfile)

# Create spectrum plot
create_single_spectrum_plot <- function(df, precursor, intensity_threshold = 0.0001, precursor_tol = 5.0, min_mz_spacing = 3) {
  # Normalize intensities
  df <- df %>% mutate(intensity = 100 * intensity / max(intensity, na.rm = TRUE))
  
  # Create label with precursor m/z
  spectrum_label <- paste0("PC 18:2 [", round(precursor, 4), "]")
  
  # Improved peak selection logic
  select_spaced_peaks <- function(data, min_spacing) {
    if (nrow(data) == 0) return(data)
    
    data <- data %>% arrange(desc(abs(intensity)))
    selected <- data[1, ]
    if (nrow(data) == 1) return(selected)
    
    for (i in 2:nrow(data)) {
      current_mz <- data$mz[i]
      
      if (nrow(selected) > 0 && all(abs(current_mz - selected$mz) >= min_spacing, na.rm = TRUE)) {
        selected <- bind_rows(selected, data[i, ])
      }
      
      if (nrow(selected) >= 30) break
    }
    return(selected)
  }
  
  # Get labels (excluding precursor and enforcing spacing)
  label_data <- df %>%
    filter(
      intensity > intensity_threshold,
      abs(mz - precursor) > precursor_tol
    ) 
  
  if (nrow(label_data) > 0) {
    label_data <- select_spaced_peaks(label_data, min_mz_spacing)
  }
  
  # Create the plot
  p <- ggplot(df, aes(x = mz, y = intensity)) +
    geom_col(width = 0.5, fill = "blue4", alpha = 0.7) +
    labs(
      title = spectrum_label,
      x = "m/z",
      y = "Relative Intensity (%)"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 12),
      axis.text.x = element_text(angle = 45, hjust = 1),
      panel.grid.major.x = element_line(color = "gray80"),
      panel.grid.minor.x = element_line(color = "gray90")
    ) +
    scale_y_continuous(breaks = seq(0, 100, by = 1), 
                       minor_breaks = seq(0, 100, by = 1), limits = c(0, 100)) +
    scale_x_continuous(
      breaks = seq(0, 1000, by = 50), 
      minor_breaks = seq(0, 1000, by = 5),
      limits = c(0, 1000)
    )
  
  # Add labels if available
  if (nrow(label_data) > 0) {
    p <- p + geom_text(
      data = label_data,
      aes(label = sprintf("%.2f", mz)),
      size = 2.6,
      vjust = -0.5,
      color = "black"
    )
  } else {
    message("No peaks met the labeling criteria.")
  }
  
  return(p)
}

# Generate and display plot
p <- create_single_spectrum_plot(avg_spectrum, precursor)
ggplotly(p, tooltip = c("x", "y"))
