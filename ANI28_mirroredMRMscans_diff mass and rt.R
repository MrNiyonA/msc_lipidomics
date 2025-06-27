file_paths <- c(
  '16_181ANI23' = "C:/Users/alecn/Desktop/msc project_lipodomics/MS data/ANI23_sn1_sn2/ANI23_16_181-ANI23_16181.mzXML",
  '16_180ANI22' = "C:/Users/alecn/Desktop/msc project_lipodomics/MS data/ANI22_db/ANI22_1618.mzXML"
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
      abs(precursorMZ - precursor) <= mz_tol  # Only include matching precursors
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

# Open files and get headers
msfile1 <- mzR::openMSfile(file_paths["16_181ANI23"])  
msfile2 <- mzR::openMSfile(file_paths["16_180ANI22"])     
header1 <- header(msfile1)
header2 <- header(msfile2)

# Select precursors separately for each condition
cat("\nSelect precursor for PC (PC 16:0/18:1):\n")
precursor_cis <- select_precursor(unique(header1$precursorMZ))

cat("\nSelect precursor for PC (PC 16:0/18:0):\n")
precursor_trans <- select_precursor(unique(header2$precursorMZ))

# Let user select retention time for each condition
cat("\nSelect retention time for cis PC (PC 16:0/18:1):\n")
scan1 <- select_rt(header1)
target_rt1 <- header1$retentionTime[header1$acquisitionNum == scan1]

cat("\nSelect retention time for trans PC (PC 16:0/18:0):\n")
scan2 <- select_rt(header2)
target_rt2 <- header2$retentionTime[header2$acquisitionNum == scan2]

# Process spectra with their respective precursors
avg_cis <- process_averaged_spectrum(msfile1, header1, target_rt1, precursor_cis)
avg_trans <- process_averaged_spectrum(msfile2, header2, target_rt2, precursor_trans)

# Close file handles
close(msfile1)
close(msfile2)

# Create comparison plot (final corrected version)
create_spaced_labeled_plot <- function(df1, df2, precursor1, precursor2, intensity_threshold = 0.00001, precursor_tol = 5.0, min_mz_spacing = 3) {
  # Normalize intensities
  df1 <- df1 %>% mutate(intensity = 100 * intensity / max(intensity, na.rm = TRUE))
  df2 <- df2 %>% mutate(intensity = -100 * intensity / max(intensity, na.rm = TRUE))
  
  # Create technique labels with precursor m/z
  cis_label <- paste0("PC 16:0/18:1 [", round(precursor1, 4), "]")
  trans_label <- paste0("PC 16:0/18:0 [", round(precursor2, 4), "]")
  
  plot_data <- bind_rows(
    mutate(df1, technique = cis_label),
    mutate(df2, technique = trans_label)
  )
  
  # Peak selection logic with spacing requirement
  select_spaced_peaks <- function(data, min_spacing) {
    if (nrow(data) == 0) return(data)
    data <- data %>% arrange(desc(abs(intensity)))
    selected <- data[1, ]
    for (i in 2:nrow(data)) {
      current_mz <- data$mz[i]
      if (all(abs(current_mz - selected$mz) >= min_spacing)) {
        selected <- bind_rows(selected, data[i, ])
      }
      if (nrow(selected) >= 30) break
    }
    return(selected)
  }
  
  # Get labels (excluding precursors and enforcing spacing)
  label_data <- plot_data %>%
    filter(
      abs(intensity) > intensity_threshold,
      abs(mz - precursor1) > precursor_tol,  # Exclude both precursors
      abs(mz - precursor2) > precursor_tol
    ) %>%
    group_by(technique) %>%
    group_modify(~ select_spaced_peaks(.x, min_mz_spacing)) %>%
    ungroup()
  
  # Create the plot
  ggplot(plot_data, aes(x = mz, y = intensity, fill = technique)) +
    geom_col(width = 0.5, position = position_dodge(width = 0.3), alpha = 0.7) +
    scale_fill_manual(values = c("blue4", "red4"),
                      labels = c(cis_label, trans_label)) +
    geom_text(
      data = label_data,
      aes(label = sprintf("%.2f", mz)),
      size = 2.6,
      position = position_dodge(width = 0.4),
      vjust = ifelse(label_data$intensity > 0, -0.5, 1.5),
      color = "black"
    ) +
    labs(
      title = " ",
      x = "m/z",
      y = "Relative Intensity (%)",
      fill = "PC Species"
    ) +
    theme(
      plot.title = element_text(hjust = 0, margin = margin(b = 20)),  # hjust=0 aligns left
      legend.position = "top",
      axis.text.x = element_text(angle = 45, hjust = 0),
      panel.grid.major.x = element_line(color = "gray80"),
      panel.grid.minor.x = element_line(color = "gray90")
    ) +
    geom_hline(yintercept = 0, color = "black") +
    scale_y_continuous(
      breaks = seq(-100, 100, by = 1),
      limits = c(-100, 100)
    ) +
    scale_x_continuous(
      breaks = seq(0, 1000, by = 50), 
      minor_breaks = seq(0, 1000, by = 5),  # Now matches limits
      limits = c(0, 1000)
    )+
    theme_minimal() +
    theme(
      legend.position = "top",
      axis.text.x = element_text(angle = 45, hjust = 1),
      panel.grid.major.x = element_line(color = "gray80"),
      panel.grid.minor.x = element_line(color = "gray90"),
      plot.title = element_text(hjust = 0.5, size = 12)
    )
}

# Generate plot
p <- create_spaced_labeled_plot(avg_cis, avg_trans, precursor_cis, precursor_trans)
ggplotly(p, tooltip = c("x", "y")) %>% 
  layout(hovermode = "x unified", legend = list(orientation = "h", x = 0.4, y = 1.1))
