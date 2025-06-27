# files were transformed from WIFF file to mzXML (with ProteoWizards MSConvert)
file_paths <- c(
  '181_16_mrm' ="C:/Users/alecn/Desktop/msc project_lipodomics/MS data/ANI23_sn1_sn2/ANI23_181_16-ANI23_18116.mzXML",
  '16_181_mrm' ="C:/Users/alecn/Desktop/msc project_lipodomics/MS data/ANI23_sn1_sn2/ANI23_16_181-ANI23_16181.mzXML",
  'cis_meoh' = "C:/Users/alecn/Desktop/msc project_lipodomics/MS data/ANI28/ANI28_10micrM_181cis-ANI28_181cis_10micM.mzXML",
  'trans_meoh' =  "C:/Users/alecn/Desktop/msc project_lipodomics/MS data/ANI28/ANI28_10micrM_181trans-ANI28_181trans_10micM.mzXML",
  'db_bl_cid' = "C:/Users/alecn/Desktop/msc project_lipodomics/MS data/ANI28/ANI28_BLdbs_CID-ANI28_rerun_BL_dbs_IDACID_10micM.mzXML",
  'db_bl_ead' = "C:/Users/alecn/Desktop/msc project_lipodomics/MS data/ANI28/ANI28_BLdbs_EAD-ANI28_rerun_BL_dbs_IDAEAD_10micM.mzXML",
  'bl_empty' = "C:/Users/alecn/Desktop/msc project_lipodomics/MS data/ANI27/ANI27_rerun_BL.mzXML",
  'ANI22' = "C:/Users/alecn/Desktop/msc project_lipodomics/MS data/ANI22_db/ANI22_12-ANI22_MRMR12.mzXML",
  'PC16180' = "C:/Users/alecn/Desktop/msc project_lipodomics/MS data/ANI22_db/ANI22_1618.mzXML",
  'BL_182_EAD' = "C:/Users/alecn/Desktop/msc project_lipodomics/MS data/ANI28/ANI28_BL_182_16181EAD-ANI28_BL_182_16181_EAD.mzXML", 
  'BL_182_CID' = "C:/Users/alecn/Desktop/msc project_lipodomics/MS data/ANI28/ANI28_BL_182_16181CID-ANI28_BL_182_16181_CID.mzXML" 
  
  
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

# Function to select retention time (instead of scan number)
select_rt <- function(header_table) {
  cat("Retention time range (minutes):\n")
  cat(sprintf("Min: %.2f, Max: %.2f\n", min(header_table$retentionTime)/60, max(header_table$retentionTime)/60))
  
  rt_input <- as.numeric(readline("Enter retention time (minutes): ")) * 60  # Convert to seconds
  closest_scan <- header_table$acquisitionNum[which.min(abs(header_table$retentionTime - rt_input))]
  
  cat(sprintf("Selected scan: %d at %.2f min\n", closest_scan, header_table$retentionTime[header_table$acquisitionNum == closest_scan]/60))
  return(closest_scan)
}

# Function to process and average spectra around a given RT (with precursor filtering)
process_averaged_spectrum <- function(msfile_handle, header_table, target_rt, precursor, time_window = 10 , mz_tol = 1.00) {
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
msfile1 <- mzR::openMSfile(file_paths["cis_meoh"])
msfile2 <- mzR::openMSfile(file_paths["db_bl_ead"])
header1 <- header(msfile1)
header2 <- header(msfile2)

# Select precursor
precursor <- select_precursor(unique(c(header1$precursorMZ, header2$precursorMZ)))

# Let user select retention time (converts to closest scan)
cat("\nSelect retention time for ead:\n")
scan1 <- select_rt(header1)
target_rt1 <- header1$retentionTime[header1$acquisitionNum == scan1]  # Get exact RT of selected scan

cat("\nSelect retention time for cid:\n")
scan2 <- select_rt(header2)
target_rt2 <- header2$retentionTime[header2$acquisitionNum == scan2]

# Process spectra (now using target_rt instead of scan number)
avg_ead <- process_averaged_spectrum(msfile1, header1, target_rt1, precursor)
avg_cid <- process_averaged_spectrum(msfile2, header2, target_rt2, precursor)

# Close file handles
close(msfile1)
close(msfile2)

# Create comparison plot
create_spaced_labeled_plot <- function(df1, df2, precursor, intensity_threshold = 0.0001, precursor_tol = 1.0, min_mz_spacing = 3) {
  # Normalize intensities to highest peak in spectrum (normally precursor peak)
  df1 <- df1 %>% mutate(intensity = 100 * intensity / max(intensity, na.rm = TRUE))
  df2 <- df2 %>% mutate(intensity = -100 * intensity / max(intensity, na.rm = TRUE))
  
  plot_data <- bind_rows(
    mutate(df1, technique = "PC 18:1 cis standard"),
    mutate(df2, technique = "PC 18:1 cis bl")
  )
  
  # Peak selection logic (unchanged from your spacing requirement)
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
  
  # Get labels (excluding precursor and enforcing spacing)
  label_data <- plot_data %>%
    filter(
      abs(intensity) > intensity_threshold,
      abs(mz - precursor) > precursor_tol
    ) %>%
    group_by(technique) %>%
    group_modify(~ select_spaced_peaks(.x, min_mz_spacing)) %>%
    ungroup()
  
  # plotting
  ggplot(plot_data, aes(x = mz, y = intensity, fill = technique)) +
    geom_col(width = 0.5, position = position_dodge(width = 0.3), alpha = 0.7) +
    scale_fill_manual(values = c("PC 18:1 cis standard" = "red4", "PC 18:1 cis bl" = "blue4")) +
    # Peak labels (now with spacing)
    geom_text(
      data = label_data,
      aes(label = sprintf("%.2f", mz)),
      size = 3.7,
      position = position_dodge(width = 5.0),
      vjust = ifelse(label_data$intensity > 0, -50, 50),
      color = "black"
    ) +
    # Title
    labs(
      title = paste("PC18:1 EAD", round(precursor, 4)),
      x = "m/z",
      y = "Relative Intensity (%)",
      fill = "Condition"
    ) +
    geom_hline(yintercept = 0, color = "black") +
    scale_y_continuous(
      breaks = seq(-100, 100, by = 1),
      minor_breaks = seq(-100, 100, by = 1),  
      limits = c(-100, 100)
    )+
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
      plot.title = element_text(hjust = 0.5)  # Center title
    ) +
    theme(
      plot.title = element_text(
        hjust = 0.1,
        margin = margin(b = 20, unit = "pt")
      )
    )
}

# Generate plot
p <- create_spaced_labeled_plot(avg_ead, avg_cid, precursor)
ggplotly(p, tooltip = c("x", "y")) %>% 
  layout(hovermode = "x unified", legend = list(orientation = "h", x = 0.4, y = 1.1))

