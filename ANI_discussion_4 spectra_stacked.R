#file_paths <- c(
#  '10 V' = "C:/Users/alecn/Desktop/msc project_lipodomics/MS data/ANI17/ANI17MRM10KE-ANI16_200RF_MRM_KE10.mzXML",
#  '15 V' = "C:/Users/alecn/Desktop/msc project_lipodomics/MS data/ANI17/ANI17MRM15KE-ANI16_200RF_MRM_KE15.mzXML",
#  '20 V' = "C:/Users/alecn/Desktop/msc project_lipodomics/MS data/ANI17/ANI17MRM20KE-ANI16_200RF__MRM_KE20.mzXML",
#  '25 V' ="C:/Users/alecn/Desktop/msc project_lipodomics/MS data/ANI17/ANI17MRM25KE-ANI16_200RF_MRM_KE25.mzXML"
#)

file_paths <- c(
  '10 ms' = "C:/Users/alecn/Desktop/msc project_lipodomics/MS data/ANI18/ANI18_10ms.mzXML",
  '20 ms' = "C:/Users/alecn/Desktop/msc project_lipodomics/MS data/ANI18/ANI18MRM20ms-ANI18_200RF_MRM_KE10.mzXML" , 
 '50 ms' = "C:/Users/alecn/Desktop/msc project_lipodomics/MS data/ANI18/ANI18MRM50ms-ANI18_200RF_MRM_ms50.mzXML"
)

library(data.table)
library(dplyr)
library(ggplot2)
library(plotly)
library(mzR)
library(ggrepel)

# Function to select precursor (unchanged)
select_precursor <- function(available_precursors) {
  cat("Available precursors:\n")
  print(sort(available_precursors))
  precursor_input <- as.numeric(readline("Enter one available precursor m/z: "))
  precursor <- available_precursors[which.min(abs(available_precursors - precursor_input))]
  return(precursor)
}

# Function to select retention time (unchanged)
select_rt <- function(header_table) {
  cat("Retention time range (minutes):\n")
  cat(sprintf("Min: %.2f, Max: %.2f\n", min(header_table$retentionTime)/60, max(header_table$retentionTime)/60))
  
  rt_input <- as.numeric(readline("Enter retention time (minutes): ")) * 60  # Convert to seconds
  closest_scan <- header_table$acquisitionNum[which.min(abs(header_table$retentionTime - rt_input))]
  
  cat(sprintf("Selected scan: %d at %.2f min\n", closest_scan, header_table$retentionTime[header_table$acquisitionNum == closest_scan]/60))
  return(closest_scan)
}

# Function to process and average spectra (unchanged)
process_averaged_spectrum <- function(msfile_handle, header_table, target_rt, precursor, time_window = 5, mz_tol = 1.00) {
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
  
  spectra_list <- lapply(scan_nums, function(scan) {
    spec <- mzR::peaks(msfile_handle, scan)
    if (nrow(spec) == 0) return(NULL)
    data.frame(mz = spec[, 1], intensity = spec[, 2], scan = scan)
  })
  
  averaged_spectrum <- bind_rows(Filter(Negate(is.null), spectra_list)) %>%
    group_by(mz) %>%
    summarise(intensity = mean(intensity, na.rm = TRUE), .groups = 'drop')
  
  return(averaged_spectrum)
}

# Select precursor (same for all spectra)
msfile <- mzR::openMSfile(file_paths[1])  # Open first file to select precursor
header <- header(msfile)
precursor <- select_precursor(unique(header$precursorMZ))
close(msfile)

# Select retention time (same for all spectra)
msfile <- mzR::openMSfile(file_paths[1])  # Open first file to select RT
header <- header(msfile)
scan <- select_rt(header)
target_rt <- header$retentionTime[header$acquisitionNum == scan]
close(msfile)

# [Previous code remains the same until the all_spectra processing]

all_spectra <- lapply(names(file_paths), function(name) {
  msfile <- mzR::openMSfile(file_paths[name])
  header <- header(msfile)
  avg_spectrum <- process_averaged_spectrum(msfile, header, target_rt, precursor)
  close(msfile)
  if (nrow(avg_spectrum) > 0) {
    # Normalize to largest peak in this spectrum
    max_intensity <- max(avg_spectrum$intensity, na.rm = TRUE)
    avg_spectrum$intensity_normalized <- (avg_spectrum$intensity / max_intensity) * 100
    avg_spectrum$sample <- name
    return(avg_spectrum)
  } else {
    return(NULL)
  }
}) %>% bind_rows()

all_spectra <- all_spectra %>%
  group_by(sample) %>%
  # Create 5 m/z windows
  mutate(mz_window = cut(mz, breaks = seq(0, 1000, by = 5 ), include.lowest = TRUE)) %>%
           # Within each window, rank peaks by intensity
           group_by(sample, mz_window) %>%
           mutate(
             intensity_rank = rank(-intensity_normalized, ties.method = "first"),
             # Label only top 2 peaks per window that meet criteria
             peak_label = ifelse(
               intensity_rank <= 1 & 
                 intensity_normalized > 0.1 &      # Now using normalized intensity (10% threshold)
                 mz > 400,# Minimum m/z
               sprintf("%.2f", mz),
               ""
             ),
             # Vertical position alternates based on rank
             label_y = intensity_normalized + 0.05  # 3% offset for labels
           ) %>%
           ungroup()
         
         # Create the plot (now using normalized intensity)
         shared_x_plot <- ggplot(all_spectra, aes(x = mz, y = intensity_normalized)) +
           geom_col(width = 0.5, fill = "blue4", alpha = 0.7) +
           geom_text(
             aes(y = label_y, label = peak_label),
             size = 2.3,
             angle = 70,
             hjust = 0,
             check_overlap = TRUE
           ) + 
           facet_wrap(
             ~sample, 
             ncol = 1, 
             scales = "free_y",
             strip.position = "left"
           ) +
           labs(
             title = paste(round(precursor, 4),"Fragmentation Patterns"),
             x = "m/z",
             y = "Relative Intensity (%)"
           ) +
           theme_minimal() +
           theme(
             plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
             axis.text.x = element_text(angle = 45, hjust = 1),
             strip.text = element_text(face = "bold", size = 10),
             panel.spacing = unit(1, "lines")
           ) +
           scale_x_continuous(
             breaks = seq(0, 1000, by = 50),
             limits = c(50, 1000)
           ) +
           scale_y_continuous(
             breaks = seq(0, 100, by = 1),
             limits = c(0, 110)  # Allow space for labels
           )
         
         
         
         # Interactive version (optional)
         ggplotly(shared_x_plot, tooltip = c("x", "y", "sample"))
         
        