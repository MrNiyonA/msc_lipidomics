file_paths <- c(
  #'Brain Extract' = "C:/Users/alecn/Desktop/msc project_lipodomics/MS data/ANI27/ANI27_rerun_BL.mzXML",
  'PC 16:0/18:1(9)' = "C:/Users/alecn/Desktop/msc project_lipodomics/MS data/ANI28/ANI28_10micrM_16181-ANI27_16181_10micM.mzXML",
  'PC 18:1(9)/16:0' = "C:/Users/alecn/Desktop/msc project_lipodomics/MS data/ANI28/ANI28_10micrM_18116-ANI28_18116_10micM.mzXML",
  'PC 18:1(9 cis)' = "C:/Users/alecn/Desktop/msc project_lipodomics/MS data/ANI28/ANI28_10micrM_181cis-ANI28_181cis_10micM.mzXML",
  'PC 18:1(9 trans)' = "C:/Users/alecn/Desktop/msc project_lipodomics/MS data/ANI28/ANI28_10micrM_181trans-ANI28_181trans_10micM.mzXML",
  'PC 18:2(9,12 cis)' = "C:/Users/alecn/Desktop/msc project_lipodomics/MS data/ANI28/ANI28_10micrM_182IDA-ANI28_282cisIDA.mzXML",
  'PC 16:0/18:0' = "C:/Users/alecn/Desktop/msc project_lipodomics/MS data/ANI28/ANI28_10micrM_16180-ANI27_16180_micM.mzXML"
)

library(ggplot2)
library(plotly)
library(mzR)
library(dplyr)

# Function to extract Base Peak Intensity (BPI)
extract_bpi <- function(file_path) {
  msfile <- mzR::openMSfile(file_path)
  header_data <- mzR::header(msfile)
  bpi_data <- data.frame(
    retentionTime = header_data$retentionTime / 60,  # Convert to minutes
    BPI = header_data$basePeakIntensity,            # Raw BPI (no normalization)
    file = basename(file_path)
  )
  close(msfile)
  return(bpi_data)
}

# Extract BPI (without normalization or filtering)
all_bpi_data <- lapply(file_paths, extract_bpi) %>%
  bind_rows(.id = "sample_name")

# Plot with raw BPI values (full retention time range)
overlaid_bpi_plot <- ggplot(all_bpi_data, aes(
  x = retentionTime, 
  y = BPI,
  color = sample_name
)) +
  geom_line(linewidth = 0.7) +
  labs(
    title = "Overlaid Base Peak Chromatograms",
    x = "Retention Time (min)",
    y = "Base Peak Intensity (cps)",
    color = "Sample"
  ) +
  theme_minimal() +
  theme(
    legend.position = "top",
    plot.title = element_text(hjust = 0.5),
    legend.text = element_text(size = 14) 
  ) +
  scale_x_continuous(
    breaks = seq(0, 20, by = 1), 
    minor_breaks = seq(0, 10, by = 0.1),  # Now matches limits
    limits = c(0, 20)
  )+
  scale_color_brewer(palette = "Set1")

# Make interactive
ggplotly(overlaid_bpi_plot) %>%
  layout(hovermode = "x unified")

