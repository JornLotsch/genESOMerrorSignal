# --- Setup: Define file paths and load libraries --------------------------------

setwd("/home/joern/Aktuell/GenerativeESOM/MouseEAE_DIB")

# Define commonly used file paths
base_path <- "/home/joern/Aktuell/GenerativeESOM/"
r_path    <- "08AnalyseProgramme/R/"
data_path <- "09Originale/"

# Load required packages
library(readxl)

# --- Data Import: Read raw and processed data -----------------------------------

# Read behavioral and lipid profile data from Excel
behavior_lipid_data <- read_excel(
  paste0(base_path, data_path, "Data Behaviorial and lipid profile of the EAE model in SJL mice and effects of FTY720-.xlsx")
)
# Clean column names for easier handling
names(behavior_lipid_data) <- make.names(names(behavior_lipid_data))
View(behavior_lipid_data)  # Inspect imported data

# Read raw and imputed lipidomics data from CSVs
lipid_raw <- read.csv(paste0(base_path, r_path, "Mouse_lipids.csv"))
View(lipid_raw)

lipid_imputed <- read.csv(paste0(base_path, r_path, "Mouse_lipids_transformed_imputed.csv"))

# --- Data Cleaning: Standardize ID columns and extract metadata -----------------

# Ensure first column is named "ID" for consistency
names(lipid_raw)[1]          <- "ID"
names(lipid_imputed)[1]      <- "ID"
names(behavior_lipid_data)[1] <- "ID"

# Select relevant metadata columns for downstream analysis
metadata_cols <- c("ID", "DRUG", "EAE", "GROUP", "AUC", "W0", "W07", "W14", "W17", "W24")
metadata <- subset(behavior_lipid_data, select = metadata_cols)

# Set row names to ID for easier subsetting
rownames(lipid_raw)     <- lipid_raw$ID
rownames(lipid_imputed) <- lipid_imputed$ID
rownames(metadata)      <- metadata$ID

# --- Data Export: Save cleaned datasets for reproducibility ---------------------

write.csv(lipid_raw, "mouse_lipidomics_data_raw.csv", row.names = FALSE)
write.csv(lipid_imputed, "mouse_lipidomics_data_transformed_imputed.csv", row.names = FALSE)
write.csv(metadata, "mouse_lipidomics_metadata.csv", row.names = FALSE)

# --- Data Transformation: Prepare for ggplot visualization ---------------------

# Combine group info with imputed data and reshape to long format for plotting
library(reshape2)
lipid_long <- melt(
  cbind.data.frame(GROUP = as.factor(metadata$GROUP), lipid_imputed),
  id.vars = c("ID", "GROUP")
)

# --- Visualization: Violin plots of lipidomics data by group and variable -------

library(ggplot2)
library(ggthemes)

plot_lipidomics <- ggplot(lipid_long, aes(x = variable, y = value, color = GROUP)) +
  geom_violin(alpha = 0.7) +
  geom_point(position = position_dodge(width = 0.9), size = 1.5) +
  facet_wrap(. ~ variable, scales = "free") +
  theme_light() +
  ggthemes::scale_color_colorblind() +
  theme(
    legend.position = c(0.93, 0.03),           # Place legend inside plot area
    legend.justification = c("right", "bottom"),
    legend.direction = "horizontal",
    strip.background = element_rect(fill = "cornsilk"), # Facet label background
    strip.text = element_text(color = "black"),         # Facet label text color
    axis.text.x = element_blank(),                      # Remove redundant x labels
    axis.ticks.x = element_blank()                      # Remove x axis ticks
  )

print(plot_lipidomics)

# Save plot to file for reporting or publication
ggsave(
  "plot_mouse_lipidomics_data_transformed_imputed.svg",
  plot_lipidomics,
  width = 14, height = 14, limitsize = FALSE
)
