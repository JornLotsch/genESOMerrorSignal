# --- Setup: Define file paths and load libraries --------------------------------

setwd("/home/joern/Aktuell/GenerativeESOM/08AnalyseProgramme/R/genESOMerrorSignal/MouseEAE_DIB")

# Define commonly used file paths
base_path <- "/home/joern/Aktuell/GenerativeESOM/"
r_path    <- "08AnalyseProgramme/R/"
data_path <- "09Originale/"

# Load required packages
library(readxl)
library(ggplot2)
library(ggthemes)

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
metadata_cols <- c("ID", "DRUG", "EAE", "GROUP")
metadata <- subset(behavior_lipid_data, select = metadata_cols)
metadata$DRUG <- metadata$DRUG - 1
metadata$DRUG[metadata$DRUG < 0] <- 0

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

plot_lipidomics <- ggplot(lipid_long, aes(x = variable, y = value, color = GROUP, fill = GROUP)) +
  geom_violin(alpha = 0.2) +
  geom_point(position = position_dodge(width = 0.9), size = 1.5) +
  facet_wrap(. ~ variable, scales = "free") +
  theme_light() +
  ggthemes::scale_color_colorblind() +
  ggthemes::scale_fill_colorblind() +
  theme(
    legend.position = c(0.93, 0.03),          
    legend.justification = c("right", "bottom"),
    legend.direction = "horizontal",
    strip.background = element_rect(fill = "cornsilk"),
    strip.text = element_text(color = "black"),      
    axis.text.x = element_blank(),                   
    axis.ticks.x = element_blank()  
  ) +
  labs(title = "Lipidomics data per treatment group", color = "GROUP", fill = "GROUP")

print(plot_lipidomics)

# Save plot to file for reporting or publication
ggsave(
  "plot_mouse_lipidomics_data_transformed_imputed.svg",
  plot_lipidomics,
  width = 14, height = 14, limitsize = FALSE
)

# --- Data set sizes and some stats -------
dim(lipid_raw)
dim(lipid_imputed)
dim(metadata)

table(metadata$GROUP)
sum(is.na(lipid_raw[,-1]))

View(metadata)
