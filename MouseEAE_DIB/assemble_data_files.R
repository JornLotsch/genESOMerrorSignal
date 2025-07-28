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
library(dplyr)

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

# --- Mapping: Set up old-to-new lipid name mapping -----------------------------

# Step 1: Define mapping as a single string (old name, then new name; new name may have spaces)
mapping_txt <- "
CEREB_AG CEREB_AG
PFC_AG PFC_AG
PFC_1AG PFC_1-AG
CEREB_2AG CEREB_2-AG
PFC_2AG PFC_2-AG
PL_1AG PL_1-AG
PL_2AG PL_2-AG
PL_AG PL_AG
CEREB_AEA CEREB_AEA
PFC_AEA PFC_AEA
CEREB_OEA CEREB_OEA
PFC_OEA PFC_OEA
CEREB_PEA CEREB_PEA
PFC_PEA PFC_PEA
PL_AEA PL_AEA
PL_OEA PL_OEA
PL_PEA PL_PEA
CEREB_C16Ceramid CEREB_Cer d18:1/16:0
HPC_C16Ceramid HPC_Cer d18:1/16:0
PFC_C16Ceramid PFC_Cer d18:1/16:0
CEREB_C18Ceramid CEREB_Cer d18:1/18:0
HPC_C18Ceramid HPC_Cer d18:1/18:0
PFC_C18Ceramid PFC_Cer d18:1/18:0
CEREB_C20Ceramid CEREB_Cer d18:1/20:0
HPC_C20Ceramid HPC_Cer d18:1/20:0
PFC_C20Ceramid PFC_Cer d18:1/20:0
PL_C16Ceramid PL_Cer d18:1/16:0
PL_C20Ceramid PL_Cer d18:1/20:0
PL_C24Ceramid PL_Cer d18:1/24:0
PL_C24_1Ceramid PL_Cer d18:1/24:1
CEREB_C18Sphinganin CEREB_Cer d18:0/18:0
HPC_C18Sphinganin HPC_Cer d18:0/18:0
PFC_C18Sphinganin PFC_Cer d18:0/18:0
PL_Sphinganin PL_SPH d18:0
PL_Sphinganin.1 PPL_S1P d18:0
PL_Sphingosin PL_SPH d18:1
PL_Sphingosin.1 PPL_S1P d18:1
PL_C16Sphinganin PL_Cer d18:0/16:0
PL_C18Sphinganin PL_Cer d18:0/18:0
PL_C24Sphinganin PL_Cer d18:0/24:0
PL_C24_1Sphinganin PL_Cer d18:0/24:1
CEREB_Sphinganin CEREB_SPH d18:0
HPC_Sphinganin HPC_SPH d18:0
PFC_Sphinganin PFC_SPH d18:0
CEREB_Sphingosin CEREB_SPH d18:1
HPC_Sphingosin HPC_SPH d18:1
PFC_Sphingosin PFC_SPH d18:1
CEREB_Sphingosin.1P CEREB_S1P d18:1
HPC_Sphingosin.1P HPC_S1P d18:1
PFC_Sphingosin.1P PFC_S1P d18:1
CEREB_LPA16_0 CEREB_LPA16:0
CEREB_LPA18_0 CEREB_LPA18:0
CEREB_LPA18_1 CEREB_LPA18:1
CEREB_LPA18_2 CEREB_LPA18:2
CEREB_LPA18_3 CEREB_LPA18:3
CEREB_LPA20_4 CEREB_LPA20:4
PL_LPA16_0 PL_LPA16:0
PL_LPA18_0 PL_LPA18:0
PL_LPA18_1 PL_LPA18:1
PL_LPA18_2 PL_LPA18:2
PL_LPA18_3 PL_LPA18:3
PL_LPA20_4 PL_LPA20:4
"

# Step 2: Parse mapping string to data frame (split each line at first space)
mapping_txt <- gsub("^[ \t\r\n]+|[ \t\r\n]+$", "", mapping_txt)  # trim blank lines
lines    <- readLines(textConnection(mapping_txt))
lines    <- lines[lines != ""]
old      <- sub(" .*", "", lines)
new      <- sub("^[^ ]+ +", "", lines)
mapping  <- data.frame(old = old, new = new, stringsAsFactors = FALSE)

# --- Utility: Rename lipid columns in a consistent way (preserve 'ID' column) ---

rename_lipid_cols <- function(df, mapping) {
  old_names <- names(df)[-1]  # All columns except 'ID'
  # Map every old name to new (if found in mapping), otherwise keep unchanged
  new_names <- sapply(
    old_names,
    function(x) {
      idx <- which(mapping$old == x)
      if (length(idx) == 1) mapping$new[idx] else x
    }
  )
  # Set column names ('ID' always first)
  names(df) <- c("ID", new_names)
  return(df)
}

# --- Apply mapping to both dataframes ------------------------------------------

lipid_raw     <- rename_lipid_cols(lipid_raw, mapping)
lipid_imputed <- rename_lipid_cols(lipid_imputed, mapping)

# --- Optional: Inspect new column names ----------------------------------------
names(lipid_raw)
names(lipid_imputed)

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
    legend.position.inside = TRUE, legend.position = c(0.93, 0.03),          
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
