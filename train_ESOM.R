#' ESOM U-matrix Training Script
#'
#' This script trains an Emergent Self-Organizing Map (ESOM) and generates U-matrix visualizations.
#' It takes processed data from a previous data preparation step and:
#' - Optionally converts data to percentages
#' - Trains an ESOM model
#' - Generates U-matrix, P-matrix, and Island visualizations
#' - Creates 2D and 3D visualizations of the ESOM maps
#' - Saves results to files

# Set working directory to current file location
setwd(dirname(rstudioapi::getSourceEditorContext()$path))

#' Load required libraries
#' @note These libraries must be installed before running the script
library(Umatrix)      # For U-matrix generation and visualization
library(caret)        # For data pre-processing
library(ggplot2)      # For visualization
library(cowplot)      # For arranging multiple plots
library(rgl)          # For 3D visualization
library(dplyr)        # For data manipulation

#' Configuration settings
#' @section File paths:
input_file <- "Mouse_lipids_processed.csv"   # Input processed data file (from previous script)
output_dir <- "umatrix_output/"              # Output directory for U-matrix files
output_prefix <- "UmxMouse"                  # Prefix for output files

#' @section ESOM parameters:
esom_lines <- 50                          # Number of lines in the ESOM grid
esom_columns <- 80                        # Number of columns in the ESOM grid
esom_toroid <- TRUE                       # Whether to use a toroid topology
esom_epochs <- 20                         # Number of training epochs
enable_plots <- TRUE                      # Whether to create and display plots
enable_percent_conversion <- FALSE        # Whether to convert data to percentage
enable_3d_visualization <- TRUE           # Whether to create 3D visualizations
enable_file_output <- TRUE                # Whether to save files to disk

#' @param random_seed Random seed for reproducibility
random_seed <- 42

#' Function to convert data to percentages
#'
#' @param data Data frame to convert
#' @return Data frame with values converted to percentages
#' @examples
#' data_percent <- toPercent(raw_data)
toPercent <- function(data) {
  # Ensure data is numeric
  data_numeric <- as.matrix(data)

  # For each row, divide by the sum of the row and multiply by 100
  data_percent <- t(apply(data_numeric, 1, function(x) {
    # Skip rows that sum to 0 to avoid division by zero
    if(sum(x) == 0) return(x)
    return((x / sum(x)) * 100)
  }))

  # Convert back to data frame and preserve row and column names
  data_percent <- as.data.frame(data_percent)
  colnames(data_percent) <- colnames(data)
  rownames(data_percent) <- rownames(data)

  return(data_percent)
}

#' Function to create and save a 3D visualization of the U-matrix
#'
#' @param umatrix U-matrix object
#' @param best_matches Best matches of the ESOM
#' @param cls Class labels
#' @param cls_colors Colors for the classes
#' @param toroid Whether the ESOM has toroid topology
#' @param imx Island mask
#' @param filename Output filename
#' @param output_dir Output directory
#' @param bm_size Size of the best match points
#' @param remove_ocean Whether to remove the "ocean" areas
#' @param show_axis Whether to show the axes
#' @param smooth_slope Whether to smooth the slopes
#' @examples
#' create_save_3d_visualization(UmxMouse, BMUmxMouse, MouseCls, rainbow(2), TRUE, ImxMouse, "Mouse_Umx3D.png", "output/")
create_save_3d_visualization <- function(umatrix, best_matches, cls, cls_colors, toroid, imx,
                                         filename, output_dir, bm_size = 0.6, remove_ocean = TRUE,
                                         show_axis = TRUE, smooth_slope = TRUE) {
  # Create 3D visualization
  Umatrix::showMatrix3D(
    Matrix = umatrix,
    BestMatches = best_matches,
    Cls = cls,
    ClsColors = cls_colors,
    Toroid = toroid,
    Imx = imx,
    BmSize = bm_size,
    RemoveOcean = remove_ocean,
    ShowAxis = show_axis,
    SmoothSlope = smooth_slope
  )

  # Create directory if it doesn't exist
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }

  # Save the 3D visualization
  full_path <- file.path(output_dir, filename)
  snapshot3d(full_path, "png")
  cat(sprintf("3D visualization saved to '%s'\n", full_path))
}

#' Function to save U-matrix related files if file output is enabled
#'
#' @param umx U-matrix object
#' @param pmx P-matrix object
#' @param imx Island mask
#' @param ustar_mx Ustar-matrix object
#' @param bm Best matches
#' @param weights ESOM weights
#' @param cls Class assignments
#' @param prefix File prefix
#' @param output_dir Output directory
#' @examples
#' save_umatrix_files(UmxMouse, PmxMouse, ImxMouse, UstarmxMouse, BMUmxMouse, WTsUmxMouse, MouseCls, "UmxMouse", "output/")
save_umatrix_files <- function(umx, pmx, imx, ustar_mx, bm, weights, cls, prefix, output_dir) {
  # Create directory if it doesn't exist
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }

  # Save files
  Umatrix::WriteBM(FileName = paste0(prefix, ".bm"), BestMatches = bm, OutDirectory = output_dir)
  Umatrix::WriteUMX(FileName = paste0(prefix, ".umx"), UMatrix = umx, OutDirectory = output_dir)
  Umatrix::WriteUMX(FileName = paste0("Pmx", prefix, ".umx"), UMatrix = pmx, OutDirectory = output_dir)
  Umatrix::WriteUMX(FileName = paste0("Umxstar", prefix, ".umx"), UMatrix = ustar_mx, OutDirectory = output_dir)
  Umatrix::WriteIMX(FileName = paste0("Imx", prefix, ".imx"), MapMask = imx, OutDirectory = output_dir)
  Umatrix::WriteWTS(FileName = paste0(prefix, ".wts"), wts = weights, OutDirectory = output_dir)
  Umatrix::WriteCLS(FileName = paste0(prefix, ".cls"), Cls = cls, OutDirectory = output_dir)

  cat(sprintf("U-matrix files saved with prefix '%s' to directory '%s'\n", prefix, output_dir))
}

#' ============================
#' Main script execution
#' ============================

# Read processed data
cat("Reading processed data...\n")
processed_data <- read.csv(input_file, row.names = 1)
cat(sprintf("Data dimensions: %d rows by %d columns\n",
            nrow(processed_data), ncol(processed_data)))

# Extract or create class labels (if not available, assume single class)
if ("Class" %in% colnames(processed_data)) {
  cat("Using 'Class' column as class labels...\n")
  data_cls <- as.integer(as.factor(processed_data$Class))
  processed_data <- processed_data[, !colnames(processed_data) %in% "Class"]
} else {
  cat("No class labels found, assuming single class...\n")
  data_cls <- rep(1, nrow(processed_data))
}

# Convert to percentage if enabled
if (enable_percent_conversion) {
  cat("Converting data to percentages...\n")
  data_for_esom <- toPercent(processed_data)
} else {
  cat("Using processed data without percentage conversion...\n")
  data_for_esom <- processed_data
}

# Train ESOM and generate U-matrix
cat("Training ESOM and generating U-matrix...\n")
set.seed(random_seed)

# Train the ESOM
umx_result <- Umatrix::iEsomTrain(as.matrix(data_for_esom), Cls = data_cls)

# Generate Island mask
cat("Generating Island mask...\n")
imx_result <- Umatrix::iUmapIsland(Umatrix = umx_result$Umatrix,
                                   BestMatches = umx_result$BestMatches,
                                   Cls = data_cls)

# Generate P-matrix and U*-matrix
cat("Generating P-matrix and U*-matrix...\n")
pmx_result <- Umatrix::iUstarmatrix(Data = data_for_esom,
                                    Weights = umx_result$Weights,
                                    Lines = esom_lines,
                                    Columns = esom_columns)

# Classify points based on U-matrix
cat("Performing classification based on U-matrix...\n")
umx_cls_result <- Umatrix::iClassification(Umatrix = umx_result$Umatrix,
                                           BestMatches = umx_result$BestMatches,
                                           Imx = imx_result$Imx)

# Create visualizations if enabled
if (enable_plots) {
  cat("Creating U-matrix visualization...\n")
  umx_plot <- Umatrix::plotMatrix(
    Matrix = umx_result$Umatrix,
    BestMatches = umx_result$BestMatches,
    Cls = data_cls,
    ClsColors = rainbow(length(unique(data_cls))),
    Toroid = esom_toroid,
    BmSize = 5,
    RemoveOcean = TRUE,
    ColorStyle = "Umatrix",
    Imx = imx_result$Imx
  ) + theme_light() +
    guides(color = "none", fill = "none") +
    labs(title = "U-matrix visualization")

  cat("Creating P-matrix visualization...\n")
  pmx_plot <- Umatrix::plotMatrix(
    Matrix = pmx_result$PMatrix,
    Toroid = esom_toroid,
    Cls = data_cls,
    ClsColors = rainbow(length(unique(data_cls))),
    BmSize = 5,
    RemoveOcean = TRUE,
    ColorStyle = "Pmatrix",
    Imx = imx_result$Imx
  ) + theme_light() +
    guides(color = "none", fill = "none") +
    labs(title = "P-matrix visualization")

  # Display plots
  cat("Displaying plots...\n")
  print(umx_plot)
  print(pmx_plot)

  # Create and save 3D visualization if enabled
  if (enable_3d_visualization) {
    cat("Creating 3D visualization...\n")
    create_save_3d_visualization(
      umatrix = umx_result$Umatrix,
      best_matches = umx_result$BestMatches,
      cls = data_cls,
      cls_colors = rainbow(length(unique(data_cls))),
      toroid = esom_toroid,
      imx = imx_result$Imx,
      filename = paste0(output_prefix, "_3D.png"),
      output_dir = output_dir
    )
  }
}

# Save U-matrix files if enabled
if (enable_file_output) {
  cat("Saving U-matrix files...\n")
  save_umatrix_files(
    umx = umx_result$Umatrix,
    pmx = pmx_result$PMatrix,
    imx = imx_result$Imx,
    ustar_mx = pmx_result$UstarMatrix,
    bm = umx_result$BestMatches,
    weights = umx_result$Weights,
    cls = umx_cls_result$Cls,
    prefix = output_prefix,
    output_dir = output_dir
  )
}

# Summary statistics
cat("\n==== Summary of ESOM Training ====\n")
cat(sprintf("Input data dimensions: %d rows x %d columns\n",
            nrow(data_for_esom), ncol(data_for_esom)))
cat(sprintf("Number of classes: %d\n", length(unique(data_cls))))
cat(sprintf("ESOM grid dimensions: %d x %d\n", esom_lines, esom_columns))
cat(sprintf("Toroid topology: %s\n", ifelse(esom_toroid, "Yes", "No")))
cat("ESOM training complete!\n")