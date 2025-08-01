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
# Handle working directory setting
tryCatch({
  if (exists("rstudioapi::getSourceEditorContext")) {
    setwd(dirname(rstudioapi::getSourceEditorContext()$path))
  }
}, error = function(e) {
  message("Unable to set working directory automatically. Please set it manually if needed.")
})

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
input_file <- "test_processed.csv"   # Input processed data file (from previous script)
class_name <- "Species"                 # Name of classes column
output_dir <- "umatrix_output/"              # Output directory for U-matrix files
output_prefix <- "Umx"                  # Prefix for output files

#' @section ESOM parameters:
esom_lines <- 50                          # Number of lines in the ESOM grid
esom_columns <- 80                        # Number of columns in the ESOM grid
esom_toroid <- TRUE                       # Whether to use a toroid topology
esom_epochs <- 20                         # Number of training epochs
enable_plots <- TRUE                      # Whether to create and display plots
enable_percent_conversion <- TRUE         # Whether to convert data to percentage
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

#' Function to create a 3D visualization of the U-matrix without automatic saving
#'
#' @param umatrix U-matrix object
#' @param best_matches Best matches of the ESOM
#' @param cls Class labels
#' @param cls_colors Colors for the classes
#' @param toroid Whether the ESOM has toroid topology
#' @param imx Island mask
#' @param bm_size Size of the best match points
#' @param remove_ocean Whether to remove the "ocean" areas
#' @param show_axis Whether to show the axes
#' @param smooth_slope Whether to smooth the slopes
#' @return Message instructing how to save the visualization
#' @examples
#' create_3d_visualization(Umx, BMUmx, Cls, rainbow(2), TRUE, Imx)
create_3d_visualization <- function(umatrix, best_matches, cls, cls_colors, toroid, imx,
                                    bm_size = 0.6, remove_ocean = TRUE,
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

  # Return instruction message
  cat("3D visualization created.\n")
  cat("Manually adjust the view to your preferred perspective using the mouse:\n")
  cat(" - Left-click and drag: Rotate\n")
  cat(" - Right-click and drag: Zoom\n")
  cat(" - Middle-click and drag: Pan\n")
  cat("\nWhen satisfied with the view, save the visualization using:\n")
  cat("snapshot3d(\"path/to/output_dir/filename.png\", \"png\")\n")

  return(invisible(NULL))
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
#' @param overwrite Logical; if TRUE, existing files will be overwritten (default FALSE)
#' @examples
#' save_umatrix_files(Umx, Pmx, Imx, Ustarmx, BMUmx, WTsUmx, Cls, "Umx", "output/", overwrite = FALSE)
save_umatrix_files <- function(
  umx, pmx, imx, ustar_mx, bm, weights, cls, prefix, output_dir, overwrite = FALSE
) {
  # Create directory if it doesn't exist
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  # Main U-matrix file path
  umx_file <- file.path(output_dir, paste0(prefix, ".umx"))
  
  # Check if file exists and handle overwrite logic
  if (file.exists(umx_file) && !overwrite) {
    cat(sprintf(
      "A U-matrix file named '%s' already exists in '%s'.\nSkipping file creation. Set overwrite = TRUE to overwrite.\n",
      paste0(prefix, ".umx"), output_dir
    ))
    return(invisible(NULL))
  }
  
  # Save files
  dbt.DataIO::WriteBM(FileName = paste0(prefix, ".bm"), BestMatches = bm, OutDirectory = output_dir)
  dbt.DataIO::WriteUMX(FileName = paste0(prefix, ".umx"), UMatrix = umx, OutDirectory = output_dir)
  dbt.DataIO::WriteUMX(FileName = paste0("Pmx", prefix, ".umx"), UMatrix = pmx, OutDirectory = output_dir)
  dbt.DataIO::WriteUMX(FileName = paste0("Umxstar", prefix, ".umx"), UMatrix = ustar_mx, OutDirectory = output_dir)
  dbt.DataIO::WriteIMX(FileName = paste0("Imx", prefix, ".imx"), MapMask = imx, OutDirectory = output_dir)
  dbt.DataIO::WriteWTS(FileName = paste0(prefix, ".wts"), wts = weights, OutDirectory = output_dir)
  dbt.DataIO::WriteCLS(FileName = paste0(prefix, ".cls"), Cls = cls, OutDirectory = output_dir)
  
  if (file.exists(umx_file) && overwrite) {
    cat(sprintf(
      "An existing U-matrix file named '%s' was overwritten in '%s'.\n",
      paste0(prefix, ".umx"), output_dir
    ))
  } else {
    cat(sprintf(
      "U-matrix files saved with prefix '%s' to directory '%s'.\n",
      prefix, output_dir
    ))
  }
}

#' Function to read U-matrix related files from disk
#'
#' @param prefix File prefix (same as used in saving)
#' @param output_dir Output directory (same as used in saving)
#' @return A list with elements: umx, pmx, imx, ustar_mx, bm, weights, cls
#' @examples
#' results <- read_umatrix_files("Umx", "output/")
read_umatrix_files <- function(prefix, output_dir) {
  # Helper to construct full file paths
  file_path <- function(name) file.path(output_dir, name)
  
  # Read files
  umx      <- dbt.DataIO::ReadUMX(file_path(paste0(prefix, ".umx")))
  pmx      <- dbt.DataIO::ReadUMX(file_path(paste0("Pmx", prefix, ".umx")))
  ustar_mx <- dbt.DataIO::ReadUMX(file_path(paste0("Umxstar", prefix, ".umx")))
  imx      <- dbt.DataIO::ReadIMX(file_path(paste0("Imx", prefix, ".imx")))
  bm       <- dbt.DataIO::ReadBM(file_path(paste0(prefix, ".bm")))
  weights  <- dbt.DataIO::ReadWTS(file_path(paste0(prefix, ".wts")))
  cls      <- dbt.DataIO::ReadCLS(file_path(paste0(prefix, ".cls")))
  
  # Return as a named list
  list(
    umx      = umx,
    pmx      = pmx,
    ustar_mx = ustar_mx,
    imx      = imx,
    bm       = bm,
    weights  = weights,
    cls      = cls
  )
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
if (class_name %in% colnames(processed_data)) {
  cat(paste0("Using ", class_name,  " column as class labels...\n"))
  data_cls <- as.integer(as.factor(processed_data[[class_name]]))
  processed_data <- processed_data[, !colnames(processed_data) %in% class_name]
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
                                    Cls = data_cls, 
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

  # Create 3D visualization if enabled
  if (enable_3d_visualization) {
    cat("Creating 3D visualization...\n")
    # Prepare the output path
    if (!dir.exists(output_dir)) {
      dir.create(output_dir, recursive = TRUE)
    }

    filename_3d <- paste0(output_prefix, "_3D.png")
    full_path_3d <- file.path(output_dir, filename_3d)

    # Create the 3D visualization
    create_3d_visualization(
      umatrix = umx_result$Umatrix,
      best_matches = umx_result$BestMatches,
      cls = data_cls,
      cls_colors = rainbow(length(unique(data_cls))),
      toroid = esom_toroid,
      imx = imx_result$Imx,
      bm_size = 0.6,
      remove_ocean = TRUE
    )

    cat("\n================================\n")
    cat("3D visualization is now open for manual adjustment.\n")
    cat("- Left-click and drag: Rotate\n")
    cat("- Right-click and drag: Zoom\n")
    cat("- Middle-click and drag: Pan\n")
    cat("\nWhen you are satisfied with the view, save it using:\n")
    cat(sprintf("snapshot3d(\"%s\", \"png\")\n", full_path_3d))
    cat("\nAfter saving, you can close the 3D window manually.\n")
    cat("================================\n")

    # NOTE: The saving should be done manually by the user after adjusting the view
    # The automatic saving has been removed from this section
  }
}

# Save the 3D visualization if enabled
# Run only when 3D visualization if enabled
# snapshot3d("umatrix_output//Umx_3D.png", "png")
# End run only when 3D visualization if enabled


# Save U-matrix files if enabled
if (enable_file_output) {
  cat("Saving U-matrix files...\n")
  save_umatrix_files(
    umx      = umx_result$Umatrix,
    pmx      = pmx_result$PMatrix,
    imx      = imx_result$Imx,
    ustar_mx = pmx_result$UstarMatrix,
    bm       = umx_result$BestMatches,
    weights  = umx_result$Weights,
    cls      = umx_cls_result$Cls,
    prefix   = output_prefix,
    output_dir = output_dir,
    overwrite = FALSE  
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


# Calculate Delaunay radius for density-based visualization
cat("Calculating Delaunay radius for density visualization...\n")
Radius_umx <- Umatrix::calculate_Delauny_radius(
  Data = as.matrix(data_for_esom),
  BestMatches = umx_result$BestMatches,
  Columns = esom_columns,
  Lines = esom_lines,
  Toroid = esom_toroid
)

# Store the radius data
RadiusData <- Radius_umx$RadiusByEM

# Summary for density radius calculation
cat(sprintf("Delaunay radius calculation complete. Min: %.2f, Max: %.2f\n",
            min(RadiusData), max(RadiusData)))

# If plotting is enabled, create a density plot of neighborhood distances
if (enable_plots) {
  cat("Creating density plot of neighborhood distances...\n")

  # Create a data frame with the neighborhood distances
  neighbor_distances_df <- data.frame(distance = Radius_umx$neighbourDistances)

  # Create the density plot with ggplot2
  density_plot <- ggplot(neighbor_distances_df, aes(x = distance)) +
    geom_density(fill = "lightblue", alpha = 0.7) +
    geom_vline(xintercept = RadiusData, linetype = "dashed", color = "salmon", size = 1) +
    annotate("text", x = RadiusData * 1.1, y = Inf, label = paste("Radius =", round(RadiusData, 2)),
             vjust = 2, hjust = 0, color = "salmon") +
    labs(title = "Density of Neighborhood Distances",
         subtitle = paste("Calculated radius by EM:", round(RadiusData, 4)),
         x = "Distance",
         y = "Density") +
    theme_light() +
    theme(plot.title = element_text(hjust = 0.5),
          plot.subtitle = element_text(hjust = 0.5))

  print(density_plot)

  # Save the plot if file output is enabled
  if (enable_file_output) {
    density_plot_file <- file.path(output_dir, paste0(output_prefix, "_density_plot.svg"))
    ggsave(density_plot_file, density_plot, width = 8, height = 6)
    cat(sprintf("Density plot saved to %s\n", density_plot_file))
  }
}

# Save the radius data if file output is enabled
if (enable_file_output) {
  cat("Saving density radius data...\n")

  # Make sure the output directory exists
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
    cat(sprintf("Created output directory: %s\n", output_dir))
  }

  radius_file <- file.path(output_dir, paste0(output_prefix, "_radius.csv"))

  # Create a proper data frame with the radius value
  radius_df <- data.frame(RadiusByEM = RadiusData)

  # Write to CSV without row names and with explicit file connection
  tryCatch({
    write.csv(radius_df, file = radius_file, row.names = FALSE)
    cat(sprintf("Density radius data saved to %s\n", radius_file))
  }, error = function(e) {
    cat(sprintf("Error saving radius data: %s\n", e$message))
  })
}
