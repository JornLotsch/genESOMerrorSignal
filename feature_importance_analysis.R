# Feature_Importance_Analysis.R
# Script for feature importance assessment using Boruta

# Handle working directory setting
tryCatch({
  if (exists("rstudioapi::getSourceEditorContext")) {
    setwd(dirname(rstudioapi::getSourceEditorContext()$path))
  }
}, error = function(e) {
  message("Unable to set working directory automatically. Please set it manually if needed.")
})

# Source the synthetic data generation script
source("generate_synthetic_data.R")

# Load required libraries
library(Boruta)
library(ggplot2)
library(reshape2)
library(caret)
library(parallel)
library(pbmcapply)
library(cowplot) # For plot_grid function

# Set up parameters
seed <- 42    # Main seed for reproducibility
nIter <- 100
list.of.seeds <- 1:nIter + seed - 1  # Creates seeds from 42 to 141
nProc <- detectCores() - 1  # Use all but one core
enable_plots <- TRUE
enable_file_output <- TRUE
output_dir <- "results"
output_prefix <- "Feature_Importance"
pfad_o <- output_dir
pfad_r <- "/"

# Configuration settings for data
input_file <- "test.csv"                # Input data file
class_name <- "Species"                 # Name of classes column
output_file <- "test_processed.csv"     # Output processed data file

# Generation parameters - control synthetic data generation
generation_multipliers <- c(1, 5)  # Can be extended to any number of multipliers
base_generation_rate <- 1          # Base rate of synthetic samples per original sample


# Helper function to check for integer(0)
is.integer0 <- function(x) {
  is.integer(x) && length(x) == 0
}

# Load the processed data
cat(sprintf("Loading data from %s...\n", output_file))
if (file.exists(output_file)) {
  data_df <- read.csv(output_file)

  # Ensure class column is a factor
  data_df[[class_name]] <- as.factor(data_df[[class_name]])

  cat(sprintf("Loaded data with %d rows and %d columns\n",
              nrow(data_df), ncol(data_df)))
  cat("Class distribution:\n")
  print(table(data_df[[class_name]]))
} else {
  stop(sprintf("Data file %s not found", output_file))
}

# Load the radius data
cat("Loading density radius data...\n")
if (file.exists(file.path(output_dir, paste0("Umx_radius.csv")))) {
  radius_df <- read.csv(file.path(output_dir, paste0("Umx_radius.csv")))
  RadiusData <- radius_df$RadiusByEM
  cat(sprintf("Loaded radius: %f\n", RadiusData))
} else {
  # If file doesn't exist, check if RadiusData exists in environment
  if (exists("RadiusData")) {
    cat(sprintf("Using radius from environment: %f\n", RadiusData))
  } else {
    stop("Radius data not found. Please run radius calculation first.")
  }
}

# Handle data splitting only if needed
needsDataSplit <- any(grepl("reduced", DataSetSizes))

# Optionally split data based on dataset types
if (needsDataSplit) {
  cat("Splitting data into training, test and validation sets...\n")
  TestDataTrainingTestValidation <- opdisDownsampling::opdisDownsampling(
    Data = within(data_df, rm(get(class_name))),
    Cls = data_df[[class_name]],
    Size = 0.8 * nrow(data_df),
    Seed = seed,
    nTrials = 10000,
    MaxCores = nProc
  )

  data_TrainingTest <- data_df[rownames(data_df) %in% TestDataTrainingTestValidation$ReducedInstances, ]
  cat("Training/Test set class distribution:\n")
  print(table(data_TrainingTest[[class_name]]))

  data_Validation <- data_df[!rownames(data_df) %in% TestDataTrainingTestValidation$ReducedInstances, ]
  cat("Validation set class distribution:\n")
  print(table(data_Validation[[class_name]]))
} else {
  cat("Skipping data splitting (not required for current dataset types)\n")
}

#################################### Evaluate variable importance ########################################################################

# Define dataset types dynamically based on generation multipliers
DataSetSizes <- c("original", "engineered_0")
for (mult in generation_multipliers) {
  DataSetSizes <- c(DataSetSizes, paste0("augmented_", mult, "_engineered"))
}
cat("Dataset types to analyze:", paste(DataSetSizes, collapse=", "), "\n")

# To store significant variables from each dataset
significant_variables <- list()

Test_VarImps <- lapply(DataSetSizes, function(DatasetNr) {
  cat(sprintf("Processing dataset type: %s\n", DatasetNr))

  # Extract generation multiplier from dataset name if applicable
  gen_multiplier <- 1
  if (grepl("augmented_", DatasetNr)) {
    gen_multiplier <- as.numeric(gsub("augmented_([0-9]+)_.*", "\\1", DatasetNr))
    cat(sprintf("Using generation multiplier: %d\n", gen_multiplier))
  }

  Imps_repeated <- pbmcapply::pbmclapply(list.of.seeds, function(x) {
    cat(sprintf("  Running with seed: %d\n", x))

    # Select or generate appropriate dataset based on DatasetNr
    if (DatasetNr == "original") {
      data_actual <- data_df
    } else if (DatasetNr == "reduced") {
      # Use the already split data
      if (!exists("data_TrainingTest")) {
        stop("'reduced' dataset requested but data splitting wasn't performed")
      }
      data_actual <- data_TrainingTest
      cat(sprintf("Using reduced dataset with %d rows\n", nrow(data_actual)))
    } else if (DatasetNr == "engineered_0") {
      set.seed(x)
      # Create engineered dataset with permuted features
      feature_cols <- setdiff(names(data_df), class_name)
      data_actual <- data_df

      # Add permuted features
      for (col in feature_cols) {
        data_actual[[paste0(col, "_permuted")]] <- sample(data_df[[col]])
      }
    } else if (grepl("augmented_.*_engineered", DatasetNr)) {
      set.seed(x)
      # Create engineered dataset with permuted features
      feature_cols <- setdiff(names(data_df), class_name)
      data_engineered <- data_df

      # Add permuted features
      for (col in feature_cols) {
        data_engineered[[paste0(col, "_permuted")]] <- sample(data_df[[col]])
      }

      set.seed(x)
      data_generated <- generate_synthetic_data(
        Data = within(data_engineered, rm(get(class_name))),
        density_radius = RadiusData,
        gen_per_data = gen_multiplier * base_generation_rate,
        Cls = data_engineered[[class_name]]
      )

      data_actual <- rbind.data.frame(
        cbind.data.frame(data_generated$original_data, class = data_generated$original_classes),
        cbind.data.frame(data_generated$generated_data, class = data_generated$generated_classes)
      )
      # Rename class column to match original
      names(data_actual)[names(data_actual) == "class"] <- class_name
    } else if (grepl("augmented_", DatasetNr)) {
      set.seed(x)
      data_generated <- generate_synthetic_data(
        Data = within(data_df, rm(get(class_name))),
        density_radius = RadiusData,
        gen_per_data = gen_multiplier * base_generation_rate,
        Cls = data_df[[class_name]]
      )

      data_actual <- rbind.data.frame(
        cbind.data.frame(data_generated$original_data, class = data_generated$original_classes),
        cbind.data.frame(data_generated$generated_data, class = data_generated$generated_classes)
      )
      # Rename class column to match original
      names(data_actual)[names(data_actual) == "class"] <- class_name
    } else if (grepl("augmented_", DatasetNr)) {
      set.seed(x)
      data_generated <- generate_synthetic_data(
        Data = within(data_df, rm(get(class_name))),
        density_radius = RadiusData,
        gen_per_data = gen_multiplier * base_generation_rate,
        Cls = data_df[[class_name]]
      )

      data_actual <- rbind.data.frame(
        cbind.data.frame(data_generated$original_data, class = data_generated$original_classes),
        cbind.data.frame(data_generated$generated_data, class = data_generated$generated_classes)
      )
      # Rename class column to match original
      names(data_actual)[names(data_actual) == "class"] <- class_name
    }

    # Split data for Boruta analysis
    set.seed(x)
    inTraining <- createDataPartition(data_actual[[class_name]], p = .67, list = FALSE)
    training <- data_actual[inTraining, ]

    # Create formula for Boruta with dynamic class column name
    boruta_formula <- as.formula(paste(class_name, "~ ."))

    # Run Boruta feature selection
    Test_Boruta_actual <- Boruta::Boruta(boruta_formula, training, pValue = 0.0001, maxRuns = 100)

    # Process Boruta results
    BorutaStats <- Boruta::attStats(Test_Boruta_actual)
    BorutaStats$Var <- rownames(BorutaStats)

    # Prepare importance data for plotting
    dfVarimp_Test_Boruta_actual_long <- reshape2::melt(Test_Boruta_actual$ImpHistory)

    return(list(
      Test_Boruta_actual = Test_Boruta_actual,
      BorutaStats = BorutaStats,
      dfVarimp_Test_Boruta_actual_long = dfVarimp_Test_Boruta_actual_long
    ))
  }, mc.cores = 0.5 * nProc)

  # Process results from all repetitions
  BorutaStats_all <- do.call(rbind.data.frame, lapply(Imps_repeated, function(x) x$BorutaStats))

  # Process permuted variables if they exist
  if (!is.integer0(grep("permuted", BorutaStats_all$Var))) {
    BorutaStats_all_permuted_notrejected <- BorutaStats_all[grep("permuted", BorutaStats_all$Var), ]
    BorutaStats_all_permuted_notrejected <- BorutaStats_all_permuted_notrejected[BorutaStats_all_permuted_notrejected$decision == "Confirmed", ]

    BorutaStats_all_notrejected <- BorutaStats_all[-grep("permuted", BorutaStats_all$Var), ]
    BorutaStats_all_notrejected <- BorutaStats_all_notrejected[BorutaStats_all_notrejected$decision == "Confirmed", ]
  } else {
    BorutaStats_all_notrejected <- BorutaStats_all
    BorutaStats_all_notrejected <- BorutaStats_all_notrejected[BorutaStats_all_notrejected$decision == "Confirmed", ]
  }

  # Count feature selection frequency
  SelectedTrue <- data.frame(table(BorutaStats_all_notrejected$Var))

  if (!is.integer0(grep("permuted", BorutaStats_all$Var))) {
    SelectedPermuted <- data.frame(table(BorutaStats_all_permuted_notrejected$Var))
    SelectedPermuted$Var2 <- gsub("_permuted", "", SelectedPermuted$Var1)
    SelectedPermuted <- SelectedPermuted[SelectedPermuted$Var2 %in% SelectedTrue$Var1, ]
  }

  # Create data frame with variable selection statistics
  if (!is.integer0(grep("permuted", BorutaStats_all$Var))) {
    dfVars <- data.frame(Var = unique(BorutaStats_all$Var[-grep("permuted", BorutaStats_all$Var)]))
  } else {
    dfVars <- data.frame(Var = unique(BorutaStats_all$Var))
  }

  rownames(dfVars) <- dfVars$Var
  dfVars$SelectedTrue <- SelectedTrue$Freq[match(dfVars$Var, SelectedTrue$Var1)]

  if (!is.integer0(grep("permuted", BorutaStats_all$Var))) {
    dfVars$SelectedPermuted <- SelectedPermuted$Freq[match(dfVars$Var, SelectedPermuted$Var2)]
  } else {
    dfVars$SelectedPermuted <- NA
  }

  dfVars[is.na(dfVars)] <- 0
  dfVars$SelectedTrueCorr <- dfVars$SelectedTrue - dfVars$SelectedPermuted
  dfVars$SelectedTrueCorrRel <- (dfVars$SelectedTrue - dfVars$SelectedPermuted) * (dfVars$SelectedTrue + dfVars$SelectedPermuted)

  # Bootstrapping to determine significance thresholds
  limtFreq <- NA

  if (!is.integer0(grep("permuted", BorutaStats_all$Var))) {
    nBootstrap <- 100000
    a <- dfVars$SelectedPermuted
    b <- dfVars$SelectedTrue

    set.seed(seed)
    a_sample <- sample(a, nBootstrap, replace = TRUE)
    set.seed(seed)
    b_sample <- sample(b, nBootstrap, replace = TRUE)

    FCs_ba_bootstrap <- b_sample - a_sample
    limtFreq <- quantile(FCs_ba_bootstrap, probs = 0.95)

    # Store variables that exceed the threshold
    important_vars <- dfVars$Var[which(dfVars$SelectedTrueCorr > limtFreq)]
    significant_variables[[DatasetNr]] <- important_vars
  } else {
    # If no permuted variables, use positive SelectedTrueCorr as significant
    important_vars <- dfVars$Var[which(dfVars$SelectedTrueCorr > 0)]
    significant_variables[[DatasetNr]] <- important_vars
  }

  # Combine importance data from all repetitions
  dfVarimp_Test_Boruta_actual_long_all <- do.call(rbind.data.frame,
                                                  lapply(Imps_repeated, function(x) x$dfVarimp_Test_Boruta_actual_long))

  # Create frequency plot showing difference (True - Permuted) with consistent color scheme
  varSelectionBar <- ggplot(data = dfVars) +
    geom_bar(
      aes(
        y = reorder(Var, SelectedTrueCorr),
        x = SelectedTrueCorr,
        fill = SelectedTrueCorr,
        color = SelectedTrueCorr > 0  # Use boolean to determine outline color
      ),
      stat = "identity"
    ) +
    # Add zero line
    geom_vline(xintercept = 0, linetype = "dashed", color = "salmon") +
    # Apply the same color gradient as the circular plot
    scale_fill_gradient2(
      "Selection Frequency",
      low = "chartreuse",     # Green for negative values
      mid = "ghostwhite",     # White for values near 0
      high = "salmon",        # Salmon to red gradient for positive values
      midpoint = 0
    ) +
    # Add color scale for outlines
    scale_color_manual(
      values = c("chartreuse3", "salmon"),
      guide = "none"  # Don't show a legend for outline colors
    ) +
    # Customize labels
    labs(
      title = "Variable selection frequency",
      x = "Times selected more than permuted copy",
      y = NULL
    ) +
    # Apply theme
    theme_light() +
    theme(
      legend.position = "bottom",
      legend.title = element_text(hjust = 0.5)
    )

  # Create circular barplot with improved visibility for values near zero
  varSelectionRadial <- ggplot(data = dfVars) +
    # Make custom panel grid
    geom_hline(
      aes(yintercept = y),
      data.frame(y = seq(0, max(dfVars$SelectedTrueCorr, na.rm = TRUE) + 5, by = 5)),
      color = "lightgrey"
    ) +
    # Add bars with colored outlines based on sign
    geom_col(
      aes(
        x = reorder(Var, SelectedTrueCorr),
        y = SelectedTrueCorr,
        fill = SelectedTrueCorr,
        color = SelectedTrueCorr > 0  # Use boolean to determine outline color
      ),
      position = "dodge2",
      show.legend = TRUE,
      alpha = 0.9
    ) +
    # Add zero line
    geom_hline(yintercept = 0, linetype = "dashed", color = "salmon", size = 1) +
    # Make it circular!
    coord_polar() +
    # Scale y axis so bars don't start in the center
    scale_y_continuous(
      limits = c(min(min(dfVars$SelectedTrueCorr, na.rm = TRUE) - 5, -5),
                 max(dfVars$SelectedTrueCorr, na.rm = TRUE) + 5),
      expand = c(0, 0)
    ) +
    # Use ghostwhite for mid values
    scale_fill_gradient2(
      "Selection Frequency",
      low = "chartreuse",     # Green for negative values
      mid = "ghostwhite",     # White for values near 0
      high = "salmon",        # Salmon to red gradient for positive values
      midpoint = 0
    ) +
    # Add color scale for outlines
    scale_color_manual(
      values = c("chartreuse3", "salmon"),
      guide = "none"  # Don't show a legend for outline colors
    ) +
    # Customize guides
    guides(
      fill = guide_colorbar(
        barwidth = 15, barheight = 0.5, title.position = "top", title.hjust = 0.5
      )
    ) +
    # Customize theme
    theme_minimal() +
    theme(
      # Remove axis ticks and text
      axis.title = element_blank(),
      axis.ticks = element_blank(),
      axis.text.y = element_blank(),
      # Use gray text for the variable names
      axis.text.x = element_text(color = "gray12", size = 9),
      # Move the legend to the bottom
      legend.position = "bottom",
      panel.grid.major.x = element_blank(),
      panel.grid.minor = element_blank()
    ) +
    # Add title
    labs(
      title = "Variable Selection Frequency",
      subtitle = "Times selected more than permuted copy"
    )

  # Add annotations for scale
  # First, calculate reasonable positions based on your data
  max_y <- max(dfVars$SelectedTrueCorr, na.rm = TRUE)
  scale_positions <- seq(0, max_y, by = 5)

  # Add scale annotations
  for(i in seq_along(scale_positions)) {
    if(scale_positions[i] > 0) {  # Only annotate positive values
      varSelectionRadial <- varSelectionRadial +
        annotate(
          "text",
          x = length(dfVars$Var) + 0.5,  # Position at the end of the variables
          y = scale_positions[i],
          label = as.character(scale_positions[i]),
          color = "gray12",
          size = 3
        )
    }
  }


  # Save individual plots if enabled
  if (enable_plots && enable_file_output) {
    # Create output directory if it doesn't exist
    if (!dir.exists(output_dir)) {
      dir.create(output_dir, recursive = TRUE)
    }

    # Save frequency plot (standard bar plot)
    freq_plot_file <- file.path(output_dir, paste0(output_prefix, "_", DatasetNr, "_frequency.svg"))
    ggsave(freq_plot_file, varSelectionBar, width = 10, height = 8)

    # Save circular frequency plot
    freq_circular_plot_file <- file.path(output_dir, paste0(output_prefix, "_", DatasetNr, "_frequency_circular.svg"))
    ggsave(freq_circular_plot_file, varSelectionRadial, width = 10, height = 10)
  }

  # Return results
  return(list(
    importance_data_long = dfVarimp_Test_Boruta_actual_long_all,
    varSelectionBar = varSelectionBar,
    varSelectionRadial = varSelectionRadial,
    limtFreq = limtFreq,
    dfVars = dfVars,
    important_vars = important_vars
  ))
})

# Assign names to results list
names(Test_VarImps) <- DataSetSizes

# Create combined plot of frequency differences
cat("Creating combined plot...\n")

# Create a list of plots for the grid
plot_list <- lapply(DataSetSizes, function(ds) {
  title_text <- if(ds == "original") {
    "Normal data: Original"
  } else if(ds == "engineered_0") {
    "Normal data: Engineered 0"
  } else if(grepl("augmented_", ds)) {
    gen_multiplier <- as.numeric(gsub("augmented_([0-9]+)_.*", "\\1", ds))
    sprintf("Normal data: Augmented %d, engineered", gen_multiplier)
  } else {
    ds
  }

  Test_VarImps[[ds]]$varSelectionBar + labs(title = title_text)
})

# Create plot grid
pTest_Varfreqs_all <- do.call(plot_grid, c(
  plot_list,
  list(labels = LETTERS[1:length(DataSetSizes)],
       nrow = 1,
       align = "h",
       axis = "tb")
))

print(pTest_Varfreqs_all)

# Save combined plot
if (enable_file_output) {
  svg_file <- file.path(pfad_o, pfad_r, "pnoEffect_Varfreqs_aug100.svg")
  ggsave(filename = svg_file, plot = pTest_Varfreqs_all, width = 20, height = 10, limitsize = FALSE)

  # Also save as PNG for easier viewing
  png_file <- file.path(pfad_o, pfad_r, "pnoEffect_Varfreqs_aug100.png")
  ggsave(filename = png_file, plot = pTest_Varfreqs_all, width = 20, height = 10, limitsize = FALSE, dpi = 300)
}

# Create table of significant variables (those that exceed the threshold)
cat("Creating significant variables table...\n")

# Get all unique variables across all datasets
all_vars <- unique(unlist(sapply(DataSetSizes, function(ds) {
  Test_VarImps[[ds]]$important_vars
})))

# Create matrix with 1/0 indicating if variable is significant in each dataset
significant_matrix <- matrix(0, nrow = length(all_vars), ncol = length(DataSetSizes))
rownames(significant_matrix) <- all_vars
colnames(significant_matrix) <- DataSetSizes

# Fill the matrix
for (i in 1:length(DataSetSizes)) {
  ds <- DataSetSizes[i]
  vars <- Test_VarImps[[ds]]$important_vars
  significant_matrix[vars, i] <- 1
}

# Convert to data frame and add count of datasets where variable is significant
significant_vars_df <- as.data.frame(significant_matrix)
significant_vars_df$SignificantCount <- rowSums(significant_vars_df)

# Sort by significance count (descending)
significant_vars_df <- significant_vars_df[order(-significant_vars_df$SignificantCount, rownames(significant_vars_df)), ]

# Display the results
cat("Significant variables by dataset:\n")
print(significant_vars_df)

# Save to CSV
if (enable_file_output) {
  significant_vars_file <- file.path(output_dir, paste0(output_prefix, "_significant_variables.csv"))
  write.csv(significant_vars_df, file = significant_vars_file, row.names = TRUE)
  cat(sprintf("Significant variables matrix saved to %s\n", significant_vars_file))
}

# Save all results
if (enable_file_output) {
  results_file <- file.path(output_dir, paste0(output_prefix, "_results.RData"))
  save(Test_VarImps, significant_vars_df, generation_multipliers, file = results_file)
  cat(sprintf("All feature importance results saved to %s\n", results_file))
}

cat("Feature importance analysis complete.\n")