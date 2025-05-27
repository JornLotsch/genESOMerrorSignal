#' Data Preparation and Imputation Script
#'
#' This script performs comprehensive data preparation for mouse lipid data including:
#' - Distribution exploration
#' - Automated transformation selection using Tukey's ladder of powers or Box-Cox
#' - Outlier detection and removal using Grubbs' test
#' - Missing value imputation using missForest
#' - Back-transformation of data to original scale

#' Check and install required packages if needed
required_packages <- c("readr", "ggplot2", "ComplexHeatmap", "missForest",
                       "forecast", "outliers", "nortest", "parallel")

for (pkg in required_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg)
  }
  library(pkg, character.only = TRUE)
}

# Check for ABCstats separately as it's on GitHub, not CRAN
if (!requireNamespace("ABCstats", quietly = TRUE)) {
  message("ABCstats package is not installed. Attempting to install from GitHub...")

  # Check if devtools is available, install if needed
  if (!requireNamespace("devtools", quietly = TRUE)) {
    message("Installing devtools package...")
    install.packages("devtools")
    if (!requireNamespace("devtools", quietly = TRUE)) {
      warning("Failed to install devtools. ABCstats installation aborted. Some functionality may be limited.")
      has_abcstats <- FALSE
    } else {
      # Install ABCstats from GitHub
      tryCatch({
        devtools::install_github("huaxuyu/ABCstats")
        if (requireNamespace("ABCstats", quietly = TRUE)) {
          has_abcstats <- TRUE
          message("ABCstats package successfully installed.")
        } else {
          warning("Failed to install ABCstats package. Some functionality may be limited.")
          has_abcstats <- FALSE
        }
      }, error = function(e) {
        warning(sprintf("Failed to install ABCstats package: %s. Some functionality may be limited.", e$message))
        has_abcstats <- FALSE
      })
    }
  } else {
    # Install ABCstats from GitHub
    tryCatch({
      devtools::install_github("huaxuyu/ABCstats")
      if (requireNamespace("ABCstats", quietly = TRUE)) {
        has_abcstats <- TRUE
        message("ABCstats package successfully installed.")
      } else {
        warning("Failed to install ABCstats package. Some functionality may be limited.")
        has_abcstats <- FALSE
      }
    }, error = function(e) {
      warning(sprintf("Failed to install ABCstats package: %s. Some functionality may be limited.", e$message))
      has_abcstats <- FALSE
    })
  }
} else {
  has_abcstats <- TRUE
}

# Handle working directory setting
tryCatch({
  if (exists("rstudioapi::getSourceEditorContext")) {
    setwd(dirname(rstudioapi::getSourceEditorContext()$path))
  }
}, error = function(e) {
  message("Unable to set working directory automatically. Please set it manually if needed.")
})
setwd("/home/joern/.Datenplatte/Joerns Dateien/Aktuell/GenerativeESOM/08AnalyseProgramme/R/genESOMerrorSignal")

#' Configuration settings
#' @section File paths:
input_file <- "test.csv"                # Input data file
class_name <- "Species"                 # Name of classes column
output_file <- "test_processed.csv"     # Output processed data file
output_file_backtransformed <- "test_processed_original_scale.csv"  # Back-transformed data
lambda_file <- "test_BoxCox_lambdas.csv"             # File to save BoxCox lambda values
heatmap_file <- "test_Processed_Data_Heatmap.svg"    # Output heatmap file

#' @section Processing options:
# Data exploration options
enable_plots <- TRUE                     # Whether to create and display plots
enable_distribution_exploration <- TRUE  # Whether to explore data distributions
enable_verbose <- TRUE                   # Whether to print detailed messages

# Transformation options
enable_transformation <- TRUE           # Whether to transform data
transformation_method <- "none"         # "auto", "none", "log10", "sqrt", "reciprocal", "boxcox", "ABCtrans_BC"

# Outlier detection options
enable_outlier_removal <- TRUE          # Whether to detect and remove outliers
outlier_method <- "grubbs"              # Method for outlier detection: "grubbs" or "boxplot"
outlier_p_threshold <- 1E-5             # P-value threshold for Grubbs test
boxplot_coef <- 1.5                     # IQR multiplier for boxplot method
max_outlier_rounds <- 100               # Maximum rounds of outlier detection

# Missing value handling
enable_imputation <- TRUE               # Whether to impute missing values
enable_case_removal <- TRUE             # Whether to remove cases with too many outliers
case_outlier_limit <- 0.5               # Maximum proportion of missing values allowed per case

#' Set up parallel processing
#' @param n_cores Number of CPU cores to use (all but one)
#' @param random_seed Random seed for reproducibility
n_cores <- parallel::detectCores() - 1
random_seed <- 42

#' Helper function for verbose output
#' @param msg Message to print
#' @param ... Additional arguments to pass to cat
verbose <- function(msg, ...) {
  if (enable_verbose) {
    cat(sprintf(msg, ...), "\n")
  }
}

#' Signed log transformation (forward)
#'
#' @param x Numeric vector to transform
#' @param m Base for logarithm (default = 2)
#' @return Transformed vector
slog2Forward <- function(x, m = 2) {
  ms <- as.character(m)
  absX <- abs(x)
  s <- sign(x)
  switch(ms,
    "0" = {
      Val <- log1p(absX)
    },
    "2" = {
      Val <- log2(absX + 1)
    },
    "10" = {
      Val <- log10(absX + 1)
    },
  {
    Val <- log((absX + 1), base = m)
  }
  )
  return(Val * s)
}

#' Signed log transformation (inverse)
#'
#' @param sLog Transformed vector
#' @param m Base for logarithm (default = 2)
#' @return Back-transformed vector
slog2Inverse <- function(sLog, m = 2) {
  ms <- as.character(m)
  s <- sign(sLog)
  switch(ms,
    "0" = {
      Val <- exp(sLog / s) - 1
    },
    "2" = {
      Val <- 2^(sLog / s) - 1
    },
    "10" = {
      Val <- 10^(sLog / s) - 1
    },
  {
    Val <- m^(sLog / s) - 1
  }
  )
  return(Val * s)
}

#' Find the appropriate Tukey transformation based on BoxCox lambda
#'
#' @param x The BoxCox lambda value
#' @return The closest value on Tukey's ladder of powers
find_tukey_trans <- function(x) {
  tukeys_lop <- c(-2, -1, -0.5, 0, 0.5, 1, 2)
  # Find the closest value instead of using findInterval
  idx <- which.min(abs(tukeys_lop - x))
  return(tukeys_lop[idx])
}

#' Apply Tukey transformation based on lambda value
#'
#' @param x Numeric vector to transform
#' @param lambda Transformation parameter from Tukey's ladder
#' @return Transformed vector
apply_tukey_trans <- function(x, lambda) {
  # Convert lambda to numeric for more robust comparison
  lambda_num <- as.numeric(lambda)

  if (is.na(lambda_num)) {
    warning("Invalid lambda value in apply_tukey_trans. Using original data.")
    return(x)
  }

  if (abs(lambda_num - (-2)) < 1e-10) return(1 / (x^2))
  if (abs(lambda_num - (-1)) < 1e-10) return(1 / x)
  if (abs(lambda_num - (-0.5)) < 1e-10) return(1 / sqrt(x))
  if (abs(lambda_num - 0) < 1e-10) return(log10(x))
  if (abs(lambda_num - 0.5) < 1e-10) return(sqrt(x))
  if (abs(lambda_num - 1) < 1e-10) return(x)
  if (abs(lambda_num - 2) < 1e-10) return(x^2)

  # Default case: return original data with warning
  warning("Unrecognized lambda value in apply_tukey_trans. Using original data.")
  return(x)
}

#' Back-transform data based on the transformation that was applied
#'
#' @param x Transformed data
#' @param lambda Transformation parameter
#' @return Back-transformed data
back_transform <- function(x, lambda) {
  # Convert lambda to numeric for more robust comparison
  lambda_num <- as.numeric(lambda)

  if (is.na(lambda_num)) {
    warning("Invalid lambda value in back_transform. Using original data.")
    return(x)
  }

  if (abs(lambda_num - (-2)) < 1e-10) return(1 / sqrt(1 / x))
  if (abs(lambda_num - (-1)) < 1e-10) return(1 / x)
  if (abs(lambda_num - (-0.5)) < 1e-10) return(1 / (x^2))
  if (abs(lambda_num - 0) < 1e-10) return(10^x)
  if (abs(lambda_num - 0.5) < 1e-10) return(x^2)
  if (abs(lambda_num - 1) < 1e-10) return(x)
  if (abs(lambda_num - 2) < 1e-10) return(sqrt(x))

  # Default case: return original data with warning
  warning("Unrecognized lambda value in back_transform. Using original data.")
  return(x)
}

#' Apply Box-Cox transformation with specific lambda
#'
#' @param x Numeric vector to transform
#' @param lambda Box-Cox lambda parameter
#' @return Transformed vector
doBoxCoxTrans <- function(x, lambda) {
  # Handle non-positive values
  if (any(x <= 0, na.rm = TRUE)) {
    warning("Box-Cox transformation requires positive values. Some values may be NA.")
    # Add small offset to make all values positive
    min_val <- min(x, na.rm = TRUE)
    if (min_val <= 0) {
      x <- x - min_val + 0.01
    }
  }

  if (abs(lambda) < 1e-10) {
    # For lambda close to zero, use logarithm
    log(x)
  } else {
    # Standard Box-Cox transformation
    (x^lambda - 1) / lambda
  }
}

#' Back-transform Box-Cox transformed data
#'
#' @param x Transformed data
#' @param lambda Box-Cox lambda parameter
#' @return Back-transformed data
back_transform_boxcox <- function(x, lambda) {
  if (abs(lambda) < 1e-10) {
    # For lambda close to zero, use exponentiation
    exp(x)
  } else {
    # Back-transform Box-Cox
    (lambda * x + 1)^(1/lambda)
  }
}

#' Apply a specified transformation to all variables in a dataset
#'
#' @param data Data frame to transform
#' @param method Transformation method to apply
#' @param Cls Optional class labels for ABC transformation
#' @return Data frame with all columns transformed using the specified method
apply_single_transformation <- function(data, method, Cls = NULL) {
  transformed_data <- data
  lambdas <- NULL  # Initialize lambdas to NULL

  # Check for non-positive values in the entire dataset for log transformations
  if (method %in% c("log", "log2", "log10")) {
    if (any(data <= 0, na.rm = TRUE)) {
      base <- switch(method,
                     "log" = exp(1),
                     "log2" = 2,
                     "log10" = 10)
      verbose("Warning: Dataset contains values <= 0, switching to signed log transformation with base %s", base)
      method <- paste0("slog", ifelse(base == exp(1), "", base))
    }
  }

  # Apply the chosen transformation to each column
  if (method == "ABCtrans_BC") {
    if (!has_abcstats) {
      warning("ABCstats package not available. Using standard BoxCox instead.")
      method <- "boxcox"
    } else {
      # Use ABCstats adaptive Box-Cox transformation
      if (is.null(Cls)) {
        Cls <- rep(1, nrow(data))  # Default to single class if not provided
      }

      # Prepare data for ABCtransform
      TableForABC <- rbind.data.frame(
        c("Group", Cls),
        cbind.data.frame(Name = names(data), t(data))
      )

      # Apply ABC transformation
      TransformedTable <- ABCstats::ABCtransform(TableForABC)

      # Extract lambdas from the transformation
      BClambdas <- as.numeric(TransformedTable$lambda[-1])

      # Fallback to forecast::BoxCox.lambda for any NA lambdas
      BClambdas_1 <- apply(data, 2, function(x) {
        tryCatch(forecast::BoxCox.lambda(x), error = function(e) NA)
      })
      BClambdas_isNa <- which(is.na(BClambdas))
      BClambdas[BClambdas_isNa] <- BClambdas_1[BClambdas_isNa]

      # Apply Box-Cox transformation with the determined lambdas
      transformed_data <- data.frame(sapply(seq_along(data), function(i) {
        tryCatch(doBoxCoxTrans(x = data[, i], lambda = BClambdas[i]),
                 error = function(e) data[, i])
      }))

      # Set column names and store lambdas
      names(transformed_data) <- names(data)
      lambdas <- BClambdas
      names(lambdas) <- names(data)
    }
  }

  # If method was ABCtrans_BC but failed, or for any other method
  if (method != "ABCtrans_BC" || !has_abcstats) {
    # Apply transformations column by column
    for (i in 1:ncol(data)) {
      col_data <- data[, i]

      # Apply the transformation
      transformed_data[, i] <- switch(method,
                                      "none" = col_data,
                                      "log" = log(col_data),
                                      "log2" = log2(col_data),
                                      "log10" = log10(col_data),
                                      "slog" = slog2Forward(col_data, m = exp(1)),
                                      "slog2" = slog2Forward(col_data, m = 2),
                                      "slog10" = slog2Forward(col_data, m = 10),
                                      "sqrt" = {
                                        if (any(col_data < 0, na.rm = TRUE)) {
                                          verbose("Warning: Column %s contains negative values, cannot apply sqrt transformation",
                                                  names(data)[i])
                                          col_data
                                        } else {
                                          sqrt(col_data)
                                        }
                                      },
                                      "reciprocal" = {
                                        if (any(col_data == 0, na.rm = TRUE)) {
                                          verbose("Warning: Column %s contains zero values, cannot apply reciprocal transformation",
                                                  names(data)[i])
                                          col_data
                                        } else {
                                          1 / col_data
                                        }
                                      },
                                      "boxcox" = {
                                        tryCatch({
                                          lambda <- forecast::BoxCox.lambda(col_data)
                                          if (is.null(lambdas)) lambdas <- numeric(ncol(data))
                                          lambdas[i] <- lambda
                                          forecast::BoxCox(col_data, lambda)
                                        }, error = function(e) {
                                          verbose("BoxCox failed for column %s: %s",
                                                  names(data)[i], e$message)
                                          col_data  # Return original if BoxCox fails
                                        })
                                      },
                                      "tukey_-2" = apply_tukey_trans(col_data, "-2"),
                                      "tukey_-1" = apply_tukey_trans(col_data, "-1"),
                                      "tukey_-0.5" = apply_tukey_trans(col_data, "-0.5"),
                                      "tukey_0" = apply_tukey_trans(col_data, "0"),
                                      "tukey_0.5" = apply_tukey_trans(col_data, "0.5"),
                                      "tukey_1" = apply_tukey_trans(col_data, "1"),
                                      "tukey_2" = apply_tukey_trans(col_data, "2"),
                                      col_data  # Default is no transformation
      )
    }

    if (method == "boxcox" && !is.null(lambdas)) {
      names(lambdas) <- names(data)
    }
  }

  # Store transformation information as attributes
  attr(transformed_data, "transformation") <- method
  if (!is.null(lambdas)) {
    attr(transformed_data, "lambdas") <- lambdas
  }

  return(transformed_data)
}

#' Back-transform data to original scale
#'
#' @param data Transformed data frame
#' @param method Transformation method that was applied
#' @param lambdas Vector of lambda values (for Box-Cox transformations)
#' @return Data frame back-transformed to original scale
back_transform_data <- function(data, method, lambdas = NULL) {
  back_transformed <- data

  if (method == "ABCtrans_BC" || method == "boxcox") {
    if (is.null(lambdas)) {
      stop("Lambda values are required for back-transforming Box-Cox data")
    }
    # Back-transform each column using the corresponding lambda
    for (i in 1:ncol(data)) {
      back_transformed[, i] <- back_transform_boxcox(data[, i], lambdas[i])
    }
  } else if (grepl("^tukey_", method)) {
    # Extract lambda from method name
    lambda <- as.numeric(gsub("tukey_", "", method))
    # Back-transform each column
    for (i in 1:ncol(data)) {
      back_transformed[, i] <- back_transform(data[, i], lambda)
    }
  } else {
    # Handle other transformations
    for (i in 1:ncol(data)) {
      back_transformed[, i] <- switch(method,
                                      "none" = data[, i],
                                      "log" = exp(data[, i]),
                                      "log2" = 2^(data[, i]),
                                      "log10" = 10^(data[, i]),
                                      "slog" = slog2Inverse(data[, i], m = exp(1)),
                                      "slog2" = slog2Inverse(data[, i], m = 2),
                                      "slog10" = slog2Inverse(data[, i], m = 10),
                                      "sqrt" = (data[, i])^2,
                                      "reciprocal" = 1 / (data[, i]),
                                      data[, i]  # Default is no transformation
      )
    }
  }
  return(back_transformed)
}

#' Apply the best Tukey transformation to each variable
#'
#' This function tests each transformation on Tukey's ladder and selects
#' the one that yields the most normal distribution according to the
#' Anderson-Darling test.
#'
#' @param data Data frame to transform
#' @return Data frame with each column optimally transformed
apply_best_tukey_transformation <- function(data) {
  # Define available transformations on Tukey's ladder
  tukey_lambdas <- c("-2", "-1", "-0.5", "0", "0.5", "1", "2")

  # Define descriptive names for each transformation
  tukey_names <- c(
    "-2" = "inverse_square",
    "-1" = "inverse",
    "-0.5" = "inverse_sqrt",
    "0" = "log10",
    "0.5" = "sqrt",
    "1" = "none",
    "2" = "square"
  )

  # Process each column
  transformed_data <- as.data.frame(matrix(NA, nrow = nrow(data), ncol = ncol(data)))
  names(transformed_data) <- names(data)
  transformation_types <- rep(NA, ncol(data))
  transformation_names <- rep(NA, ncol(data))
  variable_names <- names(data)

  for (i in 1:ncol(data)) {
    col_data <- data[, i]

    # Skip if all values are the same or there are too few non-NA values
    if (length(unique(na.omit(col_data))) <= 1 || sum(!is.na(col_data)) < 8) {
      transformed_data[, i] <- col_data
      transformation_types[i] <- "1" # No transformation
      transformation_names[i] <- "none"
      verbose("Column %s: Not enough data or variation for transformation", variable_names[i])
      next
    }

    # Test each transformation
    best_p_value <- 0
    best_transform <- "1" # Default = no transformation

    for (lambda in tukey_lambdas) {
      # Skip transformations that can't be applied to negative/zero values
      if (lambda %in% c("-2", "-1", "-0.5", "0") && any(col_data <= 0, na.rm = TRUE)) {
        next
      }

      # Transform data - handle potential errors
      tryCatch({
        transformed_col <- apply_tukey_trans(col_data, lambda)

        # Skip if transformation creates infinite or NaN values
        if (any(is.infinite(transformed_col) | is.nan(transformed_col), na.rm = TRUE)) next

        # Test normality with Anderson-Darling test
        test_result <- tryCatch(
          ad.test(scale(transformed_col, center = TRUE, scale = TRUE)),
          error = function(e) list(p.value = 0)
        )
        p_value <- test_result$p.value

        # Update if this is the best transformation so far
        if (p_value > best_p_value) {
          best_p_value <- p_value
          best_transform <- lambda
        }
      }, error = function(e) {
        # Skip this transformation if any error occurs
      })
    }

    # Apply the best transformation
    transformed_data[, i] <- apply_tukey_trans(col_data, best_transform)
    transformation_types[i] <- best_transform
    transformation_names[i] <- tukey_names[best_transform]

    verbose("Column %s: Best transformation = %s (p-value = %.4f)",
            variable_names[i], tukey_names[best_transform], best_p_value)
  }

  # Create a detailed transformation summary data frame
  transformation_summary <- data.frame(
    Variable = variable_names,
    Lambda = transformation_types,
    Transformation = transformation_names,
    stringsAsFactors = FALSE
  )

  # Store transformation information as attributes
  attr(transformed_data, "transformations") <- transformation_types
  attr(transformed_data, "transformation_names") <- transformation_names
  attr(transformed_data, "transformation_summary") <- transformation_summary

  return(transformed_data)
}

#' Explore data distribution with different transformations
#'
#' @param data Data frame to explore
#' @param classes Vector of class labels (optional)
#' @param transformation_methods Vector of transformation names to try
#' @param variables Variables to explore (NULL for all)
#' @param plot_results Whether to create and display plots
#' @param max_vars Maximum number of variables to process
#' @return Data frame with Anderson-Darling test results for each variable and transformation
explore_distribution <- function(data, classes = NULL, transformation_methods = "none",
                                 variables = NULL, plot_results = TRUE, max_vars = 20) {
  # Save original plotting parameters to restore later
  original_par <- par(no.readonly = TRUE)
  on.exit(par(original_par))

  if (is.null(variables)) {
    variables <- names(data)
  }

  # Sample only max_vars variables if there are too many
  if (length(variables) > max_vars) {
    set.seed(random_seed)
    variables <- sample(variables, max_vars)
    verbose("Too many variables, sampling %d at random for distribution exploration", max_vars)
  }

  # Create a data frame to store results
  results <- data.frame(
    Variable = character(),
    Transformation = character(),
    AD_Statistic = numeric(),
    AD_P_Value = numeric(),
    Best = character(),
    stringsAsFactors = FALSE
  )

  # If plotting is enabled, set up the plot layout
  if (plot_results) {
    par(mfrow = c(length(transformation_methods), 3))
  }

  # Process each variable
  for (i in 1:length(variables)) {
    variable_data <- subset(data, select = variables[i])[, 1]

    # Create a temporary data frame for this variable's results
    var_results <- data.frame(
      Variable = character(),
      Transformation = character(),
      AD_Statistic = numeric(),
      AD_P_Value = numeric(),
      Best = character(),
      stringsAsFactors = FALSE
    )

    for (i1 in 1:length(transformation_methods)) {
      # Apply transformation
      transformation_success <- TRUE
      ad_statistic <- NA
      ad_p_value <- NA

      if (transformation_methods[i1] == "none") {
        transformed_data <- variable_data
      } else if (transformation_methods[i1] == "log10") {
        # Check for non-positive values
        if (any(variable_data <= 0, na.rm = TRUE)) {
          transformation_success <- FALSE
          if (plot_results) {
            plot(1, 1, main = "Cannot apply log10 to non-positive values", type = "n")
          }
        } else {
          transformed_data <- log10(variable_data)
        }
      } else if (transformation_methods[i1] == "sqrt") {
        # Check for negative values
        if (any(variable_data < 0, na.rm = TRUE)) {
          transformation_success <- FALSE
          if (plot_results) {
            plot(1, 1, main = "Cannot apply sqrt to negative values", type = "n")
          }
        } else {
          transformed_data <- sqrt(variable_data)
        }
      } else if (transformation_methods[i1] == "reciprocal") {
        # Check for zero values
        if (any(variable_data == 0, na.rm = TRUE)) {
          transformation_success <- FALSE
          if (plot_results) {
            plot(1, 1, main = "Cannot apply reciprocal to zero values", type = "n")
          }
        } else {
          transformed_data <- 1 / variable_data
        }
      } else if (transformation_methods[i1] == "boxcox") {
        # Use try-catch to handle potential BoxCox errors
        tryCatch({
          # Handle non-positive values
          if (any(variable_data <= 0, na.rm = TRUE)) {
            transformation_success <- FALSE
            if (plot_results) {
              plot(1, 1, main = "Cannot apply BoxCox to non-positive values", type = "n")
            }
          } else {
            lambda <- forecast::BoxCox.lambda(variable_data, method = "guerrero", lower = -2, upper = 2)
            transformed_data <- forecast::BoxCox(variable_data, lambda)
          }
        }, error = function(e) {
          transformation_success <- FALSE
          if (plot_results) {
            plot(1, 1, main = "BoxCox transformation failed", type = "n")
          }
        })
      } else {
        transformed_data <- variable_data
      }

      # Run Anderson-Darling test if transformation was successful
      if (transformation_success && sum(!is.na(transformed_data)) >= 8) {
        ad_test_result <- tryCatch({
          test <- ad.test(scale(transformed_data, center = TRUE, scale = TRUE))
          ad_statistic <- test$statistic
          ad_p_value <- test$p.value
          test
        }, error = function(e) {
          list(statistic = NA, p.value = NA)
        })
        ad_statistic <- ad_test_result$statistic
        ad_p_value <- ad_test_result$p.value
      } else if (transformation_success && sum(!is.na(transformed_data)) < 8) {
        ad_statistic <- NA
        ad_p_value <- NA
        if (plot_results) {
          plot(1, 1, main = "Not enough data for normality test", type = "n")
        }
      }

      # Add results to the variable-specific data frame
      var_results <- rbind(var_results, data.frame(
        Variable = variables[i],
        Transformation = transformation_methods[i1],
        AD_Statistic = ad_statistic,
        AD_P_Value = ad_p_value,
        Best = "",
        stringsAsFactors = FALSE
      ))

      # Create plots if enabled and transformation was successful
      if (plot_results && transformation_success) {
        # Create histogram
        hist(transformed_data,
             prob = TRUE,
             main = paste(variables[i], "\nTransformation:", transformation_methods[i1],
                          "\nA-D p =", round(ad_p_value, 3)),
             ylab = "Density",
             xlab = "Value (transformed)"
        )

        # Create density plot
        tryCatch({
          if (length(unique(na.omit(transformed_data))) > 1) {
            plot(density(transformed_data, na.rm = TRUE),
                 main = "Density plot",
                 ylab = "Density",
                 xlab = "Value (transformed)")
          } else {
            plot(1, 1, main = "Not enough unique values for density plot", type = "n")
          }
        }, error = function(e) {
          # If density estimation fails
          plot(1, 1, main = "Density estimation failed", type = "n")
        })

        # Create QQ plot
        if (sum(!is.na(transformed_data)) >= 3) {
          qqnorm(transformed_data, main = "Q-Q Plot")
          qqline(transformed_data, col = 2)
        } else {
          plot(1, 1, main = "Not enough data for Q-Q plot", type = "n")
        }
      }
    }

    # Identify the best transformation (highest p-value) for this variable
    max_p_idx <- which.max(var_results$AD_P_Value)
    if (length(max_p_idx) > 0 && !is.na(var_results$AD_P_Value[max_p_idx])) {
      var_results$Best[max_p_idx] <- "*"
    }

    # Add this variable's results to the overall results
    results <- rbind(results, var_results)
  }

  # Return the results data frame
  return(results)
}

#' Evaluate the normality of a dataset
#'
#' @param data Data frame to evaluate
#' @return Data frame with Anderson-Darling test p-values
evaluate_normality <- function(data) {
  # Check normality with Anderson-Darling test
  p_values <- sapply(data, function(x) {
    if (sum(!is.na(x)) < 8 || length(unique(na.omit(x))) <= 1) {
      return(NA)
    }
    tryCatch(ad.test(scale(x))$p.value, error = function(e) NA)
  })

  normality_summary <- data.frame(
    Variable = names(p_values),
    P_value = p_values,
    Is_normal = p_values > 0.05
  )

  verbose("Variables with p < 0.05 (non-normal): %d out of %d (%.1f%%)",
          sum(p_values < 0.05, na.rm = TRUE),
          sum(!is.na(p_values)),
          100 * sum(p_values < 0.05, na.rm = TRUE) / sum(!is.na(p_values)))

  verbose("Median p-value: %.4f", median(p_values, na.rm = TRUE))

  return(normality_summary)
}

#' Create a comparison plot of normality before and after transformation
#'
#' @param original_normality Data frame with original normality assessment
#' @param transformed_normality Data frame with transformed normality assessment
create_normality_comparison_plot <- function(original_normality, transformed_normality) {
  # Save original plotting parameters to restore later
  original_par <- par(no.readonly = TRUE)
  on.exit(par(original_par))

  par(mfrow = c(1, 2))
  plot(original_normality$P_value, ylim = c(0, 1), main = "Original",
       ylab = "Anderson-Darling p-value")
  abline(h = 0.05, col = "red")
  abline(h = median(original_normality$P_value, na.rm = TRUE), col = "blue")

  plot(transformed_normality$P_value, ylim = c(0, 1), main = "Transformed",
       ylab = "Anderson-Darling p-value")
  abline(h = 0.05, col = "red")
  abline(h = median(transformed_normality$P_value, na.rm = TRUE), col = "blue")
}

#' Remove outliers using the BoxPlot method
#'
#' @param data Data frame from which to remove outliers
#' @param coef Coefficient for determining outliers (default = 1.5)
#' @return Data frame with outliers replaced by NA
remove_outliers_boxplot <- function(data, coef = 1.5) {
  # Clone the data
  data_for_outlier_removal <- data

  # Count initial missing values
  na_count_before <- sum(sapply(data_for_outlier_removal, function(x) sum(is.na(x))))
  total_values <- prod(dim(data_for_outlier_removal))

  # Function to identify outliers in a vector using the boxplot method
  is_boxplot_outlier <- function(x, coef = 1.5) {
    # Skip if too many NAs or not enough unique values
    if (sum(!is.na(x)) < 5 || length(unique(na.omit(x))) < 3) {
      return(rep(FALSE, length(x)))
    }

    # Calculate quartiles and IQR
    quants <- quantile(x, probs = c(0.25, 0.75), na.rm = TRUE)
    iqr <- quants[2] - quants[1]

    # Define bounds for outliers
    upper_bound <- quants[2] + coef * iqr
    lower_bound <- quants[1] - coef * iqr

    # Identify outliers
    x < lower_bound | x > upper_bound
  }

  # Track removed outliers by variable
  outliers_by_variable <- numeric(ncol(data))
  names(outliers_by_variable) <- names(data)

  # Apply the outlier detection and removal to each column
  for (i in 1:ncol(data)) {
    col_data <- data_for_outlier_removal[, i]
    outlier_mask <- is_boxplot_outlier(col_data, coef)

    # Count outliers in this column
    outliers_by_variable[i] <- sum(outlier_mask, na.rm = TRUE)

    # Replace outliers with NA
    col_data[outlier_mask] <- NA
    data_for_outlier_removal[, i] <- col_data
  }

  # Calculate statistics after removal
  na_count_after <- sum(sapply(data_for_outlier_removal, function(x) sum(is.na(x))))
  total_outliers_removed <- na_count_after - na_count_before

  verbose("BoxPlot outlier removal: %d outliers removed (%.2f%% of data)",
          total_outliers_removed, 100 * total_outliers_removed / total_values)

  # Store outlier removal info as attributes
  attr(data_for_outlier_removal, "outliers_removed") <- total_outliers_removed
  attr(data_for_outlier_removal, "outliers_by_variable") <- outliers_by_variable
  attr(data_for_outlier_removal, "na_before") <- na_count_before
  attr(data_for_outlier_removal, "na_after") <- na_count_after
  attr(data_for_outlier_removal, "method") <- "boxplot"

  return(data_for_outlier_removal)
}

#' Remove outliers using Grubbs' test
#'
#' @param data Data frame from which to remove outliers
#' @param p_threshold P-value threshold for Grubbs' test
#' @param max_rounds Maximum number of elimination rounds
#' @param random_seed Seed for reproducibility
#' @return Data frame with outliers replaced by NA
remove_outliers_grubbs <- function(data, p_threshold = 1E-5, max_rounds = 100, random_seed = 42) {
  # Clone the data
  data_for_outlier_removal <- data

  # Count initial missing values
  na_count_before <- sum(sapply(data_for_outlier_removal, function(x) sum(is.na(x))))
  total_values <- prod(dim(data_for_outlier_removal))

  # Initialize tracking variables
  parameters_with_outliers_high <- "placeholder"  # Placeholder to start the loop
  parameters_with_outliers_low <- "placeholder"   # Placeholder to start the loop
  elimination_round <- 1
  outliers_by_variable <- numeric(ncol(data))
  names(outliers_by_variable) <- names(data)

  # Set random seed for reproducibility
  set.seed(random_seed)

  # Keep track of outliers removed
  outliers_removed <- 0

  while (length(c(parameters_with_outliers_high, parameters_with_outliers_low)) > 0 &&
    elimination_round <= max_rounds) {

    if (elimination_round > 1) {
      parameters_with_outliers_high <- character(0)
      parameters_with_outliers_low <- character(0)
    }

    verbose("Outlier elimination round: %d", elimination_round)

    # Use tryCatch to handle any errors in Grubbs test
    for (col in names(data_for_outlier_removal)) {
      x <- data_for_outlier_removal[[col]]

      # Skip if too few observations or not enough variation
      if (sum(!is.na(x)) < 8 || length(unique(na.omit(x))) <= 2) {
        next
      }

      # Try to detect highest outlier
      tryCatch({
        test_result <- grubbs.test(x, type = 10, opposite = FALSE, two.sided = FALSE)
        if (test_result$p.value < p_threshold) {
          if (grepl("highest", test_result$alternative)) {
            parameters_with_outliers_high <- c(parameters_with_outliers_high, col)
          } else if (grepl("lowest", test_result$alternative)) {
            parameters_with_outliers_low <- c(parameters_with_outliers_low, col)
          }
        }
      }, error = function(e) {
        # Skip if Grubbs test fails
      })
    }

    # Remove placeholder if it's still there
    parameters_with_outliers_high <- setdiff(parameters_with_outliers_high, "placeholder")
    parameters_with_outliers_low <- setdiff(parameters_with_outliers_low, "placeholder")

    # Track how many outliers are being removed
    outliers_this_round <- 0

    if (length(parameters_with_outliers_high) > 0) {
      for (col in parameters_with_outliers_high) {
        x <- data_for_outlier_removal[[col]]
        to_replace <- which.max(x)
        if (length(to_replace) > 0 && !is.na(x[to_replace])) {
          data_for_outlier_removal[to_replace, col] <- NA
          outliers_this_round <- outliers_this_round + 1
          outliers_by_variable[col] <- outliers_by_variable[col] + 1
        }
      }
    }

    if (length(parameters_with_outliers_low) > 0) {
      for (col in parameters_with_outliers_low) {
        x <- data_for_outlier_removal[[col]]
        to_replace <- which.min(x)
        if (length(to_replace) > 0 && !is.na(x[to_replace])) {
          data_for_outlier_removal[to_replace, col] <- NA
          outliers_this_round <- outliers_this_round + 1
          outliers_by_variable[col] <- outliers_by_variable[col] + 1
        }
      }
    }

    outliers_removed <- outliers_removed + outliers_this_round
    verbose("  Outliers removed this round: %d", outliers_this_round)

    elimination_round <- elimination_round + 1

    # Break if no outliers were found
    if (outliers_this_round == 0) break
  }

  # Calculate statistics after removal
  na_count_after <- sum(sapply(data_for_outlier_removal, function(x) sum(is.na(x))))

  verbose("Grubbs outlier removal: %d outliers removed (%.2f%% of data)",
          outliers_removed, 100 * outliers_removed / total_values)

  # Store outlier removal info as attributes
  attr(data_for_outlier_removal, "outliers_removed") <- outliers_removed
  attr(data_for_outlier_removal, "outliers_by_variable") <- outliers_by_variable
  attr(data_for_outlier_removal, "na_before") <- na_count_before
  attr(data_for_outlier_removal, "na_after") <- na_count_after
  attr(data_for_outlier_removal, "method") <- "grubbs"

  return(data_for_outlier_removal)
}

#' Remove cases (rows) with too many missing values
#'
#' @param data Data frame from which to remove cases
#' @param threshold Maximum proportion of missing values allowed per case
#' @return Data frame with high-missing cases removed
remove_high_missing_cases <- function(data, threshold = 0.5) {
  # Calculate proportion of missing values for each row
  na_proportions <- rowSums(is.na(data)) / ncol(data)

  # Identify rows to keep (those with missing proportion <= threshold)
  rows_to_keep <- na_proportions <= threshold

  # Count how many rows will be removed
  rows_removed <- sum(!rows_to_keep)

  verbose("Removing %d cases with more than %.0f%% missing values",
          rows_removed, threshold * 100)

  # Return filtered data
  return(data[rows_to_keep, , drop = FALSE])
}

#' Impute missing values using missForest
#'
#' @param data Data frame with missing values to impute
#' @param n_cores Number of CPU cores to use
#' @param random_seed Random seed for reproducibility
#' @return Data frame with imputed values
impute_missing_values <- function(data, n_cores = 1, random_seed = 42) {
  # Check if there are any missing values
  if (sum(is.na(data)) == 0) {
    verbose("No missing values to impute")
    return(data)
  }

  verbose("Imputing missing values with missForest using %d cores", n_cores)

  # Set random seed for reproducibility
  set.seed(random_seed)

  # Use missForest for imputation
  imputed_data <- tryCatch({
    missForest(data, parallelize = ifelse(n_cores > 1, "forests", "no"),
               ntree = 100, maxiter = 10, verbose = enable_verbose,
               mtry = floor(sqrt(ncol(data))), replace = TRUE,
               classwt = NULL, cutoff = NULL, strata = NULL,
               sampsize = NULL, nodesize = NULL, maxnodes = NULL,
               xtrue = NULL)
  }, error = function(e) {
    warning("missForest failed: ", e$message, ". Trying simpler imputation.")
    # Fall back to simpler imputation method
    for (col in names(data)) {
      if (any(is.na(data[[col]]))) {
        if (is.numeric(data[[col]])) {
          data[[col]][is.na(data[[col]])] <- median(data[[col]], na.rm = TRUE)
        } else {
          data[[col]][is.na(data[[col]])] <- names(sort(table(data[[col]]), decreasing = TRUE))[1]
        }
      }
    }
    return(list(ximp = data, OOBerror = NA))
  })

  if (class(imputed_data) == "list" && "ximp" %in% names(imputed_data)) {
    return(imputed_data$ximp)
  } else {
    return(data)  # Return original if imputation fails
  }
}

#' Create a heatmap of the data
#'
#' @param data Data frame to visualize
#' @param file_path File path to save the heatmap
#' @param class_column Name of the column containing class labels (or NULL)
#' @param scale_data Logical indicating whether to scale the data
#' @return NULL (creates a plot)
create_heatmap <- function(data, file_path = NULL, class_column = NULL, scale_data = TRUE) {
  if (!requireNamespace("ComplexHeatmap", quietly = TRUE)) {
    warning("ComplexHeatmap package not available. Skipping heatmap creation.")
    return(NULL)
  }

  # Extract classes if class_column is provided
  classes <- NULL
  data_for_heatmap <- data

  # Check if class_column exists in data
  if (!is.null(class_column) && class_column %in% colnames(data)) {
    verbose("Using '%s' column for class annotation in heatmap", class_column)
    classes <- data[[class_column]]
    # Remove class column from data for heatmap
    data_for_heatmap <- data[, !colnames(data) %in% class_column, drop = FALSE]
  } else if (!is.null(class_column)) {
    warning("Class column '%s' not found in data. No class annotation will be used.", class_column)
  }

  # If no classes are provided, create a default class vector
  if (is.null(classes)) {
    classes <- rep(1, nrow(data))  # Default to single class if not provided
    verbose("No class information provided. Using single class for all observations.")
  }

  # Prepare data for heatmap
  heatmap_data <- as.matrix(data_for_heatmap)

  # Scale data if it's numeric and scaling is requested
  if (scale_data) {
    if (is.numeric(heatmap_data)) {
      heatmap_data <- t(scale(t(heatmap_data)))
      verbose("Data scaled for heatmap visualization")
    } else {
      warning("Data contains non-numeric columns. Skipping scaling.")
    }
  }

  # Create annotation with the classes
  row_annotation <- NULL
  class_df <- data.frame(Class = classes)
  class_colors <- list(Class = setNames(
    rainbow(length(unique(classes))),
    unique(classes)
  ))
  row_annotation <- ComplexHeatmap::rowAnnotation(
    df = class_df,
    col = class_colors
  )

  # Create the heatmap
  heatmap <- ComplexHeatmap::Heatmap(
    heatmap_data,
    name = ifelse(scale_data, "Z-score", "Data"),
    col = colorRampPalette(c("dodgerblue", "white", "chartreuse"))(100),
    show_row_names = FALSE,  # Changed to FALSE as typically there are many rows
    show_column_names = TRUE,
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    right_annotation = row_annotation
  )

  # Save to file if path is provided
  if (!is.null(file_path)) {
    svg(file_path, width = 12, height = 12)
    ComplexHeatmap::draw(heatmap)
    dev.off()
    verbose("Heatmap saved to %s", file_path)
  }

  # Display in current device
  ComplexHeatmap::draw(heatmap)

  return(NULL)
}



#' ============================
#' Main script execution
#' ============================

# Check if input file exists
if (!file.exists(input_file)) {
  stop("Input file not found: ", input_file)
}

# Read data
verbose("Reading data from %s...", input_file)
raw_data <- tryCatch(
  read.csv(input_file, row.names = 1),
  error = function(e) stop("Error reading input file: ", e$message)
)

verbose("Data dimensions: %d rows by %d columns",
        dim(raw_data)[1], dim(raw_data)[2])

# Extract class variable if it exists
class_column <- NULL
if (class_name %in% colnames(raw_data)) {
  verbose("Found '%s' column in the data. Extracting as class variable...", class_name)
  class_column <- raw_data[[class_name]]
  raw_data <- raw_data[, !colnames(raw_data) %in% class_name, drop = FALSE]
  verbose("After removing class column: %d rows by %d columns",
          dim(raw_data)[1], dim(raw_data)[2])
}

# Store original data for comparison
original_data <- raw_data

# Explore data distribution if enabled
if (enable_distribution_exploration) {
  verbose("Exploring data distributions...")
  transformation_methods <- c("none", "log10", "sqrt", "reciprocal", "boxcox")
  distribution_results <- explore_distribution(
    data = raw_data,
    classes = class_column,
    transformation_methods = transformation_methods,
    plot_results = enable_plots
  )

  # Print summary of best transformations
  best_transforms <- distribution_results[distribution_results$Best == "*", ]
  verbose("Best transformations by variable:")
  print(best_transforms[, c("Variable", "Transformation", "AD_P_Value")])
}

# Transform data if enabled
transformed_data <- raw_data
lambdas <- NULL  # Initialize lambdas for potential Box-Cox transformations

if (enable_transformation) {
  verbose("Evaluating normality of original data...")
  original_normality <- evaluate_normality(raw_data)

  # Apply transformations based on the specified method
  if (transformation_method == "auto") {
    verbose("Finding and applying optimal transformations for each variable...")
    # Apply optimal transformations
    transformed_data <- apply_best_tukey_transformation(raw_data)

    # Extract transformation information
    transformations <- attr(transformed_data, "transformations")
    transformation_names <- attr(transformed_data, "transformation_names")
    transformation_summary <- attr(transformed_data, "transformation_summary")

    # Print transformation summary
    verbose("Transformation summary:")
    print(transformation_summary)

    # Evaluate normality after transformation
    verbose("Evaluating normality after transformation...")
    transformed_normality <- evaluate_normality(transformed_data)

    # Compare normality before and after
    if (enable_plots) {
      create_normality_comparison_plot(original_normality, transformed_normality)
    }

  } else {
    # Apply a single transformation to all variables
    verbose("Applying %s transformation to all variables...", transformation_method)
    transformed_data <- apply_single_transformation(raw_data, transformation_method, class_column)

    # Extract lambdas if Box-Cox transformation was used
    if (transformation_method %in% c("boxcox", "ABCtrans_BC")) {
      lambdas <- attr(transformed_data, "lambdas")
      if (!is.null(lambdas)) {
        # Write lambdas to file
        lambda_df <- data.frame(
          Variable = names(lambdas),
          Lambda = lambdas
        )
        write.csv(lambda_df, lambda_file, row.names = FALSE)
        verbose("Box-Cox lambdas saved to %s", lambda_file)
      }
    }

    # Evaluate normality after transformation
    verbose("Evaluating normality after transformation...")
    transformed_normality <- evaluate_normality(transformed_data)

    # Compare normality before and after
    if (enable_plots) {
      create_normality_comparison_plot(original_normality, transformed_normality)
    }
  }
}

# Remove outliers if enabled
outlier_removed_data <- transformed_data

if (enable_outlier_removal) {
  verbose("Detecting and removing outliers using %s method...", outlier_method)

  if (outlier_method == "grubbs") {
    outlier_removed_data <- remove_outliers_grubbs(
      transformed_data,
      p_threshold = outlier_p_threshold,
      max_rounds = max_outlier_rounds,
      random_seed = random_seed
    )
  } else if (outlier_method == "boxplot") {
    outlier_removed_data <- remove_outliers_boxplot(
      transformed_data,
      coef = boxplot_coef
    )
  } else {
    warning("Unknown outlier method: ", outlier_method, ". Skipping outlier removal.")
  }

  # Print outlier removal summary
  if (attr(outlier_removed_data, "method") %in% c("grubbs", "boxplot")) {
    outliers_by_variable <- attr(outlier_removed_data, "outliers_by_variable")
    outliers_removed <- attr(outlier_removed_data, "outliers_removed")

    # Sort variables by number of outliers
    sorted_idx <- order(outliers_by_variable, decreasing = TRUE)
    top_outlier_vars <- head(names(outliers_by_variable)[sorted_idx], 10)
    top_outlier_counts <- head(outliers_by_variable[sorted_idx], 10)

    verbose("Top variables with outliers:")
    for (i in 1:length(top_outlier_vars)) {
      if (top_outlier_counts[i] > 0) {
        verbose("  %s: %d outliers", top_outlier_vars[i], top_outlier_counts[i])
      }
    }

    # Update class_column to match the removed rows
    if (!is.null(class_column) && outliers_removed > 0) {
      # Get the rows that were kept
      kept_rows <- attr(outlier_removed_data, "kept_rows")
      if (!is.null(kept_rows)) {
        class_column <- class_column[kept_rows]
        verbose("Updated class column to match outlier removal: %d entries", length(class_column))
      }
    }
  }

  # Remove cases with too many missing values if enabled
  if (enable_case_removal) {
    outlier_removed_data <- remove_high_missing_cases(
      outlier_removed_data,
      threshold = case_outlier_limit
    )

    # Update class_column to match the removed cases
    if (!is.null(class_column)) {
      removed_cases <- attr(outlier_removed_data, "removed_cases")
      if (!is.null(removed_cases) && length(removed_cases) > 0) {
        # Keep only cases that weren't removed
        kept_cases <- setdiff(1:length(class_column), removed_cases)
        class_column <- class_column[kept_cases]
        verbose("Updated class column after case removal: %d entries", length(class_column))
      }
    }
  }
}

# Impute missing values if enabled
processed_data <- outlier_removed_data

if (enable_imputation) {
  verbose("Imputing missing values...")
  processed_data <- impute_missing_values(
    outlier_removed_data,
    n_cores = n_cores,
    random_seed = random_seed
  )
}

# Reinsert class column if it exists
if (!is.null(class_column)) {
  verbose("Reinserting '%s' column into processed data...", class_name)
  # Verify that dimensions match
  if (nrow(processed_data) != length(class_column)) {
    warning("Dimensions mismatch between processed data (%d rows) and class column (%d entries). Class column not added.",
            nrow(processed_data), length(class_column))
  } else {
    processed_data[[class_name]] <- class_column
    verbose("Class column '%s' successfully added to processed data", class_name)
  }
}

# Save processed data
verbose("Saving processed data to %s...", output_file)
write.csv(processed_data, output_file)

# Back-transform data to original scale if needed
if (enable_transformation && transformation_method != "none") {
  verbose("Back-transforming data to original scale...")

  if (transformation_method == "auto") {
    # For auto-transformation, we need to back-transform each variable separately
    back_transformed_data <- processed_data
    transformations <- attr(transformed_data, "transformations")

    if (!is.null(transformations)) {
      for (i in 1:ncol(processed_data)) {
        col_name <- names(processed_data)[i]
        # Skip the class column if it exists
        if (col_name == class_name) next
        lambda <- transformations[i]
        back_transformed_data[, i] <- back_transform(processed_data[, i], lambda)
      }
    }
  } else {
    # For single transformation method, use the appropriate back-transform function
    # Create a copy of processed_data first
    back_transformed_data <- processed_data

    # Get columns to transform (all except class column)
    cols_to_transform <- setdiff(colnames(processed_data), class_name)

    # Apply back-transformation only to the data columns
    back_transformed_data[, cols_to_transform] <- back_transform_data(
      processed_data[, cols_to_transform, drop = FALSE],
      method = transformation_method,
      lambdas = lambdas
    )
  }

  # Save back-transformed data
  verbose("Saving back-transformed data to %s...", output_file_backtransformed)
  write.csv(back_transformed_data, output_file_backtransformed)
}

# Create heatmap of processed data if enabled
if (enable_plots) {
  verbose("Creating heatmap of processed data...")
  create_heatmap(processed_data, heatmap_file, class_column = class_name, scale_data = FALSE)
}

verbose("Data preparation and imputation complete!")