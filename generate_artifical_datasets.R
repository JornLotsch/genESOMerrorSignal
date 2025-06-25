#' Artifical Data Sets Generation Script
#'
#' This script generates artificial test datasets for evaluating oversampling methods.
#' It creates two types of datasets:
#' - A "ascending significance" dataset with clear class separation
#' - A "No effect" dataset where classes are indistinguishable
#' Both datasets are useful for benchmarking feature selection and classification algorithms.

# Handle working directory setting
tryCatch({
  if (exists("rstudioapi::getSourceEditorContext")) {
    setwd(dirname(rstudioapi::getSourceEditorContext()$path))
  }
}, error = function(e) {
  message("Unable to set working directory automatically. Please set it manually if needed.")
})

#' Load required libraries
library(stats)     # For statistical functions
library(ggplot2)   # For visualization
library(cowplot)   # For arranging multiple plots

#' Configuration settings
random_seed <- 42                # Base random seed for reproducibility
n_iter <- 100             # Number of iterations for repeated testing

#' ascending significance dataset parameters
ascending_significance_n <- 10     # Number of samples per class
ascending_significance_vars <- 50  # Number of variables to create

#' No effect dataset parameters
no_effect_cases <- 20     # Total number of samples
no_effect_vars <- 150     # Number of initial variables to create
no_effect_final <- 50     # Number of variables to select for final dataset

#' Helper functions
quantiles_95 <- function(x) {
  r <- quantile(x, probs = c(0.025, 0.25, 0.5, 0.75, 0.975))
  names(r) <- c("ymin", "lower", "middle", "upper", "ymax")
  r
}

quantiles_100 <- function(x) {
  r <- quantile(x, probs = c(0.0, 0.25, 0.5, 0.75, 1))
  names(r) <- c("ymin", "lower", "middle", "upper", "ymax")
  r
}

is.integer0 <- function(x) {
  is.integer(x) && length(x) == 0L
}

#' Function to generate dataset with ascending significance
generate_ascending_significance_data <- function(n, n_vars, seed) {
  set.seed(seed)

  # Create sequences of means for the two classes
  ma <- seq(from = 1, to = 50, length.out = n_vars)
  mb <- seq(from = 1, to = 40, length.out = n_vars)

  # Generate variables with different means between classes
  dfTest <- lapply(1:n_vars, function(i) {
    a <- rnorm(n, mean = ma[i], sd = 1)  # Class 1 samples
    b <- rnorm(n, mean = mb[i], sd = 1)  # Class 2 samples
    return(c(a, b))  # Combine samples from both classes
  })

  # Combine all variables into a data frame and add class labels
  dfTestall <- cbind.data.frame(Target = rep(1:2, each = n),
                                do.call(cbind, dfTest))
  return(dfTestall)
}

#' Function to generate dataset with no effect
generate_no_effect_data <- function(n_cases, n_vars, n_select, seed, p_threshold = 0.6) {
  # Create class labels
  classes <- rep(1:2, each = n_cases / 2)

  # Generate list of seeds for each variable
  list_of_seeds <- 1:n_vars + seed - 1

  # Generate variables with no difference between classes
  data_list <- lapply(list_of_seeds, function(i) {
    set.seed(i)
    mean_i <- sample(10:30, 1)
    SD_i <- sample(seq(1, 3, 0.01), size = 1)

    # Generate data for both classes with the same mean and SD
    data_class1 <- rnorm(n = n_cases / 2, mean = mean_i, sd = SD_i)
    set.seed(i + 1000)
    data_class2 <- rnorm(n = n_cases / 2, mean = mean_i, sd = SD_i)

    data_i <- c(data_class1, data_class2)
    return(list(data_i = data_i, sd = SD_i))
  })

  # Extract the data
  data_full <- do.call(cbind, lapply(data_list, "[[", 1))

  # Calculate p-values to ensure no significant differences
  pvals_t <- apply(data_full, 2, function(x) t.test(x ~ classes)$p.value)

  # Find non-significant variables based on p-value threshold
  nonSig <- which(pvals_t > p_threshold)

  # Randomly select n_select non-significant variables
  set.seed(seed)
  nonSig <- sample(nonSig, n_select)

  # Create final dataset with selected variables
  data_final <- data_full[, nonSig]
  colnames(data_final) <- paste0("V", 1:ncol(data_final))

  # Add class labels and return
  result_data <- cbind.data.frame(Target = classes, data_final)
  return(result_data)
}

#' Function to evaluate datasets
evaluate_dataset <- function(data, description) {
  cat(paste0("\n===== Evaluation of ", description, " =====\n"))

  # Check dimensions
  cat(sprintf("Dataset dimensions: %d observations, %d variables\n",
              nrow(data), ncol(data) - 1))

  # Calculate and display p-values for t-tests
  pvals <- apply(data[, -1], 2, function(x) t.test(x ~ data$Target)$p.value)

  cat(sprintf("Variables with p < 0.05: %d (%.1f%%)\n",
              sum(pvals < 0.05),
              100 * sum(pvals < 0.05) / length(pvals)))

  # Display quantiles of p-values
  cat("P-value quantiles:\n")
  print(quantile(pvals, probs = c(0, 0.25, 0.5, 0.75, 1)))
}

#' ============================
#' Main script execution
#' ============================

# Set random seed for reproducibility
set.seed(random_seed)

# Generate dataset with ascending significance (clear class separation)
cat("Generating ascending significance dataset...\n")
ascending_significance_data <- generate_ascending_significance_data(
  n = ascending_significance_n,
  n_vars = ascending_significance_vars,
  seed = random_seed
)

# Generate dataset with no effect (no class separation)
cat("Generating no effect dataset...\n")
no_effect_data <- generate_no_effect_data(
  n_cases = no_effect_cases,
  n_vars = no_effect_vars,
  n_select = no_effect_final,
  seed = random_seed
)

# Evaluate datasets
evaluate_dataset(ascending_significance_data, "ascending significance Dataset")
evaluate_dataset(no_effect_data, "No Effect Dataset")

# Save datasets
write.csv(ascending_significance_data, "ascending_significance_test_data.csv", row.names = FALSE)
write.csv(no_effect_data, "no_effect_test_data.csv", row.names = FALSE)

cat("\nTest data generation complete!\n")