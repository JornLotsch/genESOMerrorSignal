#' Artificial Data Sets Generation Script
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

#' Configuration settings
random_seed <- 42                # Base random seed for reproducibility

#' ascending significance dataset parameters
ascending_significance_n <- 10     # Number of samples per class
ascending_significance_vars <- 50  # Number of variables to create

#' No effect dataset parameters
no_effect_cases <- 20     # Total number of samples
no_effect_vars <- 150     # Number of initial variables to create
no_effect_final <- 50     # Number of variables to select for final dataset

#' ============================
#' Dataset Generation Functions
#' ============================

#' Generate dataset with ascending significance (clear class separation)
generate_ascending_significance_data <- function(n, n_vars, seed) {
  ma <- seq(from = 1, to = 50, length.out = n_vars)
  mb <- seq(from = 1, to = 40, length.out = n_vars)
  dfTest <- lapply(1:n_vars, function(i) {
    set.seed(seed)
    a <- rnorm(n, mean = ma[i], sd = 1)
    b <- rnorm(n, mean = mb[i], sd = 1)
    c(a, b)
  })
  dfTestall <- cbind.data.frame(Target = rep(1:2, each = n), do.call(cbind, dfTest))
  return(dfTestall)
}

#' Generate dataset with no effect (no class separation)
generate_no_effect_data <- function(n_cases, n_vars, n_select, seed, p_threshold = 0.6) {
  list_of_seeds <- 1:n_vars + seed - 1
  classes <- rep(1:2, each = n_cases / 2)
  data_noDiff_list <- lapply(list_of_seeds, function(i) {
    set.seed(i)
    mean_i <- sample(10:30, 1)
    SD_i_1 <- sample(seq(1, 3, 0.01), size = 1)
    Data_i_1 <- rnorm(n = n_cases / 2, mean = mean_i, sd = SD_i_1)
    set.seed(i + 1000)
    SD_i_2 <- SD_i_1
    Data_i_2 <- rnorm(n = n_cases / 2, mean = mean_i, sd = SD_i_2)
    Data_i <- c(Data_i_1, Data_i_2)
    return(Data_i)
  })
  data_noDiff <- do.call(cbind, data_noDiff_list)
  pvals_t <- apply(data_noDiff, 2, function(x) t.test(x ~ classes)$p.value)
  nonSig05 <- which(pvals_t > p_threshold)
  set.seed(seed)
  data_noDiff_50 <- data_noDiff[, sample(nonSig05, n_select)]
  result_data <- cbind.data.frame(Target = classes, data_noDiff_50)
  return(result_data)
}

#' Function to evaluate datasets
evaluate_dataset <- function(data, description) {
  cat(paste0("\n===== Evaluation of ", description, " =====\n"))
  cat(sprintf("Dataset dimensions: %d observations, %d variables\n",
              nrow(data), ncol(data) - 1))
  pvals <- apply(data[, -1], 2, function(x) t.test(x ~ data$Target)$p.value)
  cat(sprintf("Variables with p < 0.05: %d (%.1f%%)\n",
              sum(pvals < 0.05),
              100 * sum(pvals < 0.05) / length(pvals)))
  cat("P-value quantiles:\n")
  print(quantile(pvals, probs = c(0, 0.25, 0.5, 0.75, 1)))
}

#' ============================
#' Main script execution
#' ============================

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
