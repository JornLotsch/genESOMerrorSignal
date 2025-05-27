#' AI based data generation local helper function
#'
#' This script provides functions for generating synthetic data points around existing data.
#' The generation process uses a distance-based approach to create new samples with
#' controlled variance from the original points.

# Handle working directory setting
tryCatch({
  if (exists("rstudioapi::getSourceEditorContext")) {
    setwd(dirname(rstudioapi::getSourceEditorContext()$path))
  }
}, error = function(e) {
  message("Unable to set working directory automatically. Please set it manually if needed.")
})

#' Function to generate synthetic data samples around existing data points
#'
#' This function creates new data points in the vicinity of existing samples
#' with a controlled distribution of distances. The generated points maintain
#' the class membership of their parent samples.
#'
#' @param Data Original data matrix or data frame
#' @param density_radius Radius controlling the spread of generated points
#' @param Cls Class labels for the original data (optional)
#' @param gen_per_data Number of synthetic samples to generate per original data point
#' @return A list containing original and generated data with their respective class labels
#' @examples
#' # Generate 10 synthetic samples per original data point
#' result <- generate_synthetic_data(original_data, density_radius = 0.5,
#'                                   Cls = original_classes, gen_per_data = 10)
#' synthetic_samples <- result$generated_data
#' synthetic_classes <- result$generated_classes
generate_synthetic_data <- function(Data, density_radius, Cls = NULL, gen_per_data = 10) {
  # Distance parameters for generation
  max_distance <- 2             # Maximum distance for generated points

  # Distribution parameters for three zones of generation
  limit_ab <- 0.72              # Upper limit of zone A
  limit_bc <- 1.22              # Upper limit of zone B

  # Percentage of points to generate in each zone
  percent_a <- 0.8              # Closest zone (densest)
  percent_b <- 0.15             # Middle zone
  percent_c <- 1 - (percent_a + percent_b)  # Furthest zone (sparsest)

  # Scale radius for generation
  r <- 0.4 * density_radius

  # Get dimensions of input data
  n <- nrow(Data)               # Number of original data points
  d <- ncol(Data)               # Number of features/dimensions

  # Handle class labels
  if (is.null(Cls)) {
    Cls <- rep(1, n)            # Default to single class if not provided
  } else {
    if (length(Cls) != n) {
      stop("Unequal number of cases and class labels.")
    }
  }

  # Calculate total number of points to generate
  n_generated <- n * gen_per_data

  # Create random direction vectors for generation
  jitter <- matrix(rnorm(n_generated * d, mean = 0, sd = 1.5),
                   n_generated, d) * r

  # Normalize jitter vectors to have unit length (pure direction)
  jitter_lengths <- sqrt(rowSums(jitter ^ 2))
  jitter_indices <- which(jitter_lengths > 0)
  jitter[jitter_indices,] <- jitter[jitter_indices,] / jitter_lengths[jitter_indices]

  # Determine how many points to generate in each distance zone
  n_a <- round(n_generated * percent_a)
  n_b <- round(n_generated * percent_b)
  n_c <- n_generated - n_a - n_b

  # Generate distances following the three-zone distribution
  sigmoid_lengths <- c(
    runif(n_a, 0, limit_ab),              # Zone A (closest)
    runif(n_b, limit_ab, limit_bc),       # Zone B (middle)
    runif(n_c, limit_bc, max_distance)    # Zone C (furthest)
  ) * r

  # Create distance matrix for all dimensions
  sigmoid_matrix <- matrix(rep(sigmoid_lengths, d), ncol = d)

  # Generate new data points by adding scaled jitter to repeated original points
  generated_data <- Data[rep(1:n, gen_per_data),] + jitter * sigmoid_matrix

  # Assign class labels to generated points (same as original points)
  generated_classes <- rep(Cls, times = gen_per_data)

  # Return original and generated data with their classes
  return(list(
    original_data = Data,
    original_classes = Cls,
    generated_data = generated_data,
    generated_classes = generated_classes
  ))
}

#' Function to evaluate generated data distribution
#'
#' @param original Original data matrix
#' @param generated Generated data matrix
#' @param original_cls Original class labels
#' @param generated_cls Generated class labels
#' @return Summary statistics about the generated data
evaluate_generated_data <- function(original, generated, original_cls, generated_cls) {
  # Calculate distances between generated points and their parent points
  n_original <- nrow(original)
  n_generated <- nrow(generated)
  gen_per_data <- n_generated / n_original

  # Sample some points for evaluation (if there are too many)
  max_sample <- 1000
  if (n_generated > max_sample) {
    sample_idx <- sample(n_generated, max_sample)
    generated_sample <- generated[sample_idx, ]
    generated_cls_sample <- generated_cls[sample_idx]
    parent_indices <- ceiling(sample_idx / gen_per_data)
  } else {
    generated_sample <- generated
    generated_cls_sample <- generated_cls
    parent_indices <- ceiling(seq_len(n_generated) / gen_per_data)
  }

  # Calculate distances to parent points
  distances <- numeric(length(parent_indices))
  for (i in seq_along(parent_indices)) {
    parent_idx <- parent_indices[i]
    distances[i] <- sqrt(sum((generated_sample[i, ] - original[parent_idx, ])^2))
  }

  # Return summary statistics
  return(list(
    mean_distance = mean(distances),
    median_distance = median(distances),
    min_distance = min(distances),
    max_distance = max(distances),
    quantiles = quantile(distances, probs = c(0.25, 0.75)),
    class_distribution = table(generated_cls)
  ))
}

#' ============================
#' Simple demonstration
#' ============================

if (FALSE) {  # This code only runs if explicitly called
  # Create a simple 2D dataset for demonstration
  original_data <- matrix(c(
    0, 0,
    1, 0,
    0, 1,
    1, 1
  ), ncol = 2, byrow = TRUE)

  original_cls <- c(1, 1, 2, 2)

  # Generate synthetic data
  result <- generate_synthetic_data(
    Data = original_data,
    density_radius = 0.5,
    Cls = original_cls,
    gen_per_data = 20
  )

  # Plot original and generated data
  if (require(ggplot2)) {
    # Create data frame for plotting
    plot_data <- data.frame(
      x = c(original_data[,1], result$generated_data[,1]),
      y = c(original_data[,2], result$generated_data[,2]),
      class = factor(c(original_cls, result$generated_classes)),
      type = factor(c(rep("Original", nrow(original_data)),
                      rep("Generated", nrow(result$generated_data))))
    )

    # Create plot
    p <- ggplot(plot_data, aes(x = x, y = y, color = class, shape = type)) +
      geom_point(alpha = 0.7, size = 3) +
      scale_shape_manual(values = c(Original = 17, Generated = 16)) +
      labs(title = "Original vs. Generated Data Points",
           subtitle = paste("Generated", result$gen_per_data, "points per original point"),
           x = "Feature 1", y = "Feature 2") +
      theme_minimal()

    print(p)
  }
}