###############################################################################
# Helper Functions for Custom Boxplot Statistics
# These functions are used to define the quantiles for boxplots in the 
# variable importance analysis, ensuring that the boxes and whiskers represent
# specific percentiles (e.g., central 95% or 100% intervals).
###############################################################################

#' Calculate 95% percentile boxplot statistics
#' Returns a named vector for use with ggplot2::stat_summary (2.5%, 25%, 50%, 75%, 97.5%)
quantiles_95 <- function(x) {
  r <- quantile(x, probs = c(0.025, 0.25, 0.5, 0.75, 0.975))
  names(r) <- c("ymin", "lower", "middle", "upper", "ymax")
  r
}

#' Calculate 100% percentile boxplot statistics
#' Returns a named vector for use with ggplot2::stat_summary (0%, 25%, 50%, 75%, 100%)
quantiles_100 <- function(x) {
  r <- quantile(x, probs = c(0.0, 0.25, 0.5, 0.75, 1))
  names(r) <- c("ymin", "lower", "middle", "upper", "ymax")
  r
}

#' Helper function to check for integer(0) objects
is.integer0 <- function(x) {
  is.integer(x) && length(x) == 0L
}

###############################################################################
# Main Function: Analyze Variable Importance Across Different Data Regimes
###############################################################################
#' Analyze Variable Importance Across Different Data Regimes
#' This function assesses the importance of features (variables) for predicting a target variable,
#' using the Boruta algorithm across multiple data configurations (original, reduced, engineered, augmented).
#' It supports parallel execution and returns importance statistics and plots.
#'
#' @param data Main data frame containing features and the target variable.
#' @param class_name Name of the target column.
#' @param data_reduced Data frame with reduced features (optional).
#' @param DataSetSizes Character vector specifying data regimes to analyze.
#' @param GenPerData Number of synthetic samples to generate per data point (for augmented data).
#' @param seed Random seed for reproducibility.
#' @param nIter Number of iterations to repeat the analysis for each data regime.
#' @param max_cores Maximum number of CPU cores to use for parallel processing.
#'
#' @return A list of results for each data regime, containing importance data, plots, and selection statistics.
#'
#' @examples
#' # Example usage (not run):
#' # results <- analyze_variable_importance(data = iris, class_name = "Species", nIter = 10)
analyze_variable_importance <- function(data = Test_DataCls,
                                        class_name = "Target",
                                        data_reduced = Test_DataCls_TrainingTest,
                                        DataSetSizes = c("original", "engineered_0", "augmented_1_engineered", "augmented_5_engineered"),
                                        density_radius, 
                                        GenPerData = 1,
                                        show_varfreq_limit = TRUE,
                                        show_varimp_limit = TRUE,
                                        colorblind = TRUE,
                                        seed = 42,
                                        nIter = 100,
                                        max_cores = 10) {
  # Set up random seeds for reproducibility
  seeds <- 1:nIter + seed - 1

  # Determine the number of cores to use (min of available -1 or max_cores)
  n_cores <- min((parallel::detectCores() - 1), max_cores)

  # Initialize list to store results for each data regime
  results_by_regime <- lapply(DataSetSizes, function(data_regime) {
    # Run Boruta in parallel for each seed in this data regime
    importance_results <- pbmcapply::pbmclapply(seeds, function(iter_seed) {
      # Prepare the data according to the current regime
      if (data_regime == "original") {
        current_data <- data
      } else if (data_regime == "reduced") {
        current_data <- data_reduced
      } else if (grepl("^augmented_\\d+(_engineered)?$", data_regime)) {
        # Extract augmentation factor and check if engineered
        aug_n <- as.integer(sub("^augmented_(\\d+).*", "\\1", data_regime))
        is_engineered <- grepl("_engineered$", data_regime)
        set.seed(iter_seed)
        if (is_engineered) {
          # Combine original data with permuted features
          engineered_data <- cbind.data.frame(
            data,
            apply(data[, !(names(data) %in% class_name), drop = FALSE], 2, sample)
          )
          feature_names <- names(data[, !(names(data) %in% class_name), drop = FALSE])
          names(engineered_data) <- c("Target", feature_names, paste0(feature_names, "_permuted"))
          # Generate synthetic data
          generated_data_list <- generate_synthetic_data(
            Data = engineered_data[, !(names(engineered_data) %in% "Target"), drop = FALSE],
            density_radius = density_radius,
            gen_per_data = aug_n * GenPerData,
            Cls = engineered_data[["Target"]]
          )
          current_data <- rbind.data.frame(
            cbind.data.frame(Target = generated_data_list$original_classes, generated_data_list$original_data),
            cbind.data.frame(Target = generated_data_list$generated_classes, generated_data_list$generated_data)
          )
        } else {
          # Generate synthetic data from original features only
          generated_data_list <- generate_synthetic_data(
            Data = data[, !(names(data) %in% class_name), drop = FALSE],
            density_radius = density_radius,
            gen_per_data = aug_n * GenPerData,
            Cls = data[[class_name]]
          )
          current_data <- rbind.data.frame(
            cbind.data.frame(Target = generated_data_list$original_classes, generated_data_list$original_data),
            cbind.data.frame(Target = generated_data_list$generated_classes, generated_data_list$generated_data)
          )
        }
      } else if (data_regime == "engineered_0") {
        # Just add permuted features to original data (no augmentation)
        set.seed(iter_seed)
        current_data <- cbind.data.frame(
          data,
          apply(data[, !(names(data) %in% class_name), drop = FALSE], 2, sample)
        )
        feature_names <- names(data[, !(names(data) %in% class_name), drop = FALSE])
        names(current_data) <- c("Target", feature_names, paste0(feature_names, "_permuted"))
      }

      # Split data into training and test (not used for Boruta here, but for consistency)
      set.seed(iter_seed)
      in_training <- caret::createDataPartition(current_data$Target, p = .67, list = FALSE)
      training_data <- current_data[in_training,]

      # Run Boruta for feature selection
      boruta_result <- Boruta::Boruta(Target ~ ., training_data, pValue = 0.0001, maxRuns = 100)
      # Extract importance statistics
      boruta_stats <- Boruta::attStats(boruta_result)
      boruta_stats$Var <- rownames(boruta_stats)

      # Reshape importance history for plotting
      importance_long <- reshape2::melt(boruta_result$ImpHistory)
      # Assign color codes for plotting: 1 = dummy, 2 = true, 3 = permuted
      importance_long$ColorVar <- ifelse(importance_long$Var2 %in% names(current_data), 2, 1)
      importance_long$ColorVar[grep("permuted", importance_long$Var2)] <- 3

      return(list(
        boruta_result = boruta_result,
        boruta_stats = boruta_stats,
        importance_long = importance_long
      ))
    }, mc.cores = n_cores)

    # Aggregate Boruta statistics across all iterations
    all_boruta_stats <- do.call(rbind.data.frame, lapply(importance_results, function(x) x$boruta_stats))

    # Initialize feature importance summary list
    feature_importance <- list(
      df_features = NA,
      freq_threshold = NA,
      rel_freq_threshold = NA
    )

    # Process selection frequencies for true and permuted features
    if (!is.integer0(grep("permuted", all_boruta_stats$Var))) {
      permuted_selected <- all_boruta_stats[grep("permuted", all_boruta_stats$Var),]
      permuted_selected <- permuted_selected[permuted_selected$decision == "Confirmed",]
      true_selected <- all_boruta_stats[-grep("permuted", all_boruta_stats$Var),]
      true_selected <- true_selected[true_selected$decision == "Confirmed",]
    } else {
      true_selected <- all_boruta_stats
      true_selected <- true_selected[true_selected$decision == "Confirmed",]
    }

    # Count selection frequencies for true features
    true_selection_counts <- data.frame(table(true_selected$Var))

    # Count selection frequencies for permuted features (if any)
    if (!is.integer0(grep("permuted", all_boruta_stats$Var))) {
      permuted_selection_counts <- data.frame(table(permuted_selected$Var))
      permuted_selection_counts$Var2 <- gsub("_permuted", "", permuted_selection_counts$Var1)
      permuted_selection_counts <- permuted_selection_counts[permuted_selection_counts$Var2 %in% true_selection_counts$Var1,]
    }

    # Prepare data frame with selection frequencies
    if (!is.integer0(grep("permuted", all_boruta_stats$Var))) {
      df_features <- data.frame(Var = unique(all_boruta_stats$Var[-grep("permuted", all_boruta_stats$Var)]))
    } else {
      df_features <- data.frame(Var = unique(all_boruta_stats$Var))
    }
    rownames(df_features) <- df_features$Var
    df_features$SelectedTrue <- true_selection_counts$Freq[match(df_features$Var, true_selection_counts$Var1)]
    if (!is.integer0(grep("permuted", all_boruta_stats$Var))) {
      df_features$SelectedPermuted <- permuted_selection_counts$Freq[match(df_features$Var, permuted_selection_counts$Var2)]
    } else {
      df_features$SelectedPermuted <- NA
    }

    # Replace NAs with 0
    df_features[is.na(df_features)] <- 0
    # Calculate corrected selection frequency
    df_features$SelectedTrueCorr <- df_features$SelectedTrue - df_features$SelectedPermuted
    df_features$SelectedTrueCorrRel <- (df_features$SelectedTrue - df_features$SelectedPermuted) * (df_features$SelectedTrue + df_features$SelectedPermuted)

    # Bootstrap to determine significance thresholds (if permuted features present)
    if (!is.integer0(grep("permuted", all_boruta_stats$Var))) {
      n_bootstrap <- 100000
      set.seed(seed)
      permuted_sample <- sample(df_features$SelectedPermuted, n_bootstrap, replace = TRUE)
      set.seed(seed)
      true_sample <- sample(df_features$SelectedTrue, n_bootstrap, replace = TRUE)
      diff_bootstrap <- true_sample - permuted_sample
      freq_threshold <- quantile(diff_bootstrap, probs = 0.95)

      # Relative frequency threshold
      set.seed(seed)
      permuted_sample_rel <- sample(df_features$SelectedPermuted, n_bootstrap, replace = TRUE)
      set.seed(seed)
      true_sample_rel <- sample(df_features$SelectedTrue, n_bootstrap, replace = TRUE)
      df_bootstrap_rel <- cbind.data.frame(true_sample_rel, permuted_sample_rel)
      rel_diff_bootstrap <- apply(df_bootstrap_rel, 1, function(x)(x["permuted_sample_rel"] - x["true_sample_rel"]) * (x["permuted_sample_rel"] + x["true_sample_rel"]))
      rel_freq_threshold <- quantile(rel_diff_bootstrap, probs = 0.95)
    } else {
      freq_threshold <- NA
      rel_freq_threshold <- NA
    }

    # Prepare selection frequency vector
    select_freq <- as.vector(df_features$SelectedTrueCorr)
    names(select_freq) <- df_features$Var
    select_freq <- na.omit(select_freq)
    attributes(select_freq)$na.action <- NULL

    # Update feature importance summary
    feature_importance <- list(
      df_features = df_features,
      freq_threshold = freq_threshold,
      rel_freq_threshold = rel_freq_threshold
    )

    # Aggregate importance data for plotting
    all_importance_long <- do.call(rbind.data.frame, lapply(importance_results, function(x) x$importance_long))

    # Determine upper limit of non-importance (for plotting)
    upper_limit_non_importance <- NA
    if (!is.integer0(grep("permuted", all_importance_long$Var2))) {
      upper_limit_non_importance <- quantile(all_boruta_stats$maxImp[grep("permuted", all_boruta_stats$Var)], prob = 1)
    }

    # Plot importance distributions
    p_importance <- ggplot2::ggplot(data = all_importance_long, aes(x = reorder(Var2, value), y = value, fill = factor(ColorVar), color = factor(ColorVar))) +
      stat_summary(fun.data = quantiles_95, geom = "boxplot", alpha = 0.2, width = 0.5, position = "dodge") +
      labs(title = "Variable importances", y = "Importance [% decrease in accuracy]", x = NULL, color = "Feature class", fill = "Feature class") +
      theme_light() +
      theme(
        legend.position.inside = TRUE, legend.position = c(.2, .8), legend.direction = "vertical",
        legend.background = element_rect(colour = "transparent", fill = ggplot2::alpha("white", 0.2)),
        strip.background = element_rect(fill = "cornsilk"), strip.text = element_text(colour = "black"),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)
      )


    if (colorblind) {
      p_importance <- p_importance +
        scale_fill_manual(values = ggthemes::colorblind_pal()(8)[c(3, 2, 1)], labels = c("Dummy", "True features", "Permuted features")) +
        scale_color_manual(values = ggthemes::colorblind_pal()(8)[c(3, 2, 1)], labels = c("Dummy", "True features", "Permuted features")) +
        guides(fill = "none")
    } else {
      p_importance <- p_importance +
        scale_fill_manual(values = c("dodgerblue4", "chartreuse2", "salmon"), labels = c("Dummy", "True features", "Permuted features")) +
        scale_color_manual(values = c("dodgerblue4", "chartreuse2", "salmon"), labels = c("Dummy", "True features", "Permuted features")) +
        guides(fill = "none")
    }

    if (!is.na(upper_limit_non_importance) && show_varimp_limit) {
      p_importance <- p_importance +
        geom_hline(yintercept = upper_limit_non_importance, linetype = "dashed", color = "red") +
        annotate("text", x = .5, y = 1.05 * upper_limit_non_importance, label = "Limit of alpha error inflation", color = "red", hjust = -.5)
    }

    # Plot selection frequencies if features were selected
    if (!is.null(dim(feature_importance$df_features))) {
      p_selection_freq2 <- ggplot2::ggplot(data = df_features) +
        geom_bar(aes(y = reorder(Var, SelectedTrueCorr), x = SelectedTrue), fill = "#B37500", color = "#E69F00", alpha = 0.3, stat = "identity") +
        geom_bar(aes(y = reorder(Var, SelectedTrueCorr), x = -SelectedPermuted), fill = "#8C5C00", stat = "identity") +
        labs(title = "Variable selection frequency", x = "Times selected", y = NULL, fill = "Feature class", color = "Feature class") +
        theme_light()

      p_selection_freq <- ggplot2::ggplot(data = df_features) +
        geom_bar(aes(y = reorder(Var, SelectedTrueCorr), x = SelectedTrueCorr), fill = "#8C5C00", color = "#E69F00", alpha = 0.3, stat = "identity") +
        labs(title = "Variable selection frequency", x = "Times selected more than permuted copy", y = NULL, fill = "Feature class", color = "Feature class") +
        theme_light()

      if (!is.na(feature_importance$freq_threshold) && show_varfreq_limit) {
        p_selection_freq <- p_selection_freq + geom_vline(xintercept = feature_importance$freq_threshold, linetype = "dashed", color = "salmon")
      }
    }

    # Return results for this data regime
    return(list(
      importance_data_long = all_importance_long,
      p_importance = p_importance,
      p_selection_freq2 = p_selection_freq2,
      p_selection_freq = p_selection_freq,
      feature_importance = feature_importance
    ))
  })

  # Name the results by data regime
  names(results_by_regime) <- DataSetSizes

  return(results_by_regime)
}
