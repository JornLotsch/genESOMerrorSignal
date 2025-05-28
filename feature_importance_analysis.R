#' Feature Importance Analysis Script
#'
#' This script analyzes feature importance using the Boruta algorithm with synthetic data enhancement.
#' It takes processed data and:
#' - Optionally generates synthetic samples based on original data distribution
#' - Creates engineered variants with permuted features
#' - Runs Boruta feature selection multiple times on different data variants
#' - Calculates statistical significance of features across iterations
#' - Visualizes feature importance with bar and radial plots
#' - Exports results and visualizations to files

#### 1. Utility Functions ####

set_working_directory <- function() {
  tryCatch({
    if (exists("rstudioapi::getSourceEditorContext")) {
      setwd(dirname(rstudioapi::getSourceEditorContext()$path))
    }
  }, error = function(e) {
    message("Unable to set working directory automatically. Please set it manually if needed.")
  })
}

is_integer0 <- function(x) is.integer(x) && length(x) == 0

#### 2. Data Loading Functions ####

load_data <- function(file, class_col) {
  if (!file.exists(file)) stop(sprintf("Data file %s not found", file))
  df <- read.csv(file, row.names = 1)
  df[[class_col]] <- as.factor(df[[class_col]])
  cat(sprintf("Loaded data: %d rows, %d columns\n", nrow(df), ncol(df)))
  print(table(df[[class_col]]))
  df
}

load_radius <- function(output_dir) {
  file <- file.path(output_dir, "Umx_radius.csv")
  if (file.exists(file)) {
    radius <- read.csv(file)$RadiusByEM
    cat(sprintf("Loaded radius: %f\n", radius))
    return(radius)
  } else if (exists("RadiusData")) {
    cat(sprintf("Using radius from environment: %f\n", RadiusData))
    return(RadiusData)
  } else {
    stop("Radius data not found. Please run radius calculation first.")
  }
}

#### 3. Data Preparation ####

split_data <- function(df, class_col, seed, nProc) {
  library(opdisDownsampling)
  cat("Splitting data into train/test/validation...\n")
  set.seed(seed)
  split <- opdisDownsampling::opdisDownsampling(
    Data = within(df, rm(get(class_col))),
    Cls = df[[class_col]],
    Size = 0.8 * nrow(df),
    Seed = seed,
    nTrials = 10000,
    MaxCores = nProc
  )
  train_test <- df[rownames(df) %in% split$ReducedInstances,]
  validation <- df[!rownames(df) %in% split$ReducedInstances,]
  list(train_test = train_test, validation = validation)
}

make_engineered <- function(df, class_col, seed) {
  set.seed(seed)
  features <- setdiff(names(df), class_col)
  for (col in features) {
    df[[paste0(col, "_permuted")]] <- sample(df[[col]])
  }
  df
}

make_augmented <- function(df, class_col, radius, gen_mult, base_rate, seed, gen_func, engineered = FALSE) {
  set.seed(seed)
  if (engineered) df <- make_engineered(df, class_col, seed)
  temp <- df;
  temp[[class_col]] <- NULL
  gen <- gen_func(
    Data = temp,
    density_radius = radius,
    gen_per_data = gen_mult * base_rate,
    Cls = df[[class_col]]
  )
  out <- rbind.data.frame(
    cbind.data.frame(gen$original_data, class = gen$original_classes),
    cbind.data.frame(gen$generated_data, class = gen$generated_classes)
  )
  names(out)[names(out) == "class"] <- class_col
  out
}

get_dataset_variant <- function(type, df, class_col, seed, radius, base_rate, gen_mult, gen_func, train_test = NULL) {
  if (type == "original") return(df)
  if (type == "reduced") {
    if (is.null(train_test)) stop("'reduced' dataset requested but not available")
    return(train_test)
  }
  if (type == "engineered_0") return(make_engineered(df, class_col, seed))
  if (grepl("^augmented_\\d+_engineered$", type)) {
    return(make_augmented(df, class_col, radius, gen_mult, base_rate, seed, gen_func, engineered = TRUE))
  }
  if (grepl("^augmented_\\d+$", type)) {
    return(make_augmented(df, class_col, radius, gen_mult, base_rate, seed, gen_func, engineered = FALSE))
  }
  stop(paste("Unknown dataset type:", type))
}

#### 4. Boruta Feature Selection ####

run_boruta <- function(df, class_col, seed) {
  library(Boruta)
  library(caret)
  library(reshape2)
  set.seed(seed)
  idx <- createDataPartition(df[[class_col]], p = 0.67, list = FALSE)
  train <- df[idx,]
  formula <- as.formula(paste(class_col, "~ ."))
  bor <- Boruta(formula, train, pValue = 0.0001, maxRuns = 100)
  stats <- Boruta::attStats(bor)
  stats$Var <- rownames(stats)
  imp_long <- melt(bor$ImpHistory)
  list(boruta = bor, stats = stats, imp_long = imp_long)
}

#### 5. Importance Aggregation and Plotting ####

process_importance <- function(imps, seed) {
  library(ggplot2)
  stats_all <- do.call(rbind, lapply(imps, function(x) x$stats))
  perm_idx <- grep("permuted", stats_all$Var)
  if (!is_integer0(perm_idx)) {
    perm_notrej <- stats_all[perm_idx,]
    perm_notrej <- perm_notrej[perm_notrej$decision == "Confirmed",]
    notrej <- stats_all[-perm_idx,]
    notrej <- notrej[notrej$decision == "Confirmed",]
  } else {
    notrej <- stats_all[stats_all$decision == "Confirmed",]
  }
  sel_true <- data.frame(table(notrej$Var))
  if (!is_integer0(perm_idx)) {
    sel_perm <- data.frame(table(perm_notrej$Var))
    sel_perm$Var2 <- gsub("_permuted", "", sel_perm$Var1)
    sel_perm <- sel_perm[sel_perm$Var2 %in% sel_true$Var1,]
  }
  vars <- if (!is_integer0(perm_idx)) {
    unique(stats_all$Var[-perm_idx])
  } else unique(stats_all$Var)
  df_vars <- data.frame(Var = vars)
  rownames(df_vars) <- vars
  df_vars$SelectedTrue <- sel_true$Freq[match(df_vars$Var, sel_true$Var1)]
  if (!is_integer0(perm_idx)) {
    df_vars$SelectedPermuted <- sel_perm$Freq[match(df_vars$Var, sel_perm$Var2)]
  } else {
    df_vars$SelectedPermuted <- 0
  }
  df_vars[is.na(df_vars)] <- 0
  df_vars$SelectedTrueCorr <- df_vars$SelectedTrue - df_vars$SelectedPermuted
  limtFreq <- NA
  if (!is_integer0(perm_idx)) {
    nBootstrap <- 100000
    set.seed(seed)
    a_sample <- sample(df_vars$SelectedPermuted, nBootstrap, replace = TRUE)
    b_sample <- sample(df_vars$SelectedTrue, nBootstrap, replace = TRUE)
    FCs_ba_bootstrap <- b_sample - a_sample
    limtFreq <- quantile(FCs_ba_bootstrap, probs = 0.95)
    important_vars <- df_vars$Var[df_vars$SelectedTrueCorr > limtFreq]
  } else {
    important_vars <- df_vars$Var[df_vars$SelectedTrueCorr > 0]
  }
  list(df_vars = df_vars, important_vars = important_vars, limtFreq = limtFreq, stats_all = stats_all)
}

plot_importance <- function(df_vars, title_suffix) {
  library(ggplot2)
  # Temporarily suppress automatic plotting
  old_option <- options(ggplot2.print.object = FALSE)
  on.exit(options(old_option))

  bar <- ggplot(df_vars) +
    geom_bar(aes(y = reorder(Var, SelectedTrueCorr), x = SelectedTrueCorr, fill = SelectedTrueCorr, color = SelectedTrueCorr > 0), stat = "identity") +
    geom_vline(xintercept = 0, linetype = "dashed", color = "salmon") +
    scale_fill_gradient2("Selection Frequency", low = "salmon", mid = "ghostwhite", high = "chartreuse2", midpoint = 0) +
    scale_color_manual(values = c("chartreuse3", "salmon"), guide = "none") +
    labs(title = paste("Variable selection frequency", title_suffix), x = "Times selected more than permuted copy", y = NULL) +
    theme_light() +
    theme(legend.position = "bottom", legend.title = element_text(hjust = 0.5))
  radial <- ggplot(df_vars) +
    geom_hline(aes(yintercept = y), data.frame(y = seq(0, max(df_vars$SelectedTrueCorr, na.rm = TRUE) + 5, by = 5)), color = "lightgrey") +
    geom_col(aes(x = reorder(Var, SelectedTrueCorr), y = SelectedTrueCorr, fill = SelectedTrueCorr, color = SelectedTrueCorr > 0), position = "dodge2", alpha = 0.9) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "salmon", linewidth = 1) +
    coord_polar() +
    scale_y_continuous(limits = c(min(min(df_vars$SelectedTrueCorr, na.rm = TRUE) - 5, -5), max(df_vars$SelectedTrueCorr, na.rm = TRUE) + 5), expand = c(0, 0)) +
    scale_fill_gradient2("Selection Frequency", low = "salmon", mid = "ghostwhite", high = "chartreuse2", midpoint = 0) +
    scale_color_manual(values = c("chartreuse3", "salmon"), guide = "none") +
    guides(fill = guide_colorbar(barwidth = 15, barheight = 0.5, title.position = "top", title.hjust = 0.5)) +
    theme_minimal() +
    theme(axis.title = element_blank(), axis.ticks = element_blank(), axis.text.y = element_blank(), axis.text.x = element_text(color = "gray12", size = 9), legend.position = "bottom", panel.grid.major.x = element_blank(), panel.grid.minor = element_blank()) +
    labs(title = paste("Variable Selection Frequency", title_suffix), subtitle = "Times selected more than permuted copy")
  list(bar = bar, radial = radial)
}

plot_importance_boxplot <- function(Imps_repeated, df_vars, title_suffix, seed) {
  library(ggplot2)

  # Define custom quantile function for boxplot
  quantiles_100 <- function(x) {
    r <- quantile(x, probs = c(0, 0.25, 0.5, 0.75, 1))
    names(r) <- c("ymin", "lower", "middle", "upper", "ymax")
    r
  }

  # Calculate bootstrap thresholds
  limtFreq <- NA

  if (!is_integer0(grep("permuted", df_vars$Var))) {
    nBootstrap <- 100000
    a <- df_vars$SelectedPermuted
    b <- df_vars$SelectedTrue
    set.seed(seed)
    a_sample <- sample(a, nBootstrap, replace = TRUE)
    set.seed(seed)
    b_sample <- sample(b, nBootstrap, replace = TRUE)
    FCs_ba_bootstrap <- b_sample - a_sample
    limtFreq <- quantile(FCs_ba_bootstrap, probs = 0.95)
  }

  # Extract importance data
  imp_long_all <- do.call(rbind.data.frame, lapply(Imps_repeated, function(x) x$imp_long))

  # Add color variable for different feature types
  imp_long_all$ColorVar <- "True"
  if (!is_integer0(grep("permuted", imp_long_all$Var2))) {
    imp_long_all$ColorVar[grep("permuted", imp_long_all$Var2)] <- "Permuted"
  }
  imp_long_all$ColorVar[grep("shadow", imp_long_all$Var2)] <- "Dummy"

  # Calculate upper border of non-importance
  upperBorderOfNonImportance <- NA
  if (!is_integer0(grep("permuted", imp_long_all$Var2))) {
    stats_all <- do.call(rbind, lapply(Imps_repeated, function(x) x$stats))
    upperBorderOfNonImportance <- quantile(stats_all$maxImp[grep("permuted", stats_all$Var)], prob = 1)
  }

  # Create importance boxplot
  pVarimp_Test_actual <- ggplot(data = imp_long_all,
                                aes(x = reorder(Var2, value), y = value, fill = factor(ColorVar))) +
    stat_summary(fun.data = quantiles_100, geom = "boxplot", alpha = 0.2, width = 0.5, position = "dodge") +
    scale_fill_manual(values = c("dodgerblue4", "chartreuse2", "salmon"),
                      labels = c("Dummy", "True features", "Permuted features")) +
    labs(title = paste("Variable importances -", title_suffix),
         y = "Importance [% decrease in accuracy]",
         x = NULL,
         fill = "Feature class") +
    theme_light() +
    theme(
      legend.position = c(.2, .8),
      legend.direction = "vertical",
      legend.background = element_rect(colour = "transparent", fill = ggplot2::alpha("white", 0.2)),
      strip.background = element_rect(fill = "cornsilk"),
      strip.text = element_text(colour = "black"),
      axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)
    )

  # Add reference line for non-importance threshold if available
  if (!is.na(upperBorderOfNonImportance)) {
    pVarimp_Test_actual <- pVarimp_Test_actual +
      geom_hline(yintercept = upperBorderOfNonImportance, linetype = "dashed", color = "red") +
      annotate("text", x = 0.5, y = 1.05 * upperBorderOfNonImportance,
               label = "Limit of alpha error inflation", color = "red", hjust = -0.5)
  }

  list(
    plot = pVarimp_Test_actual,
    limtFreq = limtFreq,
    upperBorderOfNonImportance = upperBorderOfNonImportance
  )
}

#### 6. Main Pipeline ####

feature_importance_pipeline <- function(
  output_dir = "results",
  output_prefix = "Feature_Importance",
  input_file = "test_processed.csv",
  class_name = "Species",
  generation_multipliers = c(1, 5),
  base_generation_rate = 1,
  nIter = 100,
  seed = 42,
  enable_plots = TRUE,
  enable_file_output = TRUE
) {
  # Setup environment
  set_working_directory()
  source("generate_synthetic_data.R")
  library(parallel);
  library(pbmcapply);
  library(cowplot)

  nProc <- detectCores() - 1
  list.of.seeds <- seed + 0:(nIter - 1)

  # Create output directory if needed
  if (enable_file_output && !dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }

  # Load and prepare data
  data_df <- load_data(input_file, class_name)
  RadiusData <- load_radius(output_dir)

  # Configure dataset variants for analysis
  analysisDatasetVariants <- c("original", "engineered_0")
  for (mult in generation_multipliers) {
    analysisDatasetVariants <- c(analysisDatasetVariants, paste0("augmented_", mult, "_engineered"))
  }
  cat("Dataset types to analyze:", paste(analysisDatasetVariants, collapse = ", "), "\n")

  # Prepare data split if needed for "reduced" variants
  split <- if (any(grepl("reduced", analysisDatasetVariants))) {
    split_data(data_df, class_name, seed, nProc)
  } else {
    list()
  }

  # Initialize results containers
  featureImportanceResults <- list()
  significant_variables <- list()

  # Initialize plot lists with named elements
  plot_list_bars <- structure(vector("list", length(analysisDatasetVariants)), names = analysisDatasetVariants)
  plot_list_radial <- structure(vector("list", length(analysisDatasetVariants)), names = analysisDatasetVariants)
  plot_list_boxplots <- structure(vector("list", length(analysisDatasetVariants)), names = analysisDatasetVariants)

  # Process each dataset variant
  for (dataset_variant in analysisDatasetVariants) {
    cat(sprintf("Processing dataset type: %s\n", dataset_variant))

    # Extract generation multiplier for augmented datasets
    gen_multiplier <- if (grepl("augmented_", dataset_variant)) {
      as.numeric(gsub("augmented_([0-9]+)_.*", "\\1", dataset_variant))
    } else {
      1
    }

    # Run Boruta feature selection in parallel
    Imps_repeated <- pbmcapply::pbmclapply(list.of.seeds, function(s) {
      # Get appropriate dataset variant
      df_variant <- get_dataset_variant(
        dataset_variant, data_df, class_name, s, RadiusData,
        base_generation_rate, gen_multiplier, generate_synthetic_data,
        train_test = split$train_test
      )
      # Run Boruta on this dataset
      run_boruta(df_variant, class_name, s)
    }, mc.cores = nProc)

    # Process importance results
    imp_proc <- process_importance(Imps_repeated, seed)

    # Generate frequency plots
    freq_plots <- plot_importance(imp_proc$df_vars, dataset_variant)

    # Generate boxplot
    boxplot_results <- plot_importance_boxplot(Imps_repeated, imp_proc$df_vars, dataset_variant, seed)

    # Store results
    featureImportanceResults[[dataset_variant]] <- list(
      df_vars = imp_proc$df_vars,
      important_vars = imp_proc$important_vars,
      bar_plot = freq_plots$bar,
      radial_plot = freq_plots$radial,
      importance_boxplot = boxplot_results$plot,
      limtFreq = boxplot_results$limtFreq,
      upperBorderOfNonImportance = boxplot_results$upperBorderOfNonImportance
    )
    significant_variables[[dataset_variant]] <- imp_proc$important_vars

    # Store plots for combined view
    plot_list_bars[[dataset_variant]] <- freq_plots$bar
    plot_list_radial[[dataset_variant]] <- freq_plots$radial
    plot_list_boxplots[[dataset_variant]] <- boxplot_results$plot

    # # Display plots
    # if (enable_plots) {
    #   print(freq_plots$bar)
    #   print(freq_plots$radial)
    #   print(boxplot_results$plot)
    # }

    # Save individual plots if requested
    if (enable_plots && enable_file_output) {
      ggsave(
        file.path(output_dir, paste0(output_prefix, "_", dataset_variant, "_frequency.svg")),
        freq_plots$bar, width = 10, height = 8
      )
      ggsave(
        file.path(output_dir, paste0(output_prefix, "_", dataset_variant, "_frequency_circular.svg")),
        freq_plots$radial, width = 10, height = 10
      )
      ggsave(
        file.path(output_dir, paste0(output_prefix, "_", dataset_variant, "_importance_boxplot.svg")),
        boxplot_results$plot, width = 12, height = 8
      )
    }
  }

  # Create combined plots
  if (enable_plots) {
    # Extract plots from results list in correct order
    plot_list_bars <- lapply(analysisDatasetVariants, function(ds) featureImportanceResults[[ds]]$bar_plot)
    plot_list_radial <- lapply(analysisDatasetVariants, function(ds) featureImportanceResults[[ds]]$radial_plot)
    plot_list_boxplots <- lapply(analysisDatasetVariants, function(ds) featureImportanceResults[[ds]]$importance_boxplot)

    # Create combined bar plot
    p_bars <- do.call(plot_grid, c(
      plot_list_bars,
      list(labels = LETTERS[1:length(analysisDatasetVariants)], nrow = 1, align = "h", axis = "tb")
    ))

    # Create combined radial plot
    p_radial <- do.call(plot_grid, c(
      plot_list_radial,
      list(labels = LETTERS[1:length(analysisDatasetVariants)], nrow = 1, align = "h", axis = "tb")
    ))

    # Create combined boxplot
    p_boxplots <- do.call(plot_grid, c(
      plot_list_boxplots,
      list(labels = LETTERS[1:length(analysisDatasetVariants)], nrow = 1, align = "h", axis = "tb")
    ))

    # Display combined plots - force printing explicitly
    print(p_bars)
    print(p_radial)
    print(p_boxplots)

    # Save combined plots if requested
    if (enable_file_output) {
      ggsave(file.path(output_dir, "p_Varfreqs_aug_bars.svg"), p_bars,
             width = 20, height = 10, limitsize = FALSE)
      ggsave(file.path(output_dir, "p_Varfreqs_aug_radial.svg"), p_radial,
             width = 20, height = 10, limitsize = FALSE)
      ggsave(file.path(output_dir, "p_Varfreqs_aug_boxplots.svg"), p_boxplots,
             width = 20, height = 10, limitsize = FALSE)

      ggsave(file.path(output_dir, "p_Varfreqs_aug_bars.png"), p_bars,
             width = 20, height = 10, limitsize = FALSE, dpi = 300)
      ggsave(file.path(output_dir, "p_Varfreqs_aug_radial.png"), p_radial,
             width = 20, height = 10, limitsize = FALSE, dpi = 300)
      ggsave(file.path(output_dir, "p_Varfreqs_aug_boxplots.png"), p_boxplots,
             width = 20, height = 10, limitsize = FALSE, dpi = 300)
    }
  }

  # Create significance summary table
  all_vars <- unique(unlist(significant_variables))
  significant_matrix <- matrix(0, nrow = length(all_vars), ncol = length(analysisDatasetVariants))
  rownames(significant_matrix) <- all_vars
  colnames(significant_matrix) <- analysisDatasetVariants

  # Fill significance matrix
  for (i in seq_along(analysisDatasetVariants)) {
    ds <- analysisDatasetVariants[i]
    vars <- featureImportanceResults[[ds]]$important_vars
    significant_matrix[vars, i] <- 1
  }

  # Convert to data frame, add totals, and sort
  significant_vars_df <- as.data.frame(significant_matrix)
  significant_vars_df$SignificantCount <- rowSums(significant_vars_df)
  significant_vars_df <- significant_vars_df[order(-significant_vars_df$SignificantCount, rownames(significant_vars_df)),]

  # Display significance summary
  print(significant_vars_df)

  # Save results if requested
  if (enable_file_output) {
    write.csv(significant_vars_df,
              file = file.path(output_dir, paste0(output_prefix, "_significant_variables.csv")),
              row.names = TRUE)
    save(featureImportanceResults, significant_vars_df, generation_multipliers,
         file = file.path(output_dir, paste0(output_prefix, "_results.RData")))
  }

  cat("Feature importance analysis complete.\n")

  # Return results for potential further use
  invisible(list(
    feature_importance = featureImportanceResults,
    significant_variables = significant_vars_df,
    p_bars = p_bars,
    p_radial = p_radial,
    p_boxplots = p_boxplots
  ))
}

#### 7. Run the pipeline ####

# Run the pipeline with parameters for the actual dataset
test_result <- feature_importance_pipeline(
  output_dir = "results/feature_importance",
  output_prefix = "Feature_Importance_Analysis",
  input_file = "test_processed.csv",
  class_name = "Species",
  generation_multipliers = c(1, 5),
  base_generation_rate = 1,
  nIter = 100,
  seed = 123,
  enable_plots = TRUE,
  enable_file_output = TRUE
)
