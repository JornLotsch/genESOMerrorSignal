###############################################################################
# Experiment 1: Ascending Statistical Significance of Class Differences
# Author: [Your Name]
# Date: [YYYY-MM-DD]
# Description: 
#   - Loads and analyzes the artificial dataset with systematically increasing 
#     class differences for each variable.
#   - Performs t-tests, visualizes significance, splits data, and runs feature 
#     importance analysis.
###############################################################################

# ---- 1. Setup Environment ---------------------------------------------------

# Set working directory (consider using relative paths in projects)
setwd("/home/joern/Aktuell/GenerativeESOM/08AnalyseProgramme/R/genESOMerrorSignal")

# Source required helper scripts
source("generate_synthetic_data.R")      # For synthetic data generation functions
source("analyze_variable_importance.R")  # For variable importance analysis functions

# Load required libraries
library(reshape2)
library(ggplot2)
library(parallel)
library(opdisDownsampling)
library(cowplot)
library(patchwork)

# ---- 2. Load Artificial Datasets --------------------------------------------

ascending_significance_data <- read.csv("ascending_significance_test_data.csv")
# (Second dataset loading comes later)
# no_effect_data <- read.csv("no_effect_test_data.csv")

# ---- 3. Statistical Significance Analysis -----------------------------------

# Perform t-tests for each variable (excluding the target column)
t_tests_p_ascending_significance_data <- apply(
  ascending_significance_data[,-1], 2,
  function(x) t.test(x ~ ascending_significance_data$Target)$p.value
)

# Prepare data for plotting
t_tests_p_ascending_significance_data_long <- cbind.data.frame(
  variable = names(t_tests_p_ascending_significance_data), 
  value = reshape2::melt(t_tests_p_ascending_significance_data)
)
t_tests_p_ascending_significance_data_long$variable <- factor(
  t_tests_p_ascending_significance_data_long$variable,
  levels = paste0("X", 1:50)
)

# ---- 4. Visualize Significance ----------------------------------------------

barplot_t_tests_p_ascending_significance_data <- ggplot(
  data = t_tests_p_ascending_significance_data_long, 
  aes(y = variable, x = -log10(value))
) +
  geom_bar(stat = "identity", color = "cornsilk3", fill = "cornsilk2") +
  geom_vline(xintercept = -log10(0.05), color = "salmon", linetype = "dashed") +
  annotate(
    "text",
    x = -log10(0.05), 
    y = t_tests_p_ascending_significance_data_long$variable[1], # top variable
    label = "p = 0.05",
    hjust = -0.1,
    vjust = 1.5,
    color = "salmon",
    fontface = "plain",
    angle = 90
  ) +
  theme_light() +
  labs(title = "Significance of class differences")

print(barplot_t_tests_p_ascending_significance_data)

# ---- 5. Split Dataset: Training/Test/Validation -----------------------------

# Detect available cores for parallel processing (set max_cores if needed)
max_cores <- NULL
nProc <- min(parallel::detectCores() - 1, max_cores)

# Downsample for train/test/validation split
ascending_significance_data_TestValidation <- opdisDownsampling::opdisDownsampling(
  Data = within(ascending_significance_data, rm(Target)),
  Cls = ascending_significance_data$Target,
  Size = 0.8 * nrow(ascending_significance_data),
  Seed = seed,
  nTrials = 10000,
  MaxCores = nProc
)

ascending_significance_data_TrainingTest <- ascending_significance_data[
  rownames(ascending_significance_data) %in% ascending_significance_data_TestValidation$ReducedInstances, ]
table(ascending_significance_data_TrainingTest$Target)

ascending_significance_data_Validation <- ascending_significance_data[
  !rownames(ascending_significance_data) %in% ascending_significance_data_TestValidation$ReducedInstances, ]
table(ascending_significance_data_Validation$Target)

# ---- 6. Feature Importance Analysis -----------------------------------------

# Set parameters
DataSetSizes <- c("original", "engineered_0", "augmented_1_engineered", 
                  "augmented_5_engineered", "augmented_50_engineered")
GenPerData <- 1
nIter <- 100
seed <- 42
list.of.seeds <- 1:nIter + seed - 1

# Run variable importance analysis
results_analyze_variable_importance <- analyze_variable_importance(
  data = ascending_significance_data,
  class_name = "Target",
  data_reduced = ascending_significance_data_TrainingTest, 
  DataSetSizes = DataSetSizes,  
  show_varfreq_limit = FALSE,
  show_varimp_limit = TRUE,
)

# ---- 7. Plot Selection Frequency for Each Data Regime -----------------------

plot_ascending_significance_data_selection_freq <- cowplot::plot_grid(
  barplot_t_tests_p_ascending_significance_data,
  results_analyze_variable_importance$original$p_selection_freq + labs(title = "Original"),
  results_analyze_variable_importance$engineered_0$p_selection_freq + labs(title = "Engineered 0"),
  results_analyze_variable_importance$augmented_1_engineered$p_selection_freq + labs(title = "Augmented 1, engineered"),
  results_analyze_variable_importance$augmented_5_engineered$p_selection_freq + labs(title = "Augmented 5, engineered"),
  results_analyze_variable_importance$augmented_50_engineered$p_selection_freq + labs(title = "Augmented 50, engineered"),
  labels = "AUTO", 
  nrow = 1,
  align = "h", axis = "tb"
) + 
  plot_annotation(
    title = "Variable selection frequency",
    subtitle = "Dataset: ascending_significance_data"
  ) & 
  theme(
    plot.tag.position = c(0.5, 1),   # horizontally centered, at top edge
    plot.tag = element_text(size = 14, face = "bold", vjust = 0, margin = margin(b = -10))
  )


# Save the plot
ggsave(
  filename = paste0("plot_ascending_significance_data_selection_freq", ".svg"),
  plot = plot_ascending_significance_data_selection_freq, 
  width = 22, height = 10, limitsize = FALSE
)

# ---- 7. Plot variable Importance for Augmentaion 5 as an example -----------------------

plot_ascending_significance_data_var_importance <- cowplot::plot_grid(
  results_analyze_variable_importance$original$p_importance + labs(title = "Original"), 
    results_analyze_variable_importance$augmented_5_engineered$p_importance + labs(title = "Augmented 5, engineered"),
    labels = "AUTO", 
  nrow = 1,
  align = "h", axis = "tb",
  rel_widths = c(1,2)
) + 
  plot_annotation(
    title = "Variable importance",
    subtitle = "Dataset: ascending_significance_data"
  ) & 
  theme(
    plot.tag.position = c(0.5, 1),   # horizontally centered, at top edge
    plot.tag = element_text(size = 14, face = "bold", vjust = 0, margin = margin(b = -10))
  )


print(plot_ascending_significance_data_var_importance)

# Save the plot
ggsave(
  filename = paste0("plot_ascending_significance_data_var_importance", ".svg"),
  plot = plot_ascending_significance_data_var_importance, 
  width = 20, height = 8, limitsize = FALSE
)




###############################################################################
# (Code for the second dataset will follow below)
###############################################################################


# ############## Artifical data sets ##########################################
# 
# # ---- 1. Setup ----------------------------------------------------------------
# setwd("/home/joern/Aktuell/GenerativeESOM/08AnalyseProgramme/R/genESOMerrorSignal")
# 
# # Read datasets
# 
# ascending_significance_data <- read.csv("ascending_significance_test_data.csv")
# no_effect_data <- read.csv("no_effect_test_data.csv")
# 
# # Data set 1: ascending_significance_data
# # T-tests
# t_tests_p_ascending_significance_data <- 
#   apply(ascending_significance_data[,-1],2,function(x) t.test(x~ascending_significance_data$Target)$p.value)
# 
# t_tests_p_ascending_significance_data_long <- cbind.data.frame(variable = names(t_tests_p_ascending_significance_data), 
#                                                                value = reshape2::melt(t_tests_p_ascending_significance_data))
# t_tests_p_ascending_significance_data_long$variable <- factor(
#   t_tests_p_ascending_significance_data_long$variable,
#   levels = paste0("X", 1:50)
# )
# 
# barplot_t_tests_p_ascending_significance_data <- 
#   ggplot(data = t_tests_p_ascending_significance_data_long, aes(y = variable, x = -log10(value))) +
#   geom_bar(stat = "identity", color = "cornsilk3", fill = "cornsilk2") +
#   geom_vline(xintercept = -log10(0.05), color = "salmon", linetype = "dashed") +
#   annotate(
#     "text",
#     x = -log10(0.05), 
#     y = t_tests_p_ascending_significance_data_long$variable[1], # top variable
#     label = "p = 0.05",
#     hjust = -0.1, # adjust as needed to move label left/right
#     vjust = 1.5,  # adjust as needed to move label up/down
#     color = "salmon",
#     fontface = "plain",
#     angle = 90
#   ) +
#   theme_light() +
#   labs(title = "Significance of class differences")
# 
# 
# print(barplot_t_tests_p_ascending_significance_data)
# 
# 
# #################################### Split data set into training/test / validation  ########################################################################
# max_cores <- NULL
# nProc <- min(parallel::detectCores() - 1, max_cores)  
# 
# ascending_significance_data_TestValidation <- 
#   opdisDownsampling::opdisDownsampling(Data = within(ascending_significance_data, rm(Target)), Cls = ascending_significance_data$Target,
#                                                                       Size = 0.8 * nrow(ascending_significance_data), Seed = seed, nTrials = 10000, 
#                                        MaxCores = nProc)
# 
# ascending_significance_data_TrainingTest <- 
#   ascending_significance_data[rownames(ascending_significance_data) %in% ascending_significance_data_TestValidation$ReducedInstances,]
# table(ascending_significance_data_TrainingTest$Target)
# 
# ascending_significance_data_Validation <- 
#   ascending_significance_data[!rownames(ascending_significance_data) %in% ascending_significance_data_TestValidation$ReducedInstances,]
# table(ascending_significance_data_Validation$Target)
# 
# 
# #################################### Perform feature importance analysis  ########################################################################
# 
# # Parameters
# 
# DataSetSizes <- c("original", "engineered_0", "augmented_1_engineered", "augmented_5_engineered", "augmented_50_engineered")
# GenPerData <- 1
# 
# nProc <- parallel::detectCores() - 1  
# nIter <- 100
# seed <- 42
# list.of.seeds <- 1:nIter + seed - 1
# 
# 
# # Run the analysis with defaults but data names
# results_analyze_variable_importance <- analyze_variable_importance(data = ascending_significance_data,
#                                        class_name = "Target",
#                                        data_reduced = ascending_significance_data_TrainingTest, 
#                                        DataSetSizes = DataSetSizes)
# 
# # Plot selection frequency for each data regime
# plot_results_analyze_variable_importance <- cowplot::plot_grid(
#   barplot_t_tests_p_ascending_significance_data,
#   results_analyze_variable_importance$original$p_selection_freq + labs(title = "Original"),
#   results_analyze_variable_importance$engineered_0$p_selection_freq + labs(title = "Engineered 0"),
#   results_analyze_variable_importance$augmented_1_engineered$p_selection_freq + labs(title = "Augmented 1, engineered"),
#   results_analyze_variable_importance$augmented_5_engineered$p_selection_freq + labs(title = "Augmented 5, engineered"),
#   results_analyze_variable_importance$augmented_50_engineered$p_selection_freq + labs(title = "Augmented 50, engineered"),
#   # Uncomment and use the following lines if you want to include importance plots
#   # results_analyze_variable_importance$augmented_5$p_importance + labs(title = "Augmented 5"),
#   # results_analyze_variable_importance$augmented_engineered_1$p_importance + labs(title = "Augmented 1, engineered"),
#   # results_analyze_variable_importance$augmented_engineered_5$p_importance + labs(title = "Augmented 5, engineered"),
#   labels = LETTERS[1:4], # Only 4 plots in this grid; adjust if you uncomment more
#   nrow = 1,
#   align = "h", axis = "tb"
# )
# 
# print(plot_results_analyze_variable_importance)
# 
# ggsave(
#   filename = paste0("plot_results_analyze_variable_importance", ".svg"),
#   plot = plot_results_analyze_variable_importance, width = 22, height = 9, limitsize = FALSE
# )
