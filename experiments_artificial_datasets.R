###############################################################################
# Experiments 1 and 2: Ascending Statistical Significance of Class Differences
# Description:
#   Analyzes an artificial dataset with systematically increasing class differences
#   for each variable. Steps: t-tests, visualization, data splitting, and feature
#   importance analysis.
###############################################################################

# --- Setup Environment -------------------------------------------------------
setwd("/home/joern/Aktuell/GenerativeESOM/08AnalyseProgramme/R/genESOMerrorSignal")

# Source helper scripts
source("generate_synthetic_data.R")
source("analyze_variable_importance.R")
source("generate_artifical_datasets.R")

# Load required libraries
library(reshape2)
library(ggplot2)
library(parallel)
library(opdisDownsampling)
library(cowplot)
library(patchwork)

# --- Load and Scale Artificial Datasets --------------------------------------
ascending_significance_data <- read.csv("ascending_significance_data.csv")
no_effect_data <- read.csv("no_effect_data.csv")

# Scale features (excluding target column)
ascending_significance_data[, -1] <- apply(ascending_significance_data[, -1], 2, scale)
no_effect_data[, -1] <- apply(no_effect_data[, -1], 2, scale)

###############################################################################
# Analysis: Ascending Significance Dataset
###############################################################################

# --- Statistical Significance Analysis ---------------------------------------
# Perform t-tests for each variable (excluding target)
t_tests_p_ascending <- apply(
  ascending_significance_data[, -1], 2,
  function(x) t.test(x ~ ascending_significance_data$Target)$p.value
)

# Prepare data for plotting
t_tests_p_ascending_long <- data.frame(
  variable = names(t_tests_p_ascending),
  value = as.numeric(t_tests_p_ascending)
)
t_tests_p_ascending_long$variable <- factor(
  t_tests_p_ascending_long$variable, levels = paste0("X", 1:50)
)

# --- Visualize Significance --------------------------------------------------
barplot_t_tests_p_ascending <- ggplot(
  t_tests_p_ascending_long, aes(y = variable, x = -log10(value))
) +
  geom_bar(stat = "identity", color = "#8C5C00", fill = "#B37500", alpha = 0.3) +
  geom_vline(xintercept = -log10(0.05), color = "salmon", linetype = "dashed") +
  annotate(
    "text", x = -log10(0.05), y = t_tests_p_ascending_long$variable[1],
    label = "p = 0.05", hjust = -0.1, vjust = 1.5, color = "salmon",
    fontface = "plain", angle = 90
  ) +
  theme_light() +
  labs(title = "Significance of class differences")
print(barplot_t_tests_p_ascending)

# --- Split Dataset: Training/Test/Validation ---------------------------------
# Use opdisDownsampling for splitting
max_cores <- NULL
nProc <- min(parallel::detectCores() - 1, max_cores)
seed <- 42

split_asc <- opdisDownsampling::opdisDownsampling(
  Data = within(ascending_significance_data, rm(Target)),
  Cls = ascending_significance_data$Target,
  Size = 0.8 * nrow(ascending_significance_data),
  Seed = seed,
  nTrials = 10000,
  MaxCores = nProc
)
asc_TrainingTest <- ascending_significance_data[
  rownames(ascending_significance_data) %in% split_asc$ReducedInstances,]
asc_Validation <- ascending_significance_data[
  !rownames(ascending_significance_data) %in% split_asc$ReducedInstances,]

# --- Feature Importance Analysis ---------------------------------------------
DataSetSizes <- c("original", "engineered_0", "augmented_1_engineered",
                  "augmented_5_engineered", "augmented_50_engineered")
radius_gen_asc <- 3.5 #From separate assessments
results_varimp_asc <- analyze_variable_importance(
  data = ascending_significance_data,
  class_name = "Target",
  data_reduced = asc_TrainingTest,
  DataSetSizes = DataSetSizes,
  density_radius = radius_gen_asc,
  show_varfreq_limit = FALSE,
  show_varimp_limit = TRUE
)

# --- Plot Variable Importance and Selection Frequency ------------------------
plot_asc_varimp <- cowplot::plot_grid(
  results_varimp_asc$original$p_importance + labs(title = "Original"),
  results_varimp_asc$augmented_5_engineered$p_importance + labs(title = "Augmented 5, engineered"),
  labels = "AUTO", nrow = 1, align = "h", axis = "tb", rel_widths = c(1, 2)
) +
  plot_annotation(title = "Variable importance", subtitle = "Dataset: ascending_significance_data") &
  theme(plot.tag.position = c(0.5, 1), plot.tag = element_text(size = 14, face = "bold", vjust = 0))
print(plot_asc_varimp)
ggsave("plot_ascending_significance_data_var_importance.svg", plot_asc_varimp, width = 20, height = 8, limitsize = FALSE)

plot_asc_selfreq <- cowplot::plot_grid(
  barplot_t_tests_p_ascending,
  results_varimp_asc$original$p_selection_freq + labs(title = "Original"),
  results_varimp_asc$engineered_0$p_selection_freq + labs(title = "Engineered 0"),
  results_varimp_asc$augmented_1_engineered$p_selection_freq + labs(title = "Augmented 1, engineered"),
  results_varimp_asc$augmented_5_engineered$p_selection_freq + labs(title = "Augmented 5, engineered"),
  results_varimp_asc$augmented_50_engineered$p_selection_freq + labs(title = "Augmented 50, engineered"),
  labels = "AUTO", nrow = 1, align = "h", axis = "tb"
) +
  patchwork::plot_annotation(title = "Variable selection frequency", subtitle = "Dataset: ascending_significance_data") &
  theme(plot.tag.position = c(0.5, 1), plot.tag = element_text(size = 14, face = "bold", vjust = 0))
print(plot_asc_selfreq)
ggsave("plot_ascending_significance_data_selection_freq.svg", plot_asc_selfreq, width = 22, height = 10, limitsize = FALSE)

###############################################################################
# Analysis: No Effect Dataset
###############################################################################

# --- Statistical Significance Analysis ---------------------------------------
t_tests_p_no_effect <- apply(
  no_effect_data[, -1], 2,
  function(x) t.test(x ~ no_effect_data$Target)$p.value
)
t_tests_p_no_effect_long <- data.frame(
  variable = names(t_tests_p_no_effect),
  value = as.numeric(t_tests_p_no_effect)
)
t_tests_p_no_effect_long$variable <- factor(
  t_tests_p_no_effect_long$variable, levels = paste0("X", 1:50)
)

# --- Visualize Significance --------------------------------------------------
barplot_t_tests_p_no_effect <- ggplot(
  t_tests_p_no_effect_long, aes(y = variable, x = -log10(value))
) +
  geom_bar(stat = "identity", color = "#8C5C00", fill = "#B37500", alpha = 0.3) +
  geom_vline(xintercept = -log10(0.05), color = "salmon", linetype = "dashed") +
  annotate(
    "text", x = -log10(0.05), y = t_tests_p_no_effect_long$variable[1],
    label = "p = 0.05", hjust = -0.1, vjust = 1.5, color = "salmon",
    fontface = "plain", angle = 90
  ) +
  theme_light() +
  labs(title = "Significance of class differences")
print(barplot_t_tests_p_no_effect)

# --- Split Dataset: Training/Test/Validation ---------------------------------
split_noeff <- opdisDownsampling::opdisDownsampling(
  Data = within(no_effect_data, rm(Target)),
  Cls = no_effect_data$Target,
  Size = 0.8 * nrow(no_effect_data),
  Seed = seed,
  nTrials = 10000,
  MaxCores = nProc
)
noeff_TrainingTest <- no_effect_data[
  rownames(no_effect_data) %in% split_noeff$ReducedInstances,]
noeff_Validation <- no_effect_data[
  !rownames(no_effect_data) %in% split_noeff$ReducedInstances,]

# --- Feature Importance Analysis ---------------------------------------------
DataSetSizes <- c("original", "engineered_0", "augmented_1_engineered", "augmented_5_engineered")
radius_gen_noeff <- 6.6 #From separate assessments
results_varimp_noeff <- analyze_variable_importance(
  data = no_effect_data,
  class_name = "Target",
  data_reduced = noeff_TrainingTest,
  DataSetSizes = DataSetSizes,
  density_radius = radius_gen_noeff,
  show_varfreq_limit = TRUE,
  show_varimp_limit = TRUE,
  mark_sig = FALSE,
  sort_circular = TRUE
)

# --- Plot Selection Frequency ------------------------------------------------
plot_noeff_selfreq <- cowplot::plot_grid(
  barplot_t_tests_p_no_effect,
  results_varimp_noeff$original$p_selection_freq + labs(title = "Original"),
  results_varimp_noeff$engineered_0$p_selection_freq + labs(title = "Engineered 0"),
  results_varimp_noeff$augmented_1_engineered$p_selection_freq + labs(title = "Augmented 1, engineered"),
  results_varimp_noeff$augmented_5_engineered$p_selection_freq + labs(title = "Augmented 5, engineered"),
  labels = "AUTO", nrow = 1, align = "h", axis = "tb"
) +
  plot_annotation(title = "Variable selection frequency", subtitle = "Dataset: no_effect_data") &
  theme(plot.tag.position = c(0.5, 1), plot.tag = element_text(size = 14, face = "bold", vjust = 0))
print(plot_noeff_selfreq)
ggsave("plot_no_effect_data_selection_freq.svg", plot_noeff_selfreq, width = 22, height = 10, limitsize = FALSE)

# --- Plot Selection Frequency (circular) -------------------------------------

# Create the barplot panel (A)
panel_A <- cowplot::plot_grid(
  barplot_t_tests_p_no_effect,
  labels = "A"
)

# Create the 2x2 grid for circular selection frequency plots (B-E)
panel_B_to_E <- cowplot::plot_grid(
  results_varimp_noeff$original$p_selection_freq_circular + labs(title = "Original"),
  results_varimp_noeff$engineered_0$p_selection_freq_circular + labs(title = "Engineered 0"),
  results_varimp_noeff$augmented_1_engineered$p_selection_freq_circular + labs(title = "Augmented 1, engineered"),
  results_varimp_noeff$augmented_5_engineered$p_selection_freq_circular + labs(title = "Augmented 5, engineered"),
  labels = LETTERS[2:5], # B, C, D, E
  nrow = 2,
  align = "v",
  axis = "lr"
)

# Combine A and B-E panels into the final layout
plot_noeff_selfreq_circular <- cowplot::plot_grid(
  panel_A,
  panel_B_to_E,
  nrow = 1,
  align = "h",
  rel_widths = c(1, 3)
) +
  patchwork::plot_annotation(
    title = "Variable selection frequency",
    subtitle = "Dataset: no_effect_data"
  ) &
  theme(
    plot.tag.position = c(0.5, 1),
    plot.tag = element_text(size = 14, face = "bold", vjust = 0, margin = margin(b = -10))
  )

# Print and save the figure
print(plot_noeff_selfreq_circular)
ggsave("plot_no_effect_data_selection_freq_circular.svg", plot_noeff_selfreq_circular, width = 14, height = 10, limitsize = FALSE)
