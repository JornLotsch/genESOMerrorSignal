###############################################################################
# Experiment 3: Mouse Lipidomics Dataset Analysis
# Description:
#   Analyzes a real-world mouse lipidomics dataset. Steps: scaling, data splitting,
#   statistical significance analysis, feature importance analysis, and visualization.
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
library(effectsize)
library(scales)

# --- Load and Scale Mouse Lipidomics Dataset ---------------------------------
path_data <- "/home/joern/Aktuell/GenerativeESOM/08AnalyseProgramme/R/genESOMerrorSignal/MouseEAE_DIB/"

mouse_lipidomics_data_all <- read.csv(
  paste0(path_data, "mouse_lipidomics_data_transformed_imputed.csv"),
  row.names = 1
)

mouse_lipidomics_metadata <- read.csv(  paste0(path_data, "mouse_lipidomics_metadata.csv"))

mouse_lipidomics_data <- cbind.data.frame(Target = mouse_lipidomics_metadata$GROUP, mouse_lipidomics_data_all)

# Scale features (excluding target column)
mouse_lipidomics_data[,-1] <- apply(mouse_lipidomics_data[,-1], 2, scale)

# --- Split Dataset: Training/Test/Validation (Not cuurently used) ---------------------------------
max_cores <- NULL
nProc <- min(parallel::detectCores() - 1, max_cores)
seed <- 42

split_mouse_lipidomics <- opdisDownsampling::opdisDownsampling(
  Data = within(mouse_lipidomics_data, rm(Target)),
  Cls = mouse_lipidomics_data$Target,
  Size = 0.8 * nrow(mouse_lipidomics_data),
  Seed = seed,
  nTrials = 10000,
  MaxCores = nProc
)
mouse_lipidomics_TrainingTest <- mouse_lipidomics_data[
  rownames(mouse_lipidomics_data) %in% split_mouse_lipidomics$ReducedInstances,]
mouse_lipidomics_Validation <- mouse_lipidomics_data[
  !rownames(mouse_lipidomics_data) %in% split_mouse_lipidomics$ReducedInstances,]

# --- Statistical Significance Analysis ---------------------------------------
# Perform ANOVA for each variable (excluding target)
p_vals_anova <- apply(
  mouse_lipidomics_data[,-1], 2,
  function(x) summary(aov(x ~ mouse_lipidomics_data$Target))[[1]][["Pr(>F)"]][1]
)

# Identify significant variables (example list)
mouse_lipidomics_data_significant_vars <- c(
  "PL_LPA16_0", "PL_LPA20_4", "PFC_Sphingosin",
  "PFC_C18Ceramid", "PFC_Sphinganin"
)

# Prepare data for plotting
df_mouse_lipidomics_data_p_vals <- cbind.data.frame(
  p.value = p_vals_anova,
  Original_significant = 0
)
df_mouse_lipidomics_data_p_vals$Lipid <- rownames(df_mouse_lipidomics_data_p_vals)
df_mouse_lipidomics_data_p_vals$Original_significant[
  df_mouse_lipidomics_data_p_vals$Lipid %in% mouse_lipidomics_data_significant_vars
] <- 1

# --- Visualize Significance --------------------------------------------------
barplot_mouse_lipidomics_data_significant_vars <- ggplot(
  df_mouse_lipidomics_data_p_vals,
  aes(x = -log10(p.value), y = reorder(Lipid, -log10(p.value)), fill = factor(Original_significant))
) +
  geom_bar(stat = "identity", color = "#8C5C00", alpha=0.3) +
  geom_vline(xintercept = -log10(0.05), color = "salmon", linetype = "dashed") +
  annotate(
    "text", x = -log10(0.05), y = 1,
    label = "p = 0.05", hjust = -0.1, vjust = 1.5, color = "salmon",
    fontface = "plain", angle = 90
  ) +  
  scale_fill_manual(values = c("0" = "#B37500", "1" = "#B37500")) +
  theme_light() +
  theme(
    legend.position = "none",
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    axis.text.x = element_text(angle = 0, vjust = 0.5, hjust = 1)
  )

print(barplot_mouse_lipidomics_data_significant_vars)

# --- Feature Importance Analysis ---------------------------------------------
DataSetSizes <- c(
  "original", "engineered_0",
  "augmented_1_engineered", "augmented_5_engineered"
)
radius_gen_mouse_lipidomics <- 1.355932203  # From separate assessments
results_varimp_mouse_lipidomics <- analyze_variable_importance(
  data = mouse_lipidomics_data,
  class_name = "Target",
  data_reduced = mouse_lipidomics_TrainingTest,
  DataSetSizes = DataSetSizes,
  density_radius = radius_gen_mouse_lipidomics,
  show_varimp_limit = TRUE,
  mark_sig = FALSE,
  sort_circular = TRUE
)

# --- Plot Variable Selection Frequency ---------------------------------------
plot_mouse_lipidomics_data_selfreq <- cowplot::plot_grid(
  barplot_mouse_lipidomics_data_significant_vars,
  results_varimp_mouse_lipidomics$original$p_selection_freq + labs(title = "Original"),
  results_varimp_mouse_lipidomics$engineered_0$p_selection_freq + labs(title = "Engineered 0"),
  results_varimp_mouse_lipidomics$augmented_1_engineered$p_selection_freq + labs(title = "Augmented 1, engineered"),
  results_varimp_mouse_lipidomics$augmented_5_engineered$p_selection_freq + labs(title = "Augmented 5, engineered"),
  labels = "AUTO", nrow = 1, align = "h", axis = "tb"
) +
  plot_annotation(
    title = "Variable selection frequency",
    subtitle = "Dataset: mouse_lipidomics_data"
  ) &
  theme(
    plot.tag.position = c(0.5, 1),
    plot.tag = element_text(size = 14, face = "bold", vjust = 0)
  )
print(plot_mouse_lipidomics_data_selfreq)
ggsave(
  "plot_mouse_lipidomics_data_selection_freq.svg",
  plot_mouse_lipidomics_data_selfreq,
  width = 22, height = 10, limitsize = FALSE
)


# --- Correlation Analysis and Effect Size Visualization -------------
# Calculate Kendall correlations between p-values and feature importance for each DataSetSize

cor_results <- data.frame(
  DataSetSize = character(),
  Correlation = numeric(),
  p.value = numeric(),
  stringsAsFactors = FALSE
)

for (ds in DataSetSizes) {
  # Dynamically access the results object for each dataset size
  df2 <- results_varimp_mouse_lipidomics[[ds]]$feature_importance$df_features
  df1 <- df_mouse_lipidomics_data_p_vals
  merged_df <- merge(df1, df2[, c("Var", "SelectedTrueCorr")], 
                     by.x = "Lipid", by.y = "Var", all.x = TRUE)
  ct <- cor.test(merged_df$p.value, merged_df$SelectedTrueCorr, method = "kendall")
  cor_results <- rbind(
    cor_results,
    data.frame(
      DataSetSize = ds,
      Correlation = ct$estimate,
      p.value = ct$p.value
    )
  )
}

# Ensure DataSetSize is a factor with reversed levels of DataSetSizes for plotting order
cor_results$DataSetSize <- factor(
  cor_results$DataSetSize,
  levels = rev(DataSetSizes)
)

# --- Effect Size Interpretation ----------------------------------------------
# Option 1: Funder & Ozer (2019) via effectsize::interpret_r (default, recommended)
cor_results$EffectSizeLabel <- effectsize::interpret_r(cor_results$Correlation, rules = "funder2019")

# Option 2: Cohen (1988) style, more common for tau/rho (commented out)
# cor_results$EffectSizeLabel <- cut(
#   abs(cor_results$Correlation),
#   breaks = c(-Inf, 0.1, 0.3, 0.5, Inf),
#   labels = c("very small", "small", "medium", "large")
# )

# --- Visualize Correlation Effect Sizes with Annotation ----------------------
# Prepare annotation labels
cor_results$label <- sprintf("Tau = %.2f\np = %.3g", cor_results$Correlation, cor_results$p.value)

# Find the minimum correlation (for placing text just outside the bar)
# If your correlations are negative, this ensures the label is visible
cor_results$label_x <- cor_results$Correlation - 0.05 * sign(cor_results$Correlation)

p_correlations_p_versus_var_freq <- 
  ggplot(cor_results, aes(y = DataSetSize, x = Correlation, fill = EffectSizeLabel)) +
  geom_col(width = 0.7) +
  # Add annotation for Tau and p-value
  geom_text(
    aes(x = Correlation, label = label),
    hjust = ifelse(cor_results$Correlation < 0, 1.05, -0.05), 
    vjust = 0.5,
    size = 4,
    color = "ghostwhite"
  ) +
  scale_fill_manual(
    values = c(
      "tiny" = "grey80",
      "very small" = "skyblue",
      "small" = "lightblue",
      "medium" = "orange",
      "large" = "red",
      "very large" = "darkred"
    )
  ) +
  scale_x_reverse() +
  theme_light() + 
  theme(
    legend.position = c(.8, .1),
    legend.background = element_rect(fill = alpha("white", 0.5), color = NA),
    legend.box.background = element_rect(fill = alpha("white", 0.5), color = NA)
  ) +
  labs(
    title = "Correlation Effect Sizes by DataSetSize",
    y = "Data Set Size",
    x = "Kendall Correlation"
  )

print(p_correlations_p_versus_var_freq)

ggsave(
  "p_correlations_p_versus_var_freq.svg",
  p_correlations_p_versus_var_freq,
  width = 6, height = 6, limitsize = FALSE
)

# --- Annotate Variable Selection Frequency Plots with Correlation Results ----

# Prepare annotation text for each DataSetSize
cor_results$annotation <- sprintf("Tau = %.2f\np = %.3g", cor_results$Correlation, cor_results$p.value)

# Function to add annotation to a plot
add_correlation_annotation <- function(plot, label) {
  plot + annotate(
    "text",
    x = Inf, y = Inf, hjust = 1.1, vjust = 1.1,
    label = label,
    size = 4,
    fontface = "italic"
  )
}

# Add annotation to each plot
p_sel_freq_annotated <- list(
  barplot_mouse_lipidomics_data_significant_vars,
  add_correlation_annotation(
    results_varimp_mouse_lipidomics$original$p_selection_freq + labs(title = "Original"),
    cor_results$annotation[cor_results$DataSetSize == "original"]
  ),
  add_correlation_annotation(
    results_varimp_mouse_lipidomics$engineered_0$p_selection_freq + labs(title = "Engineered 0"),
    cor_results$annotation[cor_results$DataSetSize == "engineered_0"]
  ),
  add_correlation_annotation(
    results_varimp_mouse_lipidomics$augmented_1_engineered$p_selection_freq + labs(title = "Augmented 1, engineered"),
    cor_results$annotation[cor_results$DataSetSize == "augmented_1_engineered"]
  ),
  add_correlation_annotation(
    results_varimp_mouse_lipidomics$augmented_5_engineered$p_selection_freq + labs(title = "Augmented 5, engineered"),
    cor_results$annotation[cor_results$DataSetSize == "augmented_5_engineered"]
  )
)

# --- Combine plots as before, now with annotations ---------------------------
plot_mouse_lipidomics_data_selfreq_annotated <- cowplot::plot_grid(
  plotlist = p_sel_freq_annotated,
  labels = "AUTO", nrow = 1, align = "h", axis = "tb"
) +
  plot_annotation(
    title = "Variable selection frequency",
    subtitle = "Dataset: mouse_lipidomics_data"
  ) &
  theme(
    plot.tag.position = c(0.5, 1),
    plot.tag = element_text(size = 14, face = "bold", vjust = 0)
  )

print(plot_mouse_lipidomics_data_selfreq_annotated)
ggsave(
  "plot_mouse_lipidomics_data_selection_freq_annotated.svg",
  plot_mouse_lipidomics_data_selfreq_annotated,
  width = 22, height = 10, limitsize = FALSE
)
