###############################################################################
# Experiment 4: Heart Failure Clinical Records Dataset Analysis
# Description:
#   Analyzes a real-world heart failure dataset from
#   https://archive.ics.uci.edu/dataset/519/heart+failure+clinical+records
#   Steps: scaling, data splitting, statistical significance analysis,
#   feature importance analysis, and visualization.
###############################################################################

# --- Setup Environment -------------------------------------------------------
setwd("/home/joern/Aktuell/GenerativeESOM/08AnalyseProgramme/R/genESOMerrorSignal")

# Source helper scripts
source("generate_synthetic_data.R")
source("analyze_variable_importance.R")
source("generate_artifical_datasets.R")
source("read_and_prepare_data.R")

# Load required libraries
library(reshape2)
library(ggplot2)
library(parallel)
library(opdisDownsampling)
library(cowplot)
library(patchwork)
library(effectsize)
library(scales)
library(Boruta)

# --- Load and Prepare Heart Failure Dataset ----------------------------------
path_data <- "/home/joern/Aktuell/GenerativeESOM/08AnalyseProgramme/R/genESOMerrorSignal/heart+failure+clinical+records/"

heart_failure_data_all <- read.csv(
  paste0(path_data, "heart_failure_clinical_records_dataset.csv")
)

feature_vars <- c(
  "age", "anaemia", "creatinine_phosphokinase", "diabetes", "ejection_fraction",
  "high_blood_pressure", "platelets", "serum_creatinine", "serum_sodium", "sex"
)

heart_failure_data <- cbind.data.frame(
  Target = heart_failure_data_all$DEATH_EVENT,
  heart_failure_data_all[, feature_vars]
)

dim(heart_failure_data)
table(heart_failure_data$Target)
str(heart_failure_data)

# --- Explore Feature Distributions and Transformations -----------------------
class_name <- "Target"
class_column <- heart_failure_data[[class_name]]
transformation_methods <- c("none", "log10", "sqrt", "reciprocal", "boxcox")

distribution_results <- explore_distribution(
  data = heart_failure_data,
  classes = class_column,
  transformation_methods = transformation_methods,
  plot_results = TRUE
)
par(mfrow = c(1, 1))

# Print summary of best transformations
best_transforms <- distribution_results[distribution_results$Best == "*",]
verbose("Best transformations by variable:")
print(best_transforms[, c("Variable", "Transformation", "AD_P_Value")])

# --- Transform and Scale Features --------------------------------------------

heart_failure_data_transformed(data = heart_failure_data, best_transforms = best_transforms)

# Scale features (excluding target column)
heart_failure_data_transformed[, -1] <- apply(
  heart_failure_data_transformed[, -1], 2, scale
)

# --- Visualize Transformed and Scaled Data -----------------------------------
heart_failure_data_transformed_long <- reshape2::melt(
  heart_failure_data_transformed, id.vars = "Target"
)

p_heart_failure <- ggplot(
  heart_failure_data_transformed_long,
  aes(x = variable, y = value, color = as.factor(Target), fill = as.factor(Target))
) +
  geom_violin(alpha = .3) +
  geom_jitter() +
  theme_light() +
  ggthemes::scale_color_colorblind() +
  ggthemes::scale_fill_colorblind() +
  labs(
    title = "Heart Failure Clinical Records Dataset",
    fill = "Event", color = "Event"
  ) +
  theme(
    legend.position.inside = TRUE, legend.position = c(.1, .1),
    legend.background = element_rect(fill = alpha("white", 0.5)),
    legend.key = element_rect(fill = alpha("white", 0.5)),
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)
  )

# --- Feature Importance Analysis: Boruta -------------------------------------
set.seed(42)
all_boruta <- Boruta::Boruta(
  Target ~ ., data = heart_failure_data_transformed,
  pValue = 0.0001, maxRuns = 100
)

Var_decision <- data.frame(all_boruta$finalDecision)
Var_decision$Var <- rownames(Var_decision)

all_boruta_long <- reshape2::melt(all_boruta$ImpHistory)
all_boruta_long$color <- ifelse(
  all_boruta_long$Var2 %in% Var_decision$Var[Var_decision$all_boruta.finalDecision %in% c("Confirmed", "Tentative")], "Chosen",
  ifelse(
    all_boruta_long$Var2 %in% Var_decision$Var[Var_decision$all_boruta.finalDecision == "Rejected"], "Rejected", "Shadow"
  )
)

p_heart_failure_boruta <- ggplot(
  all_boruta_long,
  aes(x = reorder(Var2, value), y = value, color = color, fill = color)
) +
  geom_boxplot(alpha = .3) +
  theme_light() +
  labs(
    title = "Heart Failure Clinical Records Feature Importance",
    fill = "Decision", color = "Decision"
  ) +
  theme(
    legend.position.inside = TRUE, legend.position = c(.1, .9),
    legend.background = element_rect(fill = alpha("white", 0.5)),
    legend.key = element_rect(fill = alpha("white", 0.5)),
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)
  ) +
  scale_color_manual(values = c("chartreuse4", "salmon", "grey80")) +
  scale_fill_manual(values = c("chartreuse4", "salmon", "grey80"))

pp_heart_failure <- cowplot::plot_grid(
  p_heart_failure,
  p_heart_failure_boruta,
  labels = "AUTO",
  align = "tb"
)
print(pp_heart_failure)

# --- Split Dataset: Training/Test/Validation (Not currently used) ------------
max_cores <- NULL
nProc <- min(parallel::detectCores() - 1, max_cores)
seed <- 42

split_heart_failure <- opdisDownsampling::opdisDownsampling(
  Data = within(heart_failure_data_transformed, rm(Target)),
  Cls = heart_failure_data_transformed$Target,
  Size = 0.8 * nrow(heart_failure_data_transformed),
  Seed = seed,
  nTrials = 10000,
  MaxCores = nProc
)
heart_failure_TrainingTest <- heart_failure_data_transformed[
  rownames(heart_failure_data_transformed) %in% split_heart_failure$ReducedInstances,]
heart_failure_Validation <- heart_failure_data_transformed[
  !rownames(heart_failure_data_transformed) %in% split_heart_failure$ReducedInstances,]

# --- Statistical Significance Analysis ---------------------------------------
# Perform t-test for each variable (excluding target)
p_vals_ttest <- apply(
  heart_failure_data_transformed[, -1], 2,
  function(x) t.test(x ~ as.factor(heart_failure_data_transformed$Target))$p.value
)

# Identify significant variables (example list)
heart_failure_significant_vars <- c(
  "ejection_fraction", "serum_creatinine"
)

# Prepare data for plotting
df_heart_failure_p_vals <- cbind.data.frame(
  p.value = p_vals_ttest,
  Original_significant = 0
)
df_heart_failure_p_vals$Feature <- rownames(df_heart_failure_p_vals)
df_heart_failure_p_vals$Original_significant[
  df_heart_failure_p_vals$Feature %in% heart_failure_significant_vars
] <- 1

# --- Visualize Significance --------------------------------------------------
barplot_heart_failure_significant_vars <- ggplot(
  df_heart_failure_p_vals,
  aes(x = -log10(p.value), y = reorder(Feature, - log10(p.value)), fill = factor(Original_significant))
) +
  geom_bar(stat = "identity", color = "#8C5C00", alpha = 0.3) +
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

print(barplot_heart_failure_significant_vars)

# --- U-Matrix and Density Radius Calculation ---------------------------------
set.seed(seed)
umx_heart_failure <- Umatrix::esomTrain(
  Data = as.matrix(within(heart_failure_data_transformed, rm(Target))),
  Lines = 80, Columns = 50, Toroid = TRUE, Epochs = 100
)

radius_umx_heart_failure <- Umatrix::calculate_Delauny_radius(
  Data = as.matrix(within(heart_failure_data_transformed, rm(Target))),
  BestMatches = umx_heart_failure$BestMatches,
  Columns = 50, Lines = 80, Toroid = TRUE
)
radius_heart_failure <- radius_umx_heart_failure$RadiusByEM

# --- Feature Importance Analysis ---------------------------------------------
DataSetSizes <- c(
  "original", "engineered_0",
  "augmented_1_engineered", "augmented_5_engineered"
)

results_varimp_heart_failure <- analyze_variable_importance(
  data = heart_failure_data_transformed,
  class_name = "Target",
  data_reduced = heart_failure_TrainingTest,
  DataSetSizes = DataSetSizes,
  density_radius = radius_heart_failure,
  show_varimp_limit = TRUE,
  mark_sig = FALSE,
  sort_circular = FALSE
)

# --- Plot Variable Selection Frequency ---------------------------------------
plot_heart_failure_selfreq <- cowplot::plot_grid(
  barplot_heart_failure_significant_vars,
  results_varimp_heart_failure$original$p_selection_freq + labs(title = "Original"),
  results_varimp_heart_failure$engineered_0$p_selection_freq + labs(title = "Engineered 0"),
  results_varimp_heart_failure$augmented_1_engineered$p_selection_freq + labs(title = "Augmented 1, engineered"),
  results_varimp_heart_failure$augmented_5_engineered$p_selection_freq + labs(title = "Augmented 5, engineered"),
  labels = "AUTO", nrow = 1, align = "h", axis = "tb"
) +
  plot_annotation(
    title = "Variable selection frequency",
    subtitle = "Dataset: heart_failure_clinical_records"
  ) &
  theme(
    plot.tag.position = c(0.5, 1),
    plot.tag = element_text(size = 14, face = "bold", vjust = 0)
  )

print(plot_heart_failure_selfreq)
ggsave(
  "plot_heart_failure_clinical_records_selection_freq.svg",
  plot_heart_failure_selfreq,
  width = 22, height = 10, limitsize = FALSE
)

# --- Correlation Analysis and Effect Size Visualization ----------------------
# Calculate Kendall correlations between p-values and feature importance for each DataSetSize

cor_results_heart_failure <- data.frame(
  DataSetSize = character(),
  Correlation = numeric(),
  p.value = numeric(),
  stringsAsFactors = FALSE
)

for (ds in DataSetSizes) {
  df2 <- results_varimp_heart_failure[[ds]]$feature_importance$df_features
  df1 <- df_heart_failure_p_vals
  merged_df <- merge(df1, df2[, c("Var", "SelectedTrueCorr")],
                     by.x = "Feature", by.y = "Var", all.x = TRUE)
  ct <- cor.test(merged_df$p.value, merged_df$SelectedTrueCorr, method = "kendall")
  cor_results_heart_failure <- rbind(
    cor_results_heart_failure,
    data.frame(
      DataSetSize = ds,
      Correlation = ct$estimate,
      p.value = ct$p.value
    )
  )
}

# Ensure DataSetSize is a factor with reversed levels of DataSetSizes for plotting order
cor_results_heart_failure$DataSetSize <- factor(
  cor_results_heart_failure$DataSetSize,
  levels = rev(DataSetSizes)
)

# --- Effect Size Interpretation ----------------------------------------------
cor_results_heart_failure$EffectSizeLabel <- effectsize::interpret_r(
  cor_results_heart_failure$Correlation, rules = "funder2019"
)

# --- Visualize Correlation Effect Sizes with Annotation ----------------------
cor_results_heart_failure$label <- sprintf("Tau = %.2f\np = %.3g", cor_results_heart_failure$Correlation, cor_results_heart_failure$p.value)
cor_results_heart_failure$label_x <- cor_results_heart_failure$Correlation - 0.05 * sign(cor_results_heart_failure$Correlation)

p_correlations_p_versus_var_freq <- ggplot(
  cor_results_heart_failure, aes(y = DataSetSize, x = Correlation, fill = EffectSizeLabel)
) +
  geom_col(width = 0.7) +
  geom_text(
    aes(x = Correlation, label = label),
    hjust = ifelse(cor_results_heart_failure$Correlation < 0, 1.05, -0.05),
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
  "p_correlations_p_versus_var_freq_heart_failure.svg",
  p_correlations_p_versus_var_freq,
  width = 6, height = 6, limitsize = FALSE
)

# --- Annotate Variable Selection Frequency Plots with Correlation Results ----
cor_results_heart_failure$annotation <- sprintf("Tau = %.2f\np = %.3g", cor_results_heart_failure$Correlation, cor_results_heart_failure$p.value)

add_correlation_annotation <- function(plot, label) {
  plot + annotate(
    "text",
    x = Inf, y = Inf, hjust = 1.1, vjust = 1.1,
    label = label,
    size = 4,
    fontface = "italic"
  )
}

p_sel_freq_annotated <- list(
  barplot_heart_failure_significant_vars,
  add_correlation_annotation(
    results_varimp_heart_failure$original$p_selection_freq + labs(title = "Original"),
    cor_results_heart_failure$annotation[cor_results_heart_failure$DataSetSize == "original"]
  ),
  add_correlation_annotation(
    results_varimp_heart_failure$engineered_0$p_selection_freq + labs(title = "Engineered 0"),
    cor_results_heart_failure$annotation[cor_results_heart_failure$DataSetSize == "engineered_0"]
  ),
  add_correlation_annotation(
    results_varimp_heart_failure$augmented_1_engineered$p_selection_freq + labs(title = "Augmented 1, engineered"),
    cor_results_heart_failure$annotation[cor_results_heart_failure$DataSetSize == "augmented_1_engineered"]
  ),
  add_correlation_annotation(
    results_varimp_heart_failure$augmented_5_engineered$p_selection_freq + labs(title = "Augmented 5, engineered"),
    cor_results_heart_failure$annotation[cor_results_heart_failure$DataSetSize == "augmented_5_engineered"]
  )
)

plot_heart_failure_selfreq_annotated <- cowplot::plot_grid(
  plotlist = p_sel_freq_annotated,
  labels = "AUTO", nrow = 1, align = "h", axis = "tb"
) +
  plot_annotation(
    title = "Variable selection frequency",
    subtitle = "Dataset: heart_failure_clinical_records"
  ) &
  theme(
    plot.tag.position = c(0.5, 1),
    plot.tag = element_text(size = 14, face = "bold", vjust = 0)
  )

print(plot_heart_failure_selfreq_annotated)
ggsave(
  "plot_heart_failure_clinical_records_selection_freq_annotated.svg",
  plot_heart_failure_selfreq_annotated,
  width = 22, height = 10, limitsize = FALSE
)
