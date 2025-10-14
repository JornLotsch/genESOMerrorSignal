###############################################################################
# Experiment 5: ILPD (Indian Liver Patient Dataset) Analysis
# Description:
#   Analyzes a real-world heart failure dataset from
#   https://archive.ics.uci.edu/dataset/225/ilpd+indian+liver+patient+dataset
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
path_data <- "/home/joern/Aktuell/GenerativeESOM/08AnalyseProgramme/R/genESOMerrorSignal/ilpd+indian+liver+patient+dataset/"

Indian_liver_patient_data_all <- read.csv(
  paste0(path_data, "Indian Liver Patient Dataset (ILPD).csv")
)

feature_vars <- c(
  "Age", "Gender", "TB", "DB", "Alkphos", "Sgpt", "Sgot", "TP", "ALB", "A.G.Ratio"
)

Indian_liver_patient_data <- cbind.data.frame(
  Target = as.integer(Indian_liver_patient_data_all$Selector),
  Indian_liver_patient_data_all[, feature_vars]
)

Indian_liver_patient_data <- na.omit(Indian_liver_patient_data)
dim(Indian_liver_patient_data)
table(Indian_liver_patient_data$Target)
str(Indian_liver_patient_data)

# --- Explore Feature Distributions and Transformations -----------------------
class_name <- "Target"
class_column <- Indian_liver_patient_data[[class_name]]
transformation_methods <- c("none", "log10", "sqrt", "reciprocal", "boxcox")

distribution_results <- explore_distribution(
  data = Indian_liver_patient_data,
  classes = class_column,
  transformation_methods = transformation_methods,
  plot_results = TRUE
)
par(mfrow = c(1, 1))

# Print summary of best transformations
best_transforms <- distribution_results[distribution_results$Best == "*",]
verbose("Best transformations by variable:")
print(best_transforms[, c("Variable", "Transformation", "AD_P_Value")])
if (best_transforms$Transformation[best_transforms$Variable == "Target"] != "none") stop("Target varibale to be transformed. Please check.")

# --- Transform and Scale Features --------------------------------------------

Indian_liver_patient_data_transformed <- apply_var_wise_best_tukey_transformation(data = Indian_liver_patient_data, best_transforms = best_transforms)

# Scale features (excluding target column)
Indian_liver_patient_data_transformed[, -1] <- apply(
  Indian_liver_patient_data_transformed[, -1], 2, scale
)

# --- Visualize Transformed and Scaled Data -----------------------------------
Indian_liver_patient_data_transformed_long <- reshape2::melt(
  Indian_liver_patient_data_transformed, id.vars = "Target"
)

p_Indian_liver_patient <- ggplot(
  Indian_liver_patient_data_transformed_long,
  aes(x = variable, y = value, color = as.factor(Target), fill = as.factor(Target))
) +
  geom_violin(alpha = .3) +
  geom_jitter() +
  theme_light() +
  ggthemes::scale_color_colorblind() +
  ggthemes::scale_fill_colorblind() +
  labs(
    title = "Indian Liver Patients Dataset",
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
  Target ~ ., data = Indian_liver_patient_data_transformed,
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

p_Indian_liver_patient_boruta <- ggplot(
  all_boruta_long,
  aes(x = reorder(Var2, value), y = value, color = color, fill = color)
) +
  geom_boxplot(alpha = .3) +
  theme_light() +
  labs(
    title = "Indian Liver Patients Data Feature Importance",
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

pp_Indian_liver_patient <- cowplot::plot_grid(
  p_Indian_liver_patient,
  p_Indian_liver_patient_boruta,
  labels = "AUTO",
  align = "tb"
)
print(pp_Indian_liver_patient)

# --- Split Dataset: Training/Test/Validation (Not currently used) ------------
max_cores <- NULL
nProc <- min(parallel::detectCores() - 1, max_cores)
seed <- 42

split_Indian_liver_patient <- opdisDownsampling::opdisDownsampling(
  Data = within(Indian_liver_patient_data_transformed, rm(Target)),
  Cls = Indian_liver_patient_data_transformed$Target,
  Size = 0.8 * nrow(Indian_liver_patient_data_transformed),
  Seed = seed,
  nTrials = 10000,
  MaxCores = nProc
)
Indian_liver_patient_TrainingTest <- Indian_liver_patient_data_transformed[
  rownames(Indian_liver_patient_data_transformed) %in% split_Indian_liver_patient$ReducedInstances,]
Indian_liver_patient_Validation <- Indian_liver_patient_data_transformed[
  !rownames(Indian_liver_patient_data_transformed) %in% split_Indian_liver_patient$ReducedInstances,]

# --- Statistical Significance Analysis ---------------------------------------
# Perform t-test for each variable (excluding target)
p_vals_ttest <- apply(
  Indian_liver_patient_data_transformed[, -1], 2,
  function(x) t.test(x ~ as.factor(Indian_liver_patient_data_transformed$Target))$p.value
)

# Identify significant variables
Indian_liver_patient_significant_vars <- names(p_vals_ttest[p_vals_ttest < 0.05])

# Prepare data for plotting
df_Indian_liver_patient_p_vals <- cbind.data.frame(
  p.value = p_vals_ttest,
  Original_significant = 0
)
df_Indian_liver_patient_p_vals$Feature <- rownames(df_Indian_liver_patient_p_vals)
df_Indian_liver_patient_p_vals$Original_significant[
  df_Indian_liver_patient_p_vals$Feature %in% Indian_liver_patient_significant_vars
] <- 1

# --- Visualize Significance --------------------------------------------------
barplot_Indian_liver_patient_significant_vars <- ggplot(
  df_Indian_liver_patient_p_vals,
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

print(barplot_Indian_liver_patient_significant_vars)

# --- U-Matrix and Density Radius Calculation ---------------------------------
set.seed(seed)
umx_Indian_liver_patient <- Umatrix::esomTrain(
  Data = as.matrix(within(Indian_liver_patient_data_transformed, rm(Target))),
  Lines = 80, Columns = 50, Toroid = TRUE, Epochs = 100
)

radius_umx_Indian_liver_patient <- Umatrix::calculate_Delauny_radius(
  Data = as.matrix(within(Indian_liver_patient_data_transformed, rm(Target))),
  BestMatches = umx_Indian_liver_patient$BestMatches,
  Columns = 50, Lines = 80, Toroid = TRUE
)
radius_Indian_liver_patient <- radius_umx_Indian_liver_patient$RadiusByEM

# --- Feature Importance Analysis ---------------------------------------------
DataSetSizes <- c(
  "original", "engineered_0",
  "augmented_1_engineered", "augmented_5_engineered"
)

results_varimp_Indian_liver_patient <- analyze_variable_importance(
  data = Indian_liver_patient_data_transformed,
  class_name = "Target",
  data_reduced = Indian_liver_patient_TrainingTest,
  DataSetSizes = DataSetSizes,
  density_radius = radius_Indian_liver_patient,
  show_varimp_limit = TRUE,
  mark_sig = FALSE,
  sort_circular = FALSE
)

# --- Plot Variable Selection Frequency ---------------------------------------
plot_Indian_liver_patient_selfreq <- cowplot::plot_grid(
  barplot_Indian_liver_patient_significant_vars,
  results_varimp_Indian_liver_patient$original$p_selection_freq + labs(title = "Original"),
  results_varimp_Indian_liver_patient$engineered_0$p_selection_freq + labs(title = "Engineered 0"),
  results_varimp_Indian_liver_patient$augmented_1_engineered$p_selection_freq + labs(title = "Augmented 1, engineered"),
  results_varimp_Indian_liver_patient$augmented_5_engineered$p_selection_freq + labs(title = "Augmented 5, engineered"),
  labels = "AUTO", nrow = 1, align = "h", axis = "tb"
) +
  plot_annotation(
    title = "Variable selection frequency",
    subtitle = "Dataset: Indian_liver_patient_data"
  ) &
  theme(
    plot.tag.position = c(0.5, 1),
    plot.tag = element_text(size = 14, face = "bold", vjust = 0)
  )

print(plot_Indian_liver_patient_selfreq)
ggsave(
  "plot_Indian_liver_patient_data_selection_freq.svg",
  plot_Indian_liver_patient_selfreq,
  width = 22, height = 10, limitsize = FALSE
)

# --- Correlation Analysis and Effect Size Visualization ----------------------
# Calculate Kendall correlations between p-values and feature importance for each DataSetSize

cor_results_liver <- data.frame(
  DataSetSize = character(),
  Correlation = numeric(),
  p.value = numeric(),
  stringsAsFactors = FALSE
)

for (ds in DataSetSizes) {
  df2 <- results_varimp_Indian_liver_patient[[ds]]$feature_importance$df_features
  df1 <- df_Indian_liver_patient_p_vals
  merged_df <- merge(df1, df2[, c("Var", "SelectedTrueCorr")],
                     by.x = "Feature", by.y = "Var", all.x = TRUE)
  ct <- cor.test(merged_df$p.value, merged_df$SelectedTrueCorr, method = "kendall")
  cor_results_liver <- rbind(
    cor_results_liver,
    data.frame(
      DataSetSize = ds,
      Correlation = ct$estimate,
      p.value = ct$p.value
    )
  )
}

# Ensure DataSetSize is a factor with reversed levels of DataSetSizes for plotting order
cor_results_liver$DataSetSize <- factor(
  cor_results_liver$DataSetSize,
  levels = rev(DataSetSizes)
)

# --- Effect Size Interpretation ----------------------------------------------
cor_results_liver$EffectSizeLabel <- effectsize::interpret_r(
  cor_results_liver$Correlation, rules = "funder2019"
)

# --- Visualize Correlation Effect Sizes with Annotation ----------------------
cor_results_liver$label <- sprintf("Tau = %.2f\np = %.3g", cor_results_liver$Correlation, cor_results_liver$p.value)
cor_results_liver$label_x <- cor_results_liver$Correlation - 0.05 * sign(cor_results_liver$Correlation)

p_correlations_p_versus_var_freq <- ggplot(
  cor_results_liver, aes(y = DataSetSize, x = Correlation, fill = EffectSizeLabel)
) +
  geom_col(width = 0.7) +
  geom_text(
    aes(x = Correlation, label = label),
    hjust = ifelse(cor_results_liver$Correlation < 0, 1.05, -0.05),
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
  "p_correlations_p_versus_var_freq_Indian_liver_patient.svg",
  p_correlations_p_versus_var_freq,
  width = 6, height = 6, limitsize = FALSE
)

# --- Annotate Variable Selection Frequency Plots with Correlation Results ----
cor_results_liver$annotation <- sprintf("Tau = %.2f\np = %.3g", cor_results_liver$Correlation, cor_results_liver$p.value)

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
  barplot_Indian_liver_patient_significant_vars,
  add_correlation_annotation(
    results_varimp_Indian_liver_patient$original$p_selection_freq + labs(title = "Original"),
    cor_results_liver$annotation[cor_results_liver$DataSetSize == "original"]
  ),
  add_correlation_annotation(
    results_varimp_Indian_liver_patient$engineered_0$p_selection_freq + labs(title = "Engineered 0"),
    cor_results_liver$annotation[cor_results_liver$DataSetSize == "engineered_0"]
  ),
  add_correlation_annotation(
    results_varimp_Indian_liver_patient$augmented_1_engineered$p_selection_freq + labs(title = "Augmented 1, engineered"),
    cor_results_liver$annotation[cor_results_liver$DataSetSize == "augmented_1_engineered"]
  ),
  add_correlation_annotation(
    results_varimp_Indian_liver_patient$augmented_5_engineered$p_selection_freq + labs(title = "Augmented 5, engineered"),
    cor_results_liver$annotation[cor_results_liver$DataSetSize == "augmented_5_engineered"]
  )
)

plot_Indian_liver_patient_selfreq_annotated <- cowplot::plot_grid(
  plotlist = p_sel_freq_annotated,
  labels = "AUTO", nrow = 1, align = "h", axis = "tb"
) +
  plot_annotation(
    title = "Variable selection frequency",
    subtitle = "Dataset: Indian_liver_patient_clinical_records"
  ) &
  theme(
    plot.tag.position = c(0.5, 1),
    plot.tag = element_text(size = 14, face = "bold", vjust = 0)
  )

print(plot_Indian_liver_patient_selfreq_annotated)
ggsave(
  "plot_Indian_liver_patient_clinical_records_selection_freq_annotated.svg",
  plot_Indian_liver_patient_selfreq_annotated,
  width = 22, height = 10, limitsize = FALSE
)
