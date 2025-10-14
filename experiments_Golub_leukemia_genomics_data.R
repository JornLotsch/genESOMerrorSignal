###########################################################################################
# Experiment 8: Cancer Genetics "Golub_data" Dataset  Analysis
# Description:
#   Analyzes a real-world genomics dataset from
#   https://bioconductor.org/packages/golubEsets, reduced to 150 features 
#   according to 
#   https://rstudio-pubs-static.s3.amazonaws.com/3773_0afaead59a02436889abc68753e6c20a.html
#   Steps: scaling, data splitting, statistical significance analysis,
#   feature importance analysis, and visualization.
###########################################################################################

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
path_data <- "/home/joern/Aktuell/GenerativeESOM/08AnalyseProgramme/R/genESOMerrorSignal/Golub/"

Golub_cancer_genomics_data <- read.csv(
  paste0(path_data, "GolubData.csv")
)

dim(Golub_cancer_genomics_data)
Golub_cancer_genomics_data <- na.omit(Golub_cancer_genomics_data)
dim(Golub_cancer_genomics_data)
table(Golub_cancer_genomics_data$Target)
str(Golub_cancer_genomics_data)

Golub_cancer_genomics_data$Target <- as.integer(as.factor(Golub_cancer_genomics_data$Target))

# --- Explore Feature Distributions and Transformations -----------------------
class_name <- "Target"
class_column <- Golub_cancer_genomics_data[[class_name]]
transformation_methods <- c("none", "log10", "sqrt", "reciprocal", "boxcox")

distribution_results <- explore_distribution(
  data = Golub_cancer_genomics_data,
  classes = class_column,
  transformation_methods = transformation_methods,
  plot_results = FALSE, 
  max_vars = 200
)

# Print summary of best transformations
best_transforms <- distribution_results[distribution_results$Best == "*",]
verbose("Best transformations by variable:")
print(best_transforms[, c("Variable", "Transformation", "AD_P_Value")])
if (best_transforms$Transformation[best_transforms$Variable == "Target"] != "none") stop("Target varibale to be transformed. Please check.")

# --- Transform and Scale Features --------------------------------------------

Golub_cancer_genomics_data_transformed <- apply_var_wise_best_tukey_transformation(data = Golub_cancer_genomics_data, best_transforms = best_transforms)

# Scale features (excluding target column)
Golub_cancer_genomics_data_transformed[, -1] <- apply(
  Golub_cancer_genomics_data_transformed[, -1], 2, scale
)

# --- Visualize Transformed and Scaled Data -----------------------------------
Golub_cancer_genomics_data_transformed_long <- reshape2::melt(
  Golub_cancer_genomics_data_transformed, id.vars = "Target"
)
dim(Golub_cancer_genomics_data_transformed)

# --- Feature Importance Analysis: Boruta -------------------------------------
set.seed(42)
all_boruta <- Boruta::Boruta(
  Target ~ ., data = Golub_cancer_genomics_data_transformed,
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

p_Golub_cancer_genomics_boruta <- ggplot(
  all_boruta_long,
  aes(x = reorder(Var2, value), y = value, color = color, fill = color)
) +
  geom_boxplot(alpha = .3) +
  theme_light() +
  labs(
    title = "Golub Leukemia Genomics Feature Importance",
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

print(p_Golub_cancer_genomics_boruta)

# --- Split Dataset: Training/Test/Validation (Not currently used) ------------
max_cores <- NULL
nProc <- min(parallel::detectCores() - 1, max_cores)
seed <- 42

split_Golub_cancer_genomics <- opdisDownsampling::opdisDownsampling(
  Data = within(Golub_cancer_genomics_data_transformed, rm(Target)),
  Cls = Golub_cancer_genomics_data_transformed$Target,
  Size = 0.8 * nrow(Golub_cancer_genomics_data_transformed),
  Seed = seed,
  nTrials = 10000,
  MaxCores = nProc
)
Golub_cancer_genomics_TrainingTest <- Golub_cancer_genomics_data_transformed[
  rownames(Golub_cancer_genomics_data_transformed) %in% split_Golub_cancer_genomics$ReducedInstances,]
Golub_cancer_genomics_Validation <- Golub_cancer_genomics_data_transformed[
  !rownames(Golub_cancer_genomics_data_transformed) %in% split_Golub_cancer_genomics$ReducedInstances,]

# --- Statistical Significance Analysis ---------------------------------------
# Perform t-test for each variable (excluding target)
p_vals_ttest <- apply(
  Golub_cancer_genomics_data_transformed[, -1], 2,
  function(x) t.test(x ~ as.factor(Golub_cancer_genomics_data_transformed$Target))$p.value
)

# Identify significant variables (example list)
Golub_cancer_genomics_significant_vars <- names(p_vals_ttest[p_vals_ttest < 0.05])

# Prepare data for plotting
df_Golub_cancer_genomics_p_vals <- cbind.data.frame(
  p.value = p_vals_ttest,
  Original_significant = 0
)
df_Golub_cancer_genomics_p_vals$Feature <- rownames(df_Golub_cancer_genomics_p_vals)
df_Golub_cancer_genomics_p_vals$Original_significant[
  df_Golub_cancer_genomics_p_vals$Feature %in% Golub_cancer_genomics_significant_vars
] <- 1

# --- Visualize Significance --------------------------------------------------
barplot_Golub_cancer_genomics_significant_vars <- ggplot(
  df_Golub_cancer_genomics_p_vals,
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

print(barplot_Golub_cancer_genomics_significant_vars)

# --- U-Matrix and Density Radius Calculation ---------------------------------
set.seed(seed)
umx_Golub_cancer_genomics <- Umatrix::esomTrain(
  Data = as.matrix(within(Golub_cancer_genomics_data_transformed, rm(Target))),
  Lines = 80, Columns = 50, Toroid = TRUE, Epochs = 100
)

radius_umx_Golub_cancer_genomics <- Umatrix::calculate_Delauny_radius(
  Data = as.matrix(within(Golub_cancer_genomics_data_transformed, rm(Target))),
  BestMatches = umx_Golub_cancer_genomics$BestMatches,
  Columns = 50, Lines = 80, Toroid = TRUE
)
radius_Golub_cancer_genomics <- radius_umx_Golub_cancer_genomics$RadiusByEM

# --- Feature Importance Analysis ---------------------------------------------
DataSetSizes <- c(
  "original", "engineered_0",
  "augmented_1_engineered", "augmented_5_engineered"
)

results_varimp_Golub_cancer_genomics <- analyze_variable_importance(
  data = Golub_cancer_genomics_data_transformed,
  class_name = "Target",
  data_reduced = Golub_cancer_genomics_TrainingTest,
  DataSetSizes = DataSetSizes,
  density_radius = radius_Golub_cancer_genomics,
  show_varimp_limit = TRUE,
  mark_sig = FALSE,
  sort_circular = FALSE
)

# --- Plot Variable Selection Frequency ---------------------------------------
plot_Golub_cancer_genomics_selfreq <- cowplot::plot_grid(
  barplot_Golub_cancer_genomics_significant_vars,
  results_varimp_Golub_cancer_genomics$original$p_selection_freq + labs(title = "Original"),
  results_varimp_Golub_cancer_genomics$engineered_0$p_selection_freq + labs(title = "Engineered 0"),
  results_varimp_Golub_cancer_genomics$augmented_1_engineered$p_selection_freq + labs(title = "Augmented 1, engineered"),
  results_varimp_Golub_cancer_genomics$augmented_5_engineered$p_selection_freq + labs(title = "Augmented 5, engineered"),
  labels = "AUTO", nrow = 1, align = "h", axis = "tb"
) +
  plot_annotation(
    title = "Variable selection frequency",
    subtitle = "Dataset: Golub Cancer Genomics Dataset"
  ) &
  theme(
    plot.tag.position = c(0.5, 1),
    plot.tag = element_text(size = 14, face = "bold", vjust = 0)
  )

print(plot_Golub_cancer_genomics_selfreq)
ggsave(
  "plot_Golub_cancer_genomics_data_selection_freq.svg",
  plot_Golub_cancer_genomics_selfreq,
  width = 22, height = 10, limitsize = FALSE
)

# --- Correlation Analysis and Effect Size Visualization ----------------------
# Calculate Kendall correlations between p-values and feature importance for each DataSetSize

cor_results_Golub_cancer_genomics <- data.frame(
  DataSetSize = character(),
  Correlation = numeric(),
  p.value = numeric(),
  stringsAsFactors = FALSE
)

for (ds in DataSetSizes) {
  df2 <- results_varimp_Golub_cancer_genomics[[ds]]$feature_importance$df_features
  df1 <- df_Golub_cancer_genomics_p_vals
  merged_df <- merge(df1, df2[, c("Var", "SelectedTrueCorr")],
                     by.x = "Feature", by.y = "Var", all.x = TRUE)
  ct <- cor.test(merged_df$p.value, merged_df$SelectedTrueCorr, method = "kendall")
  cor_results_Golub_cancer_genomics <- rbind(
    cor_results_Golub_cancer_genomics,
    data.frame(
      DataSetSize = ds,
      Correlation = ct$estimate,
      p.value = ct$p.value
    )
  )
}

# Ensure DataSetSize is a factor with reversed levels of DataSetSizes for plotting order
cor_results_Golub_cancer_genomics$DataSetSize <- factor(
  cor_results_Golub_cancer_genomics$DataSetSize,
  levels = rev(DataSetSizes)
)

# --- Effect Size Interpretation ----------------------------------------------
cor_results_Golub_cancer_genomics$EffectSizeLabel <- effectsize::interpret_r(
  cor_results_Golub_cancer_genomics$Correlation, rules = "funder2019"
)

# --- Visualize Correlation Effect Sizes with Annotation ----------------------
cor_results_Golub_cancer_genomics$label <- sprintf("Tau = %.2f\np = %.3g", cor_results_Golub_cancer_genomics$Correlation, cor_results_Golub_cancer_genomics$p.value)
cor_results_Golub_cancer_genomics$label_x <- cor_results_Golub_cancer_genomics$Correlation - 0.05 * sign(cor_results_Golub_cancer_genomics$Correlation)

p_correlations_p_versus_var_freq <- ggplot(
  cor_results_Golub_cancer_genomics, aes(y = DataSetSize, x = Correlation, fill = EffectSizeLabel)
) +
  geom_col(width = 0.7) +
  geom_text(
    aes(x = Correlation, label = label),
    hjust = ifelse(cor_results_Golub_cancer_genomics$Correlation < 0, 1.05, -0.05),
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
  "p_correlations_p_versus_var_freq_Golub_cancer_genomics.svg",
  p_correlations_p_versus_var_freq,
  width = 6, height = 6, limitsize = FALSE
)

# --- Annotate Variable Selection Frequency Plots with Correlation Results ----
cor_results_Golub_cancer_genomics$annotation <- sprintf("Tau = %.2f\np = %.3g", cor_results_Golub_cancer_genomics$Correlation, cor_results_Golub_cancer_genomics$p.value)

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
  barplot_Golub_cancer_genomics_significant_vars,
  add_correlation_annotation(
    results_varimp_Golub_cancer_genomics$original$p_selection_freq + labs(title = "Original"),
    cor_results_Golub_cancer_genomics$annotation[cor_results_Golub_cancer_genomics$DataSetSize == "original"]
  ),
  add_correlation_annotation(
    results_varimp_Golub_cancer_genomics$engineered_0$p_selection_freq + labs(title = "Engineered 0"),
    cor_results_Golub_cancer_genomics$annotation[cor_results_Golub_cancer_genomics$DataSetSize == "engineered_0"]
  ),
  add_correlation_annotation(
    results_varimp_Golub_cancer_genomics$augmented_1_engineered$p_selection_freq + labs(title = "Augmented 1, engineered"),
    cor_results_Golub_cancer_genomics$annotation[cor_results_Golub_cancer_genomics$DataSetSize == "augmented_1_engineered"]
  ),
  add_correlation_annotation(
    results_varimp_Golub_cancer_genomics$augmented_5_engineered$p_selection_freq + labs(title = "Augmented 5, engineered"),
    cor_results_Golub_cancer_genomics$annotation[cor_results_Golub_cancer_genomics$DataSetSize == "augmented_5_engineered"]
  )
)

plot_Golub_cancer_genomics_selfreq_annotated <- cowplot::plot_grid(
  plotlist = p_sel_freq_annotated,
  labels = "AUTO", nrow = 1, align = "h", axis = "tb"
) +
  plot_annotation(
    title = "Variable selection frequency",
    subtitle = "Dataset: Golub_cancer_genomics"
  ) &
  theme(
    plot.tag.position = c(0.5, 1),
    plot.tag = element_text(size = 14, face = "bold", vjust = 0)
  )

print(plot_Golub_cancer_genomics_selfreq_annotated)
ggsave(
  "plot_Golub_cancer_genomics_selection_freq_annotated.svg",
  plot_Golub_cancer_genomics_selfreq_annotated,
  width = 22, height = 10, limitsize = FALSE
)
