###############################################################################
# Experiment 4: Heart failure clincal records dataset analysis
# Description:
#   Analyzes a real-world heat failure dataset from
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

# --- Load heart failure dataset ---------------------------------
path_data <- "/home/joern/Aktuell/GenerativeESOM/08AnalyseProgramme/R/genESOMerrorSignal/heart+failure+clinical+records/"

heart_failure_clinical_records_dataset_all <- read.csv(
  paste0(path_data, "heart_failure_clinical_records_dataset.csv")
)

feature_vars <-  c("age", "anaemia", "creatinine_phosphokinase", "diabetes", "ejection_fraction", 
                   "high_blood_pressure", "platelets", "serum_creatinine", "serum_sodium", "sex")

heart_failure_clinical_records_dataset <- cbind.data.frame(Target = heart_failure_clinical_records_dataset_all$DEATH_EVENT,
                                                           heart_failure_clinical_records_dataset[,feature_vars])

# --- Explore heart failure dataset ---------------------------------

class_name <- "Target"
class_column <- heart_failure_clinical_records_dataset [[class_name]]
transformation_methods <- c("none", "log10", "sqrt", "reciprocal", "boxcox")

distribution_results <- explore_distribution(
    data = heart_failure_clinical_records_dataset,
    classes = class_column,
    transformation_methods = transformation_methods,
    plot_results = TRUE
  )

# Print summary of best transformations
best_transforms <- distribution_results[distribution_results$Best == "*", ]
verbose("Best transformations by variable:")
print(best_transforms[, c("Variable", "Transformation", "AD_P_Value")])


# --- Transform and scale heart failure dataset ---------------------------------

vars <-best_transforms$Variable
methods <- best_transforms$Transformation

transform_var <- function(x, method) {
  if (method == "none") {
    return(x)
  } else if (method == "boxcox") {
    require(forecast)
    lambda <- BoxCox.lambda(x, method = "loglik")
    return(BoxCox(x, lambda))
  } else if (method == "log10") {
    return(log10(x + 1))
  } else if (method == "sqrt") {
    return(sqrt(x))
  } else if (method == "reciprocal") {
    return(1 / (x + 1e-8))
  } else {
    stop("Unknown transformation")
  }
}

transformed_list <- mapply(
  FUN = transform_var,
  x = heart_failure_clinical_records_dataset,
  method = methods,
  SIMPLIFY = FALSE
)

heart_failure_clinical_records_dataset_transformed <- as.data.frame(transformed_list)

# Scale features (excluding target column)
heart_failure_clinical_records_dataset_transformed[,-1] <- apply(heart_failure_clinical_records_dataset_transformed[,-1], 2, scale)

# Plot transformed and scaled data
heart_failure_clinical_records_dataset_transformed_long <- reshape2::melt(heart_failure_clinical_records_dataset_transformed, id.vars = "Target")
head(heart_failure_clinical_records_dataset_transformed_long)

ggplot(heart_failure_clinical_records_dataset_transformed_long, aes(x = variable, y = value, color = as.factor(Target), fill = as.factor(Target)) )+
  geom_violin(alpha = .3) +
  geom_jitter() +
  theme_light() +
  ggthemes::scale_color_colorblind() +
  ggthemes::scale_fill_colorblind() +
  labs(title = "Heart failure clincal records datase", fill = "Event", color = "Event") +
  theme(legend.position.inside = TRUE, legend.position = c(.1,.1))
  
# Check Variable importance in whole data set using simple Boruta means

all_boruta <- Boruta::Boruta(Target ~ ., data = heart_failure_clinical_records_dataset_transformed, pValue = 0.0001, maxRuns = 100)
par(mfrow=c(1,1))
plot(all_boruta, las = 3)
