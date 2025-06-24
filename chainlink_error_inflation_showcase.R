################################################################################
# Chainlink Error Inflation Showcase
#
# This script demonstrates synthetic data generation and error inflation analysis
# using the FCPS Chainlink dataset. Two generation algorithms are compared:
#   1. ESOM-based synthetic data generation
#   2. Gaussian Mixture Model (GMM)-based synthetic data generation
#
# Requirements:
#   - FCPS
#   - pbmcapply
#   - dplyr
#   - mclust
#   - reshape2
#   - ggplot2
#   - ggbeeswarm
#   - ggpubr
#   - ggh4x
#   - cowplot
################################################################################

# ---- 1. Setup ----------------------------------------------------------------
setwd("/home/joern/Aktuell/GenerativeESOM/08AnalyseProgramme/R/genESOMerrorSignal")

# Set reproducibility seed
seed <- 42

# Load required packages
library(FCPS)
library(pbmcapply)
library(dplyr)
library(mclust)
library(reshape2)
library(ggplot2)
library(ggbeeswarm)
library(ggpubr)
library(ggh4x)
library(cowplot)

# ---- 2. Data Preparation -----------------------------------------------------

# Load Chainlink data and class labels
chainlink_data <- FCPS::Chainlink$Data
chainlink_Cls  <- FCPS::Chainlink$Cls

# Basic statistical comparison for reference
apply(chainlink_data, 2, function(x) t.test(x ~ chainlink_Cls))

# Set radius for ESOM-based generation (from previous paper)
chainlink_radius <- 0.5

# Define synthetic generation ratios
generated_per_original <- c(1, 10, 50, 100)

# ---- 3. Synthetic Data Generation --------------------------------------------

# -- 3.1 ESOM-based synthetic data
chainlink_generated_ESOM <- pbmcapply::pbmclapply(
  generated_per_original,
  function(gen) {
    set.seed(seed)
    generate_synthetic_data(
      Data = chainlink_data,
      density_radius = chainlink_radius,
      Cls = chainlink_Cls,
      gen_per_data = gen
    )
  },
  mc.cores = min(parallel::detectCores(), length(generated_per_original))
)
names(chainlink_generated_ESOM) <- paste0("n_generated_", generated_per_original)
str(chainlink_generated_ESOM)

# -- 3.2 GMM-based synthetic data
chainlink_generated_GMM <- pbmcapply::pbmclapply(
  generated_per_original,
  function(gen) {
    set.seed(seed)
    gmm_fit <- Mclust(chainlink_data, G = 2)
    gmm_samples <- sim(gmm_fit$modelName, gmm_fit$parameters, gen * nrow(chainlink_data))
    gmm_data    <- gmm_samples[, 2:4]  # Feature columns
    gmm_classes <- gmm_samples[, 1]    # Class labels
    colnames(gmm_data) <- colnames(chainlink_data)
    data_list <- list(
      original_data = chainlink_data, 
      original_classes = chainlink_Cls,
      generated_data = gmm_data,
      generated_classes = gmm_classes
    )
  },
  mc.cores = min(parallel::detectCores(), length(generated_per_original))
)
names(chainlink_generated_GMM) <- paste0("n_generated_", generated_per_original)
str(chainlink_generated_GMM)

# ---- 4. Data Extraction and Wrangling ----------------------------------------

# Helper function to extract and format data from generated lists
extract_chainlink_df <- function(generated_list, algorithm_name) {
  result_list <- lapply(names(generated_list), function(name) {
    item <- generated_list[[name]]
    original_df <- as.data.frame(item$original_data) %>%
      mutate(
        list_name = name,
        type = "original",
        class = item$original_classes,
        GenerationAlgorithm = algorithm_name
      )
    generated_df <- as.data.frame(item$generated_data) %>%
      mutate(
        list_name = name,
        type = "generated",
        class = item$generated_classes,
        GenerationAlgorithm = algorithm_name
      )
    bind_rows(original_df, generated_df)
  })
  bind_rows(result_list)
}

# Extract and combine ESOM and GMM data
df_generated_ESOM <- extract_chainlink_df(chainlink_generated_ESOM, "genESOM")
df_generated_GMM  <- extract_chainlink_df(chainlink_generated_GMM,  "genGMM")
df_generated <- bind_rows(df_generated_ESOM, df_generated_GMM)

# Build analysis data frame
df_analysis <- df_generated %>%
  filter(list_name == "n_generated_1", type == "original") %>%
  mutate(Data = "Original")
for (n in generated_per_original) {
  name <- paste0("n_generated_", n)
  tmp <- df_generated %>%
    filter(list_name == name, type == "generated") %>%
    mutate(Data = name)
  df_analysis <- bind_rows(df_analysis, tmp)
}
df_analysis <- df_analysis %>%
  select(GenerationAlgorithm, Data, list_name, type, class, everything())

# Reshape for plotting
df_analysis_long <- reshape2::melt(
  df_analysis,
  id.vars = c("GenerationAlgorithm", "Data", "list_name", "type", "class")
)
desired_order <- c("Original", paste0("n_generated_", sort(generated_per_original)))
df_analysis_long$Data <- factor(df_analysis_long$Data, levels = desired_order)

write.csv(x = df_analysis, file = "df_analysis.csv")

# ---- 5. Visualization --------------------------------------------------------

# Define color palette for Data levels
data_levels  <- unique(df_analysis_long$Data)
data_palette <- setNames(
  colorRampPalette(c("cornsilk1", "cornsilk3"))(length(data_levels)),
  data_levels
)

# Generate plots for each algorithm
plots <- lapply(unique(df_analysis_long$GenerationAlgorithm), function(algo) {
  ggplot(
    subset(df_analysis_long, GenerationAlgorithm == algo),
    aes(x = as.factor(class), y = value, color = as.factor(class), fill = as.factor(class))
  ) +
    ggbeeswarm::geom_quasirandom(
      dodge.width = 0.8, shape = 21, alpha = 0.7, show.legend = FALSE
    ) +
    ggh4x::facet_grid2(
      rows = vars(Data),
      cols = vars(variable),
      strip = strip_themed(
        background_x = elem_list_rect(fill = "cornsilk", color = "cornsilk"),
        background_y = elem_list_rect(fill = data_palette, color = data_palette),
        text_x = elem_list_text(color = "black"),
        text_y = elem_list_text(color = "black")
      )
    ) +
    ggthemes::scale_color_colorblind() +
    ggthemes::scale_fill_colorblind() +
    ggpubr::stat_compare_means(
      aes(group = class), label = "p", label.y.npc = 0.9
    ) +
    theme(
      strip.text = element_text(face = "plain", size = 12, color = "black")
    ) +
    theme_light() +
    labs(
      title = paste("Chainlink:", algo),
      color = "Class", fill = "Class", x = "Class"
    ) +
    guides(color = "none", fill = "none")
})

# Combine and save the plots
combined_Chainlink_plot <- cowplot::plot_grid(plotlist = plots, ncol = 2, labels = "AUTO")
ggsave(filename = "combined_Chainlink_plot.png", combined_Chainlink_plot, width = 16, height = 16)


################################################################################
# End of script
################################################################################
