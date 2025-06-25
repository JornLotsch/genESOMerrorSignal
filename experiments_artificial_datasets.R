############## Artifical data sets ##########################################

# ---- 1. Setup ----------------------------------------------------------------
setwd("/home/joern/Aktuell/GenerativeESOM/08AnalyseProgramme/R/genESOMerrorSignal")

# Read datasets

ascending_significance_data <- read.csv("ascending_significance_test_data.csv")
no_effect_data <- read.csv("no_effect_test_data.csv")

# Data set 1: ascending_significance_data
# T-tests
t_tests_p_ascending_significance_data <- 
  apply(ascending_significance_data[,-1],2,function(x) t.test(x~ascending_significance_data$Target)$p.value)

t_tests_p_ascending_significance_data_long <- cbind.data.frame(variable = names(t_tests_p_ascending_significance_data), 
                                                               value = reshape2::melt(t_tests_p_ascending_significance_data))
t_tests_p_ascending_significance_data_long$variable <- factor(
  t_tests_p_ascending_significance_data_long$variable,
  levels = paste0("X", 1:50)
)

barplot_t_tests_p_ascending_significance_data <- 
  ggplot(data = t_tests_p_ascending_significance_data_long, aes(y = variable, x = -log10(value))) +
  geom_bar(stat = "identity", color = "cornsilk3", fill = "cornsilk2") +
  geom_vline(xintercept = -log10(0.05), color = "salmon", linetype = "dashed") +
  annotate(
    "text",
    x = -log10(0.05), 
    y = t_tests_p_ascending_significance_data_long$variable[1], # top variable
    label = "p = 0.05",
    hjust = -0.1, # adjust as needed to move label left/right
    vjust = 1.5,  # adjust as needed to move label up/down
    color = "salmon",
    fontface = "plain",
    angle = 90
  ) +
  theme_light()


print(barplot_t_tests_p_ascending_significance_data)



# Set reproducibility seed
seed <- 42


