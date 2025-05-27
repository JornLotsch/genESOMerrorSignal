# genESOM Data Processing Pipeline

## 🚧 Project Status: Under Development 🚧

This repository contains a comprehensive data preparation, analysis, and augmentation pipeline specifically designed for biomedical research data using a neural network based generative AI (genESOM). The project is currently in active development.

## Overview

This pipeline offers a robust framework for processing, transforming, analyzing, and augmenting biomedical data. It's designed to handle common challenges in omics data including outlier detection, missing value imputation, data normalization, feature importance assessment, and synthetic data generation.

## Key Features

### Data Preprocessing
- **Distribution exploration** with automated visualization
- **Intelligent transformation selection** using:
  - Tukey's ladder of powers
  - Box-Cox transformations
- **Outlier detection and removal** using:
  - Grubbs' test
  - Boxplot method
- **Missing value imputation** using missForest algorithm
- **Back-transformation** of data to original scale
- **Parallel processing** support for improved performance

### ESOM U-matrix Generation
- **Emergent Self-Organizing Map (ESOM)** training
- **U-matrix, P-matrix, and Island visualizations**
- **2D and 3D interactive visualizations**
- **Density radius calculation** for synthetic data generation
- **Percentage-based data normalization** options
- **Toroid and non-toroid topology** support
- **Detailed visual outputs** with customizable parameters

### Feature Importance Analysis
- **Boruta algorithm** for robust feature selection
- **Multiple dataset comparison** to assess feature stability across:
  - Original data
  - Engineered data with permuted features
  - Synthetic augmented data
  - Reduced dataset option for split testing
- **Statistical significance assessment** of features using bootstrapping
- **Comprehensive visualization** of feature importance metrics
- **Circular and bar plot visualizations** for easier feature comparison

### Synthetic Data Generation
- **Density-based data augmentation** using ESOM
- **Configurable synthetic data generation** with flexible multipliers
- **Quality assessment** of generated synthetic samples

## Dependencies

The pipeline requires several R packages:

### Core Data Processing
- readr, ggplot2, ComplexHeatmap
- missForest, forecast, outliers, nortest
- parallel, devtools (for installing ABCstats)
- ABCstats (installed from GitHub)

### ESOM Visualization
- Umatrix, dbt.DataIO
- rgl (for 3D visualization)
- cowplot, dplyr

### Feature Importance Analysis
- Boruta, caret, reshape2
- cowplot, pbmcapply
- opdisDownsampling (for data splitting)

## Configuration

The pipeline is highly configurable, allowing you to control:

- Data exploration options
- Transformation methods
- Outlier detection parameters
- Missing value handling
- Parallel processing settings
- ESOM grid dimensions and visualization options
- Feature importance analysis parameters
- Synthetic data generation multipliers
- Dataset types for comparative analysis

## Usage

### ESOM U-matrix Training

The ESOM U-matrix training module allows you to:

1. Generate U-matrix visualizations from your processed data
2. Create interactive 3D visualizations for better data exploration
3. Calculate density radius for synthetic data generation

<img src="./neighborhood_distances.svg">

4. Save all necessary files for further analysis

Configuration parameters at the top of the script allow you to customize:
- ESOM grid dimensions
- Visualization options
- Output file paths and formats
- Data processing options

### Feature Importance Analysis

The feature importance analysis module allows you to:

1. Process your data and identify significant features
2. Compare feature importance across different data types:
  - Original data
  - Engineered data with permuted features
  - Synthetic augmented data
  - Reduced datasets (training/test splits)
3. Generate visualizations showing feature selection frequency
4. Create a summary table of significant features across all datasets

Configuration parameters at the top of the script allow you to customize:
- Number of iterations for feature selection
- Generation multipliers for synthetic data
- Dataset types to include in the analysis
- Output file paths and formats

Example configuration:
```r
# Define dataset types to analyze
DataSetSizes <- c("original", "engineered_0", "reduced", "augmented_1_engineered", "augmented_5_engineered")
```