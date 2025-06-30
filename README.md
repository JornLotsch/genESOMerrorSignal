# Error-controlled generative AI with feature importance analysis

A pipeline implementing built-in error control mechanisms for self-organizing neural network generative AI (genESOM) through dimensionality change detection.

## Overview

This repository provides a comprehensive framework for biomedical data augmentation with automatic error control mechanisms. By leveraging Emergent Self-Organizing Maps (ESOM), the pipeline separates data structure learning from data generation, allowing for controlled dimensionality changes that can be exploited to detect error inflation in synthetic data generation.

The key innovation is the injection of a "diagnostic" signal that provides a data-based stopping point for augmentation, preserving the validity of AI-augmented datasets and preventing error inflation - a common problem in many generative AI approaches.

The repository contains relevant R and Python code.

## Core Concept

Unlike most generative models, genESOM AI separates the learning of data structure from the actual data generation process. This unique approach allows for:

1. Generation of data with altered dimensionality from the training data
2. Implementation of stopping criteria for data augmentation

## Key Components

### Data Structure Learning
- **ESOM U-matrix training** to capture the topological structure of data
- **Density radius calculation** to guide synthetic data generation

<img src="./neighborhood_distances.svg">

### Synthetic Data Generation
- **Density-based augmentation** with configurable generation rates
- **Diagnostic signal injection** for error inflation detection
- **Stopping criteria** based on feature importance stability

### Feature Importance Analysis
- **Cross-variant comparison** of feature importance across:
  - Original data
  - Engineered data with permuted features
  - Synthetic augmented data at various generation rates
- **Error inflation detection** through feature importance shifts

## Dependencies

The pipeline requires some R packages:

- Umatrix, dbt.DataIO (for ESOM training)
- Boruta, caret, reshape2 (for feature importance analysis)
- ggplot2, cowplot (for visualization)
- parallel, pbmcapply (for parallel processing)

### Optional
- missForest, forecast, outliers (for data preprocessing, if used)

## Installation

Clone this repository and install the required R packages. 

Some additional code is provided in Python. It has been tested so far in Python version 3.11.7 for Linux.


## Feature Importance Pipeline Parameters

- `output_dir`: Directory for saving results
- `output_prefix`: Prefix for output files
- `input_file`: Path to the input CSV file
- `class_name`: Name of the target class column
- `generation_multipliers`: Vector of multipliers for synthetic data generation
- `base_generation_rate`: Base rate for data generation
- `nIter`: Number of iterations for Boruta algorithm
- `seed`: Random seed for reproducibility
- `enable_plots`: Whether to generate and display plots
- `enable_file_output`: Whether to save results to files

## Outputs and Interpretation

The pipeline generates visualizations and data that help detect error inflation:

1. **Feature importance stability plots** across generation multipliers
2. **Dimensionality shift indicators** showing when generation introduces errors
3. **Bar and radial plots** showing feature selection frequency changes
4. **Comparative analysis** of feature importance across dataset variants
5. **Statistical summaries** of significant variables and their stability

## Applications

This approach is particularly valuable for:

- Biomedical data augmentation with built-in quality control
- Detecting the onset of error inflation in synthetic data generation
- Establishing data-driven stopping points for augmentation
- Preserving the validity of AI-augmented datasets
- Comparative analysis of feature importance stability

## License

This project is licensed under the GNU General Public License v3.0 (GPL-3.0)

## Citation

 Lötsch, J, & Kringel D. (2025). Generative AI with self-organizing neural networks enables implicit error inflation control via dimensionality modulation. [in preparation]
