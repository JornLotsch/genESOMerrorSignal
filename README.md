# Error-controlled generative AI with feature importance analysis

This is a pipeline that implements built-in error control mechanisms for self-organising neural network generative AI (genESOM), through dimensionality change detection.

## Overview

This repository provides a comprehensive framework for biomedical data augmentation with automatic error control mechanisms. Leveraging Emergent Self-Organising Maps (ESOM), the pipeline separates data structure learning from data generation. This allows for controlled dimensionality changes to be exploited to detect error inflation in synthetic data generation.

The key innovation is the injection of a 'diagnostic' signal that provides a data-based stopping point for augmentation, thereby preserving the validity of AI-augmented datasets and preventing error inflation, which is a common problem in many generative AI approaches.

The repository contains relevant R and Python code.

## Core concept

Unlike most generative models, genESOM AI separates the learning of the data structure from the data generation process itself. This unique approach allows for:

1. Generation of data with altered dimensionality from the training data
2. Implementation of stopping criteria for data augmentation

## Key components

### Data structure learning
- **ESOM U-matrix training** to capture the topological structure of data
- **Density radius calculation** to guide synthetic data generation

<img src="./neighborhood_distances.svg">

### Synthetic data generation
- **Density-based augmentation** with configurable generation rates
- **Diagnostic signal injection** for error inflation detection
- **Stopping criteria** based on feature importance stability

### Feature importance analysis
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


## Feature importance pipeline parameters

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

## Outputs and interpretation

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

This project is licensed under the GNU General Public Licence v3.0 (GPL-3.0).

## Citation

 Lötsch, J, Himmelspach A, Kringel D. (2025). Mitigating error propagation in generative AI via dimensionality modulation as a method for safe biomedical dataset augmentation. iScience [under revision]
