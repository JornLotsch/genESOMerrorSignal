# Error-controlled generative AI with feature importance analysis 
## (Mitigating error propagation in generative AI via dimensionality modulation as a method for safe biomedical dataset augmentation)
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

## Experiments and datasets

### Artificial datasets  
The artificial datasets are generated using the script [generate_artifical_datasets.R](https://github.com/JornLotsch/genESOMerrorSignal/blob/main/generate_artifical_datasets.R). Two key datasets include:

- **Ascending significance dataset** with systematically increasing class differences to rigorously test error inflation control.  
- **No-effect dataset** where classes show no true statistical differences, testing robustness against false positives.

### Biomedical datasets  
Evaluated with dedicated scripts tailored for each dataset:

- **Mouse lipidomics dataset** (code: `mouse_lipidomics_analysis.R`)  
  In-house preclinical lipidomics data from mice with immune modulation.

- **Heart failure clinical records dataset** (code: `heart_failure_analysis.R`)  
  Public clinical data on heart failure patients from UCI ML Repository.

- **Indian liver patient dataset** (code: `indian_liver_patient_analysis.R`)  
  Demographic and biochemical patient data from UCI ML Repository.

- **Pain thresholds sex dataset** (code: `pain_thresholds_sex_analysis.R`)  
  Quantitative pain sensitivity data assessing biological sex effects.

- **Gallstone clinical dataset** (code: `gallstone_data_analysis.R`)  
  Bioimpedance and clinical laboratory data predicting gallstone disease.

- **Golub cancer genomics dataset** (code: `golub_cancer_genomics_analysis.R`)  
  Leukemia gene expression data from Bioconductor’s golubEsets package.

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

The pipeline generates visualizations and data to detect error inflation:

1. **Feature importance stability plots** across generation multipliers  
2. **Dimensionality shift indicators** showing when synthetic data degrades signal  
3. **Bar and radial plots** depicting change in feature selection frequency  
4. **Comparative analysis** of feature importance across dataset variants  
5. **Statistical summaries** of significant variables and their stability  

## Applications

This approach is valuable for:

- Biomedical data augmentation with built-in quality control  
- Detecting the onset of error inflation in synthetic data generation  
- Establishing data-driven stopping points for augmentation  
- Preserving the validity of AI-augmented datasets  
- Comparative analysis of feature importance stability  

## License

This project is licensed under the GNU General Public Licence v3.0 (GPL-3.0).

## Citation

Lötsch, J, Himmelspach A, Kringel D. (2025). Mitigating error propagation in generative AI via dimensionality modulation as a method for safe biomedical dataset augmentation. iScience [in revision]

**Note:** This repository specifically targets the detection and control of error inflation signals during synthetic data generation. The core self-organising neural network method underpinning this approach, genESOM, was originally described in an earlier publication: Ultsch and Lötsch, 2024, *Brief Bioinformatics*, DOI: [10.1093/bib/bbae640](https://doi.org/10.1093/bib/bbae640). This work extends the original method by integrating error detection mechanisms, making it suitable for safe biomedical dataset augmentation.
