# Hi-C Fine-Scale Features Analysis

Analysis of TADs and chromatin loops in schizophrenia using Hi-C data from neuronal and non-neuronal cells.

This repository focuses on fine-scale chromatin features (topologically associating domains and chromatin loops). The main laboratory repository contains the complete multi-scale analysis from all lab members.

## Setup

Create and activate the conda environment:

```bash
conda env create -f environment_hic.yml
conda activate hic
```

Set up your data paths in a `.env` file:

```
PATH_TO_PROCESSED_MAPS=/path/to/your/processed/maps
```

## Analysis Pipeline

The analysis is organized into four main parts:

### 1. TAD Analysis - Neurons
Detection and analysis of topologically associating domains (TADs) in neuronal cells.

- **1.1** - Call TAD borders using insulation scores
- **1.2** - Average TAD border analysis
- **1.3** - Average TAD properties
- **1.4** - Differential TAD analysis between groups
- **1.5** - TAD changes with age
- **1.6** - TAD heterogeneity across samples

### 2. TAD Analysis - Non-Neurons
Same pipeline as above but for non-neuronal cells.

### 3. Chromatin Loops - Neurons
Detection and analysis of chromatin loops in neuronal cells.

- **3.1** - Call chromatin loops with custom kernels
- **3.2** - Create loop layouts and extract features
- **3.3** - Filter loops by quality
- **3.4** - Differential loop analysis
- **3.5** - Loop changes with age
- **3.6** - Loop heterogeneity across samples

### 4. Chromatin Loops - Non-Neurons
Same pipeline as above but for non-neuronal cells.

## Additional Data

The `0.additional_data/` folder contains:
- Gene annotations and promoter information
- Custom convolution kernels for loop calling
- Expected Hi-C maps

## Running the Analysis

Run notebooks in order within each directory. Most notebooks generate output files that are used by subsequent steps.

Key output directories are created automatically:
- `tads_borders_layouts/` - TAD border positions
- `loops_cooltools_data/` - Called chromatin loops
- Visualization outputs in respective `*_visualizations/` folders

## Notes

- TAD calling uses 15kb resolution with 150kb windows
- Loop calling uses multiple custom kernels to detect loops of different sizes
- All analysis compares schizophrenia (SZ) vs healthy control (HC) samples
- Both neuronal (NeuN+) and non-neuronal (NeuN-) populations are analyzed separately
