# Copilot Instructions for Phosphoproteomics PEAKS Online Analysis

## Project Overview
This R codebase analyzes phosphoproteomics data exported from PEAKS® Online, focusing on identifying, annotating, and visualizing phosphopeptides under different stress conditions. The workflow is tailored for biological research, with custom normalization, annotation, and visualization steps.

## Key Files & Structure
- `phosphoproteomics_PEAKS_Online.R`: Main analysis script. Handles data import, processing, annotation, and visualization.
- Data folders: `HS_phosph/`, `Arsenite_phosph/`, `H2O2_phosph/` contain experiment-specific CSVs.
- Output: Results and plots are written to dynamically created subfolders (e.g., `HS_phosph/output_0.05sign_2fold/`).

## Workflow Patterns
- **Experiment Selection**: Set variables (`dir_proteomics`, `condition`, `norm_opt`, colors) at the top of the script to switch between experiments.
- **Normalization**: Controlled by `norm_opt` variable. If set to `norm`, normalization is performed using protein-level data; otherwise, raw values are used.
- **Data Annotation**: Merges peptide and protein tables, prioritizes SwissProt entries, and annotates phosphosites and gene names.
- **Visualization**: Generates volcano plots, bar plots, and heatmaps using ggplot2, plotly, and pheatmap. Output filenames and directories are auto-generated based on parameters.

## Conventions & Patterns
- **Column Renaming**: Standardizes input columns to `control_1`, `control_2`, `control_3`, `exp_1`, `exp_2`, `exp_3` for consistent downstream processing.
- **Phosphosite Annotation**: Custom functions extract phosphosite positions and residue identities from peptide strings.
- **Group Filtering**: Uses gene/description patterns to categorize proteins (e.g., MT system groups) for targeted analysis and visualization.
- **Output**: Results are written as `.xlsx` and `.csv` files in output directories named by cutoff parameters and experiment type.

## External Dependencies
- R packages: `tidyverse`, `data.table`, `matrixTests`, `limma`, `ggrepel`, `plotly`, `RColorBrewer`, `pheatmap`, `gghighlight`, `ggvenn`, `qdap`, `Biostrings`, `writexl`, `plyr`.
- Data files: CSVs in experiment folders, SwissProt priority list (`Rat_SP_for_priority.txt`).

## Example: Switching Experiments
To analyze a different condition, update experiment_choice variable.

Uncomment and set as needed, then run the script.

## Tips for AI Agents
- Always update experiment variables before running analysis.
- Output directories and filenames are auto-generated; ensure required folders exist or are created by the script.
- Follow the column renaming and annotation patterns for new data sources.
- Use provided custom functions for phosphosite extraction and annotation.
- Use dplyr when appropriate, but reference all functions with library::function() to avoid namespace conflicts.

---
If any section is unclear or missing details, please specify which part needs further explanation or examples.
