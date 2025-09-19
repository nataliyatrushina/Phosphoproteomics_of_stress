# Phosphoproteomics Stress Analysis

This project provides an R-based workflow for analyzing phosphoproteomics data exported from PEAKS® Online, focusing on the identification, annotation, and visualization of phosphopeptides under different stress conditions.

## Author
Nataliya Trushina  
Email: nataliya.trushin@uni-osnabrueck.de

## Project Overview
- **Main script:** `phosphoproteomics_PEAKS_Online.R`  
  Handles data import, processing, annotation, and visualization.
- **Functions:** `phospho_functions.R`  
  Contains custom functions for annotation and processing.
- **Experiment configuration:** `experiment_config.yaml`  
  Defines all experiment conditions and parameters.
- **Batch runner:** `run_all_experiments.R`  
  Automates analysis for all experiments listed in the YAML config.

## Folder Structure
- `HS_phosph/`, `Arsenite_phosph/`, `H2O2_phosph/`: Experiment-specific data folders (CSV files).
- Output folders are auto-generated per experiment and cutoff parameters.

## Getting Started
1. **Install R and required packages:**
   - Required R packages: `tidyverse`, `data.table`, `matrixTests`, `limma`, `ggrepel`, `plotly`, `RColorBrewer`, `pheatmap`, `gghighlight`, `ggvenn`, `qdap`, `Biostrings`, `writexl`, `plyr`.
   - Install with:
     ```R
     install.packages(c('tidyverse','data.table','matrixTests','limma','ggrepel','plotly','RColorBrewer','pheatmap','gghighlight','ggvenn','qdap','Biostrings','writexl','plyr'))
     ```
2. **Configure experiments:**
   - Edit `experiment_config.yaml` to add or modify experiment conditions.
3. **Run analysis:**
   - For a single experiment: Edit the experiment choice in the main script or set the environment variable `EXPERIMENT_CHOICE`.
   - For batch analysis: Run `run_all_experiments.R` to process all experiments in the YAML config.

## Usage Examples
### Single Experiment
```R
# In R console or script
Sys.setenv(EXPERIMENT_CHOICE = "HS_phosph")
source("phosphoproteomics_PEAKS_Online.R")
```

### Batch Analysis
```R
# In R console or script
source("run_all_experiments.R")
```

## Output
- Results and plots are saved in output directories named by cutoff parameters and experiment type.
- Output includes annotated datasets and plots.

## Contact
For questions or collaboration, contact Nataliya Trushina at nataliya.trushin@uni-osnabrueck.de.


## License
This project is licensed under the MIT License. See below for details.

---

MIT License

Copyright (c) 2025 Nataliya Trushina

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
