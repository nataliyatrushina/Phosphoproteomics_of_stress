# Script to loop through all YAML experiment conditions and run the main analysis script for each

library(yaml)

# Load experiment config
config <- yaml.load_file("experiment_config.yaml")
experiment_names <- names(config$experiments)

for (exp_name in experiment_names) {
  cat("\nRunning experiment:", exp_name, "\n")
  # Set environment variable for experiment selection
  Sys.setenv(EXPERIMENT_CHOICE = exp_name)
  # Run the main analysis script
  source("phosphoproteomics_PEAKS_Online.R")
}
