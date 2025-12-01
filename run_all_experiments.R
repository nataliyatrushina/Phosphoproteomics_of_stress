# Script to loop through all YAML experiment conditions and run the main analysis script for each

library(yaml)

# Load experiment config
config <- yaml.load_file("experiment_config.yaml")
experiment_names <- names(config$experiments)
RUN_ALL <- TRUE

for (exp_name in experiment_names) {
  cat("\nRunning experiment:", exp_name, "\n")
  # Set environment variable for experiment selection
  # Sys.setenv(EXPERIMENT_CHOICE = exp_name)
  # Run the main analysis script
  source("phosphoproteomics_PEAKS_Online.R")
}

# ---------------------------------------------------------------------------- #

plots <- list()

for (exp_name in experiment_names) {
  output_dir <- paste(
    paste0(exp_name, "_phosph"),
    "/output_",
    formatC(pValue_cutoff, digits = 2, format = "f"),
    "sign_",
    FC_cutoff,
    "fold",
    sep = ""
  )
  plots[[exp_name]] <- readRDS(file.path(output_dir, paste0(norm_opt, "_volcano_plot.rds")))
}

patchwork::wrap_plots(plots)
