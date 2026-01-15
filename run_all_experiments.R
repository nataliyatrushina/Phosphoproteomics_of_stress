# Script to loop through all YAML experiment conditions and run the main analysis script for each
library(yaml)
library(patchwork)

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

plots_all <- patchwork::wrap_plots(plots)
ggsave(plots_all, filename = paste("plots_volcano.png", sep = ""), width = 10, height = 8)

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
  plots[[exp_name]] <- readRDS(file.path(output_dir, paste0(norm_opt, "_volcano_plot_cut.rds")))
}

plots_all <- patchwork::wrap_plots(plots)
ggsave(plots_all, filename = paste("plots_volcano_cut.png", sep = ""), width = 10, height = 8)

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
  plots[[exp_name]] <- readRDS(file.path(output_dir, paste0(norm_opt, "_volcano_plot_cut_selected.rds")))
}

plots_all <- patchwork::wrap_plots(plots)
# ggsave(plots_all, filename = paste("plots_volcano.svg", sep = ""), width = 10, height = 8)
ggsave(plots_all, filename = paste("plots_volcano_cut_selected.png", sep = ""), width = 10, height = 8)

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
  plots[[exp_name]] <- readRDS(file.path(output_dir, paste0(norm_opt, "_phosphoprofile.rds")))
}

plots_all <- patchwork::wrap_plots(plots, ncol = 1)
ggsave(plots_all, filename = paste("plots_selected_phosphoprofile.png", sep = ""), width = 8, height = 8)

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
  plots[[exp_name]] <- readRDS(file.path(output_dir, paste0(norm_opt, "_volcano_plot_cut_G3bp1.rds")))
}

plots_all <- patchwork::wrap_plots(plots)
ggsave(plots_all, filename = paste("plots_volcano_cut_G3bp1.png", sep = ""), width = 10, height = 8)
