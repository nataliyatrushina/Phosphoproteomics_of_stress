library(tidyverse)
# library(ggVennDiagram)
library(ggvenn)
library(qpcR)
library(gplots)
library(yaml)

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
regulation <- "upregulated"
# regulation <- "downregulated"

# ---------------------------------------------------------------------------- #

list_filenames <- base::list.files(path = "..", pattern = "Dataset_processesed_annotated_for_iGPS_and_KEA2_full.xlsx", recursive = T, full.names = TRUE)
list_filenames <- list_filenames[!grepl("archive", list_filenames)]

list_of_dat <- lapply(list_filenames, function(file) {
  readxl::read_excel(file) %>%
    dplyr::filter(Change == regulation)
})

dat_all <- dplyr::bind_rows(list_of_dat)

# ---------------------------------------------------------------------------- #

genes_keep <- readr::read_lines("../knownSG_rat_orthologs.txt")
dat_SG <- dplyr::filter(dat_all, Gene %in% genes_keep)

# ---------------------------------------------------------------------------- #

dat_wide <- dat_SG %>%
  dplyr::group_by(Condition) %>%
  dplyr::mutate(row_number_bygroup = 1:n()) %>% 
  tidyr::spread(Condition, Gene)

dat_wide <- subset(dat_wide, select = -c(row_number_bygroup))

# ---------------------------------------------------------------------------- #

config <- yaml::yaml.load_file("../experiment_config.yaml")
experiment_names <- names(config$experiments)

dat_wide_for_plot <- dat_wide %>%
  dplyr::select(experiment_names) %>%
  dplyr::distinct(across(dplyr::all_of(experiment_names)))

# ---------------------------------------------------------------------------- #

col_exp <- unlist(lapply(config$experiments, function(x) x$col_exp))

cbp_gradient <- grDevices::colorRampPalette(col_exp)
cbp <- cbp_gradient(nlevels(as.factor(dat_all$Condition)))

# IMPORTANT for decreased - mute the colors
cbp <- if (regulation == "downregulated") scales::muted(cbp) else cbp

# ---------------------------------------------------------------------------- #

# ggVennDiagram::ggVennDiagram(dat_wide_for_plot) +
#   ggplot2::scale_color_manual(values = cbp) +
#   ggplot2::scale_fill_continuous(type = "viridis", alpha = 0.5)

# ---------------------------------------------------------------------------- #

list <- as.list(dat_wide_for_plot)

for(i in 1:length(list)) {
  print(i)
  list[[i]] <- list[[i]][!is.na(list[[i]])]
}

venn_plot <- ggvenn::ggvenn(
  list,
  fill_color = cbp,
  stroke_size = 0.5,
  set_name_size = 4
) + ggplot2::theme(panel.background = ggplot2::element_rect(fill = "white"))
venn_plot
ggsave(venn_plot, filename = paste0("plot_venn_SG_", regulation, ".png"), width = 8, height = 8)
ggsave(venn_plot, filename = paste0("plot_venn_SG_", regulation, ".svg"), width = 8, height = 8)

# ---------------------------------------------------------------------------- #

ItemsList <- gplots::venn(list, show.plot = TRUE) 
lengths(attributes(ItemsList)$intersections)
# ItemsList

# dat_intersections <- dplyr::bind_rows(lapply(names(attributes(ItemsList)$intersections), function(n) {
#   data.frame(col_intersection = n, col_item = attributes(ItemsList)$intersections[[n]], stringsAsFactors = FALSE)
# }))
# dat_intersections <- as.data.frame(attributes(ItemsList)$intersections, stringsAsFactors = FALSE)
dat_intersections <- as.data.frame(lapply(attributes(ItemsList)$intersections, `length<-`, max(lengths(attributes(ItemsList)$intersections))), stringsAsFactors = FALSE)

writexl::write_xlsx(dat_intersections, paste0("Intersections_SG_genes_", regulation, ".xlsx"))

# ---------------------------------------------------------------------------- #
# ---------------------------------------------------------------------------- #
# ---------------------------------------------------------------------------- #

dat_wide <- dat_all %>%
  dplyr::group_by(Condition) %>%
  dplyr::mutate(row_number_bygroup = 1:n()) %>% 
  tidyr::spread(Condition, Gene)

dat_wide <- subset(dat_wide, select = -c(row_number_bygroup))

# ---------------------------------------------------------------------------- #

config <- yaml::yaml.load_file("../experiment_config.yaml")
experiment_names <- names(config$experiments)

dat_wide_for_plot <- dat_wide %>%
  dplyr::select(experiment_names) %>%
  dplyr::distinct(across(dplyr::all_of(experiment_names)))  # IMPORTANT for correct numbering

# ---------------------------------------------------------------------------- #

col_exp <- unlist(lapply(config$experiments, function(x) x$col_exp))

cbp_gradient <- grDevices::colorRampPalette(col_exp)
cbp <- cbp_gradient(nlevels(as.factor(dat_all$Condition)))

# IMPORTANT for decreased - mute the colors
cbp <- if (regulation == "downregulated") scales::muted(cbp) else cbp

# ---------------------------------------------------------------------------- #

# ggVennDiagram::ggVennDiagram(dat_wide_for_plot) +
#   ggplot2::scale_color_manual(values = cbp) +
#   ggplot2::scale_fill_continuous(type = "viridis", alpha = 0.5)

# ---------------------------------------------------------------------------- #

list <- as.list(dat_wide_for_plot)

for(i in 1:length(list)) {
  print(i)
  list[[i]] <- list[[i]][!is.na(list[[i]])]
}

venn_plot <- ggvenn::ggvenn(
  list,
  fill_color = cbp,
  stroke_size = 0.5,
  set_name_size = 4
) + ggplot2::theme(panel.background = ggplot2::element_rect(fill = "white"))
venn_plot
ggsave(venn_plot, filename = paste0("plot_venn_all_", regulation, ".png"), width = 8, height = 8)
# ggsave(venn_plot, filename = paste0("plot_venn_all_", regulation, ".svg"), width = 8, height = 8)

# ---------------------------------------------------------------------------- #

ItemsList <- gplots::venn(list, show.plot = TRUE) 
lengths(attributes(ItemsList)$intersections)
# ItemsList

# dat_intersections <- dplyr::bind_rows(lapply(names(attributes(ItemsList)$intersections), function(n) {
#   data.frame(col_intersection = n, col_item = attributes(ItemsList)$intersections[[n]], stringsAsFactors = FALSE)
# }))
# dat_intersections <- as.data.frame(attributes(ItemsList)$intersections, stringsAsFactors = FALSE)
dat_intersections <- as.data.frame(lapply(attributes(ItemsList)$intersections, `length<-`, max(lengths(attributes(ItemsList)$intersections))), stringsAsFactors = FALSE)

writexl::write_xlsx(dat_intersections, paste0("intersections_all_genes_", regulation, ".xlsx"))
