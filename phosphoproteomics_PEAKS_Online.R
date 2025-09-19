# title: "Phosphoproteomics script for PEAKS® Online output analysis"
# author: "Nataliya Trushina"
# date: "2025-03-22"

# ---------------------------------------------------------------------------- #
# --- Script description ----------------------------------------------------- #
# ---------------------------------------------------------------------------- #

# - Identify phosphopeptides in exported files from PEAKS.
# - Generate a detailed annotated table for these phosphopeptides.
# - Compare the presence of phosphopeptides in experimental and control samples.
# - Produce a volcano plot based on predefined thresholds.
# Required Inputs:
# - Exported PEAKS files: 'peptides' and 'proteins' tables, and "peptide-protein" for the start position.

# ---------------------------------------------------------------------------- #
# --- Libraries -------------------------------------------------------------- #
# ---------------------------------------------------------------------------- #

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

# if (!requireNamespace("BiocManager", quietly = TRUE))
# install.packages("BiocManager")
# BiocManager::install("mixOmics")
# BiocManager::install("HTSFilter")
# BiocManager::install("MSstats")
# BiocManager::install("Biostrings")

# Package incompatibility with the newest ggrepel results in the wrong plotly outputs and more.
# detach("package:ggrepel", unload = TRUE)
# remove.packages("ggrepel")
# devtools::install_version("ggrepel", version = "0.8.1", repos = "http://cran.us.r-project.org")
# library(ggrepel)

# Data manipulation and general utilities
library(tidyverse)  # Includes ggplot2, dplyr, tidyr, readr, purrr, tibble, stringr, forcats
library(data.table)

# Statistical analysis
library(matrixTests)  # row_t_welch() and row_t_equalvar() for row-wise t-tests on matrices
library(limma)  # BiocManager::install("limma")

# Visualization
library(ggrepel)
library(plotly)
library(RColorBrewer)
library(pheatmap)
library(gghighlight)    # Enhances ggplot2 plots by highlighting
library(ggvenn)         # Venn diagrams with ggplot2

# Text analysis
library(qdap)           # remotes::install_version("qdap", version = "2.4.3")

# Bioinformatics
# BiocManager::install("Biostrings")  # see BiocManager installation above
library(Biostrings)

# ---------------------------------------------------------------------------- #

# Source external phosphoproteomics functions
source("phospho_functions.R")

# ---------------------------------------------------------------------------- #
# --- Settings and thresholds for modification ------------------------------- #
# ---------------------------------------------------------------------------- #

#Cut-offs for statistics and biologically relevant difference
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# Load experiment parameters from YAML
if (!requireNamespace("yaml", quietly = TRUE)) install.packages("yaml")
library(yaml)

# Choose experiment: "HS", "H2O2", "Arsenite", "DEX"
experiment_choice <- "Arsenite" # Change this to select experiment
config <- yaml.load_file("experiment_config.yaml")
params <- config$experiments[[experiment_choice]]

dir_proteomics <- params$dir_proteomics
condition <- params$condition
norm_opt <- params$norm_opt
col_control <- params$col_control
control_name <- params$control_name
col_exp <- params$col_exp
# dir.create(condition)

pValue_cutoff <- 0.05
Sign_cutoff <- -10*log10(pValue_cutoff)
FC_cutoff <- 2

cutoff_present_sign <- Inf
cutoff_present_FC <- Inf

cols <- c("Decreased" = col_control, "NS" = "#C2C2C2", "Increased" = col_exp)
output_dir <- paste(condition, "/output_", formatC(pValue_cutoff, digits = 2, format = "f"),"sign_",as.character(FC_cutoff),"fold", sep = "")
dir.create(output_dir)
plot_title <- paste("pValue < ", as.character(pValue_cutoff),"\nFold Change > ", as.character(FC_cutoff),"\nDifference: ", condition, "/", control_name, sep = "")

# ---------------------------------------------------------------------------- #
# --- Dataset processing ----------------------------------------------------- #
# ---------------------------------------------------------------------------- #

dat_orig <- read.csv(paste(condition, "/lfq.peptides.csv", sep = ""), dec = ".")
colnames(dat_orig)

# Rename the columns to standardize the processing
print(colnames(dat_orig[6:11]))
old_names1 <- colnames(dat_orig[6:11])
new_names1 <- c("control_1", "control_2", "control_3", "exp_1", "exp_2", "exp_3")
setnames(dat_orig, old = old_names1, new = new_names1)
print(colnames(dat_orig[13:14]))
old_names2 <- colnames(dat_orig[13:14])
new_names2 <- c("group_control", "group_exp")
setnames(dat_orig, old = old_names2, new = new_names2)

dat <- dat_orig %>%
  dplyr::select(Peptide, PTM, Significance, control_1, control_2, control_3, exp_1, exp_2, exp_3, group_control, group_exp, Accession) %>%
  
  # Filter phosphorylated peptides
  dplyr::filter(grepl("Phosphorylation", PTM)) %>%
  
  dplyr::mutate(
    Phosphopeptide = Peptide,
    
    # Remove other modification tags
    Phosphopeptide = gsub("\\(\\+57\\.02\\)", "", Phosphopeptide),
    Phosphopeptide = gsub("\\(\\+15\\.99\\)", "", Phosphopeptide),
    Phosphopeptide = gsub("\\(\\+42\\.01\\)", "", Phosphopeptide),
    
    FC = group_exp / group_control,
    logFC = log2(FC),
    Significance = as.numeric(Significance),
    pvalue = 10^(-Significance / 10),
    
    Peptide_sequence = qdap::bracketX(Peptide, fix.space = FALSE)  # delete everything between brackets, check #qdap::bracketX("TSSRES(+79.97)LNQVDLDC(+57.02)AVATFDPPS(+79.97)DMESEAEDAPWSSDSLSR", fix.space = FALSE) 
  ) %>% 
  mutate(Accession = strsplit(as.character(Accession), ":")) %>% 
  unnest(Accession) %>%
  
  dplyr::mutate(Uniprot_accession = sub("\\|.*", "", Accession))

# Load the proteins table
prot <- read.csv(paste(condition, "/lfq.proteins.csv", sep = ""), dec = ".") %>%
  select(Accession, Description)

# Sometimes gene names contain "_" this could require deleting everything after their first occurrence or
# taking a delimiter for phosphosites names that never occurs in gene names
# or what is done here: they are replaced with "-" because for their identification the description has to be read anyway
#dat_proc$Gene[grepl("_", dat_proc$Gene)] # check different IDs

# Merge the tables to get full descriptions and gene names
dat_proc <- dat %>%
  left_join(prot, by = "Accession") %>%
  dplyr::mutate(
    Gene = sub(" .*", "", sub(".*GN=", "", Description)),
    Gene = gsub("_", "-", Gene),
    Full_description = Description,
    Description = sub("OS.*", "", Description)  # remove everything after OS= (first field)
  )

# ---------------------------------------------------------------------------- #

# IMPORTANT Priority is given to reviewed entries from SwissProt
reviewed <- read.csv("Rat_SP_for_priority.txt")

# ---------------------------------------------------------------------------- #

start_table <- readr::read_csv(file.path(condition, "lfq.protein-peptides.csv")) %>%  # read.csv(paste(condition, "/lfq.protein-peptides.csv", sep = ""), sep = ',', dec = ".", header = TRUE)
  dplyr::filter(grepl("Phosphorylation", PTM)) %>%
  dplyr::select(Peptide, Start) %>%
  dplyr::mutate(Peptide = sub("\\.[A-Z]", "", sub("[A-Z]\\.", "", Peptide))) %>%
  # deduplicate, some peptides appear twice because after cut a different amino acid could be predicted which were cut by previous step
  dplyr::distinct(Peptide, .keep_all = TRUE)

# ---------------------------------------------------------------------------- #

dat_proc_out <- dat_proc %>%
  dplyr::mutate(
    Database = ifelse(Uniprot_accession %in% reviewed$SwissProt, "SwissProt", "TrEMBL"),
    Database = factor(Database, levels = c("SwissProt", "TrEMBL"))
  ) %>%
  dplyr::arrange(Database) %>%
  ### VERY IMPORTANT this filtering was at first not implemented here for the REDOX paper
  # dplyr::distinct(Peptide, Significance, logFC, .keep_all = TRUE) %>%  # to remove duplicates based on the combination of the Peptide, Significance, and logFC columns. This means that only the first occurrence of each unique combination of these three values will be retained in the dataset. All other rows that have the same combination of these three values will be considered duplicates and removed. The .keep_all = TRUE argument ensures that while it checks for duplicates based on those three columns, it retains all other columns in the first instance of each unique combination.
  left_join(start_table, by = "Peptide") %>%
  # Remove peptides where the start could not be identified - very few peptides
  filter(!is.na(Start))

# ---------------------------------------------------------------------------- #

dat_proc_out$Phosphopeptide_matrix <- as.matrix(dat_proc_out[,"Phosphopeptide"])
dat_proc_out$Phosphomarks <- apply(dat_proc_out$Phosphopeptide_matrix, 1, get_list_phosphomarks)

dat_proc_out$Phosphosite_positions <- mapply(get_list_phosphosites, dat_proc_out$Phosphopeptide, dat_proc_out$Start)
# dat_proc_out$Phosphoresidues <- mapply(get_list_residues, dat_proc_out$Protein_sequence, dat_proc_out$Phosphosite_positions)

dat_proc_out$Phosphoresidues <- mapply(get_list_residues,dat_proc_out$Peptide_sequence, dat_proc_out$Phosphomarks) 
dat_proc_out$Phosphoresidues <- strsplit(dat_proc_out$Phosphoresidues, "")
dat_proc_out$Phosphoresidues <- vapply(dat_proc_out$Phosphoresidues, paste, collapse = ", ", character(1L))

dat_proc_out$Phosphosite <- mapply(paste_phosphosites,dat_proc_out$Phosphoresidues, dat_proc_out$Phosphosite_positions) 
dat_proc_out$Phosphosite <- sapply(dat_proc_out$Phosphosite, unlist)
dat_proc_out$Phosphosite<- vapply(dat_proc_out$Phosphosite, paste, collapse = ", ", character(1L))

dat_proc_out <- dat_proc_out %>%
  dplyr::mutate(
    PhosPep = Phosphopeptide,
    PhosPep = gsub("S\\(\\+79.97\\)", "pS", PhosPep),
    PhosPep = gsub("T\\(\\+79.97\\)", "pT", PhosPep),
    PhosPep = gsub("Y\\(\\+79.97\\)", "pY", PhosPep)
  )

# ---------------------------------------------------------------------------- #
# --- Normalization ---------------------------------------------------------- #
# ---------------------------------------------------------------------------- #

# Only run if the correct version is set at the beginning of the script
if (norm_opt == "norm") {
  
  norm_prot <- read.csv(paste(dir_proteomics, "/lfq.proteins.csv", sep = ""), dec = ".")  # tried 0 peptides for LFQ, did not change the numbers
  head(norm_prot)
  colnames(norm_prot)
  # Normalizing is done by the mean area
  group_control_A <- colnames(norm_prot[16])  # "*.Area"
  group_exp_A <- colnames(norm_prot[17])  # "*.Area"
  colnames(norm_prot)[which(names(norm_prot) == group_control_A)] <- "prot_control"
  colnames(norm_prot)[which(names(norm_prot) == group_exp_A)] <- "prot_exp"
  colnames(norm_prot)
  norm_prot <- norm_prot[c("Accession", "prot_control", "prot_exp")]
  
  norm_prot$Uniprot_accession <- sub("\\|.*", "", norm_prot$Accession)
  dat_proc_out_norm <- merge(x = dat_proc_out, y = norm_prot, by = c("Accession","Uniprot_accession"), all.x = TRUE, sort = TRUE)
  nrow(dat_proc_out_norm[is.na(dat_proc_out_norm$Accession),])
  dat_proc_out_norm <- dat_proc_out_norm[!is.na(dat_proc_out_norm$prot_control),]
  
  # ---------------------------------------------------------------------------- #
  # ---------------------------------------------------------------------------- #
  # ---------------------------------------------------------------------------- #
  
  print("NORMALIZING")
  x_fin <- dat_proc_out_norm
  unique_control <- filter(x_fin, (control_1 != 0 & control_2 != 0 & control_3 != 0) & (exp_1 == 0 & exp_2 == 0 & exp_3 == 0))
  #unique_control$norm_group_control <- unique_control$prot_exp*unique_control$group_control/unique_control$prot_control
  if (nrow(unique_control) > 0) unique_control$Detected <- "Unique in control"
  unique_decreased_genes <- as.data.frame(na.omit(unique(unique_control$Gene)))
  # write.table(unique_decreased_genes, paste(output_dir, "/", norm_opt, "_unique_in_control_genes.txt", sep = ""), sep = ';', row.names = FALSE, col.names = FALSE, quote = FALSE)
  unique_exp <- filter(x_fin, (control_1 == 0 & control_2 == 0 & control_3 == 0) & (exp_1 != 0 & exp_2 != 0 & exp_3 != 0))
  #unique_exp$norm_group_exp <- unique_exp$prot_control*unique_exp$group_exp/unique_exp$prot_exp
  if (nrow(unique_exp) > 0) unique_exp$Detected <- "Unique in experiment"
  unique_increased_genes <- as.data.frame(na.omit(unique(unique_exp$Gene)))
  # write.table(unique_increased_genes, paste(output_dir, "/", norm_opt, "_unique_in_exp_genes.txt", sep = ""), sep = ';', row.names = FALSE, col.names = FALSE, quote = FALSE)
  detected_in_all <- filter(x_fin, control_1 != 0 & control_2 != 0 & control_3 != 0 & exp_1 != 0 & exp_2 != 0 & exp_3 != 0) # & exp_2 != 0 & exp_3 != 0)
  detected_in_all$Detected <- "In all samples"
  
  prot_norm_factor <- detected_in_all$prot_control / detected_in_all$prot_exp
  
  # detected_in_all$FC <- (detected_in_all$group_exp / 
  #                          detected_in_all$group_control) / prot_norm_factor
  # detected_in_all$logFC <- log2(detected_in_all$FC)
  
  detected_in_all <- detected_in_all %>%
    dplyr::mutate(
      norm_control_1 = control_1 / prot_norm_factor,
      norm_control_2 = control_2 / prot_norm_factor,
      norm_control_3 = control_3 / prot_norm_factor,
      norm_exp_1 = exp_1 / prot_norm_factor,
      norm_exp_2 = exp_2 / prot_norm_factor,
      norm_exp_3 = exp_3 / prot_norm_factor,
      norm_group_control = (norm_control_1 + norm_control_2 + norm_control_3) / 3,
      norm_group_exp = (norm_exp_1 + norm_exp_2 + norm_exp_3) / 3,
      FC = norm_group_exp / norm_group_control,
      logFC = log2(FC)
    )
  
  sum_control <- data.matrix(select(detected_in_all, norm_control_1, norm_control_2, norm_control_3))
  sum_exp <- data.matrix(select(detected_in_all, norm_exp_1, norm_exp_2, norm_exp_3))
  sum_control_log <- log2(sum_control)
  sum_exp_log <- log2(sum_exp)
  
  # ---------------------------------------------------------------------------- #
  
  BT_sum <- bartlett.test(list(sum_control_log, sum_exp_log))
  if (BT_sum$p.value < 0.05) {
    print("Unequal variances")
    stat_phosph <- matrixTests::row_t_welch(sum_control_log, sum_exp_log, alternative = "two.sided", conf.level = 0.95)
  } else {
    print("Equal variances")
    stat_phosph <- matrixTests::row_t_equalvar(sum_control_log, sum_exp_log, alternative = "two.sided", mu = 0, conf.level = 0.95)
  }
  
  # ---------------------------------------------------------------------------- #
  
  # exprs <- cbind(sum_control_log, sum_exp_log)
  # group <- factor(c(rep("control", ncol(sum_control_log)), rep("experiment", ncol(sum_exp_log))))
  # design <- model.matrix(~group)
  # 
  # fit <- limma::lmFit(exprs, design)
  # fit <- limma::eBayes(fit)
  # stat_phosph <- limma::topTable(fit, coef = 2, number = Inf, sort.by = "none")  # P.Value and adj.P.Val 
  
  # ---------------------------------------------------------------------------- #
  
  detected_in_all$new_pvalue <- stat_phosph$pvalue
  # detected_in_all$new_pvalue <- stat_phosph$P.Value
  detected_in_all$adjusted_new_pvalue <- matrix(p.adjust(as.vector(as.matrix(detected_in_all$new_pvalue)), method = 'fdr'))
  # detected_in_all$adjusted_new_pvalue <- stat_phosph$adj.P.Val
  detected_in_all$Significance_new <- -10*log10(detected_in_all$new_pvalue)
  detected_in_all$Significance_new_adjusted <- -10*log10(detected_in_all$adjusted_new_pvalue)
  x_fin_3or3 <- plyr::rbind.fill(unique_control, detected_in_all, unique_exp)
  dat_final <- x_fin_3or3[c("Accession", "PhosPep", "Phosphopeptide", "Gene", "Description", "Phosphosite","FC", "logFC", "Significance_new", "Detected")]
  
} else {
  
  print("NOT NORMALIZING")
  x_fin <- dat_proc_out
  unique_control <- filter(x_fin, (control_1 != 0 & control_2 != 0 & control_3 != 0) & (exp_1 == 0 & exp_2 == 0 & exp_3 == 0))
  # unique_control$norm_group_control <- unique_control$prot_exp*unique_control$group_control/unique_control$prot_control
  unique_control$Detected <- "Unique in control"
  unique_decreased_genes <- as.data.frame(na.omit(unique(unique_control$Gene)))
  write.table(unique_decreased_genes, paste(output_dir, "/", norm_opt, "_unique_in_control_genes.txt", sep = ""), sep = ';', row.names = FALSE, col.names = FALSE, quote = FALSE)
  unique_exp <- filter(x_fin, (control_1 == 0 & control_2 == 0 & control_3 == 0) & (exp_1 != 0 & exp_2 != 0 & exp_3 != 0))
  # unique_exp$norm_group_exp <- unique_exp$prot_control*unique_exp$group_exp/unique_exp$prot_exp
  
  unique_exp$Detected <- "Unique in experiment"
  unique_increased_genes <- as.data.frame(na.omit(unique(unique_exp$Gene)))
  write.table(unique_increased_genes, paste(output_dir, "/", norm_opt, "_unique_in_exp_genes.txt", sep = ""), sep = ';', row.names = FALSE, col.names = FALSE, quote = FALSE)
  detected_in_all <- filter(x_fin, control_1 != 0 & control_2 != 0 & control_3 != 0 & exp_1 != 0 & exp_2 != 0 & exp_3 != 0) # & exp_2 != 0 & exp_3 != 0)
  
  detected_in_all$Detected <- "In all samples"
  
  x_fin_3or3 <- plyr::rbind.fill(unique_control, detected_in_all, unique_exp)
  dat_final <- x_fin_3or3[c("Accession", "PhosPep", "Phosphopeptide", "Gene", "Description", "Phosphosite","FC", "logFC", "Significance", "Detected")]
}

# ---------------------------------------------------------------------------- #
# ---------------------------------------------------------------------------- #
# ---------------------------------------------------------------------------- #

dat_final <- dat_final %>%
  dplyr::mutate(
    Significance_new = if ("Significance_new" %in% names(.)) Significance_new else Significance
    ) %>%
  dplyr::mutate(
    Change = dplyr::case_when(
      Significance_new <= Sign_cutoff ~ "NS",
      FC >= FC_cutoff | Detected == "Unique in experiment" ~ "Increased",
      FC <= 1 / FC_cutoff | Detected == "Unique in control" ~ "Decreased",
      TRUE ~ "NS"
    ),
    Change = as.factor(Change)
  ) %>% 
  mutate(Phosphosite = strsplit(as.character(Phosphosite), ", ")) %>% 
  unnest(Phosphosite) %>%
  tidyr::unite("Gene_phosphosite", c("Gene", "Phosphosite"), remove = FALSE)

dat_final %>%
  writexl::write_xlsx(file.path(output_dir, "Dataset_processesed_annotated_for_iGPS_and_KEA2_full.xlsx"))

# write.table(apply(dat_final_out[dat_final_out$Change == "Increased",], 2, as.character), paste(output_dir, "/", norm_opt, "_out_for_iGPS_and_KEA2_increased.csv", sep = ""), row.names = FALSE)
# write.table(apply(dat_final_out[dat_final_out$Change == "Decreased",], 2, as.character), paste(output_dir, "/", norm_opt, "_out_for_iGPS_and_KEA2_decreased.csv", sep = ""), row.names = FALSE)
# write.table(apply(dat_final_out,2,as.character), paste0(output_dir, "/", norm_opt, "_out_for_iGPS_and_KEA2_all.csv"), row.names = FALSE)

# ---------------------------------------------------------------------------- #
# --- Visualization ---------------------------------------------------------- #
# ---------------------------------------------------------------------------- #

# TODO fix to have better overview of cut-offs and not accidentally cut the plot on total
dat_for_plot_cut <- dat_final %>%
  dplyr::mutate(
    Significance_new = dplyr::case_when(
      Significance_new >= cutoff_present_sign ~ Inf,
      TRUE ~ Significance_new
    ),
    logFC = dplyr::case_when(
      logFC >= cutoff_present_FC ~ Inf,
      logFC <= -cutoff_present_FC ~ -Inf,
      TRUE ~ logFC
    )
  )

# ---------------------------------------------------------------------------- #

p <- ggplot(dat_final, aes(x = logFC, y = Significance_new, colour = Change,
                                    text = paste("Gene:", Gene, "\n", "Description:", Description))) +
  geom_point(size = 0.7, alpha = 0.4) +
  scale_colour_manual(values = cols) +
  geom_hline(yintercept = Sign_cutoff, colour = "gray", linetype = "dashed") +
  geom_vline(xintercept = log2(FC_cutoff), colour = "gray", linetype = "dashed") +
  geom_vline(xintercept = -log2(FC_cutoff), colour = "gray", linetype = "dashed") +
  scale_x_continuous(name = "Fold Change [log2]", breaks = seq(-10, 10, 1),
                     limits = if (is.finite(cutoff_present_FC)) c(-cutoff_present_FC, cutoff_present_FC) else NULL) +
  scale_y_continuous(name = "P Value [-10lg]", breaks = seq(0, 1000, 10)) +
  theme_bw(base_size = 10) +
  theme(
    plot.title = element_text(color = "black", size = 10),
    axis.text.x = element_text(color = "black"),
    axis.text.y = element_text(color = "black"),
    panel.grid.minor = element_blank(),
    legend.position = "top"
  )

p
ggsave(p, filename = paste(output_dir, "/", norm_opt, "_volcano_plot.png", sep = ""), width = 6, height = 4)
ggsave(p, filename = paste(output_dir, "/", norm_opt, "_volcano_plot.svg", sep = ""), width = 6, height = 4)

# ---------------------------------------------------------------------------- #

# MT-associated genes
color_MT_struct <- "#FFBD00"
color_nucl <- "#008F1D"
color_MAP <- "#390099"
color_EB <- "#FF5400"
color_stmn <- "#ff0000"
color_sever <- "#FF0054"
color_transp <- "#9E0059"

patterns_MT_struct <- c("Tubb", "Tuba")
patterns_nucl <- c("Tubg")
patterns_MAP <- c("Map1a","Map1b","Map1s","Map2\\b","Mapt","Map4\\b","Map6","Map7") #,"Pafah1b1"
patterns_EB <- c("Mapre", "Clasp")
patterns_stmn <- c("Stmn","rb3")  # stmn2 can also be scg10, stmn3 can be sclip, stmn4 - rb3
patterns_sever <- c("katanin", "fidgetin", "spastin")
patterns_transp <- c("kinesin", "dynein")

# Check where to look for patterns for
dat_MT_struct <- dplyr::filter(dat_for_plot_cut, grepl(paste(patterns_MT_struct, collapse = "|"),ignore.case = TRUE, Gene))
dat_nucl <- dplyr::filter(dat_for_plot_cut, grepl(paste(patterns_nucl, collapse = "|"),ignore.case = TRUE, Gene))
dat_MAP <- dplyr::filter(dat_for_plot_cut, grepl(paste(patterns_MAP, collapse = "|"),ignore.case = TRUE, Gene))
dat_EB <- dplyr::filter(dat_for_plot_cut, grepl(paste(patterns_EB, collapse = "|"),ignore.case = TRUE, Gene))
dat_stmn  <- dplyr::filter(dat_for_plot_cut, grepl(paste(patterns_stmn , collapse = "|"),ignore.case = TRUE, Gene))
dat_sever <- dplyr::filter(dat_for_plot_cut, grepl(paste(patterns_sever, collapse = "|"),ignore.case = TRUE, Description))
dat_transp <- dplyr::filter(dat_for_plot_cut, grepl(paste(patterns_transp, collapse = "|"),ignore.case = TRUE, Description))
dat_MT <- rbind(dat_MT_struct, dat_nucl, dat_MAP, dat_EB, dat_stmn, dat_sever, dat_transp)

# ---------------------------------------------------------------------------- #

dat_MT <- dat_final %>%
  mutate(Category = case_when(
    grepl(paste(patterns_MT_struct, collapse = "|"), Gene, ignore.case = TRUE) ~ "Structure proteins",
    grepl(paste(patterns_nucl, collapse = "|"), Gene, ignore.case = TRUE) ~ "Nucleators",
    grepl(paste(patterns_MAP, collapse = "|"), Gene, ignore.case = TRUE) ~ "MT-binding proteins",
    grepl(paste(patterns_EB, collapse = "|"), Gene, ignore.case = TRUE) ~ "End-binding proteins",
    grepl(paste(patterns_stmn, collapse = "|"), Gene, ignore.case = TRUE) ~ "Tubulin-sequestering\nproteins",
    grepl(paste(patterns_sever, collapse = "|"), Description, ignore.case = TRUE) ~ "MT-severing proteins",
    grepl(paste(patterns_transp, collapse = "|"), Description, ignore.case = TRUE) ~ "MT-transport proteins",
    TRUE ~ NA_character_
  )) %>%
  filter(!grepl("at$", Gene, ignore.case = TRUE) & !grepl("at$", Description, ignore.case = TRUE)) %>%
  filter(!is.na(Category)) %>%
  select(Category, Gene_phosphosite, logFC, Significance_new, Change) %>%
  dplyr::arrange(Gene_phosphosite, Category) %>%
  dplyr::rename(
    "Gene Phosphosite" = Gene_phosphosite,
    "Fold Change [log₂]" = logFC,
    "P Value [-10lg]" = Significance_new,
    "Change direction" = Change
  )

options <- c("Increased", "Decreased")
for (option in options) {
  dat_MT %>%
    dplyr::filter(`Change direction` == option) %>%
    writexl::write_xlsx(paste0(output_dir, "/", condition, "_MTsystem", "_", option, ".xlsx"))
}

# ---------------------------------------------------------------------------- #

.labs_sele_MT_struct <- dat_MT_struct$Gene_phosphosite[dat_MT_struct$Change != "NS"]
.labs_sele_nucl <- dat_nucl$Gene_phosphosite[dat_nucl$Change != "NS"] 
.labs_sele_MAP <- dat_MAP$Gene_phosphosite[dat_MAP$Change != "NS"] 
.labs_sele_EB <- dat_EB$Gene_phosphosite[dat_EB$Change != "NS"] 
.labs_sele_stmn <- dat_stmn$Gene_phosphosite[dat_stmn$Change != "NS"] 
.labs_sele_sever <- dat_sever$Gene_phosphosite [dat_sever$Change != "NS"] 
.labs_sele_transp <- dat_transp$Gene_phosphosite[dat_transp$Change != "NS"] 

print("1) All found in dataset, 2) DEPs")
print("Structure proteins")
dat_MT_struct$Gene
.labs_sele_MT_struct
print("Nucleators")
dat_nucl$Gene 
.labs_sele_nucl
print("MT-binding protein")
dat_MAP$Gene
.labs_sele_MAP
print("End-binding proteins")
dat_EB$Gene
.labs_sele_EB
print("Tubulin-sequestering proteins")
dat_stmn$Gene
.labs_sele_stmn
print("MT-severing proteins")
dat_sever$Gene
.labs_sele_sever
print("MT transport proteins")
dat_transp$Gene
.labs_sele_transp

# ---------------------------------------------------------------------------- #

df_nonMT <- dat_for_plot_cut

df_nonMT <- df_nonMT[!df_nonMT$Gene %in% dat_MT_struct$Gene,]
df_nonMT <- df_nonMT[!df_nonMT$Gene %in% dat_nucl$Gene,]
df_nonMT <- df_nonMT[!df_nonMT$Gene %in% dat_MAP$Gene,]
df_nonMT <- df_nonMT[!df_nonMT$Gene %in% dat_EB$Gene,]
df_nonMT <- df_nonMT[!df_nonMT$Gene %in% dat_stmn$Gene,]
df_nonMT <- df_nonMT[!df_nonMT$Gene %in% dat_sever$Gene,]
df_nonMT <- df_nonMT[!df_nonMT$Gene %in% dat_transp$Gene,]

print("Tau hits")
dat_for_plot_cut %>% dplyr::filter(, grepl("Mapt", ignore.case = TRUE, Gene))

p_MT_groups <- ggplot(df_nonMT, aes(x = logFC, y = Significance_new, colour = Change,
                                    text = paste("Gene:", Gene, "\n", "Description:", Description))) +
  geom_point(size = 0.7, alpha = 0.4) +
  scale_colour_manual(values = cols) +
  geom_hline(yintercept = Sign_cutoff, colour = "gray", linetype = "dashed") + 
  geom_vline(xintercept = log2(FC_cutoff), colour = "gray", linetype = "dashed") + 
  geom_vline(xintercept = -log2(FC_cutoff), colour = "gray", linetype = "dashed") +
  geom_point(data = dat_MT_struct, aes(fill = color_MT_struct), size = 1.5, shape = 21, alpha = 0.8) +
  ggrepel::geom_text_repel(data = dat_MT_struct[dat_MT_struct$Change != "NS",],
                           direction = "y", #nudge_x = -2, force = 20,
                           aes(label = .labs_sele_MT_struct),
                           color = color_MT_struct,
                           size = 2, max.overlaps = Inf) +
  geom_point(data = dat_nucl, aes(fill = color_nucl), color = "white", size = 1.5, shape = 21, alpha = 0.8) +
  geom_text_repel(data = dat_nucl[dat_nucl$Change != "NS",],
                  #direction = "x", nudge_x = -2, force = 20,
                  aes(label = .labs_sele_nucl),
                  color = color_nucl,
                  size = 2, max.overlaps = Inf) +
  geom_point(data = dat_MAP, aes(fill = color_MAP), color = "white", size = 1.5, shape = 21, alpha = 0.8) +
  geom_text_repel(data = dat_MAP[dat_MAP$Change != "NS",],
                  #direction = "x", #nudge_x = -2, force = 20,
                  aes(label = .labs_sele_MAP),
                  color = color_MAP,
                  size = 2, max.overlaps = Inf) +
  geom_point(data = dat_EB, aes(fill = color_EB), color = "white", size = 1.5, shape = 21, alpha = 0.8) +
  geom_text_repel(data = dat_EB[dat_EB$Change != "NS",],
                  #direction = "x", nudge_x = -2, force = 20,
                  aes(label = .labs_sele_EB),
                  color = color_EB,
                  size = 2, max.overlaps = Inf) +
  geom_point(data = dat_stmn, aes(fill = color_stmn), color = "white", size = 1.5, shape = 21, alpha = 0.8) +
  geom_text_repel(data = dat_stmn[dat_stmn$Change != "NS",],
                  direction = "x", #nudge_x = -2, force = 20,
                  aes(label = .labs_sele_stmn),
                  color = color_stmn,
                  size = 2, max.overlaps = Inf) +
  geom_point(data = dat_sever, aes(fill = color_sever), color = "white", size = 1.5, shape = 21, alpha = 0.8) +
  geom_text_repel(data = dat_sever[dat_sever$Change != "NS",],
                  direction = "x", #nudge_x = -2, force = 20,
                  aes(label = .labs_sele_sever),
                  color = color_sever,
                  size = 2, max.overlaps = Inf) +
  scale_fill_identity(name = "MT skeleton components",
                      breaks = c(color_MT_struct, color_nucl, color_MAP, color_EB, color_stmn, color_sever
                      ),
                      labels = c("Structure proteins", "Nucleators", "MT-binding proteins",
                                 "End-binding proteins", "Tubulin-sequestering\nproteins",
                                 "MT-severing proteins"
                      ),
                      guide = guide_legend(nrow = 3)) +
  scale_x_continuous(name = "Fold Change [log2]", breaks = seq(-10, 10, 1),
                     limits = if (is.finite(cutoff_present_FC)) c(-cutoff_present_FC, cutoff_present_FC) else NULL) +
  scale_y_continuous(name = "P Value [-10lg]", breaks = seq(0, 1000, 10)) +
  theme_bw(base_size = 10) +
  theme(
    plot.title = element_text(color = "black", size = 10),
    axis.text.x = element_text(color = "black"),
    axis.text.y = element_text(color = "black"),
    panel.grid.minor = element_blank(),
    legend.position = "top"
  ) +
  guides(colour = guide_legend(nrow = 3))

p_MT_groups
ggsave(p_MT_groups, filename = paste(output_dir, "/", norm_opt, "_volcano_plot_MTsystem.png", sep = ""), width = 7, height = 6)
ggsave(p_MT_groups, filename = paste(output_dir, "/", norm_opt, "_volcano_plot_MTsystem.svg", sep = ""), width = 7, height = 6)

# ---------------------------------------------------------------------------- #
# --- Bar plots for MT system groups ----------------------------------------- #
# ---------------------------------------------------------------------------- #

dat_MT_struct$Group <- "Structure proteins"
dat_nucl$Group <-  "Nucleators"
dat_MAP$Group <- "MT-binding proteins"
dat_EB$Group <- "End-binding proteins"
dat_stmn$Group  <- "Tubulin-sequestering\nproteins"
dat_sever$Group <- "MT-severing proteins"
dat_transp$Group <- "MT transport proteins"

df_all_MT <- rbind(
  dat_MT_struct[dat_MT_struct$Change == "Increased",],
  dat_nucl[dat_nucl$Change == "Increased",],
  dat_MAP[dat_MAP$Change == "Increased",],
  dat_EB[dat_EB$Change == "Increased",],
  dat_stmn[dat_stmn$Change == "Increased",], 
  dat_sever[dat_sever$Change == "Increased",]#, 
  #dat_transp[dat_transp$Change == "Increased",]
)

df_all_MT <- df_all_MT %>% 
  group_by(Group, Gene, Gene_phosphosite) %>%
  tally() %>%
  arrange(match(Group, unique(Group)), Gene_phosphosite)

#To order the conditions manually
x <- c("MT-binding proteins", 
       "End-binding proteins", "Tubulin-sequestering\nproteins",
       "Structure proteins", "Nucleators", 
       "MT-severing proteins", "MT transport proteins")

df_all_MT <- df_all_MT %>% 
  mutate(Group =  factor(Group, levels = x)) %>%
  arrange(Group) %>%
  rownames_to_column()

df_all_MT$rowname <- as.numeric(df_all_MT$rowname )
df_all_MT$Group <- as.factor(df_all_MT$Group )

write.table(df_all_MT, paste(output_dir, "/", norm_opt, "_counts_of_MT_proteins.csv", sep = ""), row.names = FALSE)

MT_cols <- c("Structure proteins" = color_MT_struct, 
             "Nucleators" = color_nucl, 
             "MT-binding proteins" = color_MAP, 
             "End-binding proteins" = color_EB, 
             "Tubulin-sequestering\nproteins" = color_stmn, 
             "MT-severing proteins" = color_sever#, 
             # "MT transport proteins" = color_transp
)

p_count_MT_groups <- ggplot(df_all_MT, aes(x = reorder(Gene, -rowname))) + 
  geom_bar(aes(fill = Group)) +
  theme_bw(base_size = 10)+
  theme(plot.title = element_text(color = "black", size = 10),
        axis.ticks = element_line(colour = "black"),
        axis.text.x = element_text(color = "black"),
        axis.text.y = element_text(color = "black"),
        #legend.text = element_text(size = 10),
        panel.grid.minor = element_blank(),
        panel.grid.major.y = element_blank(),
        legend.position = "none",
        axis.title.y = element_blank()
  ) +
  coord_flip() +
  scale_x_discrete() +
  scale_y_continuous(name = "Number of phosphosites", limits = c(0, 30), breaks = seq(0, 100, 5)) +
  scale_fill_manual(values = MT_cols)

p_count_MT_groups
ggsave(p_count_MT_groups, filename = paste(output_dir, "/", norm_opt, "_geom_bar_MTsystem_woTransp.png", sep = ""), type = "cairo", width = 2.5, height = 3)
ggsave(p_count_MT_groups, filename = paste(output_dir, "/", norm_opt, "_geom_bar_MTsystem_woTransp.svg", sep = ""), width = 2.5, height = 3)

# ---------------------------------------------------------------------------- #
# --- Bar plots for MT transport --------------------------------------------- #
# ---------------------------------------------------------------------------- #

df_all_MT <- rbind(
  dat_transp[dat_transp$Change == "Increased",]
)
write.table(df_all_MT, paste(output_dir, "/", norm_opt, "_MTsystem_Transp.csv", sep = ""), row.names = FALSE)

df_all_MT <- df_all_MT %>% 
  group_by(Group, Gene, Gene_phosphosite) %>%
  tally() %>%
  arrange(match(Group, unique(Group)), Gene_phosphosite)

# To order the conditions manually
MT_levels <- c("MT-binding proteins", "End-binding proteins", "Tubulin-sequestering\nproteins", "Structure proteins", "Nucleators", "MT-severing proteins", "MT transport proteins")

df_all_MT <- df_all_MT %>% 
  mutate(Group =  factor(Group, levels = MT_levels)) %>%
  arrange(Group) %>%
  rownames_to_column()

df_all_MT$rowname <- as.numeric(df_all_MT$rowname )
df_all_MT$Group <- as.factor(df_all_MT$Group )

write.table(df_all_MT, paste(output_dir, "/", norm_opt, "_counts_of_MT_proteins_Transp.csv", sep = ""), row.names = FALSE)

MT_cols <- c("Structure proteins" = color_MT_struct, 
             "Nucleators" = color_nucl, 
             "MT-binding proteins" = color_MAP, 
             "End-binding proteins" = color_EB, 
             "Tubulin-sequestering\nproteins" = color_stmn, 
             "MT-severing proteins" = color_sever, 
             "MT transport proteins" = color_transp
)

p_count_MT_groups <- ggplot(df_all_MT, aes(x = reorder(Gene, rowname))) +
  geom_bar(aes(fill = Group)) +
  theme_bw(base_size = 10)+
  theme(plot.title = element_text(color = "black", size = 10),
        axis.ticks = element_line(colour = "black"),
        axis.text.x = element_text(color = "black"),
        axis.text.y = element_text(color = "black"),
        #legend.text = element_text(size = 10),
        panel.grid.minor = element_blank(),
        panel.grid.major.y = element_blank(),
        legend.position = "none",
        axis.title.y = element_blank()
  ) +
  coord_flip() +
  scale_x_discrete() +
  scale_y_continuous(name = "Number of phosphosites", limits = c(0, 30), breaks = seq(0, 100, 5)) +
  scale_fill_manual(values = MT_cols)

p_count_MT_groups
ggsave(p_count_MT_groups, filename = paste(output_dir, "/", norm_opt, "_geom_bar_MTsystem_Transp.png", sep = ""), type = "cairo", width = 2, height = 2.5)
ggsave(p_count_MT_groups, filename = paste(output_dir, "/", norm_opt, "_geom_bar_MTsystem_Transp.svg", sep = ""), width = 2, height = 2.5)

# ---------------------------------------------------------------------------- #
# --- MAP1b phosphoprofile scheme -------------------------------------------- #
# ---------------------------------------------------------------------------- #

dat_MAP_increased <- dat_MAP[dat_MAP$Change == "Increased",]
df_MAP1b <- dat_MAP_increased[dat_MAP_increased$Gene == "Map1b",]

df_MAP1b$Phosphoposition <- gsub('\\D+','', df_MAP1b$Phosphosite)
df_MAP1b$Phosphoposition <- as.numeric(df_MAP1b$Phosphoposition)
df_MAP1b$Plotting <- 1
.labs_sele_MAP1b <- df_MAP1b$Phosphosite

plot_phosphosites <- ggplot(df_MAP1b, aes(x = Phosphoposition, y = Plotting)) +
  geom_point(stat = "identity", color = col_exp, size = 1) +
  geom_hline(yintercept = 0, color = "black") +
  ggtitle("Phosphoprofile of MAP1B") +
  xlab("Residue number") +
  scale_x_continuous(limits = c(0,2774), breaks = seq(0, 2800, 200)) +
  geom_text_repel(data = df_MAP1b,
                  #direction = "x", nudge_x = -2, force = 20,
                  aes(label = .labs_sele_MAP1b),
                  color = col_exp,
                  size = 2, max.overlaps = Inf) +
  geom_linerange(aes(x = Phosphoposition, ymax = Plotting, ymin = 0),
                 size = 0.2,
                 position = position_jitter(height = 0L, seed = 1L)) +
  theme_classic() +
  theme(
    plot.title = element_blank(), # element_text(color = "black", size = 12, hjust = 0.5),
    axis.ticks = element_blank(), #element_line(colour = "black"),
    axis.text.x = element_text(color = "black", size = 8),
    axis.text.y = element_blank(), #element_text(color = "black", size = 8),
    legend.position = "none",
    legend.title = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_blank()
  )

plot_phosphosites
ggsave(plot_phosphosites, filename = paste(output_dir, "/", norm_opt, "_MAP1B_phosphoprofile.png", sep = ""), type = "cairo", width = 8, height = 1)
ggsave(plot_phosphosites, filename = paste(output_dir, "/", norm_opt, "_MAP1B_phosphoprofile.svg", sep = ""), width = 8, height = 1)

# ---------------------------------------------------------------------------- #
# ---------------------------------------------------------------------------- #
# ---------------------------------------------------------------------------- #

# p_count_MT_groups <- ggplot(df_all_MT, aes(x = reorder(Gene, -rowname))) + 
#   geom_bar(aes(fill = Group)) +
#   theme_bw(base_size = 10)+
#   theme(plot.title = element_text(color = "black", size = 10),
#         axis.ticks = element_line(colour = "black"),
#         axis.text.x = element_text(color = "black"),
#         axis.text.y = element_text(color = "black"),
#         #legend.text = element_text(size = 10),
#         panel.grid.minor = element_blank(),
#         legend.position = "none"#,
#         #text = element_text(family = "Arial")
#   ) +
#   coord_flip() +
#   scale_x_discrete(name = "") +
#   scale_y_continuous(name = "Number of phosphosites]") +
#   scale_fill_manual(values = MT_cols)
# 
# ggsave(p_count_MT_groups, filename = paste(output_dir, "/", norm_opt, "_geom_bar_MTsystem.png", sep = ""), type = "cairo", width = 3, height = 5)
# ggsave(p_count_MT_groups, filename = paste(output_dir, "/", norm_opt, "_geom_bar_MTsystem.svg", sep = ""), width = 3, height = 5)

# # ##### Heatmap for samples 
# #x_g3bp1 <- dat_proc_out[grepl("G3bp1", dat_proc_out$Gene),]
# #SG_proteins <- read.csv("ab1037_rat_orthologs.txt", dec = ".", header = F)
# #x_g3bp1 <- dat_proc_out[toupper(dat_proc_out$Gene) %in% toupper(SG_proteins$V1),]
# #x_g3bp1 <- x_g3bp1[!duplicated(x_g3bp1[,c('Peptide')]),]
# x_g3bp1 <- dat_proc_out[grepl("Mapt",dat_proc_out$Gene),]
# x_g3bp1 <- unite_(x_g3bp1, "peptide_phosphosite", c("Peptide", "Phosphosite"), remove = FALSE)
# # create heatmap using pheatmap
# data <- x_g3bp1[c("peptide_phosphosite", "control_1", "control_2", "control_3", "exp_1", "exp_2", "exp_3")]
# rownames(data) <- data$peptide_phosphosite
# data <- data[c("control_1", "control_2", "control_3", "exp_1", "exp_2", "exp_3")]
# data_matrix <- as.matrix(data)
# pheatmap(data_matrix, cluster_cols = F)
# #z-score
# cal_z_score <- function(x){
#   (x - mean(x)) / sd(x)
# }
# data_norm <- t(apply(data_matrix, 1, cal_z_score))
# data_norm <- na.omit(data_norm)
# pheatmap(data_norm)
# ##### Heatmap for groups
# data <- x_g3bp1[c("Peptide", "group_control", "group_exp")]
# rownames(data) <- data$Peptide
# data <- data[c("group_control", "group_exp_A")]
# data_matrix <- as.matrix(data)
# pheatmap(data_matrix, cluster_cols = F)
# #z-score
# cal_z_score <- function(x){
#   (x - mean(x)) / sd(x)
# }
# data_norm <- t(apply(data_matrix, 1, cal_z_score))
# data_norm <- na.omit(data_norm)
# pheatmap(data_norm, cluster_cols = F)
