# Title     : Project setup
# Objective : A place to organize TCGA project ID and file prefixes
# Created by: valengo
# Created on: 09/02/21

# ID for TCGA project of interest
# The ID is used to query and download data from TCGA's server
# Find project's IDs on https://portal.gdc.cancer.gov/projects
TCGA_project <- "TCGA-CESC"

# These prefixes are used to make sure we don't mix filenames while conducting our analysis
tumor_prefix <- "_tumor_count.txt"
normal_prefix <- "_normal_count.txt"

# Tumor data path and filenames
tumor_data_path <- paste("GDCdata", "Tumor", sep = "/")
tumor_filename <-  paste0(TCGA_project, tumor_prefix)

# Normal data path and filenames
normal_data_path <- paste("GDCdata", "Normal", sep = "/")
normal_file_name <- paste0(TCGA_project, normal_prefix)

# Create a dir to export our tables
dir.create("Tables")

count_table_filename <- paste("Tables", paste0(TCGA_project,  "-TumorXNormal.tsv"), sep="/")
groups_filename <- paste("Tables", paste0(TCGA_project, "-groups.tsv"), sep = "/")

