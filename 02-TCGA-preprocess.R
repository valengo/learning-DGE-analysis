# Title     : Learning RNAseq
# Objective : To learn about RNAseq analysis
# Created by: valengo
# Created on: 05/02/21

source("utils/utils.R")
source("01-project-setup.R")

# Download and save tumor data from TCGA using functions from utils.R
download_TCGA_data(TCGA_project, "Primary Tumor", tumor_data_path)
save_TCGA_data_as_table(paste(tumor_data_path, TCGA_project, sep = "/"), tumor_filename)

# Download and save tumor data from TCGA using functions from utils.R
download_TCGA_data(TCGA_project, "Solid Tissue Normal", normal_data_path)
save_TCGA_data_as_table(paste(normal_data_path, TCGA_project, sep = "/"), normal_file_name)

# Load normal and tumor data to build a matrix with count values
tumor_data <- read.table(paste(tumor_data_path, TCGA_project, paste0(TCGA_project, tumor_prefix), sep = "/"))
normal_data <- read.table(paste(normal_data_path, TCGA_project, paste0(TCGA_project, normal_prefix), sep = "/"))

# Subset tumor_data so we can work with a similar number of samples for each group
# You can use ncol(tumor_data) and ncol(normal_data) to verify the number of samples
# CAUTION: you probably don't want just subsetting in a real context
tumor_data <- tumor_data[, -seq(from = ncol(normal_data) + 1, to = ncol(tumor_data))]

# Use the function ncol to discover how many tumor and control samples were loaded (number of columns)
# Use rep function to create a list of size ncol(tumor_data) + ncol(normal_data) to specify each sample's group
# by repeating "tumor" ncol(tumor_data) times and "normal" ncol(normal_data) times
groups <- c(rep("tumor", ncol(tumor_data)), rep("normal", ncol(normal_data)))

# Create a count matrix by merging both tables considering their row.names (gene IDs)
count_table <- merge(tumor_data, normal_data, by = "row.names", all=TRUE)

# Load magrittr package that alows us to pipe-like operators such as > (%>%)
library(magrittr)

# These meta tags are added by htseq when counting the number of reads
meta_tags <- c("__alignment_not_unique", "__ambiguous", "__no_feature", "__not_aligned", "__too_low_aQual")
# Filtering out meta tag rows
count_table <- count_table %>% dplyr::filter(!Row.names %in% meta_tags)

# Create a list with extra numbers after . on gene IDs removed
gene_ids <- lapply(strsplit(count_table$Row.names, "\\."), "[[", 1)
# Replace row.names with gene_ids and remove extra column Row.names created by merge function
count_table <- transform(count_table,
                         row.names=gene_ids, Row.names=NULL)

# Save current count_table
write.table(count_table, file = count_table_filename)
# Save groups into a file
writeLines(groups, file(groups_filename))
