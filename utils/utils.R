# Title     : Utils
# Objective : To provide utility functions to work with RNAseq data
# Created by: valengo
# Created on: 07/02/21

download_TCGA_data <- function(project, tissueType, directory) {
  # Build a query to retrieve data from TCGA for a given tissue type
  query <- TCGAbiolinks::GDCquery(project = project,
                  data.category = "Transcriptome Profiling",
                  data.type = "Gene Expression Quantification",
                  legacy = FALSE,
                  experimental.strategy = "RNA-Seq",
                  workflow.type = "HTSeq - Counts",
                  sample.type = tissueType)

  # Download data using the query and save it on a given directory
  TCGAbiolinks::GDCdownload(query, directory = directory)
}

save_TCGA_data_as_table <- function(path, filename) {
  # List count files on path
  data <- list.files(path = path, pattern = ".htseq.counts.gz", recursive = TRUE)

  # Read and merge a set of text files containing gene expression counts
  DG1 <- edgeR::readDGE(path = path, data, header = FALSE)

  # Save count data only
  write.table(DG1$counts, file = paste(path, filename, sep="/"))
}

