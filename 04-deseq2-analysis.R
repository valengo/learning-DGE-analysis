# Title     : Differential expression analysis with DESeq2
# Objective : To employ DESeq2 for differential expression analysis
# Created by: valengo
# Created on: 09/02/21

source("01-project-setup.R")

# Loading data from files created with 02-TCGA-preprocess
count_table <- read.table(count_table_filename)
groups_file_connection <- file(groups_filename)
groups <- readLines(groups_file_connection)
close(groups_file_connection)

# Create a data.frame as needed for input in DESeq2
# Rows of this data.frame correspond to columns of countData (sample's identifiers)
design <- data.frame(
  row.names = colnames(count_table),
  condition = as.factor(groups)
)

# The DESeq2 model internally corrects for library size, so transformed or normalized values such as
# counts scaled by library size should not be used as input
dds <- DESeq2::DESeqDataSetFromMatrix(countData = count_table, colData = design, design = ~condition)

# Perform differential expression analysis based on the Negative Binomial (a.k.a. Gamma-Poisson) distribution
dds <- DESeq2::DESeq(dds)

# Collect results from a DESeq analysis
dge_results <- DESeq2::results(dds)

# Results ordered by pvalue
resOrdered <- dge_results[order(dge_results$pvalue),]

# You can also filter by p-value by replacing Inf (Infinite) with your desired value
# When using Inf, all genes will be kept
resSig <- subset(resOrdered, padj < Inf)

write.table(resSig, file = paste("Tables", paste0(TCGA_project,  "-TumorXNormal-DESeq2.tsv"), sep="/"))

