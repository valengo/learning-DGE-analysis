# Title     : Differential expression analysis with EdgeR
# Objective : To employ EdgeR for differential expression analysis
# Created by: valengo
# Created on: 09/02/21

source("01-project-setup.R")

# Loading data from files created with 02-TCGA-preprocess
count_table <- read.table(count_table_filename)
groups_file_connection <- file(groups_filename)
groups <- readLines(groups_file_connection)
close(groups_file_connection)

# Employ EdgeR to perform differential expression analysis
# using our data matrix as input together with the list of groups (tumor vs normal)
dge_list <- edgeR::DGEList(count_table, group = groups)

# From docs: filterByExpr function implements the filtering strategy that was intuitively described by Chen et al (2016).
# Roughly speaking, the strategy keeps genes that have at least min.count reads in a worthwhile number samples.
dge_list <- dge_list[edgeR::filterByExpr(dge_list), , keep.lib.sizes=FALSE]

# Library size normalization
dge_list <- edgeR::calcNormFactors(dge_list)
# Maximizes the negative binomial conditional common likelihood to estimate a common dispersion value across all genes
dge_list <- edgeR::estimateCommonDisp(dge_list)
# Estimates tagwise dispersion values by an empirical Bayes method based on weighted conditional maximum likelihood
dge_list <- edgeR::estimateTagwiseDisp(dge_list)

# Classical differential expression test
# The exact test is only applicable to experiments with a single factor
# Check EdgeR's docs for other approaches considering more complex experiments
exact_test <- edgeR::exactTest(dge_list)

# Create a table of the top differentially expressed genes/tags
# Using n = Inf to return all genes instead of only the top 10
# Resulting table is sorted by p-value
top_genes <- edgeR::topTags(exact_test, n = Inf)

# Save EdgeR's results
write.table(top_genes$table, file = paste("Tables", paste0(TCGA_project,  "-TumorXNormal-EdgeR.tsv"), sep="/"))

