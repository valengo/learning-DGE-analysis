# Learning differential gene expression analysis by doing
Here you can find scripts that can give you a starting point to work with gene expression data and differential expression analysis with different R packages, such as edgeR and deseq2. In addition, there are scripts to download and preprocess TCGA RNA-seq count files hosted on GDC portal. This project was tested using R version 4.0.3.

# Project structure
## Utils
The utils folder contains utility scripts to help us to download and preprocess TCGA data for a specific project ID making use of the TCGAbiolinks package.
## XX-scripts
Considering that scripts should be used in sequence, they are named using numbers as prefixes in their filenames. Scripts are heavily commented to provide a basic explanation of what is going on in almost every line. The idea is to provide extra help to people who are not so used to coding.
### 00-dependecies.R
This project is set to use the renv package in order to create a project-specific library that is isolated from your user library. I strongly recommend that you use this approach as specified in the script, because you'll then use packages with the same versions and sources as I used when testing everything. However, as you can read in the script as well, you can also install everything manually.
### 01-project-setup.R
### 02-TCGA-preprocess.R
### 03-edger-analysis.R 
## How to get started?
Open the project using Rstudio or PyCharm with the R plugin installed. Pay attention to the recommended sequence when using the scripts (00, 01 and so on). Before running any code, try to read each comment to grasp what is going on and how you're achieving the results you're observing. 
