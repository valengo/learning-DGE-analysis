# Title     : Project's dependecies
# Objective : To unify all depedencies in a unique place
# Created by: valengo
# Created on: 09/02/21


# In this project we're using the renv package to create a project-local R dependency management
# The idea is to create a virtual environment that isolates this project from the base environment
# In other words, renv allows us to create a local library isolated from our regular user library.
# In addition, by using renv it's possible to create reproducible environments for our projects
# Because package names, versions and sources will be installed as specified in the virtual environment
# You can read more about renv on https://rstudio.github.io/renv/
if (!requireNamespace("renv", quietly = TRUE))
  install.packages("renv")

# Run this in order to install the dependencies controlled by renv and you're ready to go
renv::restore()

# Run this ONLY if you want to update the virtual environment after installing or removing a package
renv::snapshot()

# In case you DON'T WANT TO use renv
# Install packages manually
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

packages <- c("TCGAbiolinks", "edgeR", "dplyr", "DESeq2", "pheatmap")
BiocManager::install(packages, update = TRUE, ask = FALSE)


# init the virtual environment
# renv::init(bare = TRUE)