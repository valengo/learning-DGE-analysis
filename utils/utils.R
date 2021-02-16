# Title     : Utils
# Objective : Fornecer funções para facilitar o trabalho com dados de RNAseq
# Created by: valengo
# Created on: 07/02/21

download_TCGA_data <- function(project, tissueType, directory) {
  # Constrói uma query para buscar dados do TCGA de um projeto e tipo de amostra específicos (ex: Tumor)
  query <- TCGAbiolinks::GDCquery(project = project,
                  data.category = "Transcriptome Profiling",
                  data.type = "Gene Expression Quantification",
                  legacy = FALSE,
                  experimental.strategy = "RNA-Seq",
                  workflow.type = "HTSeq - Counts",
                  sample.type = tissueType)

  # Baixa os dados usando a query e salva-os no diretório passado para a função
  TCGAbiolinks::GDCdownload(query, directory = directory)
}

save_TCGA_data_as_table <- function(path, filename) {
  # Elenca os arquivos de contagem (com extensão .htseq.counts.gz) de reads no diretório passado (path).
  data <- list.files(path = path, pattern = ".htseq.counts.gz", recursive = TRUE)

  # Lê e concatena os arquivos que apresentam esses dados de contagem de expressão gênica.
  DG1 <- edgeR::readDGE(path = path, data, header = FALSE)

  # Salva a tabela de dados de contagem em um arquivo.
  write.table(DG1$counts, file = paste(path, filename, sep="/"))
}

