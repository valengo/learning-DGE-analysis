# Title     : Configurações do Projeto
# Objective : Organizar o ID do projeto TCGA de interesse e prefixos de nomes de arquivos
# Created by: valengo
# Created on: 09/02/21

# ID TCGA do projeto de interesse.
# O ID é utilizado para buscar e baixar dados do servidor do GDC.
# Encontre projetos e seus IDs em https://portal.gdc.cancer.gov/projects.
TCGA_project <- "TCGA-CESC"

# Esses prefixos são utilizados para tentar garantir que
# não vamos misturar os nomes dos arquivos nas nossas anaálises
tumor_prefix <- "_tumor_count.txt"
normal_prefix <- "_normal_count.txt"

# Nome de pasta e arquivos relacionados aos dados de tumor.
# Nós utilizamos diversas vezes as funções paste e paste0 para concatenar strings
# de forma a criar os nomes dos arquivos e seus caminhos.
GDC_data_dir <- "../GDCdata"
tumor_data_path <- paste(GDC_data_dir, "Tumor", sep = "/")
tumor_filename <-  paste0(TCGA_project, tumor_prefix)

# Nome de pasta e arquivos relacionados aos dados de normais.
normal_data_path <- paste(GDC_data_dir, "Normal", sep = "/")
normal_file_name <- paste0(TCGA_project, normal_prefix)

# Cria um diretório com o nome "Tables" para salvar dados e resultados em arquivos.
tables_dir <- "../Tables"
dir.create(tables_dir)

# Define nome para arquivo da tabela de contagem de reads por gene.
# Essa tabela será gerada durante a fase de pré-processamento dos dados.
count_table_filename <- paste(tables_dir, paste0(TCGA_project,  "-TumorXNormal.tsv"), sep="/")

# Define nome para arquivo da tabela de grupos das amostras (tumor ou normal).
# Essa tabela será gerada durante a fase de pré-processamento dos dados.
groups_filename <- paste(tables_dir, paste0(TCGA_project, "-groups.tsv"), sep = "/")

