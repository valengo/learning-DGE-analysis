# Title     : Pré-processamento dos dados
# Objective : Baixar e pré-processar dados do TCGA
# Created by: valengo
# Created on: 05/02/21

# Carrega funções criadas para facilitar o download de dados
# de amostras de tecidos tumorais e normais de um mesmo projeto TCGA.
# Você pode abrir o arquivo (utils/utils.R) e investigar a utilização das funções do pacote TCGAbiolinks.
source("utils/utils.R")

# Carrega o arquivo com as variáveis que definimos para esse projeto. Ex: o ID do projeto TCGA de interesse.
# Se você já fez isso nessa sessão ou rodou os comandos no arquivo em si, não precisa repetir.
source("01-project-setup.R")

# Baixa e salva dados de contagem de leituras por gene a partir de dados de RNAseq
# de amostras tumorais do TCGA usando funções definidas em utils.R.
download_TCGA_data(TCGA_project, "Primary Tumor", tumor_data_path)
save_TCGA_data_as_table(paste(tumor_data_path, TCGA_project, sep = "/"), tumor_filename)

# Baixa e salva dados de contagem de leituras por gene a partir de dados de RNAseq
# de amostras não tumorais do TCGA usando funções definidas em utils.R.
save_TCGA_data_as_table(paste(normal_data_path, TCGA_project, sep = "/"), normal_file_name)

# Lê os dados que foram baixados em duas tabelas (tumor e normal).
# Nessas duas tabelas cada coluna representa uma amostra e cada linha um gene.
tumor_data <- read.table(paste(tumor_data_path, TCGA_project, paste0(TCGA_project, tumor_prefix), sep = "/"))
normal_data <- read.table(paste(normal_data_path, TCGA_project, paste0(TCGA_project, normal_prefix), sep = "/"))

# O projeto TCGA que estamos usando como exemplo apresenta poucos dados de contagem de leituras por gene
# com livre acesso para amostras de tecidos normais.
# Portanto, vamos utilizar um número menor de dados tumorais do que está realmente disponível.
# A ideia é usar o mesmo número de amostras de tecidos normais e tumorais.
# Estamos utilizando a função ncol para saber o número N de colunas (amostras) normais que foram baixadas e
# criamos uma nova tabela que utiliza apenas N colunas das disponíveis na tabela de dados tumorais.
# CUIDADO: você provavelmente não quer fazer isso num contexto de análise real.
tumor_data <- tumor_data[, -seq(from = ncol(normal_data) + 1, to = ncol(tumor_data))]

# Nós vamos montar uma única tabela com todos os dados, mas antes disso vamos montar uma lista com
# a classificação de cada amostra no grupo tumor ou normal.
# Usamos a função ncol para saber quantas amostras temos de cada grupo (lembra que cada coluna é uma amostra?).
# Depois, usamos a função rep é para repetir "tumor" ncol(tumor_data) vezes e "normal" ncol(normal_data) vezes.
# Por fim, agrupamos numa lista única com a função c()
groups <- c(rep("tumor", ncol(tumor_data)), rep("normal", ncol(normal_data)))

# Agrupa as duas tabelas de amostras (tumor e normal) em uma única.
count_table <- merge(tumor_data, normal_data, by = "row.names", all=TRUE)

# Esse pacote nos proporciona utilizar operações "pipe", como por exemplo > (%>%)
library(magrittr)

# Essas meta tags são adicionadas pelo programa gtseq quando ele está contando o número de leituras por gene
meta_tags <- c("__alignment_not_unique", "__ambiguous", "__no_feature", "__not_aligned", "__too_low_aQual")
# Remove as linhas com meta tags da tabela de contagem de leituras
count_table <- count_table %>% dplyr::filter(!Row.names %in% meta_tags)

# Se você abrir o arquivo ou espiar a tabela de contagem de reads com head(), você vai notar que os IDs dos genes
# estão com números extras depois do ".". Por exemplo: ENSG00000000.00
# É recomendado que esses números sejam removidos para evitar problemas em futuras anotações.
# Cria uma lista com os números adicionais removidos dos IDs dos genes
gene_ids <- lapply(strsplit(count_table$Row.names, "\\."), "[[", 1)
# Substitui os nomes das linhas (índices) com a nova lista de IDs de genes sem os números e o "." extras
count_table <- transform(count_table,
                         row.names=gene_ids, Row.names=NULL)

# Salva a count_table em um arquivo
write.table(count_table, file = count_table_filename)
# Salva o groups em um arquivo
writeLines(groups, file(groups_filename))
