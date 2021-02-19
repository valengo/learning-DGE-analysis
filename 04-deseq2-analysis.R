# Title     : Análise de expressão diferencial com DESeq2
# Objective : Aplicar o pacote DESeq2 para fazer análise de expressão diferencial de genes
# Created by: valengo
# Created on: 09/02/21

# Carrega o arquivo com as variáveis que definimos para esse projeto. Ex: o ID do projeto TCGA de interesse.
# Se você já fez isso nessa sessão ou rodou os comandos no arquivo em si, não precisa repetir.
source("01-project-setup.R")

# Lê as tabelas que criamos e salvamos em arquivos na fase de pré-processamento.
# dos dados do TCGA.
count_table <- read.table(count_table_filename)

# Lê os dados de grupos que cada amostra na tabela count_table pertence.
# Essa lista também foi gerada na fase de pré-processamento.
groups_file_connection <- file(groups_filename)
groups <- readLines(groups_file_connection)
close(groups_file_connection)


# Cria um data.frame que é necessário como entrada para o DESeq2.
# As linhas desse data.frame correspondem aos nomes das colunas de countData,
# que são os identificadores das amostras.
design <- data.frame(
  row.names = colnames(count_table),
  condition = as.factor(groups)
)

# DESeq2 realiza as normalizações internamente.
# Portanto, não precisamos rodar algum comando para normalizar por tamanho de biblioteca, por exemplo.
dds <- DESeq2::DESeqDataSetFromMatrix(countData = count_table, colData = design, design = ~condition)

# Embora não seja necessário normalizar, podemos realizar um pré-filtro para remover genes com baixas
# contagens considerando todas as amostras. Remover linhas com poucas reads é útil porque reduz a quantidade de
# memória RAM necessária e o tempo para rodar as análises.
# Pré-filtra os dados para manter apenas linhas com pelo menos 10 reads totais. Em outras palavras,
# são mantidas para análises somente os genes que existem pelo menos 10 reads considerando todas as amostras.
keep <- rowSums(DESeq2::counts(dds)) >= 10
dds <- dds[keep,]

# Realiza análise de expressão diferencial com base na Distribuição binomial negativa (a.k.a Gamma-Poisson).
dds <- DESeq2::DESeq(dds)

# Coleta os resultados da análise.
dge_results <- DESeq2::results(dds)

# Ordena os resultados pelo valor de p.
resOrdered <- dge_results[order(dge_results$pvalue),]

# Se quiser, pode mudar o Inf (de infinito) para o valor de corte desejado para p.
# Caso contrário, todos os genes serão mantidos na tabela e serão salvos no arquivo.
resSig <- subset(resOrdered, padj < Inf)

# Salva os resultados (tabela) das análises em um arquivo.
write.table(resSig, file = paste(tables_dir, paste0(TCGA_project,  "-TumorXNormal-DESeq2.tsv"), sep="/"))
