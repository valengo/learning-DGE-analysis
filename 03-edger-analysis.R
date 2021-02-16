# Title     : Análise de expressão diferencial com edgeR
# Objective : Aplicar o pacote edgeR para fazer análise de expressão diferencial de genes
# Created by: valengo
# Created on: 09/02/21

# Carrega o arquivo com as variáveis que definimos para esse projeto. Ex: o ID do projeto TCGA de interesse.
# Se você já fez isso nessa sessão, não precisa repetir
source("01-project-setup.R")

# Lê as tabelas que criamos e salvamos em arquivos na fase de pré-processamento
# dos dados do TCGA.
count_table <- read.table(count_table_filename)

# Lê os dados de grupos que cada amostra na tabela count_table pertence.
# Essa lista também foi gerada na fase de pré-processamento.
groups_file_connection <- file(groups_filename)
groups <- readLines(groups_file_connection)
close(groups_file_connection)

# Começa a usar edgeR para análise de expressão gênica passando a tabela
# de contagem de leituras e a lista de grupos como entrada.
# Essa função basicamente serve para construir um objeto que vai armazenar
# os dados necessários para análise e outros dados intermediários que serão gerados durante esse processo.
dge_list <- edgeR::DGEList(count_table, group = groups)

# Após a inicialização da DGEList, o primeiro passo é filtrar genes que
# não são muito expressos nas amostras normais e tumorais.
# O método utilizado pela função abaixo foi descrito por Chen et al (2016).
dge_list <- dge_list[edgeR::filterByExpr(dge_list), , keep.lib.sizes=FALSE]

# Normaliza as contagens por tamanho da biblioteca (numero de leituras totais por amostra).
dge_list <- edgeR::calcNormFactors(dge_list)
# Estima um valor de dispersão comum entre todos os genes.
dge_list <- edgeR::estimateCommonDisp(dge_list)
# Estima mais valores de dispersão.
dge_list <- edgeR::estimateTagwiseDisp(dge_list)

# Realiza um teste clássico de expressão diferencial.
# Esse teste é aplicável somente para experimentos com um fator único.
# Verifique a documentação para abordagens que consideram experimentos mais complexos.
exact_test <- edgeR::exactTest(dge_list)

# Cria uma tabela com os genes que tem melhores evidências de serem diferencialmente expressos no topo.
# Ao utilizar n = Inf (de inifinito) todos os genes são retornados, ao invés dos 10 "melhores".
# A tabela resultante está ordenada pelo valor de p.
top_genes <- edgeR::topTags(exact_test, n = Inf)

# Salva a tabela de resultados em um arquivo.
write.table(top_genes$table, file = paste(tables_dir, paste0(TCGA_project,  "-TumorXNormal-EdgeR.tsv"), sep="/"))

