# Title     : Dependências de Projeto
# Objective : Unificar todas as dependências em um único lugar
# Created by: valengo
# Created on: 09/02/21

# Nesse projeto estamos utilizando o pacote renv para gerenciar as dependências do projeto de forma local.
# A ideia é criar um ambiente virtual que isole esse projeto do ambiente base que está na sua máquina.
# Em outras palavras, o renv permite que criemos uma biblioteca (library) específica para o projeto que é
# isolada da sua library de usuário(a). Todos os pacotes vão ser instalados na pasta 'renv' desse projeto.
#
# Além disso, quando usamos renv é possível criar ambientes reproduzíveis para nossos projetos.
# Porque os nomes dos pacotes, versões e fontes serão instalados como especificado no ambiente virtual (renv.lock)
# Você pode aprender mais sobre renv em https://rstudio.github.io/renv/
if (!requireNamespace("renv", quietly = TRUE))
  install.packages("renv")

# Para instalar os pacotes utilizando renv, basta rodar o comando abaixo.
renv::restore()

# Utilize o comando abaixo SOMENTE se você quiser atualizar o ambiente virtual
# depois de instalar ou remover algum pacote.
renv::snapshot()

# Caso você NÃO QUEIRA usar o renv, instale os pacotes manualmente:
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

packages <- c("TCGAbiolinks", "edgeR", "dplyr", "DESeq2", "pheatmap")
BiocManager::install(packages, update = TRUE, ask = FALSE)


# init the virtual environment
# renv::init(bare = TRUE)