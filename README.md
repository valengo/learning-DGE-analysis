# 🤓 Aprendendo análise de expressão diferencial de genes na prática
Aqui você pode encontrar scripts que podem te dar um ponto de partida para trabalhar com expressão diferencial de genes utilizando diversos pacotes R, como [edgeR](http://bioconductor.org/packages/release/bioc/html/edgeR.html) e [DESeq2](https://bioconductor.org/packages/release/bioc/html/DESeq2.html). Além disso, estão disponíveis scripts para baixar e pré-processar dados de contagem de RNAseq do [TCGA](https://www.cancer.gov/about-nci/organization/ccg/research/structural-genomics/tcga) que estão hospedados no [portal GDC](https://portal.gdc.cancer.gov). Esse projeto foi testado utilizando a versão 4.0.3 do R.

# 🎯 Recomendações gerais
Para utilizar esse projeto, você pode clonar esse repositório ou baixar algum [zip de release](https://github.com/valengo/learning-DGE-analysis/releases). Mas antes disso, recomendo que crie uma pasta e daí clone ou baixe o zip nessa pasta que acabou de criar. Assim, os dados baixados serão salvos em uma pasta "irmã" a desse repositório, assim como os dados exportados das análises. Dessa forma fica mais fácil atualizar a versão desses scripts na sua máquina quando você escolher a opção de baixar o zip de release. Basicamente, bastaria baixar novamente o zip na pasta "mãe" que contém as pastas dos dados, tabelas e desse repositório. 

# 🗂 Estrutura do projeto
## Utils
Essa pasta contem scripts de utilidades para ajudar a baixar e pré-processar dados de um projeto TCGA específico usando o pacote [TCGAbiolinks](https://bioconductor.org/packages/release/bioc/html/TCGAbiolinks.html).
## XX-scripts
Considerando que os scripts devem ser utilizados em sequência, eles são prefixados com números que indicam a ordem de execução. Os scripts são bastante comentados para tentar explicar o que cada linha de código faz. A ideia é dar uma ajuda adicional para quem ainda não tem o hábito de programar.
### 🔗 00-dependecies.R
Esse projeto foi configurado para utilizar o [pacote renv](https://rstudio.github.io/renv/articles/renv.html) para criar uma biblioteca de pacotes que é específica para esse projeto e isolada da sua biblioteca local. Recomendo fortemente que você use essa abordagem como especificada no script, porque assim você vai usar os mesmos pacotes, versões e fontes que eu utilizei. Contudo, como você pode ler no script, você pode instalar tudo manualmente também. 
### 🔗 01-project-setup.R
### 🔗 02-TCGA-preprocess.R
### 🔗 03-edger-analysis.R 
# 🚀 Mãos na massa
Abra o projeto usando [Rstudio](https://rstudio.com) ou [PyCharm](https://www.jetbrains.com/pycharm/) com o [plugin do R](https://www.jetbrains.com/help/pycharm/r-plugin-support.html) instalado. Preste atenção na ordem recomendada quando usar os scripts (00, 01 e assim por diante). Antes de rodar qualquer código, tente ler os comentários logo acima de cada linha para entender o que está acontecendo e como você está conseguindo os resultados que está observando.
