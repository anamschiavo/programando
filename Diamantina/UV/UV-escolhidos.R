#Script Iniciação Científica Sabrina Pinheiro

#REGRAS DE ENTRADA:
#*************ORDEM DAS COLUNAS: 1. Espécie | 2. Dose | 3-5. UFC | 6. Fator de Diluição **********************************
#Decimal é dado por vírgula
#Primeira linha é cabeçalho

#********** PACOTES **********
install.packages('ggplot2')   # Pacote para gráficos
library(ggplot2)
install.packages('matrixStats')   # Pacote para média e devio padrão de matrizes, no caso das gotas
library(matrixStats)
install.packages('MASS')    # Pacote necessário para a notação certa no eixo Y do gráfico
library(MASS)
install.packages('scales')    # Pacote para transformação em log10 do eixo y do gráfico
library(scales)
install.packages('dplyr')
library(dplyr)

#********** LEITURA **********
# Lê tabela em tsv, ou seja, separado por tab, vírgula como decimal, primeira linha é cabeçalho, preenche espáços vazios com NA
uvtab <- read.table("UVCtab.tsv", sep="\t", dec=",", quote="", header=TRUE, fill=TRUE)

#********** PRÉ-PROCESSAMENTO **********
ufc <- as.matrix(uvtab[,c(3,4,5)])    # Cria matriz com as contagens das 3 gotinhas
nome <- uvtab[,1]
dose <- uvtab[,2]   # Cria um vetor com doses de UVC dadas
media <- rowMeans(ufc, na.rm=TRUE)    # Tira média das gotinhas
media_cor <- media*uvtab[,6]*100    # Faz média corrigida pegando média, multiplicando pelo fator de diluição e por 100 para ser UFC/mL
desv_pad <- rowSds(ufc, na.rm=TRUE)   # Calcula o desvio padrão das três gotinhas
n <- rowSums(!is.na(ufc))   # Calcula quantas gotinhas tem (3 ou 2)
erro_pad <- desv_pad/(sqrt(n))    # Calcula erro padrão que, por definição, é desvio padrão sobre raiz quadrada do n
erro_media <- ((uvtab[,6]*100))*(erro_pad)    # Erro associado as médias corrigidas. Propagação do erro considerando normal.
media_cor_0 <- rep(c(media_cor[1], media_cor[9], media_cor[17], media_cor[25], media_cor[33], media_cor[41], media_cor[49], media_cor[57], media_cor[65], media_cor[73], media_cor[81], media_cor[89]), each=8)    # Apenas os valores docontrole repetidos para fazer a sobrevivência
erro_media_0 <- rep(c(erro_media[1], erro_media[9], erro_media[17], erro_media[25], erro_media[33], erro_media[41], erro_media[49], erro_media[57], erro_media[65], erro_media[73], erro_media[82], erro_media[89]), each=8)   # Apenas valores do erro associados aos controles para fazer o erro das sobrevivências
sobrev <- media_cor/media_cor_0   # Calcula sobrevivência (experimento sobre controle)
sobrev <- trunc(sobrev*10^5)/10^5   # Trunca valos na 4a casa decimal
erro_sobrev <- sqrt((sobrev^2)*(((erro_media/media_cor)^2)+((erro_media_0/media_cor_0)^2)))   # Calcula erro associado a sobrevivência com as fórmulas de propagação de erro
graftab <- cbind.data.frame(nome, dose, sobrev, erro_sobrev)   # Cria tabela para fazer os gráficos, com espécie, doses, sobrevivência e erro associado a sobrevivência

#********** GRÁFICOS **********
graf_sob_facet<- ggplot(graftab, aes(x=dose, y=sobrev, color=nome)) +
    geom_point(size=1)+   # Coloca as bolinhas
    geom_line(linewidth=.7)+   # Coloca as linhas
    geom_errorbar(aes(ymin=sobrev-erro_sobrev, ymax=sobrev+erro_sobrev), width=0.15, size=.5)+    # Coloca as barras de erro
    facet_wrap(~nome)+    # Deixa facetado por espécie
    theme_bw()+   # Muda a aparência para o tema certo
    xlab(bquote("Fluence"~(J/m^2)))+    # Coloca o nome do eixo x
    ylab("Survival")+    # Coloca o nome do eixo y
    theme(
      panel.grid.major = element_line(size=.5),   # Define largura da linha maior do grid interno no gráfico
      panel.grid.minor = element_line(size=.2),   # Define largura da linha menor do grid interno do gráfico
      legend.position='none')+    # Retira a legenda
    scale_y_log10(limits = c(1e-4,10),   # Transforma o eixo y em log10 com limiter de 1 a 10^-5
                  labels = trans_format("log10", math_format(10^.x)),
                   breaks = trans_breaks("log10", function(x) 10^x, n=4))

ggsave('facetado.png')    # Salva gráfico que acabamos de criar como um arquivo
