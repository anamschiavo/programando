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

#********** FUNÇÕES **********
# Função para média corrigida
calc_media_corrigida <- function(data, col_indices) {   # col_indices -> numero de gotas (ou seja, de colunas de ufc)
  ufc <- as.matrix(data[, col_indices]) #vou pegar só as colunas de ufc
  media_ufc <- rowMeans(ufc, na.rm = TRUE)
  media_corrigida <- media_ufc * data$diluicao * 100 #multiplico por 100 no final pq quero o resultado em ml
  return(media_corrigida)
}

# Função para calcular erro associado a média corrigida
calc_erro_medcor <- function(data, col_indices) {     # col_indices -> numero de gotas (ou seja, de colunas de ufc)
  ufc <- as.matrix(data[, col_indices])
  sd <- rowSds(ufc, na.rm = TRUE)
  n <- rowSums(!is.na(ufc))
  se <- sd / sqrt(n)
  erro <- (data$diluicao)*se*100      # Multiplico por 100 no final pq quero o resultado em ml
  return(erro)
}

log10_minor_break = function (...){
  function(x) {
    minx         = floor(min(log10(x), na.rm=T))-1;
    maxx         = ceiling(max(log10(x), na.rm=T))+1;
    n_major      = maxx-minx+1;
    major_breaks = seq(minx, maxx, by=1)
    minor_breaks =
      rep(log10(seq(1, 9, by=1)), times = n_major)+
      rep(major_breaks, each = 9)
    return(10^(minor_breaks))
  }
}

#********** LEITURA **********
# Lê tabela em tsv, ou seja, separado por tab, vírgula como decimal, primeira linha é cabeçalho, preenche espáços vazios com NA
nacltab <- read.table("NaCltab.tsv", sep="\t", dec=",", quote="", header=TRUE, fill=TRUE)

#********** PROCESSAMENTO **********
col_indices <- c(3, 4, 5)   # Define quais são as colunas da tabela que tem valores de UFC
mc <- calc_media_corrigida(nacltab, col_indices)
erro_mc <- calc_erro_medcor(nacltab, col_indices)

tab_media <- cbind.data.frame(nacltab, mc, erro_mc)

rep1 <- filter(tab_media, Replicata=='1')
rep2 <- filter(tab_media, Replicata=='2')
rep3 <- filter(tab_media, Replicata=='3')

media_rep <- rowMeans(cbind(rep1$mc, rep2$mc, rep2$mc), na.rm = TRUE)
erro_rep <- sqrt(rowSums(cbind((1/3*rep1$erro_mc)^2, (1/3*rep2$erro_mc)^2, (1/3*rep1$erro_mc)^2), na.rm = TRUE))

tab <- cbind.data.frame(rep1$Isolado, rep1$Conc, media_rep, erro_rep)
tab0 <- filter(tab, rep1$Conc=='0')
media0 <- rep(tab0$media_rep, each=7)
erro0 <- rep(tab0$erro_rep, each=7)
sob <- media_rep/media0
erro_sob <- sqrt((sob^2) * (((erro_rep / media_rep)^2) + ((erro0 / media0)^2)))

tabgraf <- cbind.data.frame(rep1$Isolado, rep1$Conc, sob, erro_sob)
colnames(tabgraf) <- c('isolado', 'conc', 'sob', 'erro')

#********** GRÁFICOS **********
graf_sob<- ggplot(tabgraf, aes(x=conc, y=sob)) +
    geom_point(size=.65)+   # Coloca as bolinhas
    geom_line(linewidth=.35)+   # Coloca as linhas
    geom_errorbar(aes(ymin=sob-erro, ymax=sob+erro), width=0.15, linewidth=.5)+    # Coloca as barras de erro
    facet_wrap(~isolado, ncol=4)+    # Deixa facetado por espécie
    theme_bw()+   # Muda a aparência para o tema certo
    xlab(bquote("[NaCl]"~(mol.L^-1)))+    # Coloca o nome do eixo x
    ylab(bquote("Sobrevivência"~(N/N[0])))+    # Coloca o nome do eixo y
    theme(
      panel.grid.major = element_line(size=.25),   # Define largura da linha maior do grid interno no gráfico
      panel.grid.minor = element_line(size=.1),   # Define largura da linha menor do grid interno do gráfico
      legend.position='none')+    # Retira a legenda
    scale_y_log10(limits = c(1e-6,10),   # Transforma o eixo y em log10 com limiter de 1 a 10^-5
                  labels = trans_format("log10", math_format(10^.x)),
                   breaks = trans_breaks("log10", function(x) 10^x, n=4),
                   minor_breaks=log10_minor_break())
graf_sob

ggsave('escolhidos_NaCl.png', graf_sob, device='png', unit='cm', width = 18, height = 25, dpi=300)   # Salva gráfico que acabamos de criar como um arquivo





















#********** PRÉ-PROCESSAMENTO **********
ufc <- as.matrix(nacltab[,c(3,4,5)])    # Cria matriz com as contagens das 3 gotinhas
nome <- nacltab[,1]
dose <- nacltab[,2]   # Cria um vetor com doses de UVC dadas
media <- rowMeans(ufc, na.rm=TRUE)    # Tira média das gotinhas
media_cor <- media*nacltab[,6]*100    # Faz média corrigida pegando média, multiplicando pelo fator de diluição e por 100 para ser UFC/mL
desv_pad <- rowSds(ufc, na.rm=TRUE)   # Calcula o desvio padrão das três gotinhas
n <- rowSums(!is.na(ufc))   # Calcula quantas gotinhas tem (3 ou 2)
erro_pad <- desv_pad/(sqrt(n))    # Calcula erro padrão que, por definição, é desvio padrão sobre raiz quadrada do n
erro_media <- ((nacltab[,6]*100))*(erro_pad)    # Erro associado as médias corrigidas. Propagação do erro considerando normal.
media_cor_0 <- rep(c(media_cor[1], media_cor[8], media_cor[15], media_cor[22], media_cor[29], media_cor[36], media_cor[43], media_cor[50], media_cor[57], media_cor[64], media_cor[71]), each=7)    # Apenas os valores docontrole repetidos para fazer a sobrevivência
erro_media_0 <- rep(c(erro_media[1], erro_media[8], erro_media[15], erro_media[22], erro_media[29], erro_media[36], erro_media[43], erro_media[50], erro_media[57], erro_media[64], erro_media[71]), each=7)   # Apenas valores do erro associados aos controles para fazer o erro das sobrevivências
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
    xlab("[NaCl] (mol/L)")+    # Coloca o nome do eixo x
    ylab("Survival")+    # Coloca o nome do eixo y
    theme(
      panel.grid.major = element_line(size=.5),   # Define largura da linha maior do grid interno no gráfico
      panel.grid.minor = element_line(size=.2),   # Define largura da linha menor do grid interno do gráfico
      legend.position='none')+    # Retira a legenda
    scale_y_log10(limits = c(1e-5,10),   # Transforma o eixo y em log10 com limiter de 1 a 10^-5
                  labels = trans_format("log10", math_format(10^.x)),
                   breaks = trans_breaks("log10", function(x) 10^x, n=5))

ggsave('facetado.png')    # Salva gráfico que acabamos de criar como um arquivo
