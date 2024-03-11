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
install.packages('ggtext')
library(ggtext)

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
x <- trunc(sob*10^6)/10^6
erro_sob <- sqrt((sob^2) * (((erro_rep / media_rep)^2) + ((erro0 / media0)^2)))

tabgraf <- cbind.data.frame(rep1$Isolado, rep1$Conc, x, erro_sob)
colnames(tabgraf) <- c('isolado', 'conc', 'sob', 'erro')
tabgraf$isolado <- factor(tabgraf$isolado, labels=c("D14.15.010-3", "D14.15.100-1", "D14.30.010-1", "D14.30.010-5", "D20.15.010-12", "D20.30.010-2", "D21.30.010-2", "D23.15.010-1", "D23.30.010-1", "D23.30.010-3", "D4.15.R2A-5", "D4.30.100-3", "D4.30.R2A-2", "D4.30.R2A-4", "D5.30.010-1", "D7.15.010-6", "D7.15.R2A-6", "D7.30.100-9", "D8.15.R2A-3", "italic('E. coli')*textstyle(' MG1655')"))
tabgraf$isolado <- factor(tabgraf$isolado)
#********** GRÁFICOS **********
graf_sob<- ggplot(tabgraf, aes(x=conc, y=sob)) +
    geom_point(size=.65)+   # Coloca as bolinhas
    geom_line(linewidth=.35)+   # Coloca as linhas
    geom_errorbar(aes(ymin=sob-erro, ymax=sob+erro), width=0.15, linewidth=.5)+    # Coloca as barras de erro
    facet_wrap(~isolado, ncol=4, labeller=label_parsed)+    # Deixa facetado por espécie
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
