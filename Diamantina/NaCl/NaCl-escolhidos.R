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
tab0 <- filter(tab_media, Conc==0)

media0 <- rep(tab0$mc, each=7)
erro0 <- rep(tab0$erro_mc, each=7)
sob <- mc/media0
erro_sob <- sqrt((sob^2) * (((erro_mc / mc)^2) + ((erro0 / media0)^2)))
tab_total <- cbind.data.frame(tab_media, sob, erro_sob)

rep1 <- filter(tab_total, Replicata=='1')
rep2 <- filter(tab_total, Replicata=='2')
rep3 <- filter(tab_total, Replicata=='3')

media_rep <- rowMeans(cbind(rep1$sob, rep2$sob, rep3$sob), na.rm = TRUE)
media_rep <- trunc(media_rep*10^4)/10^4
erro_rep <- sqrt(rowSums(cbind((1/3*rep1$erro_sob)^2, (1/3*rep2$erro_sob)^2, (1/3*rep3$erro_sob)^2), na.rm = TRUE))

tabgraf <- cbind.data.frame(rep1$Isolado, rep1$Conc, media_rep, erro_rep)
colnames(tabgraf) <- c('isolado', 'conc', 'sob', 'erro')
tabgraf$isolado <- factor(tabgraf$isolado, labels=c("D14.15.010-3", "D14.15.100-1", "D14.30.010-1", "D14.30.010-5", "D20.15.010-12", "D20.30.010-2", "D21.30.010-2", "D23.15.010-1", "D23.30.010-1", "D23.30.010-3", "D4.15.R2A-5", "D4.30.100-3", "D4.30.R2A-2", "D4.30.R2A-4", "D5.30.010-1", "D7.15.010-6", "D7.15.R2A-6", "D7.30.100-9", "D8.15.R2A-3", "italic('E. coli')*textstyle(' MG1655')"))
tabgraf$isolado <- factor(tabgraf$isolado)
#********** GRÁFICOS **********
graf_sob<- ggplot(tabgraf, aes(x=conc, y=sob)) +
    geom_point(size=.5)+   # Coloca as bolinhas
    geom_line(linewidth=.3)+   # Coloca as linhas
    geom_errorbar(aes(ymin=sob-erro, ymax=sob+erro), width=0.15, linewidth=.5)+    # Coloca as barras de erro
    facet_wrap(~isolado, ncol=4, labeller=label_parsed)+    # Deixa facetado por espécie
    theme_bw()+   # Muda a aparência para o tema certo
    xlab(bquote("[NaCl]"~(mol.L^-1)))+    # Coloca o nome do eixo x
    ylab(bquote("Sobrevivência"~(N/N[0])))+    # Coloca o nome do eixo y
    theme(
      panel.grid.major = element_line(size=.25),   # Define largura da linha maior do grid interno no gráfico
      panel.grid.minor = element_line(size=.1),   # Define largura da linha menor do grid interno do gráfico
      legend.position='none')+    # Retira a legenda
    scale_y_log10(limits = c(1e-4,20),   # Transforma o eixo y em log10 com limiter de 1 a 10^-5
                  labels = trans_format("log10", math_format(10^.x)),
                   breaks = trans_breaks("log10", function(x) 10^x, n=5),
                   minor_breaks=log10_minor_break())
graf_sob

write.table(tabgraf, file='sob.tsv', sep='\t', row.names=FALSE, dec=',')

ggsave('escolhidos_NaCl.png', graf_sob, device='png', unit='cm', width = 18, height = 26, dpi=300)   # Salva gráfico que acabamos de criar como um arquivo


tab_slide <- filter(tabgraf, isolado=='D14.30.010-1'|
                            isolado=='D14.30.010-5'|
                          isolado=='D20.15.010-12'|
                          isolado=='D20.30.010-2'|
                          isolado=='D23.15.010-1'|
                          isolado=='D23.30.010-1'|
                          isolado=='D4.15.R2A-5'|
                          isolado=='D5.30.010-1'|
                          isolado=='D7.15.010-6'|
                          isolado=='D7.30.100-9')

graf_slide <- ggplot(tab_slide, aes(x=conc, y=sob))+
    geom_point(size=.5)+   # Coloca as bolinhas
    geom_line(linewidth=.3)+   # Coloca as linhas
    geom_hline(yintercept=1, color='darkgrey')+
    geom_errorbar(aes(ymin=sob-erro, ymax=sob+erro), width=0.15, linewidth=.5)+    # Coloca as barras de erro
    facet_wrap(~isolado, ncol=5, labeller=label_parsed)+    # Deixa facetado por espécie
    theme_bw()+   # Muda a aparência para o tema certo
    xlab(bquote("[NaCl]"~(mol.L^-1)))+    # Coloca o nome do eixo x
    ylab(bquote("Sobrevivência"~(N/N[0])))+    # Coloca o nome do eixo y
    theme(
      panel.grid.major = element_line(linewidth=.25),   # Define largura da linha maior do grid interno no gráfico
      panel.grid.minor = element_line(linewidth=.1),   # Define largura da linha menor do grid interno do gráfico
      legend.position='none',
      axis.title = element_text(size=16),
      axis.text = element_text(size=12),
      strip.text=element_text(size=12))+    # Retira a legenda
    scale_y_log10(limits = c(1e-4,20),   # Transforma o eixo y em log10 com limiter de 1 a 10^-5
                  labels = trans_format("log10", math_format(10^.x)),
                   breaks = trans_breaks("log10", function(x) 10^x, n=5),
                   minor_breaks=log10_minor_break())
graf_slide

ggsave('escolhidos_NaCl_slide.png', graf_slide, device='png', unit='cm', width = 30, height = 15, dpi=600)   # Salva gráfico que acabamos de criar como um arquivo
