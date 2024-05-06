# Script para processamento do experimento Irradiação UVC D. radiodurans Líquido vs. Sólido
#REGRAS DE ENTRADA:
#*************ORDEM DAS COLUNAS: 1. Replicata | 2. Experimento | 3. Dose (em J/m2) | 4-6. UFC | 7. Fator de Diluição (deve chamar 'diluicao') **********************************
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
dados <- read.table("Ecoli_UVC.tsv", sep='\t', dec=',', quote='', header=TRUE, fill=TRUE)

#********** PROCESSAMENTO **********
col_indices <- c(3, 4, 5)

mc <- calc_media_corrigida(dados, col_indices)
erro_mc <- calc_erro_medcor(dados, col_indices)

tab_media <- cbind.data.frame(dados, mc, erro_mc)
tab0 <- filter(tab_media, Dose=='0')

media0 <- rep(tab0$mc, each=8)
erro0 <- rep(tab0$erro_mc, each=8)
sob <- mc/media0
sob <- trunc(sob*10^4)/10^4
erro_sob <- sqrt((sob^2) * (((erro_mc / mc)^2) + ((erro0 / media0)^2)))
tab_total <- cbind.data.frame(tab_media, sob, erro_sob)

rep1 <- filter(tab_total, Replicata=='1')
rep2 <- filter(tab_total, Replicata=='2')
rep3 <- filter(tab_total, Replicata=='3')
rep4 <- filter(tab_total, Replicata=='4')

media_rep <- rowMeans(cbind(rep1$sob, rep2$sob, rep3$sob, rep4$sob), na.rm = TRUE)
erro_rep <- sqrt(rowSums(cbind((1/4*rep1$erro_sob)^2, (1/4*rep2$erro_sob)^2, (1/4*rep3$erro_sob)^2, (1/4*rep4$erro_sob)^2), na.rm = TRUE))

tabgraf <- cbind.data.frame(rep1$Dose, media_rep, erro_rep)
colnames(tabgraf) <- c('Dose', 'sob', 'erro')

graf_sob<- ggplot(tabgraf, aes(x=Dose, y=sob)) +
    geom_point(size=.65)+   # Coloca as bolinhas
    geom_line(linewidth=.35)+   # Coloca as linhas
    geom_errorbar(aes(ymin=sob-erro, ymax=sob+erro), width=0.15, size=.5)+    # Coloca as barras de erro
    theme_bw()+   # Muda a aparência para o tema certo
    xlab(bquote("Dose"~(J.m^-2)))+    # Coloca o nome do eixo x
    ylab(bquote("Sobrevivência"~(N/N[0])))+    # Coloca o nome do eixo y
    theme(
      panel.grid.major = element_line(size=.25),   # Define largura da linha maior do grid interno no gráfico
      panel.grid.minor = element_line(size=.1),   # Define largura da linha menor do grid interno do gráfico
      legend.position='none')+    # Retira a legenda
    scale_y_log10(limits = c(1e-4,10),   # Transforma o eixo y em log10 com limiter de 1 a 10^-5
                  labels = trans_format("log10", math_format(10^.x)),
                   breaks = trans_breaks("log10", function(x) 10^x, n=4),
                   minor_breaks=log10_minor_break())
graf_sob

#********** FITTING **********
x <- c(0, 50, 100, 200, 300, 500, 700, 1000)
f <- 0
fit_ecoli <- nls(media_rep ~ (1+f*x)/(1+(9+10*f*LD10)*exp(n*log(x/LD10))),
			        start = list(LD10 = 30, n = 2))
summary(fit_ecoli)
fit_ld1_ecoli <- nls(media_rep ~ (1+f*x)/(1+(99+100*f*LD1)*exp(n*log(x/LD1))),
			        start = list(LD1 = 90, n = 2))
summary(fit_ld1_ecoli)

graf_fit_ecoli<- ggplot(tabgraf, aes(x=Dose, y=sob)) +
    geom_point(shape=15, size=1.5)+   # Coloca as bolinhas
    geom_errorbar(aes(ymin=sob-erro, ymax=sob+erro), width=0.15, size=.5)+    # Coloca as barras de erro
    geom_smooth(method='nls', data=tabgraf,
                formula = y ~ (1+f*x)/(1+(9+10*f*LD10)*exp(n*log(x/LD10))),
                se=FALSE,
                method.args = list(start = list(LD10 = 90, n = 2)), linewidth=.5)+
    theme_bw()+   # Muda a aparência para o tema certo
    geom_hline(yintercept=.1, linetype='dashed', color= 'grey30')+
    xlab(bquote("Dose"~(J.m^-2)))+    # Coloca o nome do eixo x
    ylab(bquote("Sobrevivência"~(N/N[0])))+    # Coloca o nome do eixo y
    ggtitle('E. coli')+
    theme(
      panel.grid.major = element_line(size=.25),   # Define largura da linha maior do grid interno no gráfico
      panel.grid.minor = element_line(size=.1),   # Define largura da linha menor do grid interno do gráfico
      legend.position='none',
      plot.title = element_text(hjust = 0.5, face='italic'))+
    scale_x_continuous(breaks=c(0, 100, 200, 300, 400, 500, 600, 700, 800, 900, 1000))
graf_fit_ecoli
ggsave('fit_ecoli.png', graf_fit_ecoli, device='png', unit='cm', width = 12, height = 10, dpi=300)
