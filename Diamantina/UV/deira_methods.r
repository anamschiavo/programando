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
deira <- read.table("deira.tsv", sep='\t', dec=',', quote='', header=TRUE, fill=TRUE)

#********** PROCESSAMENTO **********
col_indices <- c(4, 5, 6)

mc <- calc_media_corrigida(deira, col_indices)
erro_mc <- calc_erro_medcor(deira, col_indices)

tab_medias <- cbind.data.frame(deira, mc, erro_mc)

rep1 <- filter(tab_medias, Replicata=='2')
rep2 <- filter(tab_medias, Replicata=='3')
rep3 <- filter(tab_medias, Replicata=='4')

media_rep <- rowMeans(cbind(rep1$mc, rep2$mc, rep2$mc), na.rm = TRUE)
erro_rep <- sqrt(rowSums(cbind((1/3*rep1$erro_mc)^2, (1/3*rep2$erro_mc)^2, (1/3*rep1$erro_mc)^2), na.rm = TRUE))

tab <- cbind.data.frame(rep1$Dose, rep1$Exp, media_rep, erro_rep)
tab0 <- filter(tab, rep1$Dose=='0')
media0 <- rep(tab0$media_rep, each=10)
erro0 <- rep(tab0$erro_rep, each=10)

sob <- media_rep/media0
erro_sob <- sqrt((sob^2) * (((erro_rep / media_rep)^2) + ((erro0 / media0)^2)))

tabgraf <- cbind.data.frame(rep1$Dose, rep1$Exp, sob, erro_sob)
colnames(tabgraf) <- c('Dose', 'Exp', 'sob', 'erro')

graf<- ggplot(tabgraf, aes(x=Dose, y=sob, color=Exp, shape=Exp)) +
    geom_point(size=.75)+   # Coloca as bolinhas
    geom_line(linewidth=.25)+   # Coloca as linhas
    geom_errorbar(aes(ymin=sob-erro, ymax=sob+erro), width=0.5, linewidth=.2)+    # Coloca as barras de erro
    theme_bw()+   # Muda a aparência para o tema certo
    xlab(bquote("Dose"~(J.m^-2)))+    # Coloca o nome do eixo x
    ylab(bquote("Sobrevivência"~(N/N[0])))+    # Coloca o nome do eixo y
    scale_color_manual(values=c('springgreen3', 'black'),labels=c('Líquido', 'Sólido'))+
    scale_shape_manual(values=c(16, 15),labels=c('Líquido', 'Sólido'))+
    labs(shape='', color='')+
    theme(
      panel.grid.major = element_line(size=.3),   # Define largura da linha maior do grid interno no gráfico
      panel.grid.minor = element_line(size=.1),   # Define largura da linha menor do grid interno do gráfico
      legend.text = element_text(size=10),
      axis.title = element_text(size=10),
      axis.text = element_text(size=8))+
    scale_y_log10(limits = c(1e-5,10),   # Transforma o eixo y em log10 com limites de 1 a 10^-5
                  labels = trans_format("log10", math_format(10^.x)),
                   breaks = trans_breaks("log10", function(x) 10^x, n=5),
                   minor_breaks=log10_minor_break())+
    scale_x_continuous(breaks=c(0, 500, 1000, 1500, 2000, 2500, 3000))
graf
ggsave('irrad_deira.png', graf, device='png', unit='cm', width = 12, height = 6.5, dpi=300)

#********** FITTING **********
tab_liq <- filter(tabgraf, Exp=='L')
tab_sol <- filter(tabgraf, Exp=='S')
x <- as.numeric(tab_liq$Dose)
y_liq <- tab_liq$sob
y_sol <- tab_sol$sob

fit_liq <- nls(y_liq ~ (1+f*x)/(1+(9+10*f*LD10)*exp(n*log(x/LD10))),
			        start = list(LD10 = 1700, n = 6, f = 0.001))
summary(fit_liq)

fit_sol <- nls(y_sol ~ (1+f*x)/(1+(9+10*f*LD10)*exp(n*log(x/LD10))),
			        start = list(LD10 = 1200, n = 6, f = 0.001))
summary(fit_sol)
