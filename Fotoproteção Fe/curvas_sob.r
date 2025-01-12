#Script Fotoproteção do Ferro para S. boulardii (e Iniciação Científica de Ismael)

#REGRAS DE ENTRADA:
#************* ORDEM DAS COLUNAS: 1-6: Replicata 1 | 7-12: Replicata 2 | 13-18: Replicata 3 | 19-24: Replicata 4 **********************************
# Ordem das colunas dentro de cada replicata: 1. Dose | 2. Ferro | 3. UFC1 | 4. UFC2 | 5. UFC3 | 6. Fator
# Decimal é dado por vírgula
# Primeira linha é cabeçalho

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
install.packages("RColorBrewer")
library(RColorBrewer)
install.packages("tidyverse")
library(tidyverse)
install.packages("ggtext")
library(ggtext)
install.packages('lamW')
library(lamW)
install.packages('gridExtra')
library(gridExtra)

#********** LEITURA **********
tab <- read.table('jaja.tsv', sep="\t", dec=",", header=TRUE, fill=TRUE)

#********** PRÉ-PROCESSAMENTO **********
# Divide tabela original em tabelas de cada replicata
rep1 <- as.matrix(tab[, c(1:6)])
rep2 <- as.matrix(tab[, c(7:12)])
rep3 <- as.matrix(tab[, c(13:18)])
rep4 <- as.matrix(tab[, c(19:24)])

# Faz média corrigida pelo fator de diluição das gotinhas de cada replicata
mc1 <- rowMeans(rep1[,c(3:5)])*rep1[,6]
mc2 <- rowMeans(rep2[,c(3:5)])*rep2[,6]
mc3 <- rowMeans(rep3[,c(3:5)])*rep3[,6]
mc4 <- rowMeans(rep4[,c(3:5)])*rep4[,6]

# Faz erro padrão das gotinhas
se1 <- rowSds(rep1[,c(3:5)])/(sqrt(3))
se2 <- rowSds(rep2[,c(3:5)])/(sqrt(3))
se3 <- rowSds(rep3[,c(3:5)])/(sqrt(3))
se4 <- rowSds(rep4[,c(3:5)])/(sqrt(3))

# Calcula erro associado à média corrigida de cada replicata
erro1 <- se1*rep1[,6]
erro2 <- se2*rep2[,6]
erro3 <- se3*rep3[,6]
erro4 <- se4*rep4[,6]

# Cria tabelas com média doses, média corrigidas e erro padrão
rep_1 <- cbind.data.frame(rep1[,1:2], mc1, erro1)
colnames(rep_1) <- c('dose', 'ferro', 'media', 'erro')
rep_2 <- cbind.data.frame(rep2[,1:2], mc2, erro2)
colnames(rep_2) <- c('dose', 'ferro', 'media', 'erro')
rep_3 <- cbind.data.frame(rep3[,1:2], mc3, erro3)
colnames(rep_3) <- c('dose', 'ferro', 'media', 'erro')
rep_4 <- cbind.data.frame(rep4[,1:2], mc4, erro4)
colnames(rep_4) <- c('dose', 'ferro', 'media', 'erro')

# Separa tabelas de replicatas por concentração de Fe3+
rep1_0 <- filter(rep_1, ferro=='0')
rep1_5 <- filter(rep_1, ferro=='5')
rep1_10 <- filter(rep_1, ferro=='10')
rep1_20 <- filter(rep_1, ferro=='20')
rep1_30 <- filter(rep_1, ferro=='30')

rep2_0 <- filter(rep_2, ferro=='0')
rep2_5 <- filter(rep_2, ferro=='5')
rep2_10 <- filter(rep_2, ferro=='10')
rep2_20 <- filter(rep_2, ferro=='20')
rep2_30 <- filter(rep_2, ferro=='30')

rep3_0 <- filter(rep_3, ferro=='0')
rep3_5 <- filter(rep_3, ferro=='5')
rep3_10 <- filter(rep_3, ferro=='10')
rep3_20 <- filter(rep_3, ferro=='20')
rep3_30 <- filter(rep_3, ferro=='30')

rep4_0 <- filter(rep_4, ferro=='0')
rep4_5 <- filter(rep_4, ferro=='5')
rep4_10 <- filter(rep_4, ferro=='10')
rep4_20 <- filter(rep_4, ferro=='20')
rep4_30 <- filter(rep_4, ferro=='30')

# Cria vetores de média corrigida das doses 0 de cada concentração e replicata
# Replicata 1
m0_r1_0 <- rep(rep1_0$media[1], times=length(rep1_0$media))
m0_r1_5 <- rep(rep1_5$media[1], times=length(rep1_5$media))
m0_r1_10 <- rep(rep1_10$media[1], times=length(rep1_10$media))
m0_r1_20 <- rep(rep1_20$media[1], times=length(rep1_20$media))
m0_r1_30 <- rep(rep1_30$media[1], times=length(rep1_30$media))
# Replicata 2
m0_r2_0 <- rep(rep2_0$media[1], times=length(rep2_0$media))
m0_r2_5 <- rep(rep2_5$media[1], times=length(rep2_5$media))
m0_r2_10 <- rep(rep2_10$media[1], times=length(rep2_10$media))
m0_r2_20 <- rep(rep2_20$media[1], times=length(rep2_20$media))
m0_r2_30 <- rep(rep2_30$media[1], times=length(rep2_30$media))
# Replicata 3
m0_r3_0 <- rep(rep3_0$media[1], times=length(rep3_0$media))
m0_r3_5 <- rep(rep3_5$media[1], times=length(rep3_5$media))
m0_r3_10 <- rep(rep3_10$media[1], times=length(rep3_10$media))
m0_r3_20 <- rep(rep3_20$media[1], times=length(rep3_20$media))
m0_r3_30 <- rep(rep3_30$media[1], times=length(rep3_30$media))
# Replicata 4
m0_r4_0 <- rep(rep4_0$media[1], times=length(rep4_0$media))
m0_r4_5 <- rep(rep4_5$media[1], times=length(rep4_5$media))
m0_r4_10 <- rep(rep4_10$media[1], times=length(rep4_10$media))
m0_r4_20 <- rep(rep4_20$media[1], times=length(rep4_20$media))
m0_r4_30 <- rep(rep4_30$media[1], times=length(rep4_30$media))

# Cria vetores de erros das doses 0 de cada concentração e replicata
# Replicata 1
e0_r1_0 <- rep(rep1_0$erro[1], times=length(rep1_0$erro))
e0_r1_5 <- rep(rep1_5$erro[1], times=length(rep1_5$erro))
e0_r1_10 <- rep(rep1_10$erro[1], times=length(rep1_10$erro))
e0_r1_20 <- rep(rep1_20$erro[1], times=length(rep1_20$erro))
e0_r1_30 <- rep(rep1_30$erro[1], times=length(rep1_30$erro))
# Replicata 2
e0_r2_0 <- rep(rep2_0$erro[1], times=length(rep2_0$erro))
e0_r2_5 <- rep(rep2_5$erro[1], times=length(rep2_5$erro))
e0_r2_10 <- rep(rep2_10$erro[1], times=length(rep2_10$erro))
e0_r2_20 <- rep(rep2_20$erro[1], times=length(rep2_20$erro))
e0_r2_30 <- rep(rep2_30$erro[1], times=length(rep2_30$erro))
# Replicata 3
e0_r3_0 <- rep(rep3_0$erro[1], times=length(rep3_0$erro))
e0_r3_5 <- rep(rep3_5$erro[1], times=length(rep3_5$erro))
e0_r3_10 <- rep(rep3_10$erro[1], times=length(rep3_10$erro))
e0_r3_20 <- rep(rep3_20$erro[1], times=length(rep3_20$erro))
e0_r3_30 <- rep(rep3_30$erro[1], times=length(rep3_30$erro))
# Replicata 4
e0_r4_0 <- rep(rep4_0$erro[1], times=length(rep4_0$erro))
e0_r4_5 <- rep(rep4_5$erro[1], times=length(rep4_5$erro))
e0_r4_10 <- rep(rep4_10$erro[1], times=length(rep4_10$erro))
e0_r4_20 <- rep(rep4_20$erro[1], times=length(rep4_20$erro))
e0_r4_30 <- rep(rep4_30$erro[1], times=length(rep4_30$erro))

# Calcula sobrevivência por replicata e faz média das replicatas
sob0 <- cbind(rep1_0$media/m0_r1_0, rep2_0$media/m0_r2_0, rep3_0$media/m0_r3_0, rep4_0$media/m0_r4_0)
m_sob0 <- rowMeans(sob0, na.rm=TRUE)
sob5 <- cbind(rep1_5$media/m0_r1_5, rep2_5$media/m0_r2_5, rep3_5$media/m0_r3_5, rep4_5$media/m0_r4_5)
m_sob5 <- rowMeans(sob5, na.rm=TRUE)
sob10 <- cbind(rep1_10$media/m0_r1_10, rep2_10$media/m0_r2_10, rep3_10$media/m0_r3_10, rep4_10$media/m0_r4_10)
m_sob10 <- rowMeans(sob10, na.rm=TRUE)
sob20 <- cbind(rep1_20$media/m0_r1_20, rep2_20$media/m0_r2_20, rep3_20$media/m0_r3_20, rep4_20$media/m0_r4_20)
m_sob20 <- rowMeans(sob20, na.rm=TRUE)
sob30 <- cbind(rep1_30$media/m0_r1_30, rep2_30$media/m0_r2_30, rep3_30$media/m0_r3_30, rep4_30$media/m0_r4_30)
m_sob30 <- rowMeans(sob30, na.rm=TRUE)

# Calcula erros associados a sobreviência por replicata e a média das sobrevivências
# Concentração 0 mM
eq_sob0_1 <- ((rep1_0$media/m0_r1_0)^2)*(((rep1_0$erro/rep1_0$media)^2)+((e0_r1_0/m0_r1_0)^2))
eq_sob0_2 <- ((rep2_0$media/m0_r2_0)^2)*(((rep2_0$erro/rep2_0$media)^2)+((e0_r2_0/m0_r2_0)^2))
eq_sob0_3 <- ((rep3_0$media/m0_r3_0)^2)*(((rep3_0$erro/rep3_0$media)^2)+((e0_r3_0/m0_r3_0)^2))
eq_sob0_4 <- ((rep4_0$media/m0_r4_0)^2)*(((rep4_0$erro/rep4_0$media)^2)+((e0_r4_0/m0_r4_0)^2))
eq_sob0 <- cbind(1/16*eq_sob0_1, 1/16*eq_sob0_2, 1/16*eq_sob0_3, 1/16*eq_sob0_4)
e_sob0 <- sqrt(rowSums(eq_sob0, na.rm=TRUE))
# Concentração 5 mM
eq_sob5_1 <- ((rep1_5$media/m0_r1_5)^2)*(((rep1_5$erro/rep1_5$media)^2)+((e0_r1_5/m0_r1_5)^2))
eq_sob5_2 <- ((rep2_5$media/m0_r2_5)^2)*(((rep2_5$erro/rep2_5$media)^2)+((e0_r2_5/m0_r2_5)^2))
eq_sob5_3 <- ((rep3_5$media/m0_r3_5)^2)*(((rep3_5$erro/rep3_5$media)^2)+((e0_r3_5/m0_r3_5)^2))
eq_sob5_4 <- ((rep4_5$media/m0_r4_5)^2)*(((rep4_5$erro/rep4_5$media)^2)+((e0_r4_5/m0_r4_5)^2))
eq_sob5 <- cbind(1/16*eq_sob5_1, 1/16*eq_sob5_2, 1/16*eq_sob5_3, 1/16*eq_sob5_4)
e_sob5 <- sqrt(rowSums(eq_sob5, na.rm=TRUE))
# Concentração 10 mM
eq_sob10_1 <- ((rep1_10$media/m0_r1_10)^2)*(((rep1_10$erro/rep1_10$media)^2)+((e0_r1_10/m0_r1_10)^2))
eq_sob10_2 <- ((rep2_10$media/m0_r2_10)^2)*(((rep2_10$erro/rep2_10$media)^2)+((e0_r2_10/m0_r2_10)^2))
eq_sob10_3 <- ((rep3_10$media/m0_r3_10)^2)*(((rep3_10$erro/rep3_10$media)^2)+((e0_r3_10/m0_r3_10)^2))
eq_sob10_4 <- ((rep4_10$media/m0_r4_10)^2)*(((rep4_10$erro/rep4_10$media)^2)+((e0_r4_10/m0_r4_10)^2))
eq_sob10 <- cbind(1/16*eq_sob10_1, 1/16*eq_sob10_2, 1/16*eq_sob10_3, 1/16*eq_sob10_4)
e_sob10 <- sqrt(rowSums(eq_sob10, na.rm=TRUE))
# Concentração 20 mM
eq_sob20_1 <- ((rep1_20$media/m0_r1_20)^2)*(((rep1_20$erro/rep1_20$media)^2)+((e0_r1_20/m0_r1_20)^2))
eq_sob20_2 <- ((rep2_20$media/m0_r2_20)^2)*(((rep2_20$erro/rep2_20$media)^2)+((e0_r2_20/m0_r2_20)^2))
eq_sob20_3 <- ((rep3_20$media/m0_r3_20)^2)*(((rep3_20$erro/rep3_20$media)^2)+((e0_r3_20/m0_r3_20)^2))
eq_sob20_4 <- ((rep4_20$media/m0_r4_20)^2)*(((rep4_20$erro/rep4_20$media)^2)+((e0_r4_20/m0_r4_20)^2))
eq_sob20 <- cbind(1/16*eq_sob20_1, 1/16*eq_sob20_2, 1/16*eq_sob20_3, 1/16*eq_sob20_4)
e_sob20 <- sqrt(rowSums(eq_sob20, na.rm=TRUE))
# Concentração 30 mM
eq_sob30_1 <- ((rep1_30$media/m0_r1_30)^2)*(((rep1_30$erro/rep1_30$media)^2)+((e0_r1_30/m0_r1_30)^2))
eq_sob30_2 <- ((rep2_30$media/m0_r2_30)^2)*(((rep2_30$erro/rep2_30$media)^2)+((e0_r2_30/m0_r2_30)^2))
eq_sob30_3 <- ((rep3_30$media/m0_r3_30)^2)*(((rep3_30$erro/rep3_30$media)^2)+((e0_r3_30/m0_r3_30)^2))
eq_sob30_4 <- ((rep4_30$media/m0_r4_30)^2)*(((rep4_30$erro/rep4_30$media)^2)+((e0_r4_30/m0_r4_30)^2))
eq_sob30 <- cbind(1/16*eq_sob30_1, 1/16*eq_sob30_2, 1/16*eq_sob30_3, 1/16*eq_sob30_4)
e_sob30 <- sqrt(rowSums(eq_sob30, na.rm=TRUE))


# Calcula média das doses e seus desvios padrões
doses <- as.matrix(tab[,c(1,7,13,19)])
dose <- rowMeans(doses)
erro_dose <- rowSds(doses, na.rm=TRUE)

# Cria tabela para gráficos
erro <- c(e_sob0, e_sob5, e_sob10, e_sob20, e_sob30)
sob <- c(m_sob0, m_sob5, m_sob10, m_sob20, m_sob30)

    # Tabela pra gráfico de 30mM replicatas separadas
sob_30 <- c(rep1_30$media/m0_r1_30, rep2_30$media/m0_r2_30, rep3_30$media/m0_r3_30, rep4_30$media/m0_r4_30)
e_sob_30 <- c(sqrt(eq_sob30_1), sqrt(eq_sob30_2), sqrt(eq_sob30_3), sqrt(eq_sob30_4))
replicata <- rep(c(1, 2, 3, 4), each=length(eq_sob30_1))
dose30 <- c(rep1_30$dose, rep2_30$dose, rep3_30$dose, rep4_30$dose)
tab30 <- cbind.data.frame(dose30, sob_30, e_sob_30, replicata)
tab30$replicata <- as.factor(tab30$replicata)
graf30 <- ggplot(tab30, aes(x=dose30, y=sob_30, color=replicata))+
    geom_point(size=1.5)+
    geom_line(linewidth=.7)+
    geom_errorbar(aes(ymin=sob_30-e_sob_30, ymax=sob_30+e_sob_30), width=.15, size=.5)+
    theme_bw()+
    labs(x=bquote("Fluência"~(J/m^2)),
         y='Sobrevivência (N/N<sub>0</sub>)')+
    scale_colour_brewer(palette='Dark2', name='Replicata')+
    theme(
      panel.grid.major = element_line(size=.5),
      panel.grid.minor = element_line(size=.2),
      legend.title = element_markdown(),
      axis.title.y=element_markdown())+
    scale_y_log10(limits = c(1e-5,2),
                  labels = trans_format("log10", math_format(10^.x)),
                  breaks=trans_breaks("log10", function(x) 10^x, n=5),
                  minor_breaks=log10_minor_break())
graf30
ggsave(plot=graf30, 'S_boulardii_30_replicatas.png')




tabgraf <- cbind.data.frame(dose, tab[,2], sob, erro, erro_dose)
colnames(tabgraf) <- c('dose', 'ferro', 'sob', 'erro_sob', 'erro_dose' )
tabgraf$ferro <- as.factor(tabgraf$ferro)

#********** GRÁFICOS **********
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


junto <- ggplot(tabgraf, aes(x=dose, y=sob, color=ferro))+
    geom_point(size=1.5)+
    geom_line(linewidth=.7)+
    geom_errorbar(aes(ymin=sob-erro_sob, ymax=sob+erro_sob, xmin=dose-erro_dose, xmax=dose+erro_dose), width=.15, size=.5)+
    theme_bw()+
    labs(x=bquote("Fluência"~(J/m^2)),
         y='Sobrevivência (N/N<sub>0</sub>)')+
    scale_colour_brewer(palette='Dark2', name='[Fe<sup>3+</sup>] (mmol.L<sup>-1</sup>)')+
    theme(
      panel.grid.major = element_line(size=.5),
      panel.grid.minor = element_line(size=.2),
      legend.title = element_markdown(),
      axis.title.y=element_markdown())+
    scale_y_log10(limits = c(1e-5,2),
                  labels = trans_format("log10", math_format(10^.x)),
                  breaks=trans_breaks("log10", function(x) 10^x, n=5),
                  minor_breaks=log10_minor_break())
junto
ggsave(plot=junto, 'S_boulardii_junto_jaja.png')

junto_bw <- ggplot(tabgraf, aes(x=dose, y=sob, color=ferro, shape=ferro))+
    geom_point(size=2.5)+
    geom_line(linewidth=.7)+
    geom_errorbar(aes(ymin=sob-erro_sob, ymax=sob+erro_sob, xmin=dose-erro_dose, xmax=dose+erro_dose), width=.15, size=.5)+
    theme_bw()+
    labs(x='Fluence (J.m<sup>-2</sup>)',
         y='Viability (N/N<sub>0</sub>)')+
    scale_color_grey(start=0, end=.75, name='[Fe<sup>3+</sup>] (mmol.L<sup>-1</sup>)')+
    scale_shape_manual(values=c(15, 11, 17, 4, 16), name='[Fe<sup>3+</sup>] (mmol.L<sup>-1</sup>)')+
    theme(
      panel.grid.major = element_line(size=.5),
      panel.grid.minor = element_line(size=.2),
      legend.title = element_markdown(size=12),
      legend.text = element_text(size=10),
      axis.title=element_markdown(size=12),
      axis.text=element_text(size=10))+
    scale_y_log10(limits = c(1e-5,2),
                  labels = trans_format("log10", math_format(10^.x)),
                  breaks=trans_breaks("log10", function(x) 10^x, n=5),
                  minor_breaks=log10_minor_break())
junto_bw
ggsave(plot=junto_bw, 'S_boulardii_junto_artigo_gabriel.tiff', device='tiff', dpi=600, unit='cm', width=20, height=12)

junto_semlog <- ggplot(tabgraf, aes(x=dose, y=sob, color=ferro))+
    geom_point(size=1.5)+
    geom_line(linewidth=.7)+
    geom_errorbar(aes(ymin=sob-erro_sob, ymax=sob+erro_sob, xmin=dose-erro_dose, xmax=dose+erro_dose), width=.15, size=.5)+
    theme_bw()+
    labs(x=bquote("Fluência"~(J/m^2)),
         y='Sobrevivência (N/N<sub>0</sub>)')+
    scale_colour_brewer(palette='Dark2', name='[Fe<sup>3+</sup>] (mmol.L<sup>-1</sup>)')+
    theme(
      panel.grid.major = element_line(size=.5),
      panel.grid.minor = element_line(size=.2),
      legend.title = element_markdown(),
      axis.title.y=element_markdown())
junto_semlog
#ggsave(plot=junto, 'S_boulardii_junto_jaja_semlog.png')

tab0 <- filter(tabgraf, ferro=='0')
graf0 <- ggplot(tab0, aes(x=dose, y=sob))+
    geom_point(size=1.5)+
    geom_line(linewidth=.7)+
    geom_errorbar(aes(ymin=sob-erro_sob, ymax=sob+erro_sob), width=.15, size=.5)+
    theme_bw()+
    xlab(bquote("Fluência"~(J/m^2)))+
    ylab('Sobrevivência (N/N0)')+
    theme(
      panel.grid.major = element_line(size=.5),
      panel.grid.minor = element_line(size=.2))+
    scale_y_log10(limits = c(1e-5,5),
                  labels = trans_format("log10", math_format(10^.x)),
                   breaks=trans_breaks("log10", function(x) 10^x, n=5))
graf0
ggsave(plot=graf0, 'S_boulardii_controle_jaja.png')


#********** FITTING **********
f <- 0
# Concentração 0mM
x0 <- tab0$dose
fit0_ld10 <- nls(m_sob0 ~ (1+f*x0)/(1+(9+10*f*LD10)*exp(n*log(x0/LD10))),
			        start = list(LD10 = 200, n = 2, f=.001))
summary(fit0_ld10)

fit0_ld50 <- nls(m_sob0 ~ (1+f*x0)/(1+(2*f*LD50+1)*(x0/LD50)^n),
			start = list(LD50 = 100, n = 2, f=.001))
summary(fit0_ld50)

tabgraf0 <- cbind.data.frame(x0, m_sob0, e_sob0)
graf_fit0 <- ggplot(tabgraf0, aes(x=x0, y=m_sob0))+
    geom_smooth(method='nls',
                data = tabgraf0,
                formula = y ~ (1+f*x)/(1+(9+10*f*LD10)*exp(n*log(x/LD10))),
                se = FALSE,
                method.args = list(start = list(LD10 = 200, n = 2, f=.001)),
                color='grey')+
    geom_point()+
    geom_errorbar(aes(ymin=m_sob0-e_sob0, ymax=m_sob0+e_sob0), size=.5)+
    ggtitle('0 mM')+
    theme_bw()+
    xlab("Fluence (J/m<sup>2</sup>)")+
    ylab('Viability (N/N<sub>0</sub>)')+
    theme(axis.title = element_markdown(size=12),
          axis.text = element_text(size=10))
graf_fit0
ggsave(plot=graf_fit0, 'fit0_artigo.tiff', dpi=600, unit='cm', width=15, height=15)

# Concentração 5mM
tab5 <- filter(tabgraf, ferro=='5')
x5 <- tab5$dose
fit5_ld10 <- nls(m_sob5 ~ (1+f*x5)/(1+(9+10*f*LD10)*exp(n*log(x5/LD10))),
			        start = list(LD10 = 200, n = 2))
summary(fit5_ld10)

fit5_ld50 <- nls(m_sob5 ~ (1+f*x5)/(1+(2*f*LD50+1)*(x5/LD50)^n),
			start = list(LD50 = 200, n = 2))
summary(fit5_ld50)

tabgraf5 <- cbind.data.frame(x5, m_sob5)
graf_fit5 <- ggplot(tabgraf5, aes(x=x5, y=m_sob5))+
    geom_smooth(method='nls',
                data = tabgraf5,
                formula = y ~ (1+f*x)/(1+(9+10*f*LD10)*exp(n*log(x/LD10))),
                se = FALSE,
                method.args = list(start = list(LD10 = 200, n = 2)),
                color='grey')+
    geom_point()+
    geom_errorbar(aes(ymin=m_sob5-e_sob5, ymax=m_sob5+e_sob5), width=100)+
    ggtitle('5 mM')+
    theme_bw()+
    xlab("Fluence (J/m<sup>2</sup>)")+
    ylab('Viability (N/N<sub>0</sub>)')+
    theme(axis.title = element_markdown(size=12),
        axis.text = element_text(size=10))
graf_fit5
ggsave(plot=graf_fit5, 'fit5_artigo.tiff', dpi=600, unit='cm', width=15, height=15)

# Concentração 10mM
tab10 <- filter(tabgraf, ferro=='10')
x10 <- tab10$dose
fit10_ld10 <- nls(m_sob10 ~ (1+f*x10)/(1+(9+10*f*LD10)*exp(n*log(x10/LD10))),
			        start = list(LD10 = 700, n = 2))
summary(fit10_ld10)

fit10_ld50 <- nls(m_sob10 ~ (1+f*x10)/(1+(2*f*LD50+1)*(x10/LD50)^n),
			start = list(LD50 = 300, n = 2))
summary(fit10_ld50)

tabgraf10 <- cbind.data.frame(x10, m_sob10)
graf_fit10 <- ggplot(tabgraf10, aes(x=x10, y=m_sob10))+
    geom_smooth(method='nls',
                data = tabgraf10,
                formula = y ~ (1+f*x)/(1+(9+10*f*LD10)*exp(n*log(x/LD10))),
                se = FALSE,
                method.args = list(start = list(LD10 = 700, n = 1)),
                color='grey')+
    geom_point()+
    geom_errorbar(aes(ymin=m_sob10-e_sob10, ymax=m_sob10+e_sob10), width=200)+
    ggtitle('10 mM')+
    theme_bw()+
    xlab("Fluence (J/m<sup>2</sup>)")+
    ylab('Viability (N/N<sub>0</sub>)')+
    theme(axis.title = element_markdown(size=12),
          axis.text = element_text(size=10))
graf_fit10
ggsave(plot=graf_fit10, 'fit10_artigo.tiff', dpi=600, unit='cm', width=15, height=15)

# Concentração 20mM
tab20 <- filter(tabgraf, ferro=='20')
x20 <- tab20$dose
fit20_ld10 <- nls(m_sob20 ~ (1+f*x20)/(1+(9+10*f*LD10)*exp(n*log(x20/LD10))),
			        start = list(LD10 = 7000, n = 1))
summary(fit20_ld10)

fit20_ld50 <- nls(m_sob20 ~ (1+f*x20)/(1+(2*f*LD50+1)*(x20/LD50)^n),
			start = list(LD50 = 5000, n = 1))
summary(fit20_ld50)

tabgraf20 <- cbind.data.frame(x20, m_sob20)
graf_fit20 <- ggplot(tabgraf20, aes(x=x20, y=m_sob20))+
    geom_smooth(method='nls',
                data = tabgraf20,
                formula = y ~ (1+f*x)/(1+(9+10*f*LD10)*exp(n*log(x/LD10))),
                se = FALSE,
                method.args = list(start = list(LD10 = 7000, n = 1)),
                color='grey')+
    geom_point()+
    geom_errorbar(aes(ymin=m_sob20-e_sob20, ymax=m_sob20+e_sob20), width=400)+
    ggtitle('20 mM')+
    theme_bw()+
    xlab("Fluence (J/m<sup>2</sup>)")+
    ylab('Viability (N/N<sub>0</sub>)')+
    theme(axis.title = element_markdown(size=12),
          axis.text = element_text(size=10))
graf_fit20
ggsave(plot=graf_fit20, 'fit20_artigo.tiff', dpi=600, unit='cm', width=15, height=15)

# Concentração 30mM
tab30 <- filter(tabgraf, ferro=='30')
x30 <- tab30$dose
fit30_ld10 <- nls(m_sob30 ~ (1+f*x30)/(1+(9+10*f*LD10)*exp(n*log(x30/LD10))),
			        start = list(LD10 = 15000, n = 1))
summary(fit30_ld10)
coef(summary(fit30_ld10))

fit30_ld50 <- nls(m_sob30 ~ (1+f*x30)/(1+(2*f*LD50+1)*(x30/LD50)^n),
			start = list(LD50 = 10000, n = 1))
summary(fit30_ld50)
coef(summary(fit30_ld50))

tabgraf30 <- cbind.data.frame(x30, m_sob30)
graf_fit30 <- ggplot(tabgraf30, aes(x=x30, y=m_sob30))+
    geom_smooth(method='nls',
                data = tabgraf30,
                formula = y ~ (1+f*x)/(1+(9+10*f*LD10)*exp(n*log(x/LD10))),
                se = FALSE,
                method.args = list(start = list(LD10 = 15000, n = 1)),
                color='grey')+
    geom_point()+
    geom_errorbar(aes(ymin=m_sob30-e_sob30, ymax=m_sob30+e_sob30), width=900)+
    ggtitle('30 mM')+
    theme_bw()+
    xlab("Fluence (J/m<sup>2</sup>)")+
    ylab('Viability (N/N<sub>0</sub>)')+
    theme(axis.title = element_markdown(size=12),
          axis.text = element_text(size=10))
graf_fit30
ggsave(plot=graf_fit30, 'fit30_artigo.tiff', dpi=600, unit='cm', width=15, height=15)


#********** LD50 TEÓRICO/EXPERIMENTAL **********
# Gráfico
fun.ld50 <- function(c) ld_i*epslon*c*l*log(10)/(1-10^(-epslon*l*c))

epslon <- 2.4
l <- 0.207876
ld_i <- summary(fit0_ld50)$parameters[1,1]

ld50_medido <- c(summary(fit0_ld50)$parameters[1,1], summary(fit5_ld50)$parameters[1,1], summary(fit10_ld50)$parameters[1,1], summary(fit20_ld50)$parameters[1,1], summary(fit30_ld50)$parameters[1,1])
ld50_se <- c(summary(fit0_ld50)$parameters[1,2], summary(fit5_ld50)$parameters[1,2], summary(fit10_ld50)$parameters[1,2], summary(fit20_ld50)$parameters[1,2], summary(fit30_ld50)$parameters[1,2])
Fe <- c(0, 5, 10, 20, 30)
graf_ld50.Fe <- cbind.data.frame(ld50_medido, ld50_se, Fe)

ld50_Fe <- ggplot(graf_ld50.Fe, aes(y=ld50_medido, x=Fe))+
    geom_function(fun = fun.ld50, aes(color='darkgrey'), linewidth=.75, xlim=c(0,30), show.legend=TRUE)+
    geom_point(size=1.5, show.legend=TRUE, aes(color='black'))+
    geom_errorbar(aes(ymin=ld50_medido-ld50_se, ymax=ld50_medido+ld50_se), width=0.5, size=.5)+
    theme_bw()+
  #  xlab('[Fe<sup>3+</sup>] (mmol.L<sup> -1</sup>)')+
  #  ylab('LD<sub>50</sub> (J.m<sup> -2</sup>)')+
    labs(x='[Fe<sup>3+</sup>] (mmol.L<sup> -1</sup>)',
         y='LD<sub>50</sub> (J.m<sup> -2</sup>)', colour='')+
    theme(
      panel.grid.major = element_line(size=.5),
      panel.grid.minor = element_line(size=.2),
      legend.text = element_markdown(size=12),
      legend.position='right',
      axis.title=element_markdown(size=14),
      axis.text=element_text(size=12))+
    scale_colour_manual(labels=c('Experimental', 'Theoretical'), values=c('black', 'darkgrey'))
ld50_Fe
ggsave(plot=ld50_Fe, 'LD50_teoricoXmedido.tiff', device='tiff', dpi=600, unit='cm', width=20, height=16)

# Kolmogorov-Smirnov
ld50_teor <- c(fun.ld50(5), fun.ld50(10), fun.ld50(20), fun.ld50(30))
ks.test(ld50_teor, ld50_medido[-1])

#********** PROFUNDIDADE X CONCENTRAÇÃO Fe3+ **********
# Função isolando profundidade (resolvido por Wolframalpha Algebra)
profun_marte <- function(conc, mu, LD50){
  profundidade <- (0.00180956*LD50*mu*lambertW0(-(0.0156651*10^(-0.00680328/(LD50*mu)))/(LD50*mu)) + 0.000028347)/(conc*LD50*mu)
  return(profundidade)
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

#Gráficos
conc <- rep(c(0.001, 5, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100), times=6)
bac <- rep(c('af', 'sb'), each=36)
mu <- rep(c(0.19/3600, 0.19/(2*3600), 0.19/36000, 0.5656/3600, 0.5656/(2*3600), 0.5656/36000), each=12)
ld50.bac <- rep(c(17.6, 151.11), each=36)

prof <- profun_marte(conc, mu, ld50.bac)
tabela <- cbind.data.frame(conc, bac, mu, prof)

prof_af <- ggplot(tabela, aes(x=conc))+
    geom_function(fun = profun_marte,
        args = list(mu=0.19/3600, LD50=9), aes(color='mu'), linewidth=.7)+
    geom_function(fun = profun_marte,
        args = list(mu=0.19/(2*3600), LD50=9), aes(color='mu/2'), linewidth=.7, linetype='longdash')+
    geom_function(fun = profun_marte,
        args = list(mu=0.19/36000, LD50=9), aes(color='mu/10'), linewidth=.7, linetype='dotdash')+
    theme_bw()+
    ggtitle(expression(paste(bold('B - '), italic('A. ferrooxidans'))))+
    labs(x='[Fe<sup>3+</sup>] (mmol.L<sup> -1</sup>)',
         y='Depth (m)', colour='')+
    scale_colour_manual(labels=c(expression(paste(mu)), expression(paste(mu,'/',10)), expression(paste(mu,'/',2))), values=c('gray50', 'gray70', 'black'))+
    theme(
      plot.title = element_text(size=15),
      panel.grid.major = element_line(size=.5),
      panel.grid.minor = element_line(size=.1),
      legend.text = element_text(size=13),
      legend.position='right',
      axis.title=element_markdown(size=13),
      axis.text=element_text(size=12))+
      scale_x_log10(limits=c(0.001, 100),
                    labels = trans_format("log10", math_format(10^.x)),
                    breaks=trans_breaks("log10", function(x) 10^x, n=5),
                    minor_breaks=log10_minor_break())+
      scale_y_log10(labels = trans_format("log10", math_format(10^.x)),
                    breaks=trans_breaks("log10", function(x) 10^x, n=5),
                    minor_breaks=log10_minor_break())

prof_af
ggsave(plot=prof_af, 'profundidade_af.tiff', device='tiff', dpi=600, unit='cm', width=20, height=16)




prof_sb <- ggplot(tabela, aes(x=conc))+
#    geom_function(fun = profun_marte,
#        args = list(mu=0.5656/3600, LD50=151.11), aes(color='mu'), linewidth=.7)+
    geom_function(fun = profun_marte,
        args = list(mu=0.5656/(2*3600), LD50=151.11), aes(color='mu/2'), linewidth=.7)+
    geom_function(fun = profun_marte,
        args = list(mu=0.5656/36000, LD50=151.11), aes(color='mu/10'), linewidth=.7, linetype='dotdash')+
    theme_bw()+
    labs(x='[Fe<sup>3+</sup>] (mmol.L<sup> -1</sup>)',
         y='Depth (m)', colour='')+
    scale_colour_manual(labels=c(expression(paste(mu,'/',10)), expression(paste(mu,'/',2))), values=c('black', 'gray70'))+
    ggtitle(expression(paste(bold('A - '), italic('S. boulardii'))))+
    theme(
      plot.title = element_text(size=15),
      panel.grid.major = element_line(size=.5),
      panel.grid.minor = element_line(size=.1),
      legend.text = element_text(size=13),
      legend.position='right',
      axis.title=element_markdown(size=13),
      axis.text=element_text(size=12))+
      scale_x_log10(limits=c(0.001, 100),
                    labels = trans_format("log10", math_format(10^.x)),
                    breaks=trans_breaks("log10", function(x) 10^x, n=5),
                    minor_breaks=log10_minor_break())+
      scale_y_log10(labels = trans_format("log10", math_format(10^.x)),
                    breaks=trans_breaks("log10", function(x) 10^x, n=5),
                    minor_breaks=log10_minor_break())

prof_sb
ggsave(plot=prof_sb, 'profundidade_sb.tiff', device='tiff', dpi=600, unit='cm', width=20, height=16)

profundidade <- grid.arrange(prof_sb, prof_af, nrow=2)
ggsave(plot=profundidade, 'profundidade.tiff', device='tiff', dpi=600, unit='cm', width=20, height=30)
