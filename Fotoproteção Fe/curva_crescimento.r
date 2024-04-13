# R 4.2.0
# PARAINICIAÇÃO CIENTÍFICA DE ISMAEL - Grupo Quimiosfera
# Scrip para gráficos e análises de dados do experimento curva de crescimento de S. boulardii

# REGRAS DE ENTRADA: tabela chamada "crescimento.tsv", vírgula para decimal, tab para separação de células
# ORDEM DAS COLUNAS: 1. Ponto | 2. Replicata | 3. Horas de experimento | 4-6. Contagem de UFC | 7. Fator de correção da diluição

# ********** PACOTES E BIBLIOTECAS **********
install.packages('dplyr')
library(dplyr)
install.packages('janitor')
library(janitor)
install.packages('ggplot2')
library(ggplot2)
install.packages('viridis')
library(viridis)
install.packages('matrixStats')
library(matrixStats)
install.packages('MASS')
library(MASS)
install.packages('scales')
library(scales)
install.packages('growthrates')
library(growthrates)
install.packages("tidyverse")
library(tidyverse)
install.packages("ggtext")
library(ggtext)


# ********** LEITURA E PRÉ-PROCESSAMENTO **********
tab <- read.table('crescimento.tsv', sep='\t', dec=',', header=TRUE, fill=TRUE)

#********** PRÉ-PROCESSAMENTO **********
# Faz média corrigida pelo fator de diluição das gotinhas e para mL
mc <- rowMeans(tab[,c(4:6)])*tab$Fator*100

# Faz erro padrão das gotinhas
se <- rowSds(as.matrix(tab[,c(4:6)]))/(sqrt(3))

# Calcula erro associado a média corrigida
ec <- se*tab$Fator*100

# Junta a média corrigida e erro padrão das gotinhas na tabela original
tab <- cbind(tab, mc, ec)

# Separa a tabela em tabelas com apenas uma replicatas
rep1 <- filter(tab, Replicata=='1')
rep2 <- filter(tab, Replicata=='2')
rep3 <- filter(tab, Replicata=='3')

media <- (rep1$mc + rep2$mc + rep3$mc)/3

erros_quadrados <- cbind((1/3*rep1$ec)^2, (1/3*rep2$ec)^2, (1/3*rep3$ec)^2)
erro <- sqrt(rowSums(erros_quadrados, na.rm=TRUE))

tabgraf <- cbind.data.frame(rep1$Tempo, media, erro)
colnames(tabgraf) <- c('tempo', 'media', 'erro')

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

graf <- ggplot(tabgraf, aes(x=tempo, y=media))+
  geom_point(size=1.5)+
  geom_line(linewidth=.7)+
  geom_errorbar(aes(ymin=media-erro, ymax=media+erro), width=.5, size=.5)+
  theme_bw()+
  labs(x='Tempo (h)', y='UFC/mL')+
  scale_y_log10(#limits = c(10, NA),
                 labels = trans_format("log10", math_format(10^.x)),
                 breaks=trans_breaks("log10", function(x) 10^x, n=8),
                 minor_breaks=log10_minor_break())

# ********** MODELAGEM **********
media_log <- log10(media)
erro_log <- erro/(media*log(10))

# Modelo Baranyi-Roberts
p0_baranyi <- c(y0=media_log[1], mumax=.2, K=max(media_log), h0=1)
fit_b <- fit_growthmodel(FUN=grow_baranyi, p=p0_baranyi, time=tabgraf$tempo, y=media_log, transform='log')
param_b <- coef(fit_b)
summary(fit_b)

tab_fit_b <- grow_baranyi(tabgraf$tempo, param_b)
tabgraf_fit_b <- cbind.data.frame(tab_fit_b, media_log, erro_log)
graf_b <- ggplot(tabgraf_fit_b, aes(x=time))+
  geom_point(aes(y=media_log), size=1.5)+
  geom_line(aes(y=y))+
  geom_errorbar(aes(ymin=media_log-erro_log, ymax=media_log+erro_log), width=.5, size=.5)+
  theme_bw()+
  labs(title='Baranyi-Roberts',
    x='Tempo (h)',
    y='log<sub>10</sub>(UFC/mL)')+
  theme(axis.title.y=element_markdown())
graf_b
ggsave(plot=graf_b, 'baranyi.png')

# Modelo Gompertz
p0_gompertz <- c(y0=media_log[1], mumax=.2, K=max(media_log), lambda=1)
fit_g <- fit_growthmodel(FUN=grow_gompertz3, p=p0_gompertz, time=tabgraf$tempo, y=media_log, transform='log')
param_g <- coef(fit_g)
summary(fit_g)

tab_fit_g <- grow_gompertz3(tabgraf$tempo, param_g)
tabgraf_fit_g <- cbind.data.frame(tab_fit_g, media_log, erro_log)
graf_g <- ggplot(tabgraf_fit_g, aes(x=time))+
  geom_point(aes(y=media_log), size=1.5)+
  geom_line(aes(y=y))+
  geom_errorbar(aes(ymin=media_log-erro_log, ymax=media_log+erro_log), width=.5, size=.5)+
  theme_bw()+
  labs(title='Gompertz',
    x='Tempo (h)',
    y='log<sub>10</sub>(UFC/mL)')+
  theme(axis.title.y=element_markdown())
graf_g
ggsave(plot=graf_g, 'gompertz.png')

# Modelo Huang
p0_huang <- c(y0=media_log[1], mumax=.2, K=max(media_log), alpha=5, lambda=1)
fit_h <- fit_growthmodel(FUN=grow_huang, p=p0_huang, time=tabgraf$tempo, y=media_log, transform='log')
param_h <- coef(fit_h)
summary(fit_h)

tab_fit_h <- grow_huang(tabgraf$tempo, param_h)
tabgraf_fit_h <- cbind.data.frame(tab_fit_h, media_log, erro_log)
graf_h <- ggplot(tabgraf_fit_h, aes(x=time))+
  geom_point(aes(y=media_log), size=1.5)+
  geom_line(aes(y=y))+
  geom_errorbar(aes(ymin=media_log-erro_log, ymax=media_log+erro_log), width=.5, size=.5)+
  theme_bw()+
  labs(title='Huang',
    x='Tempo (h)',
    y='log<sub>10</sub>(UFC/mL)')+
  theme(axis.title.y=element_markdown())
graf_h
ggsave(plot=graf_h, 'huang.png')

# Modelo Richards
p0_richards <- c(y0=media_log[1], mumax=.2, K=max(media_log), beta=2)
fit_r <- fit_growthmodel(FUN=grow_richards, p=p0_richards, time=tabgraf$tempo, y=media_log, transform='log')
param_r <- coef(fit_r)
summary(fit_r)

tab_fit_r <- grow_richards(tabgraf$tempo, param_r)
tabgraf_fit_r <- cbind.data.frame(tab_fit_r, media_log, erro_log)
graf_r <- ggplot(tabgraf_fit_r, aes(x=time))+
  geom_point(aes(y=media_log), size=1.5)+
  geom_line(aes(y=y))+
  geom_errorbar(aes(ymin=media_log-erro_log, ymax=media_log+erro_log), width=.5, size=.5)+
  theme_bw()+
  labs(title='Richards',
    x='Tempo (h)',
    y='log<sub>10</sub>(UFC/mL)')+
  theme(axis.title.y=element_markdown())
graf_r
ggsave(plot=graf_r, 'richards.png')

# Modelo logistico
p0_logistico <- c(y0=media_log[1], mumax=.2, K=max(media_log))
fit_l <- fit_growthmodel(FUN=grow_logistic, p=p0_logistico, time=tabgraf$tempo, y=media_log, transform='log')
param_l <- coef(fit_l)
summary(fit_l)

tab_fit_l <- grow_logistic(tabgraf$tempo, param_l)
tabgraf_fit_l <- cbind.data.frame(tab_fit_l, media_log, erro_log)
graf_l <- ggplot(tabgraf_fit_l, aes(x=time))+
  geom_point(aes(y=media_log), size=1.5)+
  geom_line(aes(y=y))+
  geom_errorbar(aes(ymin=media_log-erro_log, ymax=media_log+erro_log), width=.5, size=.5)+
  theme_bw()+
  labs(title='Logistico',
    x='Tempo (h)',
    y='log<sub>10</sub>(UFC/mL)')+
  theme(axis.title.y=element_markdown())
graf_l
ggsave(plot=graf_l, 'logistico.png')
