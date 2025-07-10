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
#install.packages('viridis')
#library(viridis)
install.packages('matrixStats')
library(matrixStats)
install.packages('MASS')
library(MASS)
install.packages('scales')
library(scales)
install.packages("tidyverse")
library(tidyverse)
install.packages("ggtext")
library(ggtext)

# ********** FUNÇÕES **********
baranyi <- function(params, x) {
  params[1] + params[2] * (x + (1/params[2]) * log(exp(-params[2]*x) +
  exp(-params[2] * params[3]) - exp(-params[2] * (x + params[3])))) -
  log(1 + ((exp(params[2] * (x + (1/params[2]) * log(exp(-params[2]*x) +
  exp(-params[2] * params[3]) - exp(-params[2] * (x + params[3])))))-1)/
    (exp(params[4]-params[1]))))
}

yi_fitado_baranyi <- function(x, coef){   # x é vetor com valores de tempo, coef é valores dos coeficientes fitados do modelo baranyi
  y_fitado <- c()
  y0 <- coef[1]
  mmax <- coef[2]
  lambda <- coef[3]
  ymax <- coef[4]
  for(n in seq_along(x)){
    y <-  y0 + mmax * (x[n] + (1/mmax) * log(exp(-mmax*x[n]) +
                    exp(-mmax * lambda) - exp(-mmax * (x[n] + lambda)))) -
                    log(1 + ((exp(mmax * (x[n] + (1/mmax) * log(exp(-mmax*x[n]) +
                    exp(-mmax * lambda) - exp(-mmax *
                    (x[n] + lambda)))))-1)/(exp(ymax-y0))))
    y_fitado <- append(y_fitado, y, after=length(y_fitado))
  }
  return(y_fitado)
}

calc_pseudo_R2 <- function(data, teor){   #data é vetor com os dados experimentais e teor é vetor com y teórico pros mesmos valores de x
  n <- length(data)
  y_mean <- rep(mean(data), each=n)
  ss_tot <- sum((data-y_mean)^2)
  ss_res <- sum((data-teor)^2)
  pseudo_R2 <- 1-(ss_res/ss_tot)
  return(pseudo_R2)
}


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
graf

# ********** MODELAGEM **********
# Modelo Baranyi-Roberts (fitado na mão)
y <- log(media)
erro_log <- erro/media
d <- rep1$Tempo
y0 <- min(y)
fit<- nls(y ~ y0 + mmax * (d + (1/mmax) * log(exp(-mmax*d) +
            exp(-mmax * lambda) - exp(-mmax * (d + lambda)))) -
            log(1 + ((exp(mmax * (d + (1/mmax) * log(exp(-mmax*d) +
            exp(-mmax * lambda) - exp(-mmax *
            (d + lambda)))))-1)/(exp(ymax-y0)))),
            start=list(mmax=0.2, lambda=2, ymax=max(y)))
baranyi_param <- coef(fit)
print(baranyi_param)
coef <- c(y0, baranyi_param[1], baranyi_param[2], baranyi_param[3])
y_fit <- baranyi(coef, d)
tabgraf_fit <- cbind.data.frame(d, y, y_fit, erro_log)

graf_fit <- ggplot(tabgraf_fit, aes(x=d))+
  geom_line(aes(y=y_fit), linewidth=.7, color='gray55')+
  geom_point(aes(y=y), size=.93)+
  geom_errorbar(aes(ymin=y-erro_log, ymax=y+erro_log), width=.3)+
  theme_bw()+
  scale_x_continuous(breaks = seq(0, 30, 3))+
  labs(
    x='Time (hours)',
    y='ln(CFU.mL<sup>-1</sup>)')+
    theme(axis.title = element_markdown(size=12),
          axis.text = element_text(size=10),
          panel.grid.major = element_line(size=.25),
          panel.grid.minor = element_line(size=.1))
graf_fit
ggsave(plot=graf_fit, 'boulardii_growthcurve_artigo.tiff', , dpi=600, unit='cm', width=15, height=15)

summary(fit)
