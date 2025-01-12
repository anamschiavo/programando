# R 4.3.3
# PARA ARTIGO DOS 19 ISOLADOS PERCLORATO ROBERTA ALMEIDA VINCENZI - Grupo Quimiosfera

# REGRAS DE ENTRADA: tabela chamada "dados_media.xlsx"
# ORDEM DAS COLUNAS: 1. Isolado | 2. Concentração | 3. Tempo (min) | 4. Média da DO | 5. Desvio Padrão

# ********** PACOTES E BIBLIOTECAS **********
install.packages('dplyr')
library(dplyr)
install.packages('janitor')
library(janitor)
install.packages('ggplot2')
library(ggplot2)
install.packages('matrixStats')
library(matrixStats)
install.packages('MASS')
library(MASS)
install.packages('scales')
library(scales)
install.packages("gridExtra")
library(gridExtra)
install.packages("minpack.lm")
library(minpack.lm)
install.packages('writexl')
library(writexl)

# ********** FUNÇÕES *********
gompertz <- function(params, x) {
  params[1] + (params[3] * exp(-exp(-params[2] * (x - params[4]))))
}

# ********** LEITURA **********
tab <- read.table('dados_media.tsv', sep='\t', dec=',', header=TRUE, fill=TRUE)
names <- unique(tab$Isolado)
conc <- unique(tab$Concentracao)

# ********** FITTING **********
mmax <- c()
mmax_sd <- c()
C <- c()
M <- c()
M_sd <- c()
Y0 <- c()

# Loop para fit de 0 mol/L
for(n in seq_along(names)){   #passa por cada bicho
      tabfit <- filter(tab, Isolado==names[n] & Concentracao==0)
      d <- tabfit$Time
      y <- tabfit$Media_DO
      y0 <- min(y)
      fit <- nlsLM(y ~ y0 + C * exp(-exp(-mmax * (d - M))),
                  start=list(mmax=0.00001, C=(max(y)-min(y)), M=200))
      Y0 <- append(Y0, y0, after=length(Y0))
      mmax <- append(mmax, summary(fit)$parameters[1,1], after=length(mmax))
      mmax_sd <- append(mmax_sd, summary(fit)$parameters[1,2], after=length(mmax_sd))
      C <- append(C, summary(fit)$parameters[2,1], after=length(C))
      M <- append(M, summary(fit)$parameters[3,1], after=length(M))
      M_sd <- append(M_sd, summary(fit)$parameters[3,2], after=length(M_sd))
    }

# Loop para fit de 0.1 mol/L
for(n in seq_along(names)){   #passa por cada bicho
      tabfit <- filter(tab, Isolado==names[n] & Concentracao==0.1)
      d <- tabfit$Time
      y <- tabfit$Media_DO
      y0 <- min(y)
      fit <- nlsLM(y ~ y0 + C * exp(-exp(-mmax * (d - M))),
                  start=list(mmax=0.00001, C=(max(y)-min(y)), M=500))
      Y0 <- append(Y0, y0, after=length(Y0))
      mmax <- append(mmax, summary(fit)$parameters[1,1], after=length(mmax))
      mmax_sd <- append(mmax_sd, summary(fit)$parameters[1,2], after=length(mmax_sd))
      C <- append(C, summary(fit)$parameters[2,1], after=length(C))
      M <- append(M, summary(fit)$parameters[3,1], after=length(M))
      M_sd <- append(M_sd, summary(fit)$parameters[3,2], after=length(M_sd))
    }

# Fit 0.2
tab02 <- filter(tab, Isolado!='lv7.am3.b' & Isolado!='tbe4.ext1.b2' & Isolado!='tbe5.am1.d' & Concentracao==0.2)
nome02 <- unique(tab02$Isolado)
for(n in seq_along(nome02)){   #passa por cada bicho
      tabfit <- filter(tab02, Isolado==nome02[n])
      d <- tabfit$Time
      y <- tabfit$Media_DO
      y0 <- min(y)
      fit <- nlsLM(y ~ y0 + C * exp(-exp(-mmax * (d - M))),
                  start=list(mmax=0.00001, C=(max(y)-min(y)), M=1000))
      Y0 <- append(Y0, y0, after=length(Y0))
      mmax <- append(mmax, summary(fit)$parameters[1,1], after=length(mmax))
      mmax_sd <- append(mmax_sd, summary(fit)$parameters[1,2], after=length(mmax_sd))
      C <- append(C, summary(fit)$parameters[2,1], after=length(C))
      M <- append(M, summary(fit)$parameters[3,1], after=length(M))
      M_sd <- append(M_sd, summary(fit)$parameters[3,2], after=length(M_sd))
    }


# Fit tbe5.am1.e 0,3 mol/L
tab03 <- filter(tab, Concentracao==0.3 & Isolado=='tbe5.am1.e')
d <- tab03$Time
y <- tab03$Media_DO
y0 <- min(y)
fit03 <- nlsLM(y ~ y0 + C * exp(-exp(-mmax * (d - M))),
                start=list(mmax=0.00001, C=(max(y)-min(y)), M=1200))
Y0 <- append(Y0, y0, after=length(Y0))
mmax <- append(mmax, summary(fit03)$parameters[1,1], after=length(mmax))
mmax_sd <- append(mmax_sd, summary(fit03)$parameters[1,2], after=length(mmax_sd))
C <- append(C, summary(fit03)$parameters[2,1], after=length(C))
M <- append(M, summary(fit03)$parameters[3,1], after=length(M))
M_sd <- append(M_sd, summary(fit03)$parameters[3,2], after=length(M_sd))


# ********** GRÁFICOS ************
nomes <- c(names, names, nome02, 'tbe5.am1.e')
conce <- c(rep(0, times=18), rep(0.1, times=18), rep(0.2, times=15), 0.3)
tab_param <- cbind.data.frame(nomes, conce, Y0, mmax, C, M, mmax_sd, M_sd)

write_xlsx(tab_param, 'Parâmetros fit.xlsx')

# 1 be5.am1.g
tabgraf <- filter(tab, Isolado=='be5.am1.g' & Concentracao<0.3)
param_0 <- filter(tab_param, nomes=='be5.am1.g' & conce==0)
param_0 <- c(param_0[1,3], param_0[1,4], param_0[1,5], param_0[1,6])
y_0 <- gompertz(param_0, d)
tabgraf_0 <-filter(tabgraf, Concentracao==0)
tabgraf_0 <- cbind.data.frame(tabgraf_0, y_0)
graf_0 <- ggplot(tabgraf_0, aes(x=Time))+
  geom_point(size=1, aes(y=Media_DO))+
  geom_line(aes(y=y_0), color='red')+
  geom_errorbar(aes(ymin=Media_DO-sd_DO, ymax=Media_DO+sd_DO), width=.3)+
  theme_bw()+
  labs(x='Time (min)', y='Optical Density', title='be5.am1.g')
graf_0
param_1 <- filter(tab_param, nomes=='be5.am1.g' & conce==0.1)
param_1 <- c(param_1[1,3], param_1[1,4], param_1[1,5], param_1[1,6])
y_1 <- gompertz(param_1, d)
tabgraf_1 <-filter(tabgraf, Concentracao==0.1)
tabgraf_1 <- cbind.data.frame(tabgraf_1, y_1)
graf_1 <- ggplot(tabgraf_1, aes(x=Time))+
  geom_point(size=1, aes(y=Media_DO))+
  geom_line(aes(y=y_1), color='red')+
  geom_errorbar(aes(ymin=Media_DO-sd_DO, ymax=Media_DO+sd_DO), width=.3)+
  theme_bw()+
  labs(x='Time (min)', y='Optical Density')
graf_1
param_2 <- filter(tab_param, nomes=='be5.am1.g' & conce==0.2)
param_2 <- c(param_2[1,3], param_2[1,4], param_2[1,5], param_2[1,6])
y_2 <- gompertz(param_2, d)
tabgraf_2 <-filter(tabgraf, Concentracao==0.2)
tabgraf_2 <- cbind.data.frame(tabgraf_2, y_2)
graf_2 <- ggplot(tabgraf_2, aes(x=Time))+
  geom_point(size=1, aes(y=Media_DO))+
  geom_line(aes(y=y_2), color='red')+
  geom_errorbar(aes(ymin=Media_DO-sd_DO, ymax=Media_DO+sd_DO), width=.3)+
  theme_bw()+
  labs(x='Time (min)', y='Optical Density')
graf_2
grid_be5.am1.g <- grid.arrange(graf_0, graf_1, graf_2, ncol=1)
ggsave('fit_be5am1g.png', grid_be5.am1.g, device='png', unit='cm', width=18, height=27, dpi=600)


# 2 be6.ext3.d1
tabgraf <- filter(tab, Isolado=='be6.ext3.d1' & Concentracao<0.3)
param_0 <- filter(tab_param, nomes=='be6.ext3.d1' & conce==0)
param_0 <- c(param_0[1,3], param_0[1,4], param_0[1,5], param_0[1,6])
y_0 <- gompertz(param_0, d)
tabgraf_0 <-filter(tabgraf, Concentracao==0)
tabgraf_0 <- cbind.data.frame(tabgraf_0, y_0)
graf_0 <- ggplot(tabgraf_0, aes(x=Time))+
  geom_point(size=1, aes(y=Media_DO))+
  geom_line(aes(y=y_0), color='red')+
  geom_errorbar(aes(ymin=Media_DO-sd_DO, ymax=Media_DO+sd_DO), width=.3)+
  theme_bw()+
  labs(x='Time (min)', y='Optical Density', title='be6.ext3.d1')
graf_0
param_1 <- filter(tab_param, nomes=='be6.ext3.d1' & conce==0.1)
param_1 <- c(param_1[1,3], param_1[1,4], param_1[1,5], param_1[1,6])
y_1 <- gompertz(param_1, d)
tabgraf_1 <-filter(tabgraf, Concentracao==0.1)
tabgraf_1 <- cbind.data.frame(tabgraf_1, y_1)
graf_1 <- ggplot(tabgraf_1, aes(x=Time))+
  geom_point(size=1, aes(y=Media_DO))+
  geom_line(aes(y=y_1), color='red')+
  geom_errorbar(aes(ymin=Media_DO-sd_DO, ymax=Media_DO+sd_DO), width=.3)+
  theme_bw()+
  labs(x='Time (min)', y='Optical Density')
graf_1
param_2 <- filter(tab_param, nomes=='be6.ext3.d1' & conce==0.2)
param_2 <- c(param_2[1,3], param_2[1,4], param_2[1,5], param_2[1,6])
y_2 <- gompertz(param_2, d)
tabgraf_2 <-filter(tabgraf, Concentracao==0.2)
tabgraf_2 <- cbind.data.frame(tabgraf_2, y_2)
graf_2 <- ggplot(tabgraf_2, aes(x=Time))+
  geom_point(size=1, aes(y=Media_DO))+
  geom_line(aes(y=y_2), color='red')+
  geom_errorbar(aes(ymin=Media_DO-sd_DO, ymax=Media_DO+sd_DO), width=.3)+
  theme_bw()+
  labs(x='Time (min)', y='Optical Density')
graf_2
grid_be6.ext3.d1 <- grid.arrange(graf_0, graf_1, graf_2, ncol=1)
ggsave('fit_be6ext3d1.png', grid_be6.ext3.d1, device='png', unit='cm', width=18, height=27, dpi=600)


# 3 be6.ext3.d2
tabgraf <- filter(tab, Isolado=='be6.ext3.d2' & Concentracao<0.3)
param_0 <- filter(tab_param, nomes=='be6.ext3.d2' & conce==0)
param_0 <- c(param_0[1,3], param_0[1,4], param_0[1,5], param_0[1,6])
y_0 <- gompertz(param_0, d)
tabgraf_0 <-filter(tabgraf, Concentracao==0)
tabgraf_0 <- cbind.data.frame(tabgraf_0, y_0)
graf_0 <- ggplot(tabgraf_0, aes(x=Time))+
  geom_point(size=1, aes(y=Media_DO))+
  geom_line(aes(y=y_0), color='red')+
  geom_errorbar(aes(ymin=Media_DO-sd_DO, ymax=Media_DO+sd_DO), width=.3)+
  theme_bw()+
  labs(x='Time (min)', y='Optical Density', title='be6.ext3.d2')
graf_0
param_1 <- filter(tab_param, nomes=='be6.ext3.d2' & conce==0.1)
param_1 <- c(param_1[1,3], param_1[1,4], param_1[1,5], param_1[1,6])
y_1 <- gompertz(param_1, d)
tabgraf_1 <-filter(tabgraf, Concentracao==0.1)
tabgraf_1 <- cbind.data.frame(tabgraf_1, y_1)
graf_1 <- ggplot(tabgraf_1, aes(x=Time))+
  geom_point(size=1, aes(y=Media_DO))+
  geom_line(aes(y=y_1), color='red')+
  geom_errorbar(aes(ymin=Media_DO-sd_DO, ymax=Media_DO+sd_DO), width=.3)+
  theme_bw()+
  labs(x='Time (min)', y='Optical Density')
graf_1
param_2 <- filter(tab_param, nomes=='be6.ext3.d2' & conce==0.2)
param_2 <- c(param_2[1,3], param_2[1,4], param_2[1,5], param_2[1,6])
y_2 <- gompertz(param_2, d)
tabgraf_2 <-filter(tabgraf, Concentracao==0.2)
tabgraf_2 <- cbind.data.frame(tabgraf_2, y_2)
graf_2 <- ggplot(tabgraf_2, aes(x=Time))+
  geom_point(size=1, aes(y=Media_DO))+
  geom_line(aes(y=y_2), color='red')+
  geom_errorbar(aes(ymin=Media_DO-sd_DO, ymax=Media_DO+sd_DO), width=.3)+
  theme_bw()+
  labs(x='Time (min)', y='Optical Density')
graf_2
grid_be6.ext3.d2 <- grid.arrange(graf_0, graf_1, graf_2, ncol=1)
ggsave('fit_be6ext3d2.png', grid_be6.ext3.d2, device='png', unit='cm', width=18, height=27, dpi=600)


# 4 colb
tabgraf <- filter(tab, Isolado=='colb' & Concentracao<0.3)
param_0 <- filter(tab_param, nomes=='colb' & conce==0)
param_0 <- c(param_0[1,3], param_0[1,4], param_0[1,5], param_0[1,6])
y_0 <- gompertz(param_0, d)
tabgraf_0 <-filter(tabgraf, Concentracao==0)
tabgraf_0 <- cbind.data.frame(tabgraf_0, y_0)
graf_0 <- ggplot(tabgraf_0, aes(x=Time))+
  geom_point(size=1, aes(y=Media_DO))+
  geom_line(aes(y=y_0), color='red')+
  geom_errorbar(aes(ymin=Media_DO-sd_DO, ymax=Media_DO+sd_DO), width=.3)+
  theme_bw()+
  labs(x='Time (min)', y='Optical Density', title='colb')
graf_0
param_1 <- filter(tab_param, nomes=='colb' & conce==0.1)
param_1 <- c(param_1[1,3], param_1[1,4], param_1[1,5], param_1[1,6])
y_1 <- gompertz(param_1, d)
tabgraf_1 <-filter(tabgraf, Concentracao==0.1)
tabgraf_1 <- cbind.data.frame(tabgraf_1, y_1)
graf_1 <- ggplot(tabgraf_1, aes(x=Time))+
  geom_point(size=1, aes(y=Media_DO))+
  geom_line(aes(y=y_1), color='red')+
  geom_errorbar(aes(ymin=Media_DO-sd_DO, ymax=Media_DO+sd_DO), width=.3)+
  theme_bw()+
  labs(x='Time (min)', y='Optical Density')
graf_1
param_2 <- filter(tab_param, nomes=='colb' & conce==0.2)
param_2 <- c(param_2[1,3], param_2[1,4], param_2[1,5], param_2[1,6])
y_2 <- gompertz(param_2, d)
tabgraf_2 <-filter(tabgraf, Concentracao==0.2)
tabgraf_2 <- cbind.data.frame(tabgraf_2, y_2)
graf_2 <- ggplot(tabgraf_2, aes(x=Time))+
  geom_point(size=1, aes(y=Media_DO))+
  geom_line(aes(y=y_2), color='red')+
  geom_errorbar(aes(ymin=Media_DO-sd_DO, ymax=Media_DO+sd_DO), width=.3)+
  theme_bw()+
  labs(x='Time (min)', y='Optical Density')
graf_2
grid_colb <- grid.arrange(graf_0, graf_1, graf_2, ncol=1)
ggsave('fit_colb.png', grid_colb, device='png', unit='cm', width=18, height=27, dpi=600)


# 5 lv7.am1.l
tabgraf <- filter(tab, Isolado=='lv7.am1.l' & Concentracao<0.3)
param_0 <- filter(tab_param, nomes=='lv7.am1.l' & conce==0)
param_0 <- c(param_0[1,3], param_0[1,4], param_0[1,5], param_0[1,6])
y_0 <- gompertz(param_0, d)
tabgraf_0 <-filter(tabgraf, Concentracao==0)
tabgraf_0 <- cbind.data.frame(tabgraf_0, y_0)
graf_0 <- ggplot(tabgraf_0, aes(x=Time))+
  geom_point(size=1, aes(y=Media_DO))+
  geom_line(aes(y=y_0), color='red')+
  geom_errorbar(aes(ymin=Media_DO-sd_DO, ymax=Media_DO+sd_DO), width=.3)+
  theme_bw()+
  labs(x='Time (min)', y='Optical Density', title='lv7.am1.l')
graf_0
param_1 <- filter(tab_param, nomes=='lv7.am1.l' & conce==0.1)
param_1 <- c(param_1[1,3], param_1[1,4], param_1[1,5], param_1[1,6])
y_1 <- gompertz(param_1, d)
tabgraf_1 <-filter(tabgraf, Concentracao==0.1)
tabgraf_1 <- cbind.data.frame(tabgraf_1, y_1)
graf_1 <- ggplot(tabgraf_1, aes(x=Time))+
  geom_point(size=1, aes(y=Media_DO))+
  geom_line(aes(y=y_1), color='red')+
  geom_errorbar(aes(ymin=Media_DO-sd_DO, ymax=Media_DO+sd_DO), width=.3)+
  theme_bw()+
  labs(x='Time (min)', y='Optical Density')
graf_1
param_2 <- filter(tab_param, nomes=='lv7.am1.l' & conce==0.2)
param_2 <- c(param_2[1,3], param_2[1,4], param_2[1,5], param_2[1,6])
y_2 <- gompertz(param_2, d)
tabgraf_2 <-filter(tabgraf, Concentracao==0.2)
tabgraf_2 <- cbind.data.frame(tabgraf_2, y_2)
graf_2 <- ggplot(tabgraf_2, aes(x=Time))+
  geom_point(size=1, aes(y=Media_DO))+
  geom_line(aes(y=y_2), color='red')+
  geom_errorbar(aes(ymin=Media_DO-sd_DO, ymax=Media_DO+sd_DO), width=.3)+
  theme_bw()+
  labs(x='Time (min)', y='Optical Density')
graf_2
grid_lv7.am1.l <- grid.arrange(graf_0, graf_1, graf_2, ncol=1)
ggsave('fit_lv7am1l.png', grid_lv7.am1.l, device='png', unit='cm', width=18, height=27, dpi=600)


# 6 lv7.am3.b
tabgraf <- filter(tab, Isolado=='lv7.am3.b' & Concentracao<0.3)
param_0 <- filter(tab_param, nomes=='lv7.am3.b' & conce==0)
param_0 <- c(param_0[1,3], param_0[1,4], param_0[1,5], param_0[1,6])
y_0 <- gompertz(param_0, d)
tabgraf_0 <-filter(tabgraf, Concentracao==0)
tabgraf_0 <- cbind.data.frame(tabgraf_0, y_0)
graf_0 <- ggplot(tabgraf_0, aes(x=Time))+
  geom_point(size=1, aes(y=Media_DO))+
  geom_line(aes(y=y_0), color='red')+
  geom_errorbar(aes(ymin=Media_DO-sd_DO, ymax=Media_DO+sd_DO), width=.3)+
  theme_bw()+
  labs(x='Time (min)', y='Optical Density', title='lv7.am3.b')
graf_0
param_1 <- filter(tab_param, nomes=='lv7.am3.b' & conce==0.1)
param_1 <- c(param_1[1,3], param_1[1,4], param_1[1,5], param_1[1,6])
y_1 <- gompertz(param_1, d)
tabgraf_1 <-filter(tabgraf, Concentracao==0.1)
tabgraf_1 <- cbind.data.frame(tabgraf_1, y_1)
graf_1 <- ggplot(tabgraf_1, aes(x=Time))+
  geom_point(size=1, aes(y=Media_DO))+
  geom_line(aes(y=y_1), color='red')+
  geom_errorbar(aes(ymin=Media_DO-sd_DO, ymax=Media_DO+sd_DO), width=.3)+
  theme_bw()+
  labs(x='Time (min)', y='Optical Density')
graf_1
grid_lv7.am3.b <- grid.arrange(graf_0, graf_1, ncol=1)
ggsave('fit_lv7am3b.png', grid_lv7.am3.b, device='png', unit='cm', width=18, height=27, dpi=600)


# 7 p4d
tabgraf <- filter(tab, Isolado=='p4d' & Concentracao<0.3)
param_0 <- filter(tab_param, nomes=='p4d' & conce==0)
param_0 <- c(param_0[1,3], param_0[1,4], param_0[1,5], param_0[1,6])
y_0 <- gompertz(param_0, d)
tabgraf_0 <-filter(tabgraf, Concentracao==0)
tabgraf_0 <- cbind.data.frame(tabgraf_0, y_0)
graf_0 <- ggplot(tabgraf_0, aes(x=Time))+
  geom_point(size=1, aes(y=Media_DO))+
  geom_line(aes(y=y_0), color='red')+
  geom_errorbar(aes(ymin=Media_DO-sd_DO, ymax=Media_DO+sd_DO), width=.3)+
  theme_bw()+
  labs(x='Time (min)', y='Optical Density', title='p4d')
graf_0
param_1 <- filter(tab_param, nomes=='p4d' & conce==0.1)
param_1 <- c(param_1[1,3], param_1[1,4], param_1[1,5], param_1[1,6])
y_1 <- gompertz(param_1, d)
tabgraf_1 <-filter(tabgraf, Concentracao==0.1)
tabgraf_1 <- cbind.data.frame(tabgraf_1, y_1)
graf_1 <- ggplot(tabgraf_1, aes(x=Time))+
  geom_point(size=1, aes(y=Media_DO))+
  geom_line(aes(y=y_1), color='red')+
  geom_errorbar(aes(ymin=Media_DO-sd_DO, ymax=Media_DO+sd_DO), width=.3)+
  theme_bw()+
  labs(x='Time (min)', y='Optical Density')
graf_1
param_2 <- filter(tab_param, nomes=='p4d' & conce==0.2)
param_2 <- c(param_2[1,3], param_2[1,4], param_2[1,5], param_2[1,6])
y_2 <- gompertz(param_2, d)
tabgraf_2 <-filter(tabgraf, Concentracao==0.2)
tabgraf_2 <- cbind.data.frame(tabgraf_2, y_2)
graf_2 <- ggplot(tabgraf_2, aes(x=Time))+
  geom_point(size=1, aes(y=Media_DO))+
  geom_line(aes(y=y_2), color='red')+
  geom_errorbar(aes(ymin=Media_DO-sd_DO, ymax=Media_DO+sd_DO), width=.3)+
  theme_bw()+
  labs(x='Time (min)', y='Optical Density')
graf_2
grid_p4d <- grid.arrange(graf_0, graf_1, graf_2, ncol=1)
ggsave('fit_p4d.png', grid_p4d, device='png', unit='cm', width=18, height=27, dpi=600)


# 8 tbe2.ext2.g
tabgraf <- filter(tab, Isolado=='tbe2.ext2.g' & Concentracao<0.3)
param_0 <- filter(tab_param, nomes=='tbe2.ext2.g' & conce==0)
param_0 <- c(param_0[1,3], param_0[1,4], param_0[1,5], param_0[1,6])
y_0 <- gompertz(param_0, d)
tabgraf_0 <-filter(tabgraf, Concentracao==0)
tabgraf_0 <- cbind.data.frame(tabgraf_0, y_0)
graf_0 <- ggplot(tabgraf_0, aes(x=Time))+
  geom_point(size=1, aes(y=Media_DO))+
  geom_line(aes(y=y_0), color='red')+
  geom_errorbar(aes(ymin=Media_DO-sd_DO, ymax=Media_DO+sd_DO), width=.3)+
  theme_bw()+
  labs(x='Time (min)', y='Optical Density', title='tbe2.ext2.g')
graf_0
param_1 <- filter(tab_param, nomes=='tbe2.ext2.g' & conce==0.1)
param_1 <- c(param_1[1,3], param_1[1,4], param_1[1,5], param_1[1,6])
y_1 <- gompertz(param_1, d)
tabgraf_1 <-filter(tabgraf, Concentracao==0.1)
tabgraf_1 <- cbind.data.frame(tabgraf_1, y_1)
graf_1 <- ggplot(tabgraf_1, aes(x=Time))+
  geom_point(size=1, aes(y=Media_DO))+
  geom_line(aes(y=y_1), color='red')+
  geom_errorbar(aes(ymin=Media_DO-sd_DO, ymax=Media_DO+sd_DO), width=.3)+
  theme_bw()+
  labs(x='Time (min)', y='Optical Density')
graf_1
param_2 <- filter(tab_param, nomes=='tbe2.ext2.g' & conce==0.2)
param_2 <- c(param_2[1,3], param_2[1,4], param_2[1,5], param_2[1,6])
y_2 <- gompertz(param_2, d)
tabgraf_2 <-filter(tabgraf, Concentracao==0.2)
tabgraf_2 <- cbind.data.frame(tabgraf_2, y_2)
graf_2 <- ggplot(tabgraf_2, aes(x=Time))+
  geom_point(size=1, aes(y=Media_DO))+
  geom_line(aes(y=y_2), color='red')+
  geom_errorbar(aes(ymin=Media_DO-sd_DO, ymax=Media_DO+sd_DO), width=.3)+
  theme_bw()+
  labs(x='Time (min)', y='Optical Density')
graf_2
grid_tbe2.ext2.g <- grid.arrange(graf_0, graf_1, graf_2, ncol=1)
ggsave('fit_tbe2ext2g.png', grid_tbe2.ext2.g, device='png', unit='cm', width=18, height=27, dpi=600)


# 9 tbe2.sob.d
tabgraf <- filter(tab, Isolado=='tbe2.sob.d' & Concentracao<0.3)
param_0 <- filter(tab_param, nomes=='tbe2.sob.d' & conce==0)
param_0 <- c(param_0[1,3], param_0[1,4], param_0[1,5], param_0[1,6])
y_0 <- gompertz(param_0, d)
tabgraf_0 <-filter(tabgraf, Concentracao==0)
tabgraf_0 <- cbind.data.frame(tabgraf_0, y_0)
graf_0 <- ggplot(tabgraf_0, aes(x=Time))+
  geom_point(size=1, aes(y=Media_DO))+
  geom_line(aes(y=y_0), color='red')+
  geom_errorbar(aes(ymin=Media_DO-sd_DO, ymax=Media_DO+sd_DO), width=.3)+
  theme_bw()+
  labs(x='Time (min)', y='Optical Density', title='tbe2.sob.d')
graf_0
param_1 <- filter(tab_param, nomes=='tbe2.sob.d' & conce==0.1)
param_1 <- c(param_1[1,3], param_1[1,4], param_1[1,5], param_1[1,6])
y_1 <- gompertz(param_1, d)
tabgraf_1 <-filter(tabgraf, Concentracao==0.1)
tabgraf_1 <- cbind.data.frame(tabgraf_1, y_1)
graf_1 <- ggplot(tabgraf_1, aes(x=Time))+
  geom_point(size=1, aes(y=Media_DO))+
  geom_line(aes(y=y_1), color='red')+
  geom_errorbar(aes(ymin=Media_DO-sd_DO, ymax=Media_DO+sd_DO), width=.3)+
  theme_bw()+
  labs(x='Time (min)', y='Optical Density')
graf_1
param_2 <- filter(tab_param, nomes=='tbe2.sob.d' & conce==0.2)
param_2 <- c(param_2[1,3], param_2[1,4], param_2[1,5], param_2[1,6])
y_2 <- gompertz(param_2, d)
tabgraf_2 <-filter(tabgraf, Concentracao==0.2)
tabgraf_2 <- cbind.data.frame(tabgraf_2, y_2)
graf_2 <- ggplot(tabgraf_2, aes(x=Time))+
  geom_point(size=1, aes(y=Media_DO))+
  geom_line(aes(y=y_2), color='red')+
  geom_errorbar(aes(ymin=Media_DO-sd_DO, ymax=Media_DO+sd_DO), width=.3)+
  theme_bw()+
  labs(x='Time (min)', y='Optical Density')
graf_2
grid_tbe2.sob.d <- grid.arrange(graf_0, graf_1, graf_2, ncol=1)
ggsave('fit_tbe2sobd.png', grid_tbe2.sob.d, device='png', unit='cm', width=18, height=27, dpi=600)


# 10 tbe4.am2.d
tabgraf <- filter(tab, Isolado=='tbe4.am2.d' & Concentracao<0.3)
param_0 <- filter(tab_param, nomes=='tbe4.am2.d' & conce==0)
param_0 <- c(param_0[1,3], param_0[1,4], param_0[1,5], param_0[1,6])
y_0 <- gompertz(param_0, d)
tabgraf_0 <-filter(tabgraf, Concentracao==0)
tabgraf_0 <- cbind.data.frame(tabgraf_0, y_0)
graf_0 <- ggplot(tabgraf_0, aes(x=Time))+
  geom_point(size=1, aes(y=Media_DO))+
  geom_line(aes(y=y_0), color='red')+
  geom_errorbar(aes(ymin=Media_DO-sd_DO, ymax=Media_DO+sd_DO), width=.3)+
  theme_bw()+
  labs(x='Time (min)', y='Optical Density', title='tbe4.am2.d')
graf_0
param_1 <- filter(tab_param, nomes=='tbe4.am2.d' & conce==0.1)
param_1 <- c(param_1[1,3], param_1[1,4], param_1[1,5], param_1[1,6])
y_1 <- gompertz(param_1, d)
tabgraf_1 <-filter(tabgraf, Concentracao==0.1)
tabgraf_1 <- cbind.data.frame(tabgraf_1, y_1)
graf_1 <- ggplot(tabgraf_1, aes(x=Time))+
  geom_point(size=1, aes(y=Media_DO))+
  geom_line(aes(y=y_1), color='red')+
  geom_errorbar(aes(ymin=Media_DO-sd_DO, ymax=Media_DO+sd_DO), width=.3)+
  theme_bw()+
  labs(x='Time (min)', y='Optical Density')
graf_1
param_2 <- filter(tab_param, nomes=='tbe4.am2.d' & conce==0.2)
param_2 <- c(param_2[1,3], param_2[1,4], param_2[1,5], param_2[1,6])
y_2 <- gompertz(param_2, d)
tabgraf_2 <-filter(tabgraf, Concentracao==0.2)
tabgraf_2 <- cbind.data.frame(tabgraf_2, y_2)
graf_2 <- ggplot(tabgraf_2, aes(x=Time))+
  geom_point(size=1, aes(y=Media_DO))+
  geom_line(aes(y=y_2), color='red')+
  geom_errorbar(aes(ymin=Media_DO-sd_DO, ymax=Media_DO+sd_DO), width=.3)+
  theme_bw()+
  labs(x='Time (min)', y='Optical Density')
graf_2
grid_tbe4.am2.d <- grid.arrange(graf_0, graf_1, graf_2, ncol=1)
ggsave('fit_tbe4am2d.png', grid_tbe4.am2.d, device='png', unit='cm', width=18, height=27, dpi=600)


# 11 tbe4.ext1.b2
tabgraf <- filter(tab, Isolado=='tbe4.ext1.b2' & Concentracao<0.3)
param_0 <- filter(tab_param, nomes=='tbe4.ext1.b2' & conce==0)
param_0 <- c(param_0[1,3], param_0[1,4], param_0[1,5], param_0[1,6])
y_0 <- gompertz(param_0, d)
tabgraf_0 <-filter(tabgraf, Concentracao==0)
tabgraf_0 <- cbind.data.frame(tabgraf_0, y_0)
graf_0 <- ggplot(tabgraf_0, aes(x=Time))+
  geom_point(size=1, aes(y=Media_DO))+
  geom_line(aes(y=y_0), color='red')+
  geom_errorbar(aes(ymin=Media_DO-sd_DO, ymax=Media_DO+sd_DO), width=.3)+
  theme_bw()+
  labs(x='Time (min)', y='Optical Density', title='tbe4.ext1.b2')
graf_0
param_1 <- filter(tab_param, nomes=='tbe4.ext1.b2' & conce==0.1)
param_1 <- c(param_1[1,3], param_1[1,4], param_1[1,5], param_1[1,6])
y_1 <- gompertz(param_1, d)
tabgraf_1 <-filter(tabgraf, Concentracao==0.1)
tabgraf_1 <- cbind.data.frame(tabgraf_1, y_1)
graf_1 <- ggplot(tabgraf_1, aes(x=Time))+
  geom_point(size=1, aes(y=Media_DO))+
  geom_line(aes(y=y_1), color='red')+
  geom_errorbar(aes(ymin=Media_DO-sd_DO, ymax=Media_DO+sd_DO), width=.3)+
  theme_bw()+
  labs(x='Time (min)', y='Optical Density')
graf_1
grid_tbe4.ext1.b2 <- grid.arrange(graf_0, graf_1, ncol=1)
ggsave('fit_tbe4ext1b2.png', grid_tbe4.ext1.b2, device='png', unit='cm', width=18, height=27, dpi=600)


# 12 tbe5.am1.d
tabgraf <- filter(tab, Isolado=='tbe5.am1.d' & Concentracao<0.3)
param_0 <- filter(tab_param, nomes=='tbe5.am1.d' & conce==0)
param_0 <- c(param_0[1,3], param_0[1,4], param_0[1,5], param_0[1,6])
y_0 <- gompertz(param_0, d)
tabgraf_0 <-filter(tabgraf, Concentracao==0)
tabgraf_0 <- cbind.data.frame(tabgraf_0, y_0)
graf_0 <- ggplot(tabgraf_0, aes(x=Time))+
  geom_point(size=1, aes(y=Media_DO))+
  geom_line(aes(y=y_0), color='red')+
  geom_errorbar(aes(ymin=Media_DO-sd_DO, ymax=Media_DO+sd_DO), width=.3)+
  theme_bw()+
  labs(x='Time (min)', y='Optical Density', title='tbe5.am1.d')
graf_0
param_1 <- filter(tab_param, nomes=='tbe5.am1.d' & conce==0.1)
param_1 <- c(param_1[1,3], param_1[1,4], param_1[1,5], param_1[1,6])
y_1 <- gompertz(param_1, d)
tabgraf_1 <-filter(tabgraf, Concentracao==0.1)
tabgraf_1 <- cbind.data.frame(tabgraf_1, y_1)
graf_1 <- ggplot(tabgraf_1, aes(x=Time))+
  geom_point(size=1, aes(y=Media_DO))+
  geom_line(aes(y=y_1), color='red')+
  geom_errorbar(aes(ymin=Media_DO-sd_DO, ymax=Media_DO+sd_DO), width=.3)+
  theme_bw()+
  labs(x='Time (min)', y='Optical Density')
graf_1
grid_tbe5.am1.d <- grid.arrange(graf_0, graf_1, ncol=1)
ggsave('fit_tbe5am1d.png', grid_tbe5.am1.d, device='png', unit='cm', width=18, height=27, dpi=600)


# 13 tbe5.am1.e
tabgraf <- filter(tab, Isolado=='tbe5.am1.e' & Concentracao<=0.3)
param_0 <- filter(tab_param, nomes=='tbe5.am1.e' & conce==0)
param_0 <- c(param_0[1,3], param_0[1,4], param_0[1,5], param_0[1,6])
y_0 <- gompertz(param_0, d)
tabgraf_0 <-filter(tabgraf, Concentracao==0)
tabgraf_0 <- cbind.data.frame(tabgraf_0, y_0)
graf_0 <- ggplot(tabgraf_0, aes(x=Time))+
  geom_point(size=1, aes(y=Media_DO))+
  geom_line(aes(y=y_0), color='red')+
  geom_errorbar(aes(ymin=Media_DO-sd_DO, ymax=Media_DO+sd_DO), width=.3)+
  theme_bw()+
  labs(x='Time (min)', y='Optical Density', title='tbe5.am1.e')
graf_0
param_1 <- filter(tab_param, nomes=='tbe5.am1.e' & conce==0.1)
param_1 <- c(param_1[1,3], param_1[1,4], param_1[1,5], param_1[1,6])
y_1 <- gompertz(param_1, d)
tabgraf_1 <-filter(tabgraf, Concentracao==0.1)
tabgraf_1 <- cbind.data.frame(tabgraf_1, y_1)
graf_1 <- ggplot(tabgraf_1, aes(x=Time))+
  geom_point(size=1, aes(y=Media_DO))+
  geom_line(aes(y=y_1), color='red')+
  geom_errorbar(aes(ymin=Media_DO-sd_DO, ymax=Media_DO+sd_DO), width=.3)+
  theme_bw()+
  labs(x='Time (min)', y='Optical Density')
graf_1
param_2 <- filter(tab_param, nomes=='tbe5.am1.e' & conce==0.2)
param_2 <- c(param_2[1,3], param_2[1,4], param_2[1,5], param_2[1,6])
y_2 <- gompertz(param_2, d)
tabgraf_2 <-filter(tabgraf, Concentracao==0.2)
tabgraf_2 <- cbind.data.frame(tabgraf_2, y_2)
graf_2 <- ggplot(tabgraf_2, aes(x=Time))+
  geom_point(size=1, aes(y=Media_DO))+
  geom_line(aes(y=y_2), color='red')+
  geom_errorbar(aes(ymin=Media_DO-sd_DO, ymax=Media_DO+sd_DO), width=.3)+
  theme_bw()+
  labs(x='Time (min)', y='Optical Density')
graf_2
param_3 <- filter(tab_param, nomes=='tbe5.am1.e' & conce==0.3)
param_3 <- c(param_3[1,3], param_3[1,4], param_3[1,5], param_3[1,6])
y_3 <- gompertz(param_3, d)
tabgraf_3 <-filter(tabgraf, Concentracao==0.3)
tabgraf_3 <- cbind.data.frame(tabgraf_3, y_3)
graf_3 <- ggplot(tabgraf_3, aes(x=Time))+
  geom_point(size=1, aes(y=Media_DO))+
  geom_line(aes(y=y_3), color='red')+
  geom_errorbar(aes(ymin=Media_DO-sd_DO, ymax=Media_DO+sd_DO), width=.3)+
  theme_bw()+
  labs(x='Time (min)', y='Optical Density')
graf_3
grid_tbe5.am1.e <- grid.arrange(graf_0, graf_1, graf_2, graf_3, ncol=1)
ggsave('fit_tbe5am1e.png', grid_tbe5.am1.e, device='png', unit='cm', width=18, height=27, dpi=600)


# 14 tlv6.am1.j
tabgraf <- filter(tab, Isolado=='tlv6.am1.j' & Concentracao<0.3)
param_0 <- filter(tab_param, nomes=='tlv6.am1.j' & conce==0)
param_0 <- c(param_0[1,3], param_0[1,4], param_0[1,5], param_0[1,6])
y_0 <- gompertz(param_0, d)
tabgraf_0 <-filter(tabgraf, Concentracao==0)
tabgraf_0 <- cbind.data.frame(tabgraf_0, y_0)
graf_0 <- ggplot(tabgraf_0, aes(x=Time))+
  geom_point(size=1, aes(y=Media_DO))+
  geom_line(aes(y=y_0), color='red')+
  geom_errorbar(aes(ymin=Media_DO-sd_DO, ymax=Media_DO+sd_DO), width=.3)+
  theme_bw()+
  labs(x='Time (min)', y='Optical Density', title='tlv6.am1.j')
graf_0
param_1 <- filter(tab_param, nomes=='tlv6.am1.j' & conce==0.1)
param_1 <- c(param_1[1,3], param_1[1,4], param_1[1,5], param_1[1,6])
y_1 <- gompertz(param_1, d)
tabgraf_1 <-filter(tabgraf, Concentracao==0.1)
tabgraf_1 <- cbind.data.frame(tabgraf_1, y_1)
graf_1 <- ggplot(tabgraf_1, aes(x=Time))+
  geom_point(size=1, aes(y=Media_DO))+
  geom_line(aes(y=y_1), color='red')+
  geom_errorbar(aes(ymin=Media_DO-sd_DO, ymax=Media_DO+sd_DO), width=.3)+
  theme_bw()+
  labs(x='Time (min)', y='Optical Density')
graf_1
param_2 <- filter(tab_param, nomes=='tlv6.am1.j' & conce==0.2)
param_2 <- c(param_2[1,3], param_2[1,4], param_2[1,5], param_2[1,6])
y_2 <- gompertz(param_2, d)
tabgraf_2 <-filter(tabgraf, Concentracao==0.2)
tabgraf_2 <- cbind.data.frame(tabgraf_2, y_2)
graf_2 <- ggplot(tabgraf_2, aes(x=Time))+
  geom_point(size=1, aes(y=Media_DO))+
  geom_line(aes(y=y_2), color='red')+
  geom_errorbar(aes(ymin=Media_DO-sd_DO, ymax=Media_DO+sd_DO), width=.3)+
  theme_bw()+
  labs(x='Time (min)', y='Optical Density')
graf_2
grid_tlv6.am1.j <- grid.arrange(graf_0, graf_1, graf_2, ncol=1)
ggsave('fit_tlv6am1j.png', grid_tlv6.am1.j, device='png', unit='cm', width=18, height=27, dpi=600)


# 15 tlv7.am1.c
tabgraf <- filter(tab, Isolado=='tlv7.am1.c' & Concentracao<0.3)
param_0 <- filter(tab_param, nomes=='tlv7.am1.c' & conce==0)
param_0 <- c(param_0[1,3], param_0[1,4], param_0[1,5], param_0[1,6])
y_0 <- gompertz(param_0, d)
tabgraf_0 <-filter(tabgraf, Concentracao==0)
tabgraf_0 <- cbind.data.frame(tabgraf_0, y_0)
graf_0 <- ggplot(tabgraf_0, aes(x=Time))+
  geom_point(size=1, aes(y=Media_DO))+
  geom_line(aes(y=y_0), color='red')+
  geom_errorbar(aes(ymin=Media_DO-sd_DO, ymax=Media_DO+sd_DO), width=.3)+
  theme_bw()+
  labs(x='Time (min)', y='Optical Density', title='tlv7.am1.c')
graf_0
param_1 <- filter(tab_param, nomes=='tlv7.am1.c' & conce==0.1)
param_1 <- c(param_1[1,3], param_1[1,4], param_1[1,5], param_1[1,6])
y_1 <- gompertz(param_1, d)
tabgraf_1 <-filter(tabgraf, Concentracao==0.1)
tabgraf_1 <- cbind.data.frame(tabgraf_1, y_1)
graf_1 <- ggplot(tabgraf_1, aes(x=Time))+
  geom_point(size=1, aes(y=Media_DO))+
  geom_line(aes(y=y_1), color='red')+
  geom_errorbar(aes(ymin=Media_DO-sd_DO, ymax=Media_DO+sd_DO), width=.3)+
  theme_bw()+
  labs(x='Time (min)', y='Optical Density')
graf_1
param_2 <- filter(tab_param, nomes=='tlv7.am1.c' & conce==0.2)
param_2 <- c(param_2[1,3], param_2[1,4], param_2[1,5], param_2[1,6])
y_2 <- gompertz(param_2, d)
tabgraf_2 <-filter(tabgraf, Concentracao==0.2)
tabgraf_2 <- cbind.data.frame(tabgraf_2, y_2)
graf_2 <- ggplot(tabgraf_2, aes(x=Time))+
  geom_point(size=1, aes(y=Media_DO))+
  geom_line(aes(y=y_2), color='red')+
  geom_errorbar(aes(ymin=Media_DO-sd_DO, ymax=Media_DO+sd_DO), width=.3)+
  theme_bw()+
  labs(x='Time (min)', y='Optical Density')
graf_2
grid_tlv7.am1.c <- grid.arrange(graf_0, graf_1, graf_2, ncol=1)
ggsave('fit_tlv7am1c.png', grid_tlv7.am1.c, device='png', unit='cm', width=18, height=27, dpi=600)


# 16 tlv7.am1.e
tabgraf <- filter(tab, Isolado=='tlv7.am1.e' & Concentracao<0.3)
param_0 <- filter(tab_param, nomes=='tlv7.am1.e' & conce==0)
param_0 <- c(param_0[1,3], param_0[1,4], param_0[1,5], param_0[1,6])
y_0 <- gompertz(param_0, d)
tabgraf_0 <-filter(tabgraf, Concentracao==0)
tabgraf_0 <- cbind.data.frame(tabgraf_0, y_0)
graf_0 <- ggplot(tabgraf_0, aes(x=Time))+
  geom_point(size=1, aes(y=Media_DO))+
  geom_line(aes(y=y_0), color='red')+
  geom_errorbar(aes(ymin=Media_DO-sd_DO, ymax=Media_DO+sd_DO), width=.3)+
  theme_bw()+
  labs(x='Time (min)', y='Optical Density', title='tlv7.am1.e')
graf_0
param_1 <- filter(tab_param, nomes=='tlv7.am1.e' & conce==0.1)
param_1 <- c(param_1[1,3], param_1[1,4], param_1[1,5], param_1[1,6])
y_1 <- gompertz(param_1, d)
tabgraf_1 <-filter(tabgraf, Concentracao==0.1)
tabgraf_1 <- cbind.data.frame(tabgraf_1, y_1)
graf_1 <- ggplot(tabgraf_1, aes(x=Time))+
  geom_point(size=1, aes(y=Media_DO))+
  geom_line(aes(y=y_1), color='red')+
  geom_errorbar(aes(ymin=Media_DO-sd_DO, ymax=Media_DO+sd_DO), width=.3)+
  theme_bw()+
  labs(x='Time (min)', y='Optical Density')
graf_1
param_2 <- filter(tab_param, nomes=='tlv7.am1.e' & conce==0.2)
param_2 <- c(param_2[1,3], param_2[1,4], param_2[1,5], param_2[1,6])
y_2 <- gompertz(param_2, d)
tabgraf_2 <-filter(tabgraf, Concentracao==0.2)
tabgraf_2 <- cbind.data.frame(tabgraf_2, y_2)
graf_2 <- ggplot(tabgraf_2, aes(x=Time))+
  geom_point(size=1, aes(y=Media_DO))+
  geom_line(aes(y=y_2), color='red')+
  geom_errorbar(aes(ymin=Media_DO-sd_DO, ymax=Media_DO+sd_DO), width=.3)+
  theme_bw()+
  labs(x='Time (min)', y='Optical Density')
graf_2
grid_tlv7.am1.e <- grid.arrange(graf_0, graf_1, graf_2, ncol=1)
ggsave('fit_tlv7am1e.png', grid_tlv7.am1.e, device='png', unit='cm', width=18, height=27, dpi=600)


# 17 tlv7.am4.h
tabgraf <- filter(tab, Isolado=='tlv7.am4.h' & Concentracao<0.3)
param_0 <- filter(tab_param, nomes=='tlv7.am4.h' & conce==0)
param_0 <- c(param_0[1,3], param_0[1,4], param_0[1,5], param_0[1,6])
y_0 <- gompertz(param_0, d)
tabgraf_0 <-filter(tabgraf, Concentracao==0)
tabgraf_0 <- cbind.data.frame(tabgraf_0, y_0)
graf_0 <- ggplot(tabgraf_0, aes(x=Time))+
  geom_point(size=1, aes(y=Media_DO))+
  geom_line(aes(y=y_0), color='red')+
  geom_errorbar(aes(ymin=Media_DO-sd_DO, ymax=Media_DO+sd_DO), width=.3)+
  theme_bw()+
  labs(x='Time (min)', y='Optical Density', title='tlv7.am4.h')
graf_0
param_1 <- filter(tab_param, nomes=='tlv7.am4.h' & conce==0.1)
param_1 <- c(param_1[1,3], param_1[1,4], param_1[1,5], param_1[1,6])
y_1 <- gompertz(param_1, d)
tabgraf_1 <-filter(tabgraf, Concentracao==0.1)
tabgraf_1 <- cbind.data.frame(tabgraf_1, y_1)
graf_1 <- ggplot(tabgraf_1, aes(x=Time))+
  geom_point(size=1, aes(y=Media_DO))+
  geom_line(aes(y=y_1), color='red')+
  geom_errorbar(aes(ymin=Media_DO-sd_DO, ymax=Media_DO+sd_DO), width=.3)+
  theme_bw()+
  labs(x='Time (min)', y='Optical Density')
graf_1
param_2 <- filter(tab_param, nomes=='tlv7.am4.h' & conce==0.2)
param_2 <- c(param_2[1,3], param_2[1,4], param_2[1,5], param_2[1,6])
y_2 <- gompertz(param_2, d)
tabgraf_2 <-filter(tabgraf, Concentracao==0.2)
tabgraf_2 <- cbind.data.frame(tabgraf_2, y_2)
graf_2 <- ggplot(tabgraf_2, aes(x=Time))+
  geom_point(size=1, aes(y=Media_DO))+
  geom_line(aes(y=y_2), color='red')+
  geom_errorbar(aes(ymin=Media_DO-sd_DO, ymax=Media_DO+sd_DO), width=.3)+
  theme_bw()+
  labs(x='Time (min)', y='Optical Density')
graf_2
grid_tlv7.am4.h <- grid.arrange(graf_0, graf_1, graf_2, ncol=1)
ggsave('fit_tlv7am4h.png', grid_tlv7.am4.h, device='png', unit='cm', width=18, height=27, dpi=600)


# 18 tlv7.ext2.k
tabgraf <- filter(tab, Isolado=='tlv7.ext2.k' & Concentracao<0.3)
param_0 <- filter(tab_param, nomes=='tlv7.ext2.k' & conce==0)
param_0 <- c(param_0[1,3], param_0[1,4], param_0[1,5], param_0[1,6])
y_0 <- gompertz(param_0, d)
tabgraf_0 <-filter(tabgraf, Concentracao==0)
tabgraf_0 <- cbind.data.frame(tabgraf_0, y_0)
graf_0 <- ggplot(tabgraf_0, aes(x=Time))+
  geom_point(size=1, aes(y=Media_DO))+
  geom_line(aes(y=y_0), color='red')+
  geom_errorbar(aes(ymin=Media_DO-sd_DO, ymax=Media_DO+sd_DO), width=.3)+
  theme_bw()+
  labs(x='Time (min)', y='Optical Density', title='tlv7.ext2.k')
graf_0
param_1 <- filter(tab_param, nomes=='tlv7.ext2.k' & conce==0.1)
param_1 <- c(param_1[1,3], param_1[1,4], param_1[1,5], param_1[1,6])
y_1 <- gompertz(param_1, d)
tabgraf_1 <-filter(tabgraf, Concentracao==0.1)
tabgraf_1 <- cbind.data.frame(tabgraf_1, y_1)
graf_1 <- ggplot(tabgraf_1, aes(x=Time))+
  geom_point(size=1, aes(y=Media_DO))+
  geom_line(aes(y=y_1), color='red')+
  geom_errorbar(aes(ymin=Media_DO-sd_DO, ymax=Media_DO+sd_DO), width=.3)+
  theme_bw()+
  labs(x='Time (min)', y='Optical Density')
graf_1
param_2 <- filter(tab_param, nomes=='tlv7.ext2.k' & conce==0.2)
param_2 <- c(param_2[1,3], param_2[1,4], param_2[1,5], param_2[1,6])
y_2 <- gompertz(param_2, d)
tabgraf_2 <-filter(tabgraf, Concentracao==0.2)
tabgraf_2 <- cbind.data.frame(tabgraf_2, y_2)
graf_2 <- ggplot(tabgraf_2, aes(x=Time))+
  geom_point(size=1, aes(y=Media_DO))+
  geom_line(aes(y=y_2), color='red')+
  geom_errorbar(aes(ymin=Media_DO-sd_DO, ymax=Media_DO+sd_DO), width=.3)+
  theme_bw()+
  labs(x='Time (min)', y='Optical Density')
graf_2
grid_tlv7.ext2.k <- grid.arrange(graf_0, graf_1, graf_2, ncol=1)
ggsave('fit_tlv7ext2k.png', grid_tlv7.ext2.k, device='png', unit='cm', width=18, height=27, dpi=600)
