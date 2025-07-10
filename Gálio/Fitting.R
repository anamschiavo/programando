# ********** PACOTES E BIBLIOTECAS **********
install.packages('dplyr')
library(dplyr)
install.packages('janitor')
library(janitor)
install.packages('ggplot2')
library(ggplot2)
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

baranyi <- function(params, x) {
  params[1] + params[2] * (x + (1/params[2]) * log(exp(-params[2]*x) +
  exp(-params[2] * params[3]) - exp(-params[2] * (x + params[3])))) -
  log(1 + ((exp(params[2] * (x + (1/params[2]) * log(exp(-params[2]*x) +
  exp(-params[2] * params[3]) - exp(-params[2] * (x + params[3])))))-1)/
    (exp(params[4]-params[1]))))
}


# Com green rust
# ********** FITTING BARANYI *********
dados <- read.table('430nm.tsv', sep='\t', dec=',', header=TRUE, fill=TRUE)
d <- as.numeric(dados$Time)
y <- as.numeric(dados$DO)
sd <- as.numeric(dados$sd)
y0 <- min(y)
baranyi_430 <- nlsLM(y ~ y0 + mmax * (d + (1/mmax) * log(exp(-mmax*d) +
                exp(-mmax * lambda) - exp(-mmax * (d + lambda)))) -
                log(1 + ((exp(mmax * (d + (1/mmax) * log(exp(-mmax*d) +
                exp(-mmax * lambda) - exp(-mmax *
                (d + lambda)))))-1)/(exp(ymax-y0)))),
                start=list(mmax=0.3, ymax=max(y), lambda=60))
summary(baranyi_430)

param_baranyi <- coef(baranyi_430)
coef_baranyi <- c(y0, param_baranyi[1], param_baranyi[3], param_baranyi[2]) # y0, mmax, lambda, ymax
y_fit_baranyi <- baranyi(coef_baranyi, d)

tab_baranyi <- cbind.data.frame(d, y, y_fit_baranyi, sd)
graf_baranyi <- ggplot(tab_baranyi, aes(x=d))+
  geom_point(size=.85, aes(y=y))+
  geom_line(aes(y=y_fit_baranyi))+
  geom_errorbar(aes(ymin=y-sd, ymax=y+sd), width=.3)+
  theme_bw()+
  labs(
    x='Tempo (h)',
    y=expression(paste('DO'[430])))
graf_baranyi

# ********** FITTING GOMPERRTZ *********

gompertz_430 <- nlsLM(y ~ A + C * exp(-exp(-B * (d - M))),
        start=list(A=max(y), B=0.08, C=max(y)-min(y), M=100))
summary(gompertz_430)

param_gompertz <- coef(gompertz_430)
y_fit_gompertz <- gompertz(param_gompertz, d)

tab_gompertz <- cbind.data.frame(d, y, y_fit_gompertz, sd)
graf_gompertz <- ggplot(tab_gompertz, aes(x=d))+
  geom_point(size=.85, aes(y=y))+
  geom_line(aes(y=y_fit_gompertz))+
  geom_errorbar(aes(ymin=y-sd, ymax=y+sd), width=.3)+
  theme_bw()+
  labs(
    x='Tempo (h)',
    y=expression(paste('DO'[430])))
graf_gompertz

# Sem green rust
# ********** FITTING BARANYI *********
dados <- read.table('430nm.tsv', sep='\t', dec=',', header=TRUE, fill=TRUE)
dados$Time <- as.numeric(dados$Time)
dado <- filter(dados, Time>60)
d <- dado$Time
y <- as.numeric(dado$DO)
sd <- as.numeric(dado$sd)
y0 <- min(y)
baranyi_430 <- nlsLM(y ~ y0 + mmax * (d + (1/mmax) * log(exp(-mmax*d) +
                exp(-mmax * lambda) - exp(-mmax * (d + lambda)))) -
                log(1 + ((exp(mmax * (d + (1/mmax) * log(exp(-mmax*d) +
                exp(-mmax * lambda) - exp(-mmax *
                (d + lambda)))))-1)/(exp(ymax-y0)))),
                start=list(mmax=0.03, ymax=max(y), lambda=2))
summary(baranyi_430)

param_baranyi <- coef(baranyi_430)
coef_baranyi <- c(y0, param_baranyi[1], param_baranyi[3], param_baranyi[2]) # y0, mmax, lambda, ymax
y_fit_baranyi <- baranyi(coef_baranyi, d)

tab_baranyi <- cbind.data.frame(d, y, y_fit_baranyi, sd)
graf_baranyi <- ggplot(tab_baranyi, aes(x=d))+
  geom_point(size=.85, aes(y=y))+
  geom_line(aes(y=y_fit_baranyi))+
  geom_errorbar(aes(ymin=y-sd, ymax=y+sd), width=.3)+
  theme_bw()+
  labs(
    x='Tempo (h)',
    y=expression(paste('DO'[430])))
graf_baranyi

# ********** FITTING GOMPERRTZ *********

gompertz_430 <- nlsLM(y ~ A + C * exp(-exp(-B * (d - M))),
        start=list(A=max(y), B=0.08, C=max(y)-min(y), M=100))
summary(gompertz_430)

param_gompertz <- coef(gompertz_430)
y_fit_gompertz <- gompertz(param_gompertz, d)

tab_gompertz <- cbind.data.frame(d, y, y_fit_gompertz, sd)
graf_gompertz <- ggplot(tab_gompertz, aes(x=d))+
  geom_point(size=.85, aes(y=y))+
  geom_line(aes(y=y_fit_gompertz))+
  geom_errorbar(aes(ymin=y-sd, ymax=y+sd), width=.3)+
  theme_bw()+
  labs(
    x='Tempo (h)',
    y=expression(paste('DO'[430])))
graf_gompertz
