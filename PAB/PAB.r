# R 4.3.3
# PARA MESTRADO DE PEDRO ABILHERA - Grupo Quimiosfera

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
install.packages('viridis')
library(viridis)


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



# ********** CURVAS DE CRESCIMENTO *********
curvas <- read.table('curvas_crescimento.tsv', sep='\t', dec=',', header=TRUE, fill=TRUE)

# Curva lâmpada fluorescente
fluo <- filter(curvas, Experimento=='fluo')
fluoA <- filter(fluo, Erlen=='A')
fluoB <- filter(fluo, Erlen=='B')
fluoC <- filter(fluo, Erlen=='C')
fluoD <- filter(fluo, Erlen=='D')
fluoE <- filter(fluo, Erlen=='E')
fluo_mean <- (1/5*fluoA$Avg.DO)+(1/5*fluoB$Avg.DO)+(1/5*fluoC$Avg.DO)+(1/5*fluoD$Avg.DO)+(1/5*fluoE$Avg.DO)
fluo_std <- sqrt((1/25*(fluoA$DO.Std)^2)+(1/25*(fluoB$DO.Std)^2)+(1/25*(fluoC$DO.Std)^2)+(1/25*(fluoD$DO.Std)^2)+(1/25*(fluoE$DO.Std)^2))
d <-fluoA$Dia
y0 <- min(fluo_mean)
y <- fluo_mean
baranyi_fluo <- nlsLM(y ~ y0 + mmax * (d + (1/mmax) * log(exp(-mmax*d) +
                exp(-mmax * lambda) - exp(-mmax * (d + lambda)))) -
                log(1 + ((exp(mmax * (d + (1/mmax) * log(exp(-mmax*d) +
                exp(-mmax * lambda) - exp(-mmax *
                (d + lambda)))))-1)/(exp(ymax-y0)))),
                start=list(mmax=0.1, ymax=max(fluo_mean), lambda=2))
Y0 <- c(y0)

param_fluo <- coef(baranyi_fluo)
coef_fluo <- c(y0, param_fluo[1], param_fluo[3], param_fluo[2]) # y0, mmax, lambda, ymax
y_fit_fluo <- baranyi(coef_fluo, d)

fluo_graf <- cbind.data.frame(d, y, y_fit_fluo, fluo_std, fluoA$Experimento)
graf_fluo <- ggplot(fluo_graf, aes(x=d))+
  geom_point(size=.85, aes(y=y))+
  geom_line(aes(y=y_fit_fluo))+
  geom_errorbar(aes(ymin=y-fluo_std, ymax=y+fluo_std), width=.3)+
  theme_bw()+
  labs(
    x='Tempo (dia)',
    y=expression(paste('DO'[750])))
graf_fluo
ggsave('curva_fluo.png', graf_fluo, device='png', unit='cm', width=15, height=12, dpi=600)

# Curva buffer
buf <- filter(curvas, Experimento=='buffer')
bufA <- filter(buf, Erlen=='A')
bufB <- filter(buf, Erlen=='B')
bufC <- filter(buf, Erlen=='C')
buf_mean <- (1/3*bufA$Avg.DO)+(1/3*bufB$Avg.DO)+(1/3*bufC$Avg.DO)
buf_std <- sqrt((1/9*(bufA$DO.Std)^2)+(1/9*(bufB$DO.Std)^2)+(1/9*(bufC$DO.Std)^2))
d <-bufA$Dia
y0 <- min(buf_mean)
y <- buf_mean
baranyi_buf <- nlsLM(y ~ y0 + mmax * (d + (1/mmax) * log(exp(-mmax*d) +
                exp(-mmax * lambda) - exp(-mmax * (d + lambda)))) -
                log(1 + ((exp(mmax * (d + (1/mmax) * log(exp(-mmax*d) +
                exp(-mmax * lambda) - exp(-mmax *
                (d + lambda)))))-1)/(exp(ymax-y0)))),
                start=list(mmax=0.1, ymax=max(buf_mean), lambda=2))
Y0 <- append(Y0, y0, after=length(Y0))

param_buf <- coef(baranyi_buf)
coef_buf <- c(y0, param_buf[1], param_buf[3], param_buf[2]) # y0, mmax, lambda, ymax
y_fit_buf <- baranyi(coef_buf, d)

buf_graf <- cbind.data.frame(d, y, y_fit_buf, buf_std, bufA$Experimento)
graf_buf <- ggplot(buf_graf, aes(x=d))+
  geom_point(size=.7, aes(y=y))+
  geom_line(aes(y=y_fit_buf))+
  geom_errorbar(aes(ymin=y-buf_std, ymax=y+buf_std), width=.3)+
  theme_bw()+
  labs(
    x='Tempo (dia)',
    y=expression(paste('DO'[750])))
graf_buf
ggsave('curva_buffer.png', graf_buf, device='png', unit='cm', width=15, height=12, dpi=600)

# Curva led
led <- filter(curvas, Experimento=='led')
ledA <- filter(led, Erlen=='A')
ledB <- filter(led, Erlen=='B')
ledC <- filter(led, Erlen=='C')
led_mean <- (1/3*ledA$Avg.DO)+(1/3*ledB$Avg.DO)+(1/3*ledC$Avg.DO)
led_std <- sqrt((1/9*(ledA$DO.Std)^2)+(1/9*(ledB$DO.Std)^2)+(1/9*(ledC$DO.Std)^2))
y0 <- min(led_mean)
y <- led_mean

baranyi_led <- nlsLM(y ~ y0 + mmax * (d + (1/mmax) * log(exp(-mmax*d) +
                exp(-mmax * lambda) - exp(-mmax * (d + lambda)))) -
                log(1 + ((exp(mmax * (d + (1/mmax) * log(exp(-mmax*d) +
                exp(-mmax * lambda) - exp(-mmax *
                (d + lambda)))))-1)/(exp(ymax-y0)))),
                start=list(mmax=0.1, ymax=max(buf_mean), lambda=2))
Y0 <- append(Y0, y0, after=length(Y0))

param_led <- coef(baranyi_led)
coef_led <- c(y0, param_led[1], param_led[3], param_led[2]) # y0, mmax, lambda, ymax
y_fit_led <- baranyi(coef_led, d)

led_graf <- cbind.data.frame(d, y, y_fit_led, led_std, ledA$Experimento)
graf_led <- ggplot(led_graf, aes(x=d))+
  geom_point(size=.7, aes(y=y))+
  geom_line(aes(y=y_fit_led))+
  geom_errorbar(aes(ymin=y-led_std, ymax=y+led_std), width=.3)+
  theme_bw()+
  labs(x='Tempo (dia)',
    y=expression(paste('DO'[750])))
graf_led
ggsave('curva_led.png', graf_led, device='png', unit='cm', width=15, height=12, dpi=600)

# Gráfico conjunto das 3 curvas
colnames(fluo_graf) <- c('d', 'y', 'y_fit', 'std', 'exper')
colnames(buf_graf) <- c('d', 'y', 'y_fit', 'std', 'exper')
colnames(led_graf) <- c('d', 'y', 'y_fit', 'std', 'exper')
tabgraf <- rbind.data.frame(fluo_graf, buf_graf, led_graf)
graf_curvas <- ggplot(tabgraf, aes(x=d))+
  geom_point(aes(y=y, color=exper, shape=exper))+
  geom_line(aes(y=y_fit, color=exper))+
  geom_errorbar(aes(ymin=y-std, ymax=y+std, color=exper), width=.3)+
  theme_bw()+
  labs(x='Tempo (dia)', y=expression(paste('DO'[750])), color='Experimento', shape='Experimento')+
  scale_color_manual(values=c('black', 'gold', 'darkorange'), labels=c('Lâmpada Fluorescente', 'Buffer', 'LED'))+
  scale_shape_manual(values=c(15, 16, 17), labels=c('Lâmpada Fluorescente', 'Buffer', 'LED'))
graf_curvas
ggsave('curvas_crescimento_fit.png', graf_curvas, device='png', unit='cm', width=18, height=12, dpi=600)

# Tabela com os parâmetros do fittings
mmax <- c(summary(baranyi_fluo)$parameters[1,1], summary(baranyi_buf)$parameters[1,1], summary(baranyi_led)$parameters[1,1])
mmax_sd <- c(summary(baranyi_fluo)$parameters[1,2], summary(baranyi_buf)$parameters[1,2], summary(baranyi_led)$parameters[1,2])
lambda <- c(summary(baranyi_fluo)$parameters[3,1], summary(baranyi_buf)$parameters[3,1], summary(baranyi_led)$parameters[3,1])
lambda_sd <- c(summary(baranyi_fluo)$parameters[3,2], summary(baranyi_buf)$parameters[3,2], summary(baranyi_led)$parameters[3,2])
ymax <- c(summary(baranyi_fluo)$parameters[2,1], summary(baranyi_buf)$parameters[2,1], summary(baranyi_led)$parameters[2,1])
ymax_sd <- c(summary(baranyi_fluo)$parameters[2,2], summary(baranyi_buf)$parameters[2,2], summary(baranyi_led)$parameters[2,2])
expmt <- c('fluo', 'buffer', 'led')

params <- cbind.data.frame(expmt, Y0, mmax, mmax_sd, lambda, lambda_sd, ymax, ymax_sd)
write_xlsx(params, 'Parâmetros Baranyi.xlsx')




# ********** EXTRAÇÕES DE CLOROFILA *********
# Triplicata analítica
ext <- read.table('DO-clorofila_barras.tsv', sep='\t', dec=',', header=TRUE, fill=TRUE)
ext$Dia <- as.factor(ext$Dia)
clrf_bar <- ggplot(ext, aes(x=Dia, y=Chl, fill=Extração))+
  geom_bar(position="dodge", stat="identity", color='black')+
  geom_errorbar(aes(ymin=Chl-Chl.std, ymax=Chl+Chl.std), width=.2,position=position_dodge(.9))+
  theme_bw()+
  labs(x='Tempo (dia)', y='Clorofila (\u00b5g/mL)')+
  scale_fill_manual(values=c('gray90', 'gray50', 'gray25'))
clrf_bar
ggsave('clorofila_analitica.png', clrf_bar, device='png', unit='cm', width=15, height=12, dpi=600)


chl.dia <- read.table('clorofila_analitica_box.tsv', sep='\t', dec=',', header=TRUE, fill=TRUE)
dia <- rep(c(6, 12, 18, 24, 30), each=9)
chl <- c(chl.dia$Dia.6, chl.dia$Dia.12, chl.dia$Dia.18, chl.dia$Dia.24, chl.dia$Dia.30)
chl_dia <- cbind.data.frame(dia, chl)
chl_dia$dia <- as.factor(chl_dia$dia)
chl_box <- ggplot(chl_dia, aes(x=dia, y=chl))+
  geom_boxplot(fill='gray90')+
  theme_bw()+
  labs(x='Tempo (dia)', y='Clorofila (\u00b5g/mL)')
chl_box
ggsave('clorofila_analitica_boxplot.png', chl_box, device='png', unit='cm', width=10, height=15, dpi=600)

# Triplicata biológica
DO_chl <- read.table('DO-clorofila.tsv', sep='\t', dec=',', header=TRUE, fill=TRUE)

rl <- lm(Chl.avg ~ DO.avg+0, data=DO_chl)
summary(rl)
summary(rl)$r.squared

graf_dc <- ggplot(DO_chl, aes(x=DO.avg, y=Chl.avg))+
  geom_point(aes(x=DO.avg, y=Chl.avg), size=.75)+
  geom_errorbar(aes(xmin=DO.avg-DO.std, xmax=DO.avg+DO.std), width=0)+
  geom_errorbar(aes( ymin=Chl.avg-Chl.std, ymax=Chl.avg+Chl.std), width=0)+
  geom_line(aes(x=DO.avg,
                y=predict(rl, newdata=DO_chl)))+
  geom_label(label=expression(paste('y = 3,85192x    R'^2~'= 0,9950101')),
              x=1.25, y=7.5)+
  theme_bw()+
  labs(x=expression(paste('DO'[750])), y='Clorofila (\u00b5g/mL)')
graf_dc
ggsave('clorofila_DO.png', graf_dc, device='png', unit='cm', width=18, height=12, dpi=600)


# ********** PESO SECO *********
dw_do <- read.table('PesoSeco_DO.tsv', sep='\t', dec=',', header=TRUE, fill=TRUE)

rl_dw <- lm(DW.avg ~ DO.avg+0, data=dw_chl)
summary(rl_dw)
summary(rl_dw)$r.squared

graf_dw <- ggplot(dw_chl, aes(x=DO.avg, y=DW.avg))+
  geom_point(aes(x=DO.avg, y=DW.avg), size=.75)+
  geom_errorbar(aes(xmin=DO.avg-DO.std, xmax=DO.avg+DO.std), width=0)+
  geom_errorbar(aes( ymin=DW.avg-DW.std, ymax=DW.avg+DW.std), width=0)+
  geom_line(aes(x=DO.avg,
                y=predict(rl_dw, newdata=dw_chl)))+
  geom_label(label=expression(paste('y = 0,477984x     R'^2~'= 0,9941634')),
              x=1, y=1250)+
  theme_bw()+
  labs(x=expression(paste('DO'[750~nm])), y='Peso Seco (g/L)')
graf_dw
ggsave('DW_DO.png', graf_dw, device='png', unit='cm', width=18, height=12, dpi=600)



# ********** ASTROCAM *********
astro <- read.table('astrocam.tsv', sep='\t', dec=',', header=TRUE, fill=TRUE)

ctrl <- filter(astro, Erlen=='Controle')
ctrl1 <- filter(ctrl, Replicata=='1')
ctrl2 <- filter(ctrl, Replicata=='2')
ctrl3 <- filter(ctrl, Replicata=='3')
do.avg_c <- (1/3*ctrl1$DO.avg)+(1/3*ctrl2$DO.avg)+(1/3*ctrl3$DO.avg)
do.std_c <- sqrt((1/9*(ctrl1$DO.std)^2)+(1/9*(ctrl2$DO.std)^2)+(1/9*(ctrl3$DO.std)^2))
chl.avg_c <- (1/3*ctrl1$Chl.avg)+(1/3*ctrl2$Chl.avg)+(1/3*ctrl3$Chl.avg)
chl.std_c <- sqrt((1/9*(ctrl1$Chl.std)^2)+(1/9*(ctrl2$Chl.std)^2)+(1/9*(ctrl3$Chl.std)^2))

d14 <- filter(astro, Erlen=='14d')
d14_1 <- filter(d14, Replicata=='1')
d14_2 <- filter(d14, Replicata=='2')
d14_3 <- filter(d14, Replicata=='3')
do.avg_14 <- (1/3*d14_1$DO.avg)+(1/3*d14_2$DO.avg)+(1/3*d14_3$DO.avg)
do.std_14 <- sqrt((1/9*(d14_1$DO.std)^2)+(1/9*(d14_2$DO.std)^2)+(1/9*(d14_3$DO.std)^2))
chl.avg_14 <- (1/3*d14_1$Chl.avg)+(1/3*d14_2$Chl.avg)+(1/3*d14_3$Chl.avg)
chl.std_14 <- sqrt((1/9*(d14_1$Chl.std)^2)+(1/9*(d14_2$Chl.std)^2)+(1/9*(d14_3$Chl.std)^2))

d28 <- filter(astro, Erlen=='28d')
d28_1 <- filter(d28, Replicata=='1')
d28_2 <- filter(d28, Replicata=='2')
d28_3 <- filter(d28, Replicata=='3')
do.avg_28 <- (1/3*d28_1$DO.avg)+(1/3*d28_2$DO.avg)+(1/3*d28_3$DO.avg)
do.std_28 <- sqrt((1/9*(d28_1$DO.std)^2)+(1/9*(d28_2$DO.std)^2)+(1/9*(d28_3$DO.std)^2))
chl.avg_28 <- (1/3*d28_1$Chl.avg)+(1/3*d28_2$Chl.avg)+(1/3*d28_3$Chl.avg)
chl.std_28 <- sqrt((1/9*(d28_1$Chl.std)^2)+(1/9*(d28_2$Chl.std)^2)+(1/9*(d28_3$Chl.std)^2))

expmt <- c('Controle', 'Controle', 'Controle', 'Astrocam 14 dias', 'Astrocam 14 dias','Astrocam 28 dias','Astrocam 28 dias' )
tempo <- c('0 dia', '14 dias', '28 dias', '0 dia', '14 dias', '0 dia', '28 dias')
do.avg <- c(do.avg_c, do.avg_14, do.avg_28)
do.std <- c(do.std_c, do.std_14, do.std_28)
chl.avg <- c(chl.avg_c, chl.avg_14, chl.avg_28)
chl.std <- c(chl.std_c, chl.std_14, chl.std_28)

tabastro <- cbind.data.frame(expmt, tempo, do.avg, do.std, chl.avg, chl.std)
tabastro$expmt <- factor(tabastro$expmt, levels=c('Controle', 'Astrocam 14 dias','Astrocam 28 dias'))
write_xlsx(tabastro, 'Astrocam.xlsx')

graf_astro_do <- ggplot(tabastro, aes(x=tempo, y=do.avg, shape=expmt, color=expmt))+
  geom_point()+
  geom_line()+
  geom_errorbar(aes(ymin=do.avg-do.std, ymax=do.avg+do.std), width=0.3)+
  theme_bw()+
  labs(x='Tempo (dia)', y=expression(paste('DO'[750])), shape='Experimento', color='Experimento')+
  scale_shape_manual(values=c(15, 16, 17), labels=c('14 dias', '28 dias', 'Controle'))+
  scale_color_manual(values=c('black', 'gold', 'darkorange'), labels=c('14 dias', '28 dias', 'Controle'))+
  xlim(0, 30)
graf_astro_do
ggsave('astrocam_do.png', graf_astro_do, device='png', unit='cm', width=18, height=12, dpi=600)

bar_astro_do <- ggplot(tabastro, aes(x=expmt, y=do.avg, fill=tempo))+
  geom_bar(position='dodge', stat='identity', color='black')+
  geom_errorbar(aes(ymin=do.avg-do.std, ymax=do.avg+do.std), width=.2,position=position_dodge(.9))+
  theme_bw()+
  labs(x='', y=expression(paste('DO'[750])), fill='')+
  scale_fill_manual(values=c('gray90', 'gray50', 'gray25'))
bar_astro_do
ggsave('astrocam_do_barras.png', bar_astro_do, device='png', unit='cm', width=18, height=12, dpi=600)

graf_astro_chl <- ggplot(tabastro, aes(x=tempo, y=chl.avg, shape=expmt, color=expmt))+
  geom_point()+
  geom_line()+
  geom_errorbar(aes(ymin=chl.avg-chl.std, ymax=chl.avg+chl.std), width=0.3)+
  theme_bw()+
  labs(x='Tempo (dia)', y='Clorofila (\u00b5g/mL)', shape='Experimento', color='Experimento')+
  scale_shape_manual(values=c(15, 16, 17), labels=c('14 dias', '28 dias', 'Controle'))+
  scale_color_manual(values=c('black', 'gold', 'darkorange'), labels=c('14 dias', '28 dias', 'Controle'))+
  xlim(0, 30)
graf_astro_chl
ggsave('astrocam_chl.png', graf_astro_chl, device='png', unit='cm', width=18, height=12, dpi=600)

bar_astro_chl <- ggplot(tabastro, aes(x=expmt, y=chl.avg, fill=tempo))+
  geom_bar(position='dodge', stat='identity', color='black')+
  geom_errorbar(aes(ymin=chl.avg-chl.std, ymax=chl.avg+chl.std), width=.2,position=position_dodge(.9))+
  theme_bw()+
  labs(x='', y='Clorofila (\u00b5g/mL)', fill='')+
  scale_fill_manual(values=c('gray90', 'gray50', 'gray25'))
bar_astro_chl
ggsave('astrocam_chl_barras.png', bar_astro_chl, device='png', unit='cm', width=18, height=12, dpi=600)


# ********** LIVE/DEAD *********
# Versão 1
livedead1 <- read.table('LiveDead1.tsv', sep='\t', dec=',', header=TRUE, fill=TRUE)

total1 <- livedead1$Live.avg+livedead1$Dead.avg
total1_sd <- sqrt((1/4*(livedead1$Live.std)^2)+(1/4*(livedead1$Dead.std)^2))
live1 <- livedead1$Live.avg/total1
live1_sd <- sqrt((live1^2)*(((livedead1$Live.std^2)/(livedead1$Live.avg^2))+((total1_sd^2)/(total1^2))))
live1 <- 100*live1
live1_sd <- 100*live1_sd

tab_ld1 <- cbind.data.frame(livedead1$Experimento, live1, live1_sd)
colnames(tab_ld1) <- c('dia', 'live', 'livesd')
tab_ld1$dia <- as.factor(tab_ld1$dia)
tab_ld1$dia <- factor(tab_ld1$dia, levels=c('Controle', '14 dias','28 dias'))

graf_ld1 <- ggplot(tab_ld1, aes(x=dia, y=live))+
  geom_bar(stat='identity', width=.5, color='black', fill='gray80')+
  geom_errorbar(aes(ymin=live-livesd, ymax=live+livesd), width=.3)+
  theme_bw()+
  labs(x='', y='Viabilidade (%)')
graf_ld1
ggsave('livedead_v1.png', graf_ld1, device='png', unit='cm', width=15, height=12, dpi=600)


# Versão 2
livedead2 <- read.table('LiveDead2.tsv', sep='\t', dec=',', header=TRUE, fill=TRUE)

total2 <- livedead2$Live.avg+livedead2$Dead.avg
total2_sd <- sqrt((1/4*(livedead2$Live.std)^2)+(1/4*(livedead2$Dead.std)^2))
live2 <- livedead2$Live.avg/total2
live2_sd <- sqrt((live2^2)*(((livedead2$Live.std^2)/(livedead2$Live.avg^2))+((total2_sd^2)/(total2^2))))
live2 <- 100*live2
live2_sd <- 100*live2_sd

tab_ld2 <- cbind.data.frame(livedead2$Experimento, live2, live2_sd)
colnames(tab_ld2) <- c('dia', 'live', 'livesd')
tab_ld2$dia <- as.factor(tab_ld2$dia)
tab_ld2$dia <- factor(tab_ld2$dia, levels=c('Controle', '14 dias','28 dias', '120 dias'))

graf_ld2 <- ggplot(tab_ld2, aes(x=dia, y=live))+
  geom_bar(stat='identity', width=.5, color='black', fill='gray80')+
  geom_errorbar(aes(ymin=live-livesd, ymax=live+livesd), width=.3)+
  theme_bw()+
  labs(x='', y='Viabilidade (%)')
graf_ld2
ggsave('livedead_v2.png', graf_ld2, device='png', unit='cm', width=15, height=12, dpi=600)



# ********** CRESCIMENTO PÓSD-DESSECAÇÃO *********
dessec <- read.table('crescimento_dessecação.tsv', sep='\t', dec=',', header=TRUE, fill=TRUE)

des_A <- filter(dessec, Erlen=='A')
des_B <- filter(dessec, Erlen=='B')
des_C <- filter(dessec, Erlen=='C')
do <- (1/3*des_A$Avg.DO)+(1/3*des_B$Avg.DO)+(1/3*des_C$Avg.DO)
do_sd <- sqrt((1/9*(des_A$DO.Std^2))+(1/9*(des_B$DO.Std^2))+(1/9*(des_C$DO.Std^2)))

tab_des <- cbind.data.frame(des_A$Experimento, des_A$Dia, do, do_sd)
colnames(tab_des) <- c('exper', 'tempo', 'do', 'do_sd')

graf_des <- ggplot(tab_des, aes(x=tempo, y=do, color=exper))+
  geom_point(aes(shape=exper))+
  geom_line()+
  geom_errorbar(aes(ymin=do-do_sd, ymax=do+do_sd), width=.5)+
  theme_bw()+
  labs(x='Tempo (dia)', y=expression(paste('DO'[750])), color='', shape='')+
  scale_color_manual(values=c('black', 'gold', 'darkorange', 'blue'))+
  scale_shape_manual(values=c(15, 16, 17, 18))
graf_des
ggsave('curvas_crescimento_dessecação.png', graf_des, device='png', unit='cm', width=18, height=12, dpi=600)

# Fit Controle
des_c <- filter(tab_des, exper=='Controle')
d <- des_c$tempo
y0 <- min(des_c$do)
y <- des_c$do

baranyi_c <- nlsLM(y ~ y0 + mmax * (d + (1/mmax) * log(exp(-mmax*d) +
                exp(-mmax * lambda) - exp(-mmax * (d + lambda)))) -
                log(1 + ((exp(mmax * (d + (1/mmax) * log(exp(-mmax*d) +
                exp(-mmax * lambda) - exp(-mmax *
                (d + lambda)))))-1)/(exp(ymax-y0)))),
                start=list(mmax=0.1, ymax=max(y), lambda=2))
summary(baranyi_c)
Y0 <- c(y0)
param_c <- coef(baranyi_c)
coef_c <- c(y0, param_c[1], param_c[3], param_c[2]) # y0, mmax, lambda, ymax
y_fit_c <- baranyi(coef_c, d)

# Fit 14 dias
des_14 <- filter(tab_des, exper=='14 dias')
d <- des_14$tempo
y0 <- min(des_14$do)
y <- des_14$do

baranyi_14 <- nlsLM(y ~ y0 + mmax * (d + (1/mmax) * log(exp(-mmax*d) +
                exp(-mmax * lambda) - exp(-mmax * (d + lambda)))) -
                log(1 + ((exp(mmax * (d + (1/mmax) * log(exp(-mmax*d) +
                exp(-mmax * lambda) - exp(-mmax *
                (d + lambda)))))-1)/(exp(ymax-y0)))),
                start=list(mmax=0.1, ymax=max(y), lambda=2))
summary(baranyi_14)
Y0 <- append(Y0, y0, after=length(Y0))
param_14 <- coef(baranyi_14)
coef_14 <- c(y0, param_14[1], param_14[3], param_14[2]) # y0, mmax, lambda, ymax
y_fit_14 <- baranyi(coef_14, d)

# Fit 28 dias
des_28 <- filter(tab_des, exper=='28 dias')
d <- des_28$tempo
y0 <- min(des_28$do)
y <- des_28$do

baranyi_28 <- nlsLM(y ~ y0 + mmax * (d + (1/mmax) * log(exp(-mmax*d) +
                exp(-mmax * lambda) - exp(-mmax * (d + lambda)))) -
                log(1 + ((exp(mmax * (d + (1/mmax) * log(exp(-mmax*d) +
                exp(-mmax * lambda) - exp(-mmax *
                (d + lambda)))))-1)/(exp(ymax-y0)))),
                start=list(mmax=0.1, ymax=max(y), lambda=2))
summary(baranyi_28)
Y0 <- append(Y0, y0, after=length(Y0))
param_28 <- coef(baranyi_28)
coef_28 <- c(y0, param_28[1], param_28[3], param_28[2]) # y0, mmax, lambda, ymax
y_fit_28 <- baranyi(coef_28, d)

# Fit 120 dias
des_120 <- filter(tab_des, exper=='120 dias')
d <- des_120$tempo
y0 <- min(des_120$do)
y <- des_120$do

baranyi_120 <- nlsLM(y ~ y0 + mmax * (d + (1/mmax) * log(exp(-mmax*d) +
                exp(-mmax * lambda) - exp(-mmax * (d + lambda)))) -
                log(1 + ((exp(mmax * (d + (1/mmax) * log(exp(-mmax*d) +
                exp(-mmax * lambda) - exp(-mmax *
                (d + lambda)))))-1)/(exp(ymax-y0)))),
                start=list(mmax=0.1, ymax=max(y), lambda=2))
summary(baranyi_120)
Y0 <- append(Y0, y0, after=length(Y0))
param_120 <- coef(baranyi_120)
coef_120 <- c(y0, param_120[1], param_120[3], param_120[2]) # y0, mmax, lambda, ymax
y_fit_120 <- baranyi(coef_120, d)


y_fit <- c(y_fit_c, y_fit_14, y_fit_28, y_fit_120)
tab_des <- cbind.data.frame(tab_des, y_fit)

graf_des_fit <- ggplot(tab_des, aes(x=tempo, y=do, color=exper))+
  geom_point(aes(shape=exper))+
  geom_line(aes(y=y_fit))+
  geom_errorbar(aes(ymin=do-do_sd, ymax=do+do_sd), width=.5)+
  theme_bw()+
  labs(x='Tempo (dia)', y=expression(paste('DO'[750])), color='', shape='')+
  scale_color_manual(values=c('black', 'gold', 'darkorange', 'blue'))+
  scale_shape_manual(values=c(15, 16, 17, 18))
graf_des_fit
ggsave('curvas_crescimento_dessecação_fitado.png', graf_des_fit, device='png', unit='cm', width=18, height=12, dpi=600)

mmax <- c(summary(baranyi_c)$parameters[1,1], summary(baranyi_14)$parameters[1,1], summary(baranyi_28)$parameters[1,1], summary(baranyi_120)$parameters[1,1])
mmax_sd <- c(summary(baranyi_c)$parameters[1,2], summary(baranyi_14)$parameters[1,2], summary(baranyi_28)$parameters[1,2], summary(baranyi_120)$parameters[1,2])
lambda <- c(summary(baranyi_c)$parameters[3,1], summary(baranyi_14)$parameters[3,1], summary(baranyi_28)$parameters[3,1], summary(baranyi_120)$parameters[3,1])
lambda_sd <- c(summary(baranyi_c)$parameters[3,2], summary(baranyi_14)$parameters[3,2], summary(baranyi_28)$parameters[3,2], summary(baranyi_120)$parameters[3,2])
ymax <- c(summary(baranyi_c)$parameters[2,1], summary(baranyi_14)$parameters[2,1], summary(baranyi_28)$parameters[2,1], summary(baranyi_120)$parameters[2,1])
ymax_sd <- c(summary(baranyi_c)$parameters[2,2], summary(baranyi_14)$parameters[2,2], summary(baranyi_28)$parameters[2,2], summary(baranyi_120)$parameters[2,2])
expmt <- c('Controle', '14 dias', '28 dias', '120 dias')

params <- cbind.data.frame(expmt, Y0, mmax, mmax_sd, lambda, lambda_sd, ymax, ymax_sd)
write_xlsx(params, 'Parâmetros Baranyi Dessecação.xlsx')




# ********** AUTOFLUORESCÊNCIA DA CLOROFILA *********
# Versão 1
autof <- read.table('autofluorescencia.tsv', sep='\t', dec=',', header=TRUE, fill=TRUE)

autof1 <- filter(autof, Experimento=='1')
autof2 <- filter(autof, Experimento=='2')
autof3 <- filter(autof, Experimento=='3')
autof4 <- filter(autof, Experimento=='4')
autof5 <- filter(autof, Experimento=='5')
autof6 <- filter(autof, Experimento=='6')

autof14_avg <- (1/3*autof1$Autof.avg)+(1/3*autof2$Autof.avg)+(1/3*autof3$Autof.avg)
autof14_sd <- sqrt(((1/3*autof1$Autof.std)^2)+((1/3*autof2$Autof.std)^2)+((1/3*autof3$Autof.std)^2))
autof28_avg <- (1/3*autof5$Autof.avg)+(1/3*autof6$Autof.avg)+(1/3*autof4$Autof.avg)
autof28_sd <- sqrt(((1/3*autof4$Autof.std)^2)+((1/3*autof5$Autof.std)^2)+((1/3*autof6$Autof.std)^2))
autof0_avg <- (1/2*autof14_avg[1])+(1/2*autof28_avg[1])
autof0_sd <- sqrt(((1/2*autof14_sd[1])^2)+((1/2*autof28_sd[1])^2))

autof_avg <- c(autof0_avg, autof14_avg[-1], autof28_avg[-1])
autof_avg.p <- (autof_avg/autof0_avg)*100
autof_sd <- c(autof0_sd, autof14_sd[-1], autof28_sd[-1])
autof_sd.p <- sqrt(((autof_sd/autof_avg)^2)+((autof0_sd/autof0_avg)^2))*100

dia <- c(0, 14, 14, 14, 28, 28, 28)
tempo <- c('Controle', '0', '40', '180', '0', '40', '180')
tab_autof <- cbind.data.frame(dia, tempo, autof_avg.p, autof_sd.p)
tab_autof$tempo <- factor(tab_autof$tempo, levels=c('Controle', '0', '40', '180'))
colnames(tab_autof) <- c('dia', 'tempo', 'autof', 'sd')
tab_autof$dia <- as.factor(tab_autof$dia)
tab_autof$tempo <- as.factor(tab_autof$tempo)

graf_autof <- ggplot(tab_autof, aes(x=dia, y=autof, fill=tempo))+
  geom_bar(position='dodge', stat='identity', color='black')+
  geom_errorbar(aes(ymin=autof-sd, ymax=autof+sd), width=.2,position=position_dodge(.9))+
  theme_bw()+
  labs(x='Tempo (dia)', y='Autofluorescência da clorofila (%)', fill='Tempo (min)')+
  scale_fill_manual(values=c('gray95', 'gray75', 'gray50', 'gray25'))
graf_autof

ggsave('autofluorescencia_v1.png', graf_autof, device='png', unit='cm', width=18, height=12, dpi=600)


# Versão 2
autof2 <- read.table('autofluorescencia2.tsv', sep='\t', dec=',', header=TRUE, fill=TRUE)

autof2_1 <- filter(autof2, Experimento=='1')
autof2_2 <- filter(autof2, Experimento=='2')
autof2_3 <- filter(autof2, Experimento=='3')
autof2_4 <- filter(autof2, Experimento=='4')
autof2_5 <- filter(autof2, Experimento=='5')
autof2_6 <- filter(autof2, Experimento=='6')
autof2_7 <- filter(autof2, Experimento=='7')
autof2_8 <- filter(autof2, Experimento=='8')
autof2_9 <- filter(autof2, Experimento=='9')

autof14_avg <- (1/3*autof2_1$Autof.avg)+(1/3*autof2_2$Autof.avg)+(1/3*autof2_3$Autof.avg)
autof14_sd <- sqrt(((1/3*autof2_1$Autof.std)^2)+((1/3*autof2_2$Autof.std)^2)+((1/3*autof2_3$Autof.std)^2))
autof28_avg <- (1/3*autof2_4$Autof.avg)+(1/3*autof2_5$Autof.avg)+(1/3*autof2_6$Autof.avg)
autof28_sd <- sqrt(((1/3*autof2_4$Autof.std)^2)+((1/3*autof2_5$Autof.std)^2)+((1/3*autof2_6$Autof.std)^2))
autof120_avg <- (1/3*autof2_7$Autof.avg)+(1/3*autof2_8$Autof.avg)+(1/3*autof2_9$Autof.avg)
autof120_sd <- sqrt(((1/3*autof2_7$Autof.std)^2)+((1/3*autof2_8$Autof.std)^2)+((1/3*autof2_9$Autof.std)^2))

autof0_avg <- (1/3*autof14_avg[1])+(1/3*autof28_avg[1])+(1/3*autof120_avg[1])
autof0_sd <- sqrt(((1/3*autof14_sd[1])^2)+((1/3*autof28_sd[1])^2)+((1/3*autof120_sd[1])^2))

autof_avg <- c(autof0_avg, autof14_avg[-1], autof28_avg[-1], autof120_avg[-1])
autof_avg.p <- (autof_avg/autof0_avg)*100

autof_sd <- c(autof0_sd, autof14_sd[-1], autof28_sd[-1], autof120_sd[-1])
autof_sd.p <- sqrt(((autof_sd/autof_avg)^2)+((autof0_sd/autof0_avg)^2))*100

dia <- c(0, 14, 14, 28, 28, 120, 120)
tempo <- c(0, 40, 180, 40, 180, 40, 180)
tab_autof <- cbind.data.frame(dia, tempo, autof_avg.p, autof_sd.p)
colnames(tab_autof) <- c('dia', 'tempo', 'autof', 'sd')
tab_autof$dia <- as.factor(tab_autof$dia)
tab_autof$tempo <- as.factor(tab_autof$tempo)

graf_autof2 <- ggplot(tab_autof, aes(x=dia, y=autof, fill=tempo))+
  geom_bar(position='dodge', stat='identity', color='black')+
  geom_errorbar(aes(ymin=autof-sd, ymax=autof+sd), width=.2,position=position_dodge(.9))+
  theme_bw()+
  labs(x='Tempo (dia)', y='Autofluorescência da clorofila (%)', fill='Tempo (min)')+
  scale_fill_manual(values=c('gray90', 'gray50', 'gray25'))
graf_autof2

ggsave('autofluorescencia2.png', graf_autof2, device='png', unit='cm', width=18, height=12, dpi=600)


# ********** REGOLITO *********
regtab <- read.table('regolito.tsv', sep='\t', dec=',', header=TRUE, fill=TRUE)
reg1 <- filter(regtab, Replicata=='1')
reg2 <- filter(regtab, Replicata=='2')
reg3 <- filter(regtab, Replicata=='3')
dia <- reg1$Dia
erlen <- reg1$Erlen
chl.avg <- (1/3*reg1$Chl.avg)+(1/3*reg2$Chl.avg)+(1/3*reg3$Chl.avg)
chl.sd <- sqrt(((1/3*reg1$Chl.std)^2)+((1/3*reg2$Chl.std)^2)+((1/3*reg3$Chl.std)^2))
reg <- cbind.data.frame(erlen, dia, chl.avg, chl.sd)

grafreg <- ggplot(reg, aes(x=dia, y=chl.avg, color=erlen))+
  geom_line(linewidth=.5)+
  geom_point(aes(shape=erlen))+
  geom_errorbar(aes(ymin=chl.avg-chl.sd, ymax=chl.avg+chl.sd), width=.5)+
  theme_bw()+
  labs(x='Tempo (dia)', y='Clorofila (\u00b5g/mL)', color='', shape='')+
  scale_color_manual(values=c('black', 'gold', 'darkorange', 'blue', 'dimgray'))+
  scale_shape_manual(values=c(15, 16, 17, 18, 11))
grafreg
ggsave('regolito.png', grafreg, device='png', unit='cm', width=18, height=12, dpi=600)

reg$dia <- as.factor(reg$dia)
exper <- c('BG-11', 'BG-11', 'BG-11', 'BG-11 + MGS-1', 'BG-11 + MGS-1', 'BG-11 + MGS-1', 'BG-11 + MGS-1', 'Água + MGS-1', 'Água + MGS-1', 'Água + MGS-1', 'Água + MGS-1')
dias <- c('0 dia', '14 dias', '28 dias', '0 dia', '14 dias', '0 dia', '28 dias', '0 dia', '14 dias', '0 dia', '28 dias')
reg <- cbind.data.frame(reg, exper, dias)
reg$exper <- factor(reg$exper, levels=c('BG-11', 'BG-11 + MGS-1', 'Água + MGS-1'))
grafreg_bar <- ggplot(reg, aes(x=exper, y=chl.avg, fill=dias))+
  geom_bar(position='dodge', stat='identity', color='black')+
  geom_errorbar(aes(ymin=chl.avg-chl.sd, ymax=chl.avg+chl.sd), width=.2,position=position_dodge(.9))+
  theme_bw()+
  labs(x='', y='Clorofila (\u00b5g/mL)', fill='')+
  scale_fill_manual(values=c('gray90', 'gray50', 'gray25'))
grafreg_bar
ggsave('regolito_barras.png', grafreg_bar, device='png', unit='cm', width=18, height=12, dpi=600)

write_xlsx(reg, 'Regolito.xlsx')

# ********** PERCLORATO *********
perctab <- read.table('perclorato.tsv', sep='\t', dec=',', header=TRUE, fill=TRUE)
Y0 <-c()

#Fit controle
perc_0 <- filter(perctab, Experimento=='BG-11')
d <- perc_0$Dia
y0 <- min(perc_0$Avg.DO)
y <- perc_0$Avg.DO

baranyi_0 <- nlsLM(y ~ y0 + mmax * (d + (1/mmax) * log(exp(-mmax*d) +
                exp(-mmax * lambda) - exp(-mmax * (d + lambda)))) -
                log(1 + ((exp(mmax * (d + (1/mmax) * log(exp(-mmax*d) +
                exp(-mmax * lambda) - exp(-mmax *
                (d + lambda)))))-1)/(exp(ymax-y0)))),
                start=list(mmax=0.1, ymax=max(y), lambda=2))
summary(baranyi_0)
Y0 <- c(y0)
param_0 <- coef(baranyi_0)
coef_0 <- c(y0, param_0[1], param_0[3], param_0[2]) # y0, mmax, lambda, ymax
y_fit_0 <- baranyi(coef_0, d)

#Fit 0,05%
perc_005 <- filter(perctab, Experimento=='0,05%')
d <- perc_005$Dia
y0 <- min(perc_005$Avg.DO)
y <- perc_005$Avg.DO

baranyi_005 <- nlsLM(y ~ y0 + mmax * (d + (1/mmax) * log(exp(-mmax*d) +
                exp(-mmax * lambda) - exp(-mmax * (d + lambda)))) -
                log(1 + ((exp(mmax * (d + (1/mmax) * log(exp(-mmax*d) +
                exp(-mmax * lambda) - exp(-mmax *
                (d + lambda)))))-1)/(exp(ymax-y0)))),
                start=list(mmax=0.1, ymax=max(y), lambda=2))
summary(baranyi_005)
Y0 <- append(Y0, y0, after=length(Y0))
param_005 <- coef(baranyi_005)
coef_005 <- c(y0, param_005[1], param_005[3], param_005[2]) # y0, mmax, lambda, ymax
y_fit_005 <- baranyi(coef_005, d)

#Fit 0,20%
perc_02 <- filter(perctab, Experimento=='0,20%')
d <- perc_02$Dia
y0 <- min(perc_02$Avg.DO)
y <- perc_02$Avg.DO

baranyi_02 <- nlsLM(y ~ y0 + mmax * (d + (1/mmax) * log(exp(-mmax*d) +
                exp(-mmax * lambda) - exp(-mmax * (d + lambda)))) -
                log(1 + ((exp(mmax * (d + (1/mmax) * log(exp(-mmax*d) +
                exp(-mmax * lambda) - exp(-mmax *
                (d + lambda)))))-1)/(exp(ymax-y0)))),
                start=list(mmax=0.1, ymax=max(y), lambda=2))
summary(baranyi_02)
Y0 <- append(Y0, y0, after=length(Y0))
param_02 <- coef(baranyi_02)
coef_02 <- c(y0, param_02[1], param_02[3], param_02[2]) # y0, mmax, lambda, ymax
y_fit_02 <- baranyi(coef_02, d)

#Fit 0,50%
perc_05 <- filter(perctab, Experimento=='0,50%')
d <- perc_05$Dia
y0 <- min(perc_05$Avg.DO)
y <- perc_05$Avg.DO

baranyi_05 <- nlsLM(y ~ y0 + mmax * (d + (1/mmax) * log(exp(-mmax*d) +
                exp(-mmax * lambda) - exp(-mmax * (d + lambda)))) -
                log(1 + ((exp(mmax * (d + (1/mmax) * log(exp(-mmax*d) +
                exp(-mmax * lambda) - exp(-mmax *
                (d + lambda)))))-1)/(exp(ymax-y0)))),
                start=list(mmax=0.1, ymax=max(y), lambda=2))
summary(baranyi_05)
Y0 <- append(Y0, y0, after=length(Y0))
param_05 <- coef(baranyi_05)
coef_05 <- c(y0, param_05[1], param_05[3], param_05[2]) # y0, mmax, lambda, ymax
y_fit_05 <- baranyi(coef_05, d)

#Fit 1%
perc_1 <- filter(perctab, Experimento=='1%')
d <- perc_1$Dia
y0 <- min(perc_1$Avg.DO)
y <- perc_1$Avg.DO

baranyi_1 <- nlsLM(y ~ y0 + mmax * (d + (1/mmax) * log(exp(-mmax*d) +
                exp(-mmax * lambda) - exp(-mmax * (d + lambda)))) -
                log(1 + ((exp(mmax * (d + (1/mmax) * log(exp(-mmax*d) +
                exp(-mmax * lambda) - exp(-mmax *
                (d + lambda)))))-1)/(exp(ymax-y0)))),
                start=list(mmax=0.1, ymax=max(y), lambda=2))
summary(baranyi_1)
Y0 <- append(Y0, y0, after=length(Y0))
param_1 <- coef(baranyi_1)
coef_1 <- c(y0, param_1[1], param_1[3], param_1[2]) # y0, mmax, lambda, ymax
y_fit_1 <- baranyi(coef_1, d)

#Fit 5%
perc_5 <- filter(perctab, Experimento=='5%')
d <- perc_5$Dia
y0 <- min(perc_5$Avg.DO)
y <- perc_5$Avg.DO
ymax <- max(perc_5$Avg.DO)

baranyi_5 <- nlsLM(y ~ y0 + mmax * (d + (1/mmax) * log(exp(-mmax*d) +
                exp(-mmax * lambda) - exp(-mmax * (d + lambda)))) -
                log(1 + ((exp(mmax * (d + (1/mmax) * log(exp(-mmax*d) +
                exp(-mmax * lambda) - exp(-mmax *
                (d + lambda)))))-1)/(exp(ymax-y0)))),
                start=list(mmax=0.1, lambda=2))
summary(baranyi_5)
Y0 <- append(Y0, y0, after=length(Y0))
param_5 <- coef(baranyi_5)
coef_5 <- c(y0, param_5[1], param_5[3], ymax) # y0, mmax, lambda, ymax
y_fit_5 <- baranyi(coef_5, d)

# Gráfico
y.fit <- c(y_fit_0, y_fit_005, y_fit_02, y_fit_05, y_fit_1, y_fit_5)
exper <- rep(c('0,00%', '0,05%', '0,20%', '0,50%', '1,00%', '5,00%'), each=11)
perc <- cbind.data.frame(perctab, y.fit, exper)


grafperc <- ggplot(perc, aes(x=Dia, y=Avg.DO, color=exper))+
  geom_point(aes(shape=exper))+
  geom_line(aes(y=y.fit))+
  geom_errorbar(aes(ymin=Avg.DO-DO.Std, ymax=Avg.DO+DO.Std), width=.5)+
  theme_bw()+
  labs(x='Tempo (dia)', y=expression(paste('DO'[750])), color='[Perclorato]', shape='[Perclorato]')+
  scale_color_manual(values=c('black', 'gold', 'darkorange', 'blue', 'dimgray', 'tomato'))+
  scale_shape_manual(values=c(15, 16, 17, 18, 11, 7))
grafperc
ggsave('curvas_perclorato_fitado.png', grafperc, device='png', unit='cm', width=18, height=12, dpi=600)

mmax <- c(summary(baranyi_0)$parameters[1,1], summary(baranyi_005)$parameters[1,1], summary(baranyi_02)$parameters[1,1], summary(baranyi_05)$parameters[1,1], summary(baranyi_1)$parameters[1,1])
mmax_sd <- c(summary(baranyi_0)$parameters[1,2], summary(baranyi_005)$parameters[1,2], summary(baranyi_02)$parameters[1,2], summary(baranyi_05)$parameters[1,2], summary(baranyi_1)$parameters[1,2])
lambda <- c(summary(baranyi_0)$parameters[3,1], summary(baranyi_005)$parameters[3,1], summary(baranyi_02)$parameters[3,1], summary(baranyi_05)$parameters[3,1], summary(baranyi_1)$parameters[3,1])
lambda_sd <- c(summary(baranyi_0)$parameters[3,2], summary(baranyi_005)$parameters[3,2], summary(baranyi_02)$parameters[3,2], summary(baranyi_05)$parameters[3,2], summary(baranyi_1)$parameters[3,2])
ymax <- c(summary(baranyi_0)$parameters[2,1], summary(baranyi_005)$parameters[2,1], summary(baranyi_02)$parameters[2,1], summary(baranyi_05)$parameters[2,1], summary(baranyi_1)$parameters[2,1])
ymax_sd <- c(summary(baranyi_0)$parameters[2,2], summary(baranyi_005)$parameters[2,2], summary(baranyi_02)$parameters[2,2], summary(baranyi_05)$parameters[2,2], summary(baranyi_1)$parameters[2,2])
expm <- c('0,00%', '0,05%', '0,20%', '0,50%', '1,00%')

paramsperc <- cbind.data.frame(expm, Y0[-6], mmax, mmax_sd, lambda, lambda_sd, ymax, ymax_sd)

write_xlsx(paramsperc, 'Parâmetros Baranyi Perclorato.xlsx')


# ********** SEM NITROGÊNIO *********
ntab <- read.table('semn.tsv', sep='\t', dec=',', header=TRUE, fill=TRUE)
n.a <- filter(ntab, Erlen=='A')
n.b <- filter(ntab, Erlen=='B')
n.c <- filter(ntab, Erlen=='C')
dia <- n.a$Dia
exper <- n.a$Experimento
do.avg <- (1/3*n.a$Avg.DO)+(1/3*n.b$Avg.DO)+(1/3*n.c$Avg.DO)
do.sd <- sqrt(((1/3*n.a$DO.Std)^2)+((1/3*n.b$DO.Std)^2)+((1/3*n.c$DO.Std)^2))
tabn <- cbind.data.frame(exper, dia, do.avg, do.sd)

Y0 <-c()

#Fit 1x
n1x <- filter(tabn, exper=='1x')
d <- n1x$dia
y0 <- min(n1x$do.avg)
y <- n1x$do.avg

baranyi_1x <- nlsLM(y ~ y0 + mmax * (d + (1/mmax) * log(exp(-mmax*d) +
                exp(-mmax * lambda) - exp(-mmax * (d + lambda)))) -
                log(1 + ((exp(mmax * (d + (1/mmax) * log(exp(-mmax*d) +
                exp(-mmax * lambda) - exp(-mmax *
                (d + lambda)))))-1)/(exp(ymax-y0)))),
                start=list(mmax=0.1, ymax=max(y), lambda=2))
summary(baranyi_1x)
Y0 <- append(Y0, y0, after=length(Y0))
param_1x <- coef(baranyi_1x)
coef_1x <- c(y0, param_1x[1], param_1x[3], param_1x[2]) # y0, mmax, lambda, ymax
y_fit_1x <- baranyi(coef_1x, d)

#Fit 05x
n05x <- filter(tabn, exper=='0,5x')
d <- n05x$dia
y0 <- min(n05x$do.avg)
y <- n05x$do.avg

baranyi_05x <- nlsLM(y ~ y0 + mmax * (d + (1/mmax) * log(exp(-mmax*d) +
                exp(-mmax * lambda) - exp(-mmax * (d + lambda)))) -
                log(1 + ((exp(mmax * (d + (1/mmax) * log(exp(-mmax*d) +
                exp(-mmax * lambda) - exp(-mmax *
                (d + lambda)))))-1)/(exp(ymax-y0)))),
                start=list(mmax=0.1, ymax=max(y), lambda=2))
summary(baranyi_05x)
Y0 <- append(Y0, y0, after=length(Y0))
param_05x <- coef(baranyi_05x)
coef_05x <- c(y0, param_05x[1], param_05x[3], param_05x[2]) # y0, mmax, lambda, ymax
y_fit_05x <- baranyi(coef_05x, d)

#Fit 0x
n0x <- filter(tabn, exper=='0x')
d <- n0x$dia
y0 <- min(n0x$do.avg)
y <- n0x$do.avg

baranyi_0x <- nlsLM(y ~ y0 + mmax * (d + (1/mmax) * log(exp(-mmax*d) +
                exp(-mmax * lambda) - exp(-mmax * (d + lambda)))) -
                log(1 + ((exp(mmax * (d + (1/mmax) * log(exp(-mmax*d) +
                exp(-mmax * lambda) - exp(-mmax *
                (d + lambda)))))-1)/(exp(ymax-y0)))),
                start=list(mmax=0.1, ymax=max(y), lambda=2))
summary(baranyi_0x)
Y0 <- append(Y0, y0, after=length(Y0))
param_0x <- coef(baranyi_0x)
coef_0x <- c(y0, param_0x[1], param_0x[3], param_0x[2]) # y0, mmax, lambda, ymax
y_fit_0x <- baranyi(coef_0x, d)

# Gráfico
y.fit <- c(y_fit_1x, y_fit_05x, y_fit_0x)
tabn <- cbind.data.frame(tabn, y.fit)
tabn$exper <- factor(tabn$exper, levels=c('1x', '0,5x', '0x'))
grafn <- ggplot(tabn, aes(x=dia, y=do.avg, color=exper))+
  geom_point(aes(shape=exper))+
  geom_line(aes(y=y.fit))+
  geom_errorbar(aes(ymin=do.avg-do.sd, ymax=do.avg+do.sd), width=.5)+
  theme_bw()+
  labs(x='Tempo (dia)', y=expression(paste('DO'[750])), color=expression(paste('[NaNO'[3]~']')), shape=expression(paste('[NaNO'[3]~']')))+
  scale_color_manual(values=c('black', 'gold', 'darkorange'))+
  scale_shape_manual(values=c(15, 16, 17))
grafn
ggsave('curvas_semn_fitado.png', grafn, device='png', unit='cm', width=18, height=12, dpi=600)

mmax <- c(summary(baranyi_1x)$parameters[1,1], summary(baranyi_05x)$parameters[1,1], summary(baranyi_0x)$parameters[1,1])
mmax_sd <- c(summary(baranyi_1x)$parameters[1,2], summary(baranyi_05x)$parameters[1,2], summary(baranyi_0x)$parameters[1,2])
lambda <- c(summary(baranyi_1x)$parameters[3,1], summary(baranyi_05x)$parameters[3,1], summary(baranyi_0x)$parameters[3,1])
lambda_sd <- c(summary(baranyi_1x)$parameters[3,2], summary(baranyi_05x)$parameters[3,2], summary(baranyi_0x)$parameters[3,2])
ymax <- c(summary(baranyi_1x)$parameters[2,1], summary(baranyi_05x)$parameters[2,1], summary(baranyi_0x)$parameters[2,1])
ymax_sd <- c(summary(baranyi_1x)$parameters[2,2], summary(baranyi_05x)$parameters[2,2], summary(baranyi_0x)$parameters[2,2])
expm <- c('1x', '0,5x', '0x')

paramsN <- cbind.data.frame(expm, Y0, mmax, mmax_sd, lambda, lambda_sd, ymax, ymax_sd)
write_xlsx(paramsN, 'Parâmetros Baranyi Sem N.xlsx')
