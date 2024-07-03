# R 4.2.0
# PARA DOUTORADO DE ANA PAULA MUCHE SCHIAVO - Grupo Quimiosfera
# Scrip para gráficos e análises de dados do experimento E. coli extremófila

# REGRAS DE ENTRADA: tabela chamada "cepa_replicata.tsv", vírgula para decimal, tab para separação de células
# ORDEM DAS COLUNAS: 1. Ponto | 2. Concentração de NaCl | 3. Horas de experimento | 4-6. Contagem de UFC | 7. Fator de correção da diluição | 8. Média corrigida

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
#install.packages('growthrates')
#library(growthrates)

# ********** FUNÇÕES **********
baranyi <- function(params, x) {
  params[1] + params[2] * (x + (1/params[2]) * log(exp(-params[2]*x) +
  exp(-params[2] * params[3]) - exp(-params[2] * (x + params[3])))) -
  log(1 + ((exp(params[2] * (x + (1/params[2]) * log(exp(-params[2]*x) +
  exp(-params[2] * params[3]) - exp(-params[2] * (x + params[3])))))-1)/
    (exp(params[4]-params[1]))))
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

calc_mse <- function(data, teor){   #data é vetor com os dados experimentais e teor é vetor com y teórico pros mesmos valores de x
  n <- length(data)
  mse <- (1/n)*(sum((data-teor)^2))
  return(mse)
}

# ********** LEITURA E PRÉ-PROCESSAMENTO **********

# Lê tabela em TSV, vírgula como decimal e primeira linha como cabeçalho
ms1655_1 <- read.table('MS1655_1.tsv', sep='\t', dec=',', header=TRUE, fill=TRUE)
ms1655_2 <- read.table('MS1655_2.tsv', sep='\t', dec=',', header=TRUE, fill=TRUE)
ms1655_3 <- read.table('MS1655_3.tsv', sep='\t', dec=',', header=TRUE, fill=TRUE)

ms1655_1$Cor.Mean <- ms1655_1$Cor.Mean*100    # Corrige a média corrigida pela diluição para a densidade de UFC seja por mL
ms1655_2$Cor.Mean <- ms1655_2$Cor.Mean*100
ms1655_3$Cor.Mean <- ms1655_3$Cor.Mean*100

ms1655_ufc_1 <- as.matrix(ms1655_1[,c(4,5,6)])    # Cria matriz só com com as 3 colunas de contagem das gotas
ms1655_ufc_2 <- as.matrix(ms1655_2[,c(4,5,6)])
ms1655_ufc_3 <- as.matrix(ms1655_3[,c(4,5,6)])

media_ms1655_1 <- rowMeans(ms1655_ufc_1, na.rm=TRUE)    # Faz média das 3 colunas com valores contados de UFC
media_ms1655_2 <- rowMeans(ms1655_ufc_2, na.rm=TRUE)
media_ms1655_3 <- rowMeans(ms1655_ufc_3, na.rm=TRUE)

sd_ms1655_1 <- rowSds(ms1655_ufc_1, na.rm=TRUE)   # Calcula desvio padrão entre as gostas de uma replicata
sd_ms1655_2 <- rowSds(ms1655_ufc_2, na.rm=TRUE)
sd_ms1655_3 <- rowSds(ms1655_ufc_3, na.rm=TRUE)

n_ms1655_1 <- rowSums(!is.na(ms1655_ufc_1))   # Calcula número de gostas contadas. Para isso ele conta quantos valores não NA existem naquelas 3 colunas
n_ms1655_2 <- rowSums(!is.na(ms1655_ufc_2))
n_ms1655_3 <- rowSums(!is.na(ms1655_ufc_3))

se_ms1655_1 <- sd_ms1655_1/(sqrt(n_ms1655_1))   # Calcula erro padrão (por definição: desvio padrão dividido por raiz quadrada do n)
se_ms1655_2 <- sd_ms1655_2/(sqrt(n_ms1655_2))
se_ms1655_3 <- sd_ms1655_3/(sqrt(n_ms1655_3))

erro_ms1655_1 <- ((ms1655_1$Dilution*100))*(se_ms1655_1)   # Calcula erro associado a média corrigida das três gotinhas. Propagação do erro padrão considerando Normal.
erro_ms1655_2 <- ((ms1655_2$Dilution*100))*(se_ms1655_2)
erro_ms1655_3 <- ((ms1655_3$Dilution*100))*(se_ms1655_3)



# Para gráfico total (dados reais, média entre  replicatas)
medias_ms1655 <- cbind(ms1655_1$Cor.Mean, ms1655_2$Cor.Mean, ms1655_3$Cor.Mean)   # Junta as médias corrigas de cada replicata numa tabela, cada replicata como uma coluna
media_total_ms1655 <- rowMeans(medias_ms1655, na.rm=TRUE)    # Cria vetor da média entre as médias corrigidas de cada replicata
erro_quadrado_ms1655 <- cbind((1/3*erro_ms1655_1)^2, (1/3*erro_ms1655_2)^2, (1/3*erro_ms1655_3)^2)    # Cria tabela com os erros padrões de cada replicata ao quadrado
erro_total_ms1655 <- sqrt(rowSums(erro_quadrado_ms1655, na.rm=TRUE))   # Cria vetor com erro associado a tirar a média entre replicatas
tempo_total_ms1655 <- rep(0:24, each=10)   #Cria vetor com tempo em horas para facilitar a gaficação
erro_ln <- erro_total_ms1655/media_total_ms1655

conc_total_ms1655 <- rep(c('0 mol/L', '0.25 mol/L', '0.5 mol/L', '0.75 mol/L', '1 mol/L', '1.25 mol/L', '1.5 mol/L', '2 mol/L', '3 mol/L', '4 mol/L'), times=25)   # Cria vetor com concentraçãoes para conseguir colocar cada concentração em uma cor

#ms1655_1$Conc. <- factor(ms1655_1$Conc., levels = c('0 mol/L', '0.5 mol/L', '1 mol/L', '1.25 mol/L', '1.5 mol/L', '2 mol/L', '3 mol/L', '4 mol/L'))
ms1655_graf <- cbind.data.frame(tempo_total_ms1655, conc_total_ms1655, media_total_ms1655, erro_total_ms1655)
ln_ufc <- log(media_total_ms1655)
mg1655_graf_ln <- cbind.data.frame(tempo_total_ms1655, conc_total_ms1655, ln_ufc, erro_ln)

graf_total <- ggplot(ms1655_graf, aes(x=tempo_total_ms1655, y=media_total_ms1655, color=conc_total_ms1655)) +
    geom_point(size=1.5, aes(shape=conc_total_ms1655))+
    geom_line(linewidth=.7)+
    geom_errorbar(aes(ymin=media_total_ms1655-erro_total_ms1655, ymax=media_total_ms1655+erro_total_ms1655), width=0.15, size=.5)+
    theme_bw()+
    labs(x='Time (h)', y=expression(CFU.mL^{-1}), color='', shape='')+
    theme(
      panel.grid.major = element_line(size=.5),
      panel.grid.minor = element_line(size=.2))+
    #  legend.title = element_text(size=12, family='mono'),
    #  legend.text = element_text(size = 10, family='sans'),
    #  axis.text = element_text(colour = 'black', size=10, family='sans'),
    #  axis.title = element_text(size=12, family='sans'))+
    scale_y_log10(#limits = c(10, NA),
                   labels = trans_format("log10", math_format(10^.x)),
                   breaks=trans_breaks("log10", function(x) 10^x, n=8),
                   minor_breaks=log10_minor_break()
                 )+
    scale_shape_manual(values=c(15,4,17,13,19,10,11,12,18,16))+
    scale_color_viridis(discrete=TRUE, option='turbo')
graf_total
ggsave('mg1655_continuo.png')


graf_cinza <- ggplot(ms1655_graf, aes(x=tempo_total_ms1655, y=media_total_ms1655, color=conc_total_ms1655)) +
    geom_point(size=2, aes(shape=conc_total_ms1655))+
    geom_line(linewidth=.4)+
    geom_errorbar(aes(ymin=media_total_ms1655-erro_total_ms1655, ymax=media_total_ms1655+erro_total_ms1655), width=0.15, size=.5)+
    theme_bw()+
    labs(x='Time (h)', y=expression(CFU.mL^{-1}), color='', shape='')+
    theme(
      panel.grid.major = element_line(size=.3),
      panel.grid.minor = element_line(size=.1))+
    #  legend.title = element_text(size=12, family='mono'),
    #  legend.text = element_text(size = 10, family='sans'),
    #  axis.text = element_text(colour = 'black', size=10, family='sans'),
    #  axis.title = element_text(size=12, family='sans'))+
    scale_y_log10(#limits = c(10, NA),
                   labels = trans_format("log10", math_format(10^.x)),
                   breaks=trans_breaks("log10", function(x) 10^x, n=8),
                   minor_breaks=log10_minor_break()
                 )+
    scale_shape_manual(values=c(15,19,17,23,4,8,14,3,19,15))+
    scale_color_manual(values=c('gray0', 'gray50', 'gray10', 'gray55', 'gray20', 'gray60', 'gray30', 'gray65', 'gray40', 'gray70'))
    #scale_color_grey(start=0, end=.8)
graf_cinza
ggsave('mg1655_junto_cinza.png')


facet_name <- as_labeller(c('0 mol/L'='0.00~mol.L^-1','0.25 mol/L'='0.25~mol.L^-1', '0.5 mol/L'='0.50~mol.L^-1', '0.75 mol/L'='0.75~mol.L^-1', '1 mol/L'='1.00~mol.L^-1', '1.25 mol/L'='1.25~mol.L^-1', '1.5 mol/L'='1.50~mol.L^-1', '2 mol/L'='2.00~mol.L^-1', '3 mol/L'='3.00~mol.L^-1', '4 mol/L'='4.00~mol.L^-1'), default=label_parsed)

graf_facet <- ggplot(ms1655_graf, aes(x=tempo_total_ms1655, y=media_total_ms1655, color=conc_total_ms1655)) +
    geom_point(size=.8)+
    geom_line(linewidth=.4)+
    geom_errorbar(aes(ymin=media_total_ms1655-erro_total_ms1655, ymax=media_total_ms1655+erro_total_ms1655), width=0.15, size=.5)+
    facet_wrap(~conc_total_ms1655, labeller=facet_name)+
    theme_bw()+
    labs(x='Tempo (h)', y=expression(UFC.mL^{-1}))+
    theme(
      panel.grid.major = element_line(size=.2),
      panel.grid.minor = element_line(size=.1),
      legend.position = "none")+
    scale_y_log10(#limits = c(10, NA),
                   labels = trans_format("log10", math_format(10^.x)),
                   breaks=trans_breaks("log10", function(x) 10^x, n=8),
                   minor_breaks=log10_minor_break())+
    scale_color_manual(values=c("#70846B", "#70846B", "#70846B", "#EC9C61", "#EC9C61", "#EC9C61", "#973E3A", "#973E3A", "#973E3A", "#973E3A"))
graf_facet
ggsave('mg1655_facet_tese.png', graf_facet, device='png', unit='cm', width=20, height=16)

graf_facet_cinza <- ggplot(ms1655_graf, aes(x=tempo_total_ms1655, y=media_total_ms1655)) +
    geom_point(size=.5, shape=1)+
    geom_line(linewidth=.2)+
    geom_errorbar(aes(ymin=media_total_ms1655-erro_total_ms1655, ymax=media_total_ms1655+erro_total_ms1655), width=0.15, size=.5)+
    facet_wrap(~conc_total_ms1655, labeller=facet_name, ncol=2)+
    theme_bw()+
    labs(x='Time (h)', y=expression(CFU.mL^{-1}))+
    theme(
      panel.grid.major = element_line(size=.2),
      panel.grid.minor = element_line(size=.1),
      legend.position = "none")+
    scale_y_log10(limits = c(100, 10000000000),
                   labels = trans_format("log10", math_format(10^.x)),
                   breaks=trans_breaks("log10", function(x) 10^x, n=5),
                   minor_breaks=log10_minor_break())
graf_facet_cinza
ggsave('mg1655_facet_cinza.png', graf_facet_cinza, device='png', unit='cm', width=10, height=20)

graf_facet_cinza_ln <- ggplot(mg1655_graf_ln, aes(x=tempo_total_ms1655, y=ln_ufc)) +
    geom_point(size=.5, shape=1)+
    geom_line(linewidth=.2)+
    geom_errorbar(aes(ymin=ln_ufc-erro_ln, ymax=ln_ufc+erro_ln), width=0.15, size=.5)+
    facet_wrap(~conc_total_ms1655, labeller=facet_name, ncol=2)+
    theme_bw()+
    labs(x='Time (h)', y=expression(ln~(CFU.mL^{-1})))+
    theme(
      panel.grid.major = element_line(linewidth=.2),
      panel.grid.minor = element_line(linewidth=.1),
      legend.position = "none")
graf_facet_cinza_ln
ggsave('mg1655_facet_cinza_artigo.png', graf_facet_cinza_ln, device='png', unit='cm', width=10, height=20, dpi=600)



# ********** MODELAGEM **********
colnames(ms1655_graf) <- c('tempo', 'conc', 'media', 'erro')    # Renomeia as colunas para ficar mais fácil a manipulação

# 0M fitado manualmente (Baranyi)
tab_0 <- filter(ms1655_graf, conc=='0 mol/L')   # Pega informações apenas da concentração desejada
d <- tab_0$tempo
erro_0 <- (tab_0$erro)/(tab_0$media)
y <- log(tab_0$media)   # Transforma média em log10(média). Aparentemente os modelos fitam melhor com a transformação
lambda <- 1.5
y0 <- min(y)
fit_0M <- nls(y ~ y0 + mmax * (d + (1/mmax) * log(exp(-mmax*d) +
                exp(-mmax * lambda) - exp(-mmax * (d + lambda)))) -
                log(1 + ((exp(mmax * (d + (1/mmax) * log(exp(-mmax*d) +
                exp(-mmax * lambda) - exp(-mmax *
                (d + lambda)))))-1)/(exp(ymax-y0)))),
                start=list(mmax=0.5, ymax=max(y)))
baranyi_param_0M <- coef(fit_0M)
print(baranyi_param_0M)
coef_0 <- c(min(y), baranyi_param_0M[1], 1.5, baranyi_param_0M[2])
y_0M <- baranyi(coef_0, d)
tab0_graf <- cbind.data.frame(d, y, y_0M, erro_0)
graf_0 <- ggplot(tab0_graf, aes(x=d))+
  geom_point(aes(y=y))+
  geom_line(aes(y=y_0M))+
  geom_errorbar(aes(ymin=y-erro_0, ymax=y+erro_0), width=.3)+
  theme_bw()+
  ggtitle('Baranyi - 0M (lag fixada 2)')+
  labs(
    x='Time (hours)',
    y=expression(paste('ln'*'(UFC.mL'^'-1'*')')))
graf_0
summary(fit_0M)
ggsave('graf_0.png', graf_0)


# 0,25M fitado manualmente (Baranyi)
tab_025 <- filter(ms1655_graf, conc=='0.25 mol/L')   # Pega informações apenas da concentração desejada
d <- tab_025$tempo
erro_025 <- (tab_025$erro)/(tab_025$media)
y <- log(tab_025$media)   # Transforma média em log10(média). Aparentemente os modelos fitam melhor com a transformação
lambda <- 1.5
y0 <- min(y)
fit_025M <- nls(y ~ y0 + mmax * (d + (1/mmax) * log(exp(-mmax*d) +
                exp(-mmax * lambda) - exp(-mmax * (d + lambda)))) -
                log(1 + ((exp(mmax * (d + (1/mmax) * log(exp(-mmax*d) +
                exp(-mmax * lambda) - exp(-mmax *
                (d + lambda)))))-1)/(exp(ymax-y0)))),
                start=list(mmax=0.6, ymax=max(y)))
baranyi_param_025M <- coef(fit_025M)
print(baranyi_param_025M)
coef_025 <- c(min(y), baranyi_param_025M[1], 1.5, baranyi_param_025M[2])
y_025M <- baranyi(coef_025, d)
tab025_graf <- cbind.data.frame(d, y, y_025M, erro_025)
graf_025 <- ggplot(tab025_graf, aes(x=d))+
  geom_point(aes(y=y))+
  geom_line(aes(y=y_025M))+
  geom_errorbar(aes(ymin=y-erro_025, ymax=y+erro_025), width=.3)+
  theme_bw()+
  ggtitle('Baranyi - 0,25M (lag fixada)')+
  labs(
    x='Time (hours)',
    y=expression(paste('ln'*'(UFC.mL'^'-1'*')')))
graf_025
summary(fit_025M)
ggsave('graf_025.png', graf_025)

# 0,5M fitado manualmente (Baranyi)
tab_05 <- filter(ms1655_graf, conc=='0.5 mol/L')   # Pega informações apenas da concentração desejada
d <- tab_05$tempo
erro_05 <- (tab_05$erro)/(tab_05$media)
y <- log(tab_05$media)   # Transforma média em log10(média). Aparentemente os modelos fitam melhor com a transformação
lambda <- 2
y0 <- min(y)
fit_05M <- nls(y ~ y0 + mmax * (d + (1/mmax) * log(exp(-mmax*d) +
                exp(-mmax * lambda) - exp(-mmax * (d + lambda)))) -
                log(1 + ((exp(mmax * (d + (1/mmax) * log(exp(-mmax*d) +
                exp(-mmax * lambda) - exp(-mmax *
                (d + lambda)))))-1)/(exp(ymax-y0)))),
                start=list( mmax=1.5, ymax=max(y)))
baranyi_param_05M <- coef(fit_05M)
print(baranyi_param_05M)
coef_05 <- c(min(y), baranyi_param_05M[1], 2 , baranyi_param_05M[2])
y_05M <- baranyi(coef_05, d)
tab05_graf <- cbind.data.frame(d, y, y_05M, erro_05)
graf_05 <- ggplot(tab05_graf, aes(x=d))+
  geom_point(aes(y=y))+
  geom_line(aes(y=y_05M))+
  geom_errorbar(aes(ymin=y-erro_05, ymax=y+erro_05), width=.3)+
  theme_bw()+
  ggtitle('Baranyi - 0,5M (lag fixada)')+
  labs(
    x='Time (hours)',
    y=expression(paste('ln'*'(UFC.mL'^'-1'*')')))
graf_05
summary(fit_05M)
ggsave('graf_05.png', graf_05)

# 0,75M fitado manualmente (Baranyi)
tab_075 <- filter(ms1655_graf, conc=='0.75 mol/L')   # Pega informações apenas da concentração desejada
d <- tab_075$tempo
erro_075 <- (tab_075$erro)/(tab_075$media)
y <- log(tab_075$media)   # Transforma média em log10(média). Aparentemente os modelos fitam melhor com a transformação
y0 <- min(y)
fit_075M <- nls(y ~ y0 + mmax * (d + (1/mmax) * log(exp(-mmax*d) +
                exp(-mmax * lambda) - exp(-mmax * (d + lambda)))) -
                log(1 + ((exp(mmax * (d + (1/mmax) * log(exp(-mmax*d) +
                exp(-mmax * lambda) - exp(-mmax *
                (d + lambda)))))-1)/(exp(ymax-y0)))),
                start=list(mmax=0.7, lambda=4, ymax=max(y)))
baranyi_param_075M <- coef(fit_075M)
print(baranyi_param_075M)
coef_075 <- c(min(y), baranyi_param_075M[1], baranyi_param_075M[2], baranyi_param_075M[3])
y_075M <- baranyi(coef_075, d)
tab075_graf <- cbind.data.frame(d, y, y_075M, erro_075)
graf_075 <- ggplot(tab075_graf, aes(x=d))+
  geom_point(aes(y=y))+
  geom_line(aes(y=y_075M))+
  geom_errorbar(aes(ymin=y-erro_075, ymax=y+erro_075), width=.3)+
  theme_bw()+
  ggtitle('Baranyi - 0,75M')+
  labs(
    x='Time (hours)',
    y=expression(paste('ln'*'(UFC.mL'^'-1'*')')))
graf_075
summary(fit_075M)
ggsave('graf_075.png', graf_075)

# 1M fitado manualmente (Baranyi)
tab_1 <- filter(ms1655_graf, conc=='1 mol/L')   # Pega informações apenas da concentração desejada
d <- tab_1$tempo
erro_1 <- (tab_1$erro)/(tab_1$media)
y <- log(tab_1$media)   # Transforma média em log10(média). Aparentemente os modelos fitam melhor com a transformação
y0 <- min(y)
fit_1M <- nls(y ~ y0 + mmax * (d + (1/mmax) * log(exp(-mmax*d) +
                exp(-mmax * lambda) - exp(-mmax * (d + lambda)))) -
                log(1 + ((exp(mmax * (d + (1/mmax) * log(exp(-mmax*d) +
                exp(-mmax * lambda) - exp(-mmax *
                (d + lambda)))))-1)/(exp(ymax-y0)))),
                start=list(mmax=0.7, lambda=8, ymax=max(y)))
baranyi_param_1M <- coef(fit_1M)
print(baranyi_param_1M)
coef_1 <- c(min(y), baranyi_param_1M[1], baranyi_param_1M[2], baranyi_param_1M[3])
y_1M <- baranyi(coef_1, d)
tab1_graf <- cbind.data.frame(d, y, y_1M, erro_1)
graf_1 <- ggplot(tab1_graf, aes(x=d))+
  geom_point(aes(y=y))+
  geom_line(aes(y=y_1M))+
  geom_errorbar(aes(ymin=y-erro_1, ymax=y+erro_1), width=.3)+
  theme_bw()+
  ggtitle('Baranyi - 1M')+
  labs(
    x='Time (hours)',
    y=expression(paste('ln'*'(UFC.mL'^'-1'*')')))
graf_1
summary(fit_1M)
ggsave('graf_1.png', graf_1)

# 1,25M fitado manualmente (Baranyi)
tab_125 <- filter(ms1655_graf, conc=='1.25 mol/L')   # Pega informações apenas da concentração desejada
d <- tab_125$tempo
erro_125 <- (tab_125$erro)/(tab_125$media)
y <- log(tab_125$media)   # Transforma média em log10(média). Aparentemente os modelos fitam melhor com a transformação
y0 <- min(y)
fit_125M <- nls(y ~ y0 + mmax * (d + (1/mmax) * log(exp(-mmax*d) +
                exp(-mmax * lambda) - exp(-mmax * (d + lambda)))) -
                log(1 + ((exp(mmax * (d + (1/mmax) * log(exp(-mmax*d) +
                exp(-mmax * lambda) - exp(-mmax *
                (d + lambda)))))-1)/(exp(ymax-y0)))),
                start=list(mmax=0.4, lambda=10, ymax=max(y)))
baranyi_param_125M <- coef(fit_125M)
print(baranyi_param_125M)
coef_125 <- c(y0, baranyi_param_125M[1], baranyi_param_125M[2], baranyi_param_125M[3])
y_125M <- baranyi(coef_125, d)
tab125_graf <- cbind.data.frame(d, y, y_125M, erro_125)
graf_125 <- ggplot(tab125_graf, aes(x=d))+
  geom_point(aes(y=y))+
  geom_line(aes(y=y_125M))+
  geom_errorbar(aes(ymin=y-erro_125, ymax=y+erro_125), width=.3)+
  theme_bw()+
  ggtitle('Baranyi - 1,25M')+
  labs(
    x='Time (hours)',
    y=expression(paste('ln'*'(UFC.mL'^'-1'*')')))
graf_125
summary(fit_125M)
ggsave('graf_125.png', graf_125)




# Gráfico de velocidades de Crescimento
mumax <- c(summary(fit_0M)$parameters[1,1], summary(fit_025M)$parameters[1,1], summary(fit_05M)$parameters[1,1], summary(fit_075M)$parameters[1,1], summary(fit_1M)$parameters[1,1], summary(fit_125M)$parameters[1,1])
conc <- c(0, 0.25, 0.5, 0.75, 1, 1.25)
erro_mu <- c(summary(fit_0M)$parameters[1,2], summary(fit_025M)$parameters[1,2], summary(fit_05M)$parameters[1,2], summary(fit_075M)$parameters[1,2], summary(fit_1M)$parameters[1,2], summary(fit_125M)$parameters[1,2])
tab_mu <- cbind.data.frame(mumax, erro_mu, conc)
graf_mu <- ggplot(tab_mu, aes(x=conc, y=mumax))+
    geom_point()+
    geom_line()+
    geom_errorbar(aes(ymin=mumax-erro_mu, ymax=mumax+erro_mu), width=.05)+
    theme_bw()+
    labs(
      x=expression(paste('[NaCl] (mol.L)'^'-1')),
      y=expression(paste(mu~(h^-1))))+
    scale_x_continuous(breaks=c(0, 0.25, 0.5, 0.75, 1, 1.25))+
    theme(
      panel.grid.major = element_line(size=.25),   # Define largura da linha maior do grid interno no gráfico
      panel.grid.minor = element_line(size=.1),   # Define largura da linha menor do grid interno do gráfico
      axis.title.y = element_text(size=12))
graf_mu
ggsave('mu.png', graf_mu, device='png', unit='cm', width=7, height=7, dpi=600)

write.table(tab_mu, file='mumax_tab.tsv', sep='\t', row.names=FALSE, dec=',')


#********** PSEUDO-R2 **********
y_baranyi_0 <- yi_fitado_baranyi(d, coef_0)
y_baranyi_025 <- yi_fitado_baranyi(d, coef_025)
y_baranyi_05 <- yi_fitado_baranyi(d, coef_05)
y_baranyi_075 <- yi_fitado_baranyi(d, coef_075)
y_baranyi_1 <- yi_fitado_baranyi(d, coef_1)
y_baranyi_125 <- yi_fitado_baranyi(d, coef_125)

y_0 <- log(tab_0$media)
y_025 <- log(tab_025$media)
y_05 <- log(tab_05$media)
y_075 <- log(tab_075$media)
y_1 <- log(tab_1$media)
y_125 <- log(tab_125$media)

pseudoR2_0M <- calc_pseudo_R2(y_0, y_baranyi_0)
pseudoR2_025M <- calc_pseudo_R2(y_025, y_baranyi_025)
pseudoR2_05M <- calc_pseudo_R2(y_05, y_baranyi_05)
pseudoR2_075M <- calc_pseudo_R2(y_075, y_baranyi_075)
pseudoR2_1M <- calc_pseudo_R2(y_1, y_baranyi_1)
pseudoR2_125M <- calc_pseudo_R2(y_125, y_baranyi_125)

mse_0M <- calc_mse(y_0, y_baranyi_0)
mse_025M <- calc_mse(y_025, y_baranyi_025)
mse_05M <- calc_mse(y_05, y_baranyi_05)
mse_075M <- calc_mse(y_075, y_baranyi_075)
mse_1M <- calc_mse(y_1, y_baranyi_1)
mse_125M <- calc_mse(y_125, y_baranyi_125)
