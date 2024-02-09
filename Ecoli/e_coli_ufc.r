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
install.packages('growthrates')
library(growthrates)


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
#erro_junto_ms1655 <- cbind(erro_ms1655_1, erro_ms1655_2, erro_ms1655_3)
#erro_total_ms1655 <- rowMeans(erro_junto_ms1655, na.rm=TRUE)
erro_quadrado_ms1655 <- cbind((1/3*erro_ms1655_1)^2, (1/3*erro_ms1655_2)^2, (1/3*erro_ms1655_3)^2)    # Cria tabela com os erros padrões de cada replicata ao quadrado
erro_total_ms1655 <- sqrt(rowSums(erro_quadrado_ms1655, na.rm=TRUE))   # Cria vetor com erro associado a tirar a média entre replicatas
tempo_total_ms1655 <- rep(0:24, each=10)   #Cria vetor com tempo em horas para facilitar a gaficação
conc_total_ms1655 <- rep(c('0 mol/L', '0.25 mol/L', '0.5 mol/L', '0.75 mol/L', '1 mol/L', '1.25 mol/L', '1.5 mol/L', '2 mol/L', '3 mol/L', '4 mol/L'), times=25)   # Cria vetor com concentraçãoes para conseguir colocar cada concentração em uma cor

#ms1655_1$Conc. <- factor(ms1655_1$Conc., levels = c('0 mol/L', '0.5 mol/L', '1 mol/L', '1.25 mol/L', '1.5 mol/L', '2 mol/L', '3 mol/L', '4 mol/L'))
ms1655_graf <- cbind.data.frame(tempo_total_ms1655, conc_total_ms1655, media_total_ms1655, erro_total_ms1655)

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

graf_total<- ggplot(ms1655_graf, aes(x=tempo_total_ms1655, y=media_total_ms1655, color=conc_total_ms1655)) +
    geom_point(size=1.5, aes(shape=conc_total_ms1655))+
    geom_line(size=.7)+
    geom_errorbar(aes(ymin=media_total_ms1655-erro_total_ms1655, ymax=media_total_ms1655+erro_total_ms1655), width=0.15, size=.5)+
    theme_bw()+
    labs(x='Time (h)', y='CFU/mL', color='', shape='')+
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

ggsave('mg1655_continuo.png')



# ********** MODELAGEM **********
# Criando tabelas para modelagem usando pacote growthrates
colnames(ms1655_graf) <- c('tempo', 'conc', 'media', 'erro')    # Renomeia as colunas para ficar mais fácil a manipulação

# 0M fitado com Gompertz de 3 parâmetros (com lag)
tab_0 <- filter(ms1655_graf, conc=='0 mol/L')   # Pega informações apenas da concentração desejada
tab0 <- tab_0[,-2]   # Retira coluna de concentração e erro
tab0$erro <- (tab0$erro)/(tab0$media*log(10))   # Faz propagação de erro associado ao log10(media)
tab0$media <- log10(tab0$media)   # Transforma média em log10(média). Aparentemente os modelos fitam melhor com a transformação


p0 <- c(y0=tab0[1,2], mumax=0.2, K=max(tab0$media), lambda=0)   # Chutes iniciais dos parâmetros do modelo de Gompertz (y0=valor de log10(UFC) do tempo 0, mumax=velocidade máxima de crescimento, K=carrying capacity, lambda=tempo de lag)
fit0 <- fit_growthmodel(FUN=grow_gompertz3, p=p0, tab0$tempo, tab0$media, transform="log", which=c('y0','mumax', 'K'))    # Faz o fit do modelo usando os dados experimentais e chutes iniciais de parâmetros. Nesse caso fixa lag=0
param_0 <- coef(fit0)   # Cria vetor com parâmetros fitados
summary(fit0)   # Mostra todos os valores de parâmetros, erro padrão associado, t value, significância, graus de liberdade e correlação entre os parâmetros

tab0_graf <- grow_gompertz3(tab0$tempo, param_0)    # Cria curva contínua do Gompertz de 3 parâmetros com os parâmetros fitados
tab0_graf <- cbind.data.frame(tab0_graf, tab0$media, tab0$erro)    # Cria tabela com valores fitados e dados experimentais
graf_0 <- ggplot(tab0_graf, aes(x=time))+   # Faz gráfico de linha (fit) e pontos (dados experimentais) sobrepostos
    geom_point(aes(y=tab0$media))+
    geom_line(aes(y=y))+
    geom_errorbar(aes(ymin=tab0$media-tab0$erro, ymax=tab0$media+tab0$erro), width=.3)+
    theme_bw()+
    ggtitle('0M - Gompertz3 (lag fixada)')+
    labs(
      x='Time (hours)',
      y=expression(paste('log'['10']*'(UFC.mL'^'-1'*')')))
graf_0
ggsave('0M-Gompertz3.png')    # Salva gráfico do fit e dados originais


# 0,25M fitado com Gompertz de 3 parâmetros (fixando lag=0)
tab_025 <- filter(ms1655_graf, conc=='0.25 mol/L')
tab025 <- tab_025[,-2]
tab025$erro <- (tab025$erro)/(tab025$media*log(10))
tab025$media <- log10(tab025$media)

p025 <- c(y0=tab025[1,2], mumax=0.2, K=max(tab025$media) , lambda=0)
fit025 <- fit_growthmodel(FUN=grow_gompertz3, p=p025, tab025$tempo, tab025$media, transform="log", which=c('y0', 'mumax', 'K'))
param_025 <- coef(fit025)
summary(fit025)

tab025_graf <- grow_gompertz3(tab025$tempo, param_025)
tab025_graf <- cbind.data.frame(tab025_graf, tab025$media, tab025$erro)
graf_025 <- ggplot(tab025_graf, aes(x=time))+   # Faz gráfico de linha (fit) e pontos (dados experimentais) sobrepostos
    geom_point(aes(y=tab025$media))+
    geom_line(aes(y=y))+
    geom_errorbar(aes(ymin=tab025$media-tab025$erro, ymax=tab025$media+tab025$erro), width=.3)+
    theme_bw()+
    ggtitle('0,25M - Gompertz3 (lag fixada)')+
    labs(
      x='Time (hours)',
      y=expression(paste('log'['10']*'(UFC.mL'^'-1'*')')))
graf_025
ggsave('0,25M-Gompertz3.png')


# 0,5M fitado com Gompertz de 3 parâmetros (fixando lag=0)
tab_05 <- filter(ms1655_graf, conc=='0.5 mol/L')
tab05 <- tab_05[,-2]
tab05$erro <- (tab05$erro)/(tab05$media*log(10))
tab05$media <- log10(tab05$media)

p05 <- c(y0=tab05[1,2], mumax=0.2, K=max(tab05$media), lambda=0)
fit05 <- fit_growthmodel(FUN=grow_gompertz3, p=p05, tab05$tempo, tab05$media, transform="log", which=c('y0', 'mumax', 'K'))
param_05 <- coef(fit05)
summary(fit05)

tab05_graf <- grow_gompertz3(tab05$tempo, param_05)
tab05_graf <- cbind.data.frame(tab05_graf, tab05$media, tab05$erro)
graf_05 <- ggplot(tab05_graf, aes(x=time))+   # Faz gráfico de linha (fit) e pontos (dados experimentais) sobrepostos
    geom_point(aes(y=tab05$media))+
    geom_line(aes(y=y))+
    geom_errorbar(aes(ymin=tab05$media-tab05$erro, ymax=tab05$media+tab05$erro), width=.3)+
    theme_bw()+
    ggtitle('0,5M - Gompertz3 (lag fixada)')+
    labs(
      x='Time (hours)',
      y=expression(paste('log'['10']*'(UFC.mL'^'-1'*')')))
graf_05
ggsave('0,5M-Gompertz3.png')



# 0,75M fitado com Gompertz de 3 parâmetros (sem fixar lag)
tab_075 <- filter(ms1655_graf, conc=='0.75 mol/L')
tab075 <- tab_075[,-2]
tab075$erro <- (tab075$erro)/(tab075$media*log(10))
tab075$media <- log10(tab075$media)

p075 <- c(y0=tab075[1,2], mumax=0.2, K=max(tab075$media), lambda=0)
fit075 <- fit_growthmodel(FUN=grow_gompertz3, p=p075, tab075$tempo, tab075$media, transform="log")
param_075 <- coef(fit075)
summary(fit075)

tab075_graf <- grow_gompertz3(tab075$tempo, param_075)
tab075_graf <- cbind.data.frame(tab075_graf, tab075$media, tab075$erro)
graf_075 <- ggplot(tab075_graf, aes(x=time))+   # Faz gráfico de linha (fit) e pontos (dados experimentais) sobrepostos
    geom_point(aes(y=tab075$media))+
    #geom_line(aes(y=y))+
    geom_errorbar(aes(ymin=tab075$media-tab075$erro, ymax=tab075$media+tab075$erro), width=.3)+
    theme_bw()+
    ggtitle('0,75M - Gompertz3')+
    labs(
      x='Time (hours)',
      y=expression(paste('log'['10']*'(UFC.mL'^'-1'*')')))
graf_075
ggsave('0,75M-Gompertz3_pontos.png')


# 1M fitado com Gompertz de 3 parâmetros (sem fixar lag)
tab_1 <- filter(ms1655_graf, conc=='1 mol/L')
tab1 <- tab_1[,-2]
tab1$erro <- (tab1$erro)/(tab1$media*log(10))
tab1$media <- log10(tab1$media)

p1 <- c(y0=tab1[1,2], mumax=0.2, K=max(tab1$media), lambda=5)
fit1 <- fit_growthmodel(FUN=grow_gompertz3, p=p1, tab1$tempo, tab1$media, transform="log")
param_1 <- coef(fit1)
summary(fit1)

tab1_graf <- grow_gompertz3(tab1$tempo, param_1)
tab1_graf <- cbind.data.frame(tab1_graf, tab1$media, tab1$erro)
graf_1 <- ggplot(tab1_graf, aes(x=time))+   # Faz gráfico de linha (fit) e pontos (dados experimentais) sobrepostos
    geom_point(aes(y=tab1$media))+
    geom_line(aes(y=y))+
    geom_errorbar(aes(ymin=tab1$media-tab1$erro, ymax=tab1$media+tab1$erro), width=.3)+
    theme_bw()+
    ggtitle('1M - Gompertz3')+
    labs(
      x='Time (hours)',
      y=expression(paste('log'['10']*'(UFC.mL'^'-1'*')')))
graf_1
ggsave('1M-Gompertz3.png')


# 1,25M fitado com Gompertz de 3 parâmetros (sem fixar lag)
tab_125 <- filter(ms1655_graf, conc=='1.25 mol/L')
tab125 <- tab_125[,-2]
tab125$erro <- (tab125$erro)/(tab125$media*log(10))
tab125$media <- log10(tab125$media)

p125 <- c(y0=tab125[1,2], mumax=0.2, K=max(tab125$media), lambda=5)
fit125 <- fit_growthmodel(FUN=grow_gompertz3, p=p125, tab125$tempo, tab125$media, transform="log")
param_125 <- coef(fit125)
summary(fit125)

tab125_graf <- grow_gompertz3(tab125$tempo, param_125)
tab125_graf <- cbind.data.frame(tab125_graf, tab125$media, tab125$erro)
graf_125 <- ggplot(tab125_graf, aes(x=time))+   # Faz gráfico de linha (fit) e pontos (dados experimentais) sobrepostos
    geom_point(aes(y=tab125$media))+
    geom_line(aes(y=y))+
    geom_errorbar(aes(ymin=tab125$media-tab125$erro, ymax=tab125$media+tab125$erro), width=.3)+
    theme_bw()+
    ggtitle('1,25M - Gompertz3')+
    labs(
      x='Time (hours)',
      y=expression(paste('log'['10']*'(UFC.mL'^'-1'*')')))
graf_125
ggsave('1,25M-Gompertz3.png')


# Gráfico de velocidades de Crescimento
mumax <- c(param_0[2], param_025[2], param_05[2], param_075[2], param_1[2], param_125[2])
conc <- c(0, 0.25, 0.5, 0.75, 1, 1.25)
erro_mu <- c(0.005265, 0.006495, 0.005802, 0.004528, 0.003579, 0.00111)
tab_mu <- cbind.data.frame(mumax, erro_mu, conc)
graf_mu <- ggplot(tab_mu, aes(x=conc, y=mumax))+
    geom_point()+
    geom_line()+
    geom_errorbar(aes(ymin=mumax-erro_mu, ymax=mumax+erro_mu), width=.01)+
    theme_bw()+
    labs(
      x=expression(paste('[NaCl] (mol.L)'^'-1')),
      y=expression(paste(mu['max']))
    )
graf_mu
ggsave('mu.png')
