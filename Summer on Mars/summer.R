# R 4.2.0
# PARA EXPERIMENTO SUMMER ON MARS - Grupo Quimiosfera
# Scrip para gráficos e análises de dados do experimento Summer on Mars, Ciclo de congelamento e descongelamento

# REGRAS DE ENTRADA: vírgula para decimal, tab para separação de células
# ORDEM DAS COLUNAS: 1. Ponto | 2. Concentração de NaCl | 3. Horas de experimento | 4-7. Contagem de UFC | 8. Fator de correção da diluição | 9. Média corrigida

# ********** PACOTES E BIBLIOTECAS **********
#install.packages('dplyr')
library(dplyr)
#install.packages('janitor')
library(janitor)
#install.packages('ggplot2')
library(ggplot2)
#install.packages('viridis')
library(viridis)
#install.packages('matrixStats')
library(matrixStats)
#install.packages('MASS')
library(MASS)
#install.packages('scales')
library(scales)
#install.packages('babar')
library(babar)


setwd("Documents/USP-Pós/Projeto Summer on Mars/redacao cientifica")

source("/home/rvincenzi/Documents/USP-Pós/Projeto Summer on Mars/redacao cientifica/summer.R")

# ********** LEITURA E PRÉ-PROCESSAMENTO **********

# Lê tabela em TSV, vírgula como decimal e primeira linha como cabeçalho
rep1 <- read.table('rep1.tsv', sep='\t', dec=',', header=TRUE, fill=TRUE)
rep2 <- read.table('rep2.tsv', sep='\t', dec=',', header=TRUE, fill=TRUE)
rep3 <- read.table('rep3.tsv', sep='\t', dec=',', header=TRUE, fill=TRUE)


ufc_1 <- as.matrix(rep1[,c(3,4,5,6)])    # Cria matriz só com com as 3 colunas de contagem das gotas
ufc_2 <- as.matrix(rep2[,c(3,4,5,6)])
ufc_3 <- as.matrix(rep3[,c(3,4,5,6)])

media_1 <- rowMeans(ufc_1, na.rm=TRUE)    # Faz média das 4 colunas com valores contados de UFC
media_2 <- rowMeans(ufc_2, na.rm=TRUE)
media_3 <- rowMeans(ufc_3, na.rm=TRUE)

media_cor_1 <- ((rep1$Dil.*100))*(media_1)   # Faz média corrigida já para densidade ser por mL
media_cor_2 <- ((rep2$Dil.*100))*(media_2)
media_cor_3 <- ((rep3$Dil.*100))*(media_3)

c_1 <- rep(c(media_cor_1[1], media_cor_1[2], media_cor_1[3]), times=9)    # Faz vetor com repetição com a quantidade de pontos (ciclos) dos 3 valores do controle
c_2 <- rep(c(media_cor_2[1], media_cor_2[2], media_cor_2[3]), times=9)
c_3 <- rep(c(media_cor_3[1], media_cor_3[2], media_cor_3[3]), times=9)

sob_1 <- media_cor_1/c_1    # Calcula a sobrevivência (experimento/controle)
sob_2 <- media_cor_2/c_2
sob_3 <- media_cor_3/c_3

sd_1 <- rowSds(ufc_1, na.rm=TRUE)   # Calcula desvio padrão entre as gostas de uma replicata
sd_2 <- rowSds(ufc_2, na.rm=TRUE)
sd_3 <- rowSds(ufc_3, na.rm=TRUE)

n_1 <- rowSums(!is.na(ufc_1))   # Calcula número de gostas contadas. Para isso ele conta quantos valores não NA existem naquelas 3 colunas
n_2 <- rowSums(!is.na(ufc_2))
n_3 <- rowSums(!is.na(ufc_3))

se_1 <- sd_1/(sqrt(n_1))   # Calcula erro padrão (por definição: desvio padrão dividido por raiz quadrada do n)
se_2 <- sd_2/(sqrt(n_2))
se_3 <- sd_3/(sqrt(n_3))

erro_1 <- ((rep1$Dil*100))*(se_1)   # Calcula erro associado a média corrigida das três gotinhas. Propagação do erro padrão considerando Normal.
erro_2 <- ((rep2$Dil*100))*(se_2)
erro_3 <- ((rep3$Dil*100))*(se_3)

erro_c_1 <- rep(c(erro_1[1], erro_1[2], erro_1[3]), times=9)    # Cria vetor apenas com o erro dos controles
erro_c_2 <- rep(c(erro_2[1], erro_2[2], erro_2[3]), times=9)
erro_c_3 <- rep(c(erro_3[1], erro_3[2], erro_3[3]), times=9)

erro_p_1 <- sqrt((sob_1^2)*(((erro_1/media_cor_1)^2)+((erro_c_1/c_1)^2)))   # Faz vetor com erro das frações de sobrevivência
erro_p_2 <- sqrt((sob_2^2)*(((erro_2/media_cor_2)^2)+((erro_c_2/c_2)^2)))
erro_p_3 <- sqrt((sob_3^2)*(((erro_3/media_cor_3)^2)+((erro_c_3/c_3)^2)))




# Para gráfico total (dados reais, média entre replicatas)
sobrevivencias <- cbind(sob_1, sob_2, sob_3)   # Junta as sobrevivências de cada replicata numa tabela, cada replicata como uma coluna
sob_total <- rowMeans(sobrevivencias, na.rm=TRUE)    # Cria vetor da média entre as médias corrigidas de cada replicata
erro_quadrado <- cbind((1/3*erro_p_1)^2, (1/3*erro_p_2)^2, (1/3*erro_p_3)^2)    # Cria tabela com os erros padrões de cada replicata ao quadrado
erro_total <- sqrt(rowSums(erro_quadrado, na.rm=TRUE))   # Cria vetor com erro associado a tirar a média entre replicatas
ciclos <- rep(0:8, each=3)   #Cria vetor com tempo em horas para facilitar a gaficação
bichos <- rep(c('ColB', 'P4D', 'E. coli'), times=9)   # Cria vetor com repetição dos nomes dos bichos pra conseguir separar por bicho no gráfico

sob_graf <- cbind.data.frame(bichos, ciclos, sob_total, erro_total)   #Junta todos os dados necessários pro gráfico num único data frame

write.table(sob_graf, file='sob.tsv', quote=FALSE, sep='\t', col.names = NA)


graf_sob<- ggplot(sob_graf, aes(x=ciclos, y=sob_total, color=bichos)) +
    geom_point(size=1)+
    geom_line(size=.7)+
    geom_errorbar(aes(ymin=sob_total-erro_total, ymax=sob_total+erro_total), width=0.15, size=.5)+
    theme_bw()+
    labs(x='Cycles', y='Survival', color='')+
    theme(
      panel.grid.major = element_line(size=.5),
      panel.grid.minor = element_line(size=.2)) +
    scale_color_manual(values = c("gray51", "firebrick1", "black"))
    #  legend.title = element_text(size=12, family='mono'),
    #  legend.text = element_text(size = 10, family='sans'),
    #  axis.text = element_text(colour = 'black', size=10, family='sans'),
    #  axis.title = element_text(size=12, family='sans'))+

ggsave('sobrevivencia.png', width = 5, height = 4)



setwd("/home/rvincenzi/Documents/USP-Pós/Projeto Summer on Mars/Summer on Mars/Curva de crescimento 1")


# Para pacote babar - MODELAGEM
halo1_c <- read.table('halo1_c.tsv', sep='\t', dec=',', header=TRUE, fill=TRUE)
halo2_c <- read.table('halo2_c.tsv', sep='\t', dec=',', header=TRUE, fill=TRUE)
halo1 <- read.table('halo1.tsv', sep='\t', dec=',', header=TRUE, fill=TRUE)
halo2 <- read.table('halo2.tsv', sep='\t', dec=',', header=TRUE, fill=TRUE)
colb1_c <- read.table('colb1_c.tsv', sep='\t', dec=',', header=TRUE, fill=TRUE)
colb2_c <- read.table('colb2_c.tsv', sep='\t', dec=',', header=TRUE, fill=TRUE)
colb1 <- read.table('colb1.tsv', sep='\t', dec=',', header=TRUE, fill=TRUE)
colb2 <- read.table('colb2.tsv', sep='\t', dec=',', header=TRUE, fill=TRUE)

#Criação de tabelas para ficar no formato certo de input - Tempo (nome da coluna deve ser 'time') | log10 da contagem de CFU (nome da coluna deve ser 'logc')


#h1_c <- data.frame(time=halo1_c[,4], logc=log10(halo1_c[,11]))
#h2_c <- data.frame(time=halo2_c[,4], logc=log10(halo2_c[,11]))
#h1 <- data.frame(time=halo1[,4], logc=log10(halo1[,10]))
#h2 <- data.frame(time=halo2[,4], logc=log10(halo2[,10]))
#c1_c <- data.frame(time=colb1_c[,4], logc=log10(colb1_c[,11]))
#c2_c <- data.frame(time=colb2_c[,4], logc=log10(colb2_c[,11]))
#c1 <- data.frame(time=colb1[,4], logc=log10(colb1[,10]))
#c2 <- data.frame(time=colb2[,4], logc=log10(colb2[,10]))


h_c <- cbind(halo1_c[,11], halo2_c[,11])
halo_media_c <- rowMeans(h_c)

h <- cbind(halo1[,10], halo2[,10])
halo_media <- rowMeans(h)

c_c <- cbind(colb1_c[,11], colb2_c[,11])
colb_media_c <- rowMeans(c_c)

c <- cbind(halo1[,10], halo2[,10])
colb_media <- rowMeans(h)

halo_c <- data.frame(time=halo1_c[,4], logc=log10(halo_media_c))
halo <- data.frame(time=halo1[,4], logc=log10(halo_media))
colb_c <- data.frame(time=colb2_c[,4], logc=log10(colb_media_c))
colb <- data.frame(time=colb1[,4], logc=log10(colb_media))




#Fit com bayesfit

#fit_h1_c <- Bayesfit(h1_c, model='logistic')
#fit_h2_c <- Bayesfit(h2_c, model='logistic')
#fit_h1 <- Bayesfit(h1, model='logistic')
#fit_h2 <- Bayesfit(h2, model='logistic')
#fit_c1_c <- Bayesfit(c1_c, model='logistic')
#fit_c2_c <- Bayesfit(c2_c, model='logistic')
#fit_c2 <- Bayesfit(c2, model='logistic')
#fit_c1 <- Bayesfit(c1, model='logistic')

fit_h_c <- Bayesfit(halo_c, model='logistic')
fit_h <- Bayesfit(halo, model='logistic')
fit_c_c <- Bayesfit(colb_c, model='logistic')
fit_c <- Bayesfit(colb, model='logistic')



# Comparação halomonas controle e experimento
h1_halo <- Bayescompare(halo_c, halo, hyp='H1', model='logistic')		# Hipótese 1: curvas são replicatas
h2_halo <- Bayescompare(halo_c, halo, hyp='H2', model='logistic')		# Hipótese 2: cuvas tem mesma taxa de crescimento
h3_halo <- Bayescompare(halo_c, halo, hyp='H3', model='logistic')		# Hipótese 3: todos os parâmetros das curvas são diferentes

le_h1_halo <- h1_halo$logevidence
le_h2_halo <- h2_halo$logevidence
le_h3_halo <- h3_halo$logevidence

# Calcula 2lnB, sendo B op fator de Bayes, valores podem ser interpretados de acordo com a Escala de Jeffrey's
lbf_12_halo <-2*log(exp(le_h1_halo)/exp(le_h2_halo))	 # Hipótese 1 contra 2, Halomonas controle versus experimento
lbf_13_halo <-2*log(exp(le_h1_halo)/exp(le_h3_halo))    # Hipótese 1 contra  3, Halomonas controle versus experimento

#[1] -4.362194

# Comparação ColB controle e experimento
h1_colb <- Bayescompare(colb_c, colb, hyp='H1', model='logistic')		# Hipótese 1: curvas são replicatas
h2_colb <- Bayescompare(colb_c, colb, hyp='H2', model='logistic')		# Hipótese 2: cuvas tem mesma taxa de crescimento
h3_colb <- Bayescompare(colb_c, colb, hyp='H3', model='logistic')		# Hipótese 3: todos os parâmetros das curvas são diferentes

le_h1_colb <- h1_colb$logevidence
le_h2_colb <- h2_colb$logevidence
le_h3_colb <- h3_colb$logevidence

# Calcula 2lnB, sendo B op fator de Bayes, valores podem ser interpretados de acordo com a Escala de Jeffrey's
lbf_12_colb <-2*log(exp(le_h1_colb)/exp(le_h2_colb))	 # Hipótese 1 contra 2, ColB controle versus experimento


# Fit das curvas individuais com modelo logístico
fit_h1_c <- Bayesfit(h1_c, model='logistic')
fit_h2_c <- Bayesfit(h2_c, model='logistic')
fit_m_1M <- Bayesfit(t_m_1M, model='logistic')
fit_m_125M <- Bayesfit(t_m_125M, model='logistic')

par_0 <- fit_m_0M$means
var_0 <- fit_m_0M$vars
par_05 <- fit_m_05M$means
var_05 <- fit_m_05M$vars
par_1 <- fit_m_1M$means
var_1 <- fit_m_1M$vars
par_125 <- fit_m_125M$means
var_125 <- fit_m_125M$vars

mu_conc <- c(0, 0.5, 1, 1.25)		# Cria vetor com concentrações de NaCl em que há crescimento
mu_max <- c(par_0[3], par_05[3], par_1[3], par_125[3])
mu_sd <- c(sqrt(var_0[3]), sqrt(var_05[3]), sqrt(var_1[3]), sqrt(var_125[3]))
mu <- data.frame(conc=mu_conc, mu=mu_max, sd=mu_sd)

#bayes <- exp(le_h1)/exp(le_h3)
#lbf <- 2*log(bayes)

# Gráficos dados do babar
#Gráfico de mu x concentração de NaCl
graf_mu<- ggplot(mu, aes(x=conc, y=mu)) +
    geom_point(size=2)+
    geom_line(size=.7)+
    geom_errorbar(aes(ymin=mu-sd, ymax=mu+sd), width=0.01, size=.5)+
    theme_bw()+
    labs(x='NaCl concentration (mol/L)', y='Growth Rate')+
    theme(
      panel.grid.major = element_line(size=.5),
      panel.grid.minor = element_line(size=.2))

ggsave('graf_mu.png')

#Gráfico modelo e dados reais
