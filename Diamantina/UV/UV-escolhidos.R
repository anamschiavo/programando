#Script Doutorado Ana Paula Muche Schiavo

#REGRAS DE ENTRADA:
#*************ORDEM DAS COLUNAS: 1. Espécie | 2. Dose | 3-5. UFC | 6. Fator de Diluição **********************************
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
# Lê tabela em tsv, ou seja, separado por tab, vírgula como decimal, primeira linha é cabeçalho, preenche espáços vazios com NA
uvtab <- read.table("UVC_escolhidos.tsv", sep="\t", dec=",", quote="", header=TRUE, fill=TRUE)

#********** PROCESSAMENTO **********
col_indices <- c(3, 4, 5)   # Define quais são as colunas da tabela que tem valores de UFC
mc <- calc_media_corrigida(uvtab, col_indices)
erro_mc <- calc_erro_medcor(uvtab, col_indices)

tab_media <- cbind.data.frame(uvtab, mc, erro_mc)

rep1 <- filter(tab_media, Replicata=='1')
rep2 <- filter(tab_media, Replicata=='2')
rep3 <- filter(tab_media, Replicata=='3')

media_rep <- rowMeans(cbind(rep1$mc, rep2$mc, rep3$mc), na.rm = TRUE)
erro_rep <- sqrt(rowSums(cbind((1/3*rep1$erro_mc)^2, (1/3*rep2$erro_mc)^2, (1/3*rep3$erro_mc)^2), na.rm = TRUE))

tab <- cbind.data.frame(rep1$Isolado, rep1$Dose, media_rep, erro_rep)
tab0 <- filter(tab, rep1$Dose=='0')
media0 <- rep(tab0$media_rep, each=8)
erro0 <- rep(tab0$erro_rep, each=8)
sob <- media_rep/media0
erro_sob <- sqrt((sob^2) * (((erro_rep / media_rep)^2) + ((erro0 / media0)^2)))

tabgraf <- cbind.data.frame(rep1$Isolado, rep1$Dose, sob, erro_sob)
colnames(tabgraf) <- c('isolado', 'Dose', 'sob', 'erro')

#********** GRÁFICOS **********
graf_sob<- ggplot(tabgraf, aes(x=Dose, y=sob)) +
    geom_point(size=.65)+   # Coloca as bolinhas
    geom_line(linewidth=.35)+   # Coloca as linhas
    geom_errorbar(aes(ymin=sob-erro, ymax=sob+erro), width=0.15, size=.5)+    # Coloca as barras de erro
    facet_wrap(~isolado, ncol=3)+    # Deixa facetado por espécie
    theme_bw()+   # Muda a aparência para o tema certo
    xlab(bquote("Dose"~(J.m^-2)))+    # Coloca o nome do eixo x
    ylab(bquote("Sobrevivência"~(N/N[0])))+    # Coloca o nome do eixo y
    theme(
      panel.grid.major = element_line(size=.25),   # Define largura da linha maior do grid interno no gráfico
      panel.grid.minor = element_line(size=.1),   # Define largura da linha menor do grid interno do gráfico
      legend.position='none')+    # Retira a legenda
    scale_y_log10(limits = c(1e-4,10),   # Transforma o eixo y em log10 com limiter de 1 a 10^-5
                  labels = trans_format("log10", math_format(10^.x)),
                   breaks = trans_breaks("log10", function(x) 10^x, n=4),
                   minor_breaks=log10_minor_break())
graf_sob

ggsave('escolhidos_UV.png', graf_sob, device='png', unit='cm', width = 15, height = 20, dpi=300)   # Salva gráfico que acabamos de criar como um arquivo

#********** FITTING **********
x <- c(0, 50, 100, 200, 300, 500, 700, 1000)
ld10 <- c()
se_ld10 <- c()
isolado <- c()

# Isolado D16a.30.010-2
d16a_30_010_2 <- filter(tabgraf, isolado=='D16a.30.010-2')
y_16a_010_2 <- d16a_30_010_2$sob
f <- 0
fit_16a_010_2 <- nls(y_16a_010_2 ~ (1+f*x)/(1+(9+10*f*LD10)*exp(n*log(x/LD10))),
			        start = list(LD10 = 100, n = 2))
summary(fit_16a_010_2)
ld10 <- append(ld10, summary(fit_16a_010_2)$parameters[1,1] , after=length(ld10))
se_ld10 <- append(se_ld10, summary(fit_16a_010_2)$parameters[1,2] , after=length(se_ld10))
isolado <- append(isolado, 'D16a.30.010-2', after=length(isolado))

# Isolado D16a.30.010-7
d16a_30_010_7 <- filter(tabgraf, isolado=='D16a.30.010-7')
y_16a_010_7 <- d16a_30_010_7$sob
fit_16a_010_7 <- nls(y_16a_010_7 ~ (1+f*x)/(1+(9+10*f*LD10)*exp(n*log(x/LD10))),
			        start = list(LD10 = 120, n = 2))
summary(fit_16a_010_7)
ld10 <- append(ld10, summary(fit_16a_010_7)$parameters[1,1] , after=length(ld10))
se_ld10 <- append(se_ld10, summary(fit_16a_010_7)$parameters[1,2] , after=length(se_ld10))
isolado <- append(isolado, 'D16a.30.010-7', after=length(isolado))

# Isolado D17b.30.010-2
d17b_30_010_2 <- filter(tabgraf, isolado=='D17b.30.010-2')
y_17b_010_2 <- d17b_30_010_2$sob
fit_17b_010_2 <- nls(y_17b_010_2 ~ (1+f*x)/(1+(9+10*f*LD10)*exp(n*log(x/LD10))),
			        start = list(LD10 = 200, n = 2))
summary(fit_17b_010_2)
ld10 <- append(ld10, summary(fit_17b_010_2)$parameters[1,1] , after=length(ld10))
se_ld10 <- append(se_ld10, summary(fit_17b_010_2)$parameters[1,2] , after=length(se_ld10))
isolado <- append(isolado, 'D17b.30.010-2', after=length(isolado))

# Isolado D17b.30.010-4
d17b_30_010_4 <- filter(tabgraf, isolado=='D17b.30.010-4')
y_17b_010_4 <- d17b_30_010_4$sob
fit_17b_010_4 <- nls(y_17b_010_4 ~ (1+f*x)/(1+(9+10*f*LD10)*exp(n*log(x/LD10))),
			        start = list(LD10 = 200, n = 2))
summary(fit_17b_010_4)
ld10 <- append(ld10, summary(fit_17b_010_4)$parameters[1,1] , after=length(ld10))
se_ld10 <- append(se_ld10, summary(fit_17b_010_4)$parameters[1,2] , after=length(se_ld10))
isolado <- append(isolado, 'D17b.30.010-4', after=length(isolado))

# Isolado D17b.30.010-4.2
d17b_30_010_42 <- filter(tabgraf, isolado=='D17b.30.010-4.2')
y_17b_010_42 <- d17b_30_010_42$sob
fit_17b_010_42 <- nls(y_17b_010_42 ~ (1+f*x)/(1+(9+10*f*LD10)*exp(n*log(x/LD10))),
			        start = list(LD10 = 200, n = 2))
summary(fit_17b_010_42)
ld10 <- append(ld10, summary(fit_17b_010_42)$parameters[1,1] , after=length(ld10))
se_ld10 <- append(se_ld10, summary(fit_17b_010_42)$parameters[1,2] , after=length(se_ld10))
isolado <- append(isolado, 'D17b.30.010-4.2', after=length(isolado))

# Isolado D17b.30.010-5
d17b_30_010_5 <- filter(tabgraf, isolado=='D17b.30.010-5')
y_17b_010_5 <- d17b_30_010_5$sob
fit_17b_010_5 <- nls(y_17b_010_5 ~ (1+f*x)/(1+(9+10*f*LD10)*exp(n*log(x/LD10))),
			        start = list(LD10 = 150, n = 2))
summary(fit_17b_010_5)
ld10 <- append(ld10, summary(fit_17b_010_5)$parameters[1,1] , after=length(ld10))
se_ld10 <- append(se_ld10, summary(fit_17b_010_5)$parameters[1,2] , after=length(se_ld10))
isolado <- append(isolado, 'D17b.30.010-5', after=length(isolado))

# Isolado D17b.30.R2A-3
d17b_30_R2A_3 <- filter(tabgraf, isolado=='D17b.30.R2A-3')
y_17b_R2A_3 <- d17b_30_R2A_3$sob
fit_17b_R2A_3 <- nls(y_17b_R2A_3 ~ (1+f*x)/(1+(9+10*f*LD10)*exp(n*log(x/LD10))),
			        start = list(LD10 = 100, n = 2))
summary(fit_17b_R2A_3)
ld10 <- append(ld10, summary(fit_17b_R2A_3)$parameters[1,1] , after=length(ld10))
se_ld10 <- append(se_ld10, summary(fit_17b_R2A_3)$parameters[1,2] , after=length(se_ld10))
isolado <- append(isolado, 'D17b.30.R2A-3', after=length(isolado))

# Isolado D17b.30.R2A-4.2
d17b_30_R2A_42 <- filter(tabgraf, isolado=='D17b.30.R2A-4.2')
y_17b_R2A_42 <- d17b_30_R2A_42$sob
fit_17b_R2A_42 <- nls(y_17b_R2A_42 ~ (1+f*x)/(1+(9+10*f*LD10)*exp(n*log(x/LD10))),
			        start = list(LD10 = 250, n = 2))
summary(fit_17b_R2A_42)
ld10 <- append(ld10, summary(fit_17b_R2A_42)$parameters[1,1] , after=length(ld10))
se_ld10 <- append(se_ld10, summary(fit_17b_R2A_42)$parameters[1,2] , after=length(se_ld10))
isolado <- append(isolado, 'D17b.30.R2A-4.2', after=length(isolado))

# Isolado D1c.30.010-2
d1c_30_010_2 <- filter(tabgraf, isolado=='D1c.30.010-2')
y_1c_010_2 <- d1c_30_010_2$sob
fit_1c_010_2 <- nls(y_1c_010_2 ~ (1+f*x)/(1+(9+10*f*LD10)*exp(n*log(x/LD10))),
			        start = list(LD10 = 100, n = 2))
summary(fit_1c_010_2)
ld10 <- append(ld10, summary(fit_1c_010_2)$parameters[1,1] , after=length(ld10))
se_ld10 <- append(se_ld10, summary(fit_1c_010_2)$parameters[1,2] , after=length(se_ld10))
isolado <- append(isolado, 'D1c.30.010-2', after=length(isolado))

# Isolado D21.30.010-1
d21_30_010_1 <- filter(tabgraf, isolado=='D21.30.010-1')
y_21_010_1 <- d21_30_010_1$sob
fit_21_010_1 <- nls(y_21_010_1 ~ (1+f*x)/(1+(9+10*f*LD10)*exp(n*log(x/LD10))),
			        start = list(LD10 = 100, n = 2))
summary(fit_21_010_1)
ld10 <- append(ld10, summary(fit_21_010_1)$parameters[1,1] , after=length(ld10))
se_ld10 <- append(se_ld10, summary(fit_21_010_1)$parameters[1,2] , after=length(se_ld10))
isolado <- append(isolado, 'D21.30.010-1', after=length(isolado))

# #### Isolado D4.30.R2A-3
d4_30_R2A_3 <- filter(tabgraf, isolado=='D4.30.R2A-3')
y_4_R2A_3 <- d4_30_R2A_3$sob
fit_4_R2A_3 <- nls(y_4_R2A_3 ~ (1+f*x)/(1+(9+10*f*LD10)*exp(n*log(x/LD10))),
			        start = list(LD10 = 300, n =  3))
summary(fit_4_R2A_3)
ld10 <- append(ld10, summary(fit_4_R2A_3)$parameters[1,1] , after=length(ld10))
se_ld10 <- append(se_ld10, summary(fit_4_R2A_3)$parameters[1,2] , after=length(se_ld10))
isolado <- append(isolado, 'D4.30.R2A-3', after=length(isolado))

graf_fit<- ggplot(d4_30_R2A_3, aes(x=Dose, y=sob)) +
    geom_point(shape=15, size=2)+   # Coloca as bolinhas
    geom_errorbar(aes(ymin=sob-erro, ymax=sob+erro), width=0.15, size=.5)+    # Coloca as barras de erro
    geom_smooth(method='nls', data=d4_30_R2A_3,
                formula = y ~ (1+f*x)/(1+(9+10*f*LD10)*exp(n*log(x/LD10))),
                se=FALSE,
                method.args = list(start = list(LD10 = 300, n = 3)))+
    theme_bw()+   # Muda a aparência para o tema certo
    xlab(bquote("Dose"~(J.m^-2)))+    # Coloca o nome do eixo x
    ylab(bquote("Sobrevivência"~(N/N[0])))+    # Coloca o nome do eixo y
    theme(
      panel.grid.major = element_line(size=.25),   # Define largura da linha maior do grid interno no gráfico
      panel.grid.minor = element_line(size=.1),   # Define largura da linha menor do grid interno do gráfico
      legend.position='none')
graf_fit




# Isolado D5.30.010-5/2
d5_30_010_52 <- filter(tabgraf, isolado=='D5.30.010-5/2')
y_5_010_52 <- d5_30_010_52$sob
fit_5_010_52 <- nls(y_5_010_52 ~ (1+f*x)/(1+(9+10*f*LD10)*exp(n*log(x/LD10))),
			        start = list(LD10 = 170, n = 2))
summary(fit_5_010_52)
ld10 <- append(ld10, summary(fit_5_010_52)$parameters[1,1] , after=length(ld10))
se_ld10 <- append(se_ld10, summary(fit_5_010_52)$parameters[1,2] , after=length(se_ld10))
isolado <- append(isolado, 'D5.30.010-5/2', after=length(isolado))
