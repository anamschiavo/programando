#-------MANIPULAÇÃO DAS TABELAS DO PLATE READER E PLOTS-----------
#------------ Micrococcus luteus?D7.30.100-9 em NaCl-----------------------

# Carregar pacotes necessários
install.packages("readxl")
library(readxl)
library(dplyr)
install.packages("openxlsx")
library(openxlsx)
install.packages('matrixStats')
library(matrixStats)

# PRÉ-PROCESSAMENTO
# Carregar os dados brutos
raw <- read_excel("/Users/ana_m/OneDrive/Desktop/USP et al/Limbo/Micrococcus/19isolados_23-nov-2024 18-43-03.xlsx", sheet = "Plate 1 - Raw Data")

# Separa coluna de Tempo
tempo <- seq(0, 24, 1/3)

# Criar tabelas das concentrações
# 0M NaCl
branco_0M <- raw$G2
corrig_0M <- cbind(raw$B2-branco_0M, raw$C2-branco_0M, raw$D2-branco_0M, raw$E2-branco_0M, raw$F2-branco_0M)
mean_0M <- rowMeans(corrig_0M)
sd_0M <- rowSds(corrig_0M)
tab0 <- cbind(tempo, mean_0M, sd_0M)

# 0.25M NaCl
branco_025M <- raw$G3
corrig_025M <- cbind(raw$B3-branco_025M, raw$C3-branco_025M, raw$D3-branco_025M, raw$E3-branco_025M, raw$F3-branco_025M)
mean_025M <- rowMeans(corrig_025M)
sd_025M <- rowSds(corrig_025M)
tab025 <- cbind(tempo, mean_025M, sd_025M)

# 0.5M NaCl
branco_05M <- raw$G4
corrig_05M <- cbind(raw$B4-branco_05M, raw$C4-branco_05M, raw$D4-branco_05M, raw$E4-branco_05M, raw$F4-branco_05M)
mean_05M <- rowMeans(corrig_05M)
sd_05M <- rowSds(corrig_05M)
tab05 <- cbind(tempo, mean_05M, sd_05M)

# 0.75M NaCl
branco_075M <- raw$G5
corrig_075M <- cbind(raw$B5-branco_075M, raw$C5-branco_075M, raw$D5-branco_075M, raw$E5-branco_075M, raw$F5-branco_075M)
mean_075M <- rowMeans(corrig_075M)
sd_075M <- rowSds(corrig_075M)
tab075 <- cbind(tempo, mean_075M, sd_075M)

# 1M NaCl
branco_1M <- raw$G6
corrig_1M <- cbind(raw$B6-branco_1M, raw$C6-branco_1M, raw$D6-branco_1M, raw$E6-branco_1M, raw$F6-branco_1M)
mean_1M <- rowMeans(corrig_1M)
sd_1M <- rowSds(corrig_1M)
tab1 <- cbind(tempo, mean_1M, sd_1M)

# 1.25M NaCl
branco_125M <- raw$G7
corrig_125M <- cbind(raw$B7-branco_125M, raw$C7-branco_125M, raw$D7-branco_125M, raw$E7-branco_125M, raw$F7-branco_125M)
mean_125M <- rowMeans(corrig_125M)
sd_125M <- rowSds(corrig_125M)
tab125 <- cbind(tempo, mean_125M, sd_125M)

# 1.5M NaCl
branco_15M <- raw$G8
corrig_15M <- cbind(raw$B8-branco_15M, raw$C8-branco_15M, raw$D8-branco_15M, raw$E8-branco_15M, raw$F8-branco_15M)
mean_15M <- rowMeans(corrig_15M)
sd_15M <- rowSds(corrig_15M)
tab15 <- cbind(tempo, mean_15M, sd_15M)

# 2M NaCl
branco_2M <- raw$G9
corrig_2M <- cbind(raw$B9-branco_2M, raw$C9-branco_2M, raw$D9-branco_2M, raw$E9-branco_2M, raw$F9-branco_2M)
mean_2M <- rowMeans(corrig_2M)
sd_2M <- rowSds(corrig_2M)
tab2 <- cbind(tempo, mean_2M, sd_2M)

# 3M NaCl
branco_3M <- raw$G10
corrig_3M <- cbind(raw$B10-branco_3M, raw$C10-branco_3M, raw$D10-branco_3M, raw$E10-branco_3M, raw$F10-branco_3M)
mean_3M <- rowMeans(corrig_3M)
sd_3M <- rowSds(corrig_3M)
tab3 <- cbind(tempo, mean_3M, sd_3M)

# 4M NaCl
branco_4M <- raw$G11
corrig_4M <- cbind(raw$B11-branco_4M, raw$C11-branco_4M, raw$D11-branco_4M, raw$E11-branco_4M, raw$F11-branco_4M)
mean_4M <- rowMeans(corrig_4M)
sd_4M <- rowSds(corrig_4M)
tab4 <- cbind(tempo, mean_4M, sd_4M)

#-----GRÁFICO---------------

library(ggplot2)
library(tidyr)
library(dplyr)

graf_0 <- ggplot(tab0, aes(x=tempo, y=mean_0M))+
        geom_point(size=1.5)+
        geom_errorbar(aes(ymin=mean_0M-sd_0M, ymax=mean_0M+sd_0M))

graf_025 <- ggplot(tab025, aes(x=tempo, y=mean_025M))+
        geom_point(size=1.5)+
        geom_errorbar(aes(ymin=mean_025M-sd_025M, ymax=mean_025M+sd_025M))

graf_05 <- ggplot(tab05, aes(x=tempo, y=mean_05M))+
        geom_point(size=1.5)+
        geom_errorbar(aes(ymin=mean_05M-sd_05M, ymax=mean_05M+sd_05M))

graf_075 <- ggplot(tab075, aes(x=tempo, y=mean_075M))+
        geom_point(size=1.5)+
        geom_errorbar(aes(ymin=mean_075M-sd_075M, ymax=mean_075M+sd_075M))

graf_1 <- ggplot(tab1, aes(x=tempo, y=mean_1M))+
        geom_point(size=1.5)+
        geom_errorbar(aes(ymin=mean_1M-sd_1M, ymax=mean_1M+sd_1M))

graf_125 <- ggplot(tab125, aes(x=tempo, y=mean_125M))+
        geom_point(size=1.5)+
        geom_errorbar(aes(ymin=mean_125M-sd_125M, ymax=mean_125M+sd_125M))

graf_15 <- ggplot(tab15, aes(x=tempo, y=mean_15M))+
        geom_point(size=1.5)+
        geom_errorbar(aes(ymin=mean_15M-sd_15M, ymax=mean_15M+sd_15M))

graf_2 <- ggplot(tab2, aes(x=tempo, y=mean_2M))+
        geom_point(size=1.5)+
        geom_errorbar(aes(ymin=mean_2M-sd_2M, ymax=mean_2M+sd_2M))

graf_3 <- ggplot(tab3, aes(x=tempo, y=mean_3M))+
        geom_point(size=1.5)+
        geom_errorbar(aes(ymin=mean_3M-sd_3M, ymax=mean_3M+sd_3M))

graf_4 <- ggplot(tab4, aes(x=tempo, y=mean_4M))+
        geom_point(size=1.5)+
        geom_errorbar(aes(ymin=mean_4M-sd_4M, ymax=mean_4M+sd_4M))

plot_medias <- ggsave("/Users/isabe/Downloads/nepalensis/Curva_nacl_DO/colB_22_7/plot_medias.png", plot = last_plot(), width = 8, height = 6, dpi = 300, bg = "white")
