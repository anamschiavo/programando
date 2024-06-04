#Código em R para analizar levantamento de dados sobre halófilos. Foram procurados papers com a termo "halophilic", coloquei na tabela nome da espécie, domínio (representado por A=Archea, B=Bacteria, E=Eukarya),
#mínimo da concentração ótima, máximo da concentração ótima, mínimo da concentração suportada, máximo da concentração suportada, como ela é classificada (extreme, slight,moderate)
#meio de cultura, referência do paper. A ideia é ver a distribuição: quantos bixos sobrevivem/tem ótimo em cera quantidade.
#Escrito em JULHO/2019 para doutorado


#**********PACOTES E BIBLIOTECAS*********
install.packages("ggplot2") #Linha necessária apenas se não tiver pacote instalado
library(ggplot2)
install.packages("gridExtra")
library(gridExtra)

#lê as duas tabelas, coloca em objetos no R
sup <- read.table("halo_sup.tsv", sep="", dec=".", header=TRUE, na.strings="NA") #lê tabela de concentração suportada
opt <- read.table("halo_opt.tsv", sep="", dec=".", header=TRUE, na.strings="NA") #lê tabela de concentração ótima

halo_sup <- na.omit(sup) #tira linhas com valores nulos
halo_opt <- na.omit(opt)


#*********** SUPORTADO **********

dom_sup <- as.vector(halo_sup[,2]) #cria um vetor de dominios correspondentes
dom_sup <- as.character(dom_sup) #transforam em vetor de character, pra depois na tabela ficar no formato certo

class_sup <- as.vector(halo_sup[,5]) #Cria um vetor de classificação dada pelo autor
class_sup <- as.character(class_sup) #transforam em vetor de character, pra depois na tabela ficar no formato certo

mean_sup <- c() #Vetor no qual será colocado a média entre concentrações suportadas mínimas e máximas na mesma ordem da tabela

for (n in seq_along(halo_sup[,4])){    #pega cada bicho
    media <- (halo_sup[n,3]+halo_sup[n,4])/2    #faz a média entre o mínimo e o máximo, tornando um único valor
    mean_sup <-append(mean_sup, media, after=length(mean_sup))  #coloca valor da média em um vetor na mesma ordem da tabela
    }

mean_sup <- as.numeric(mean_sup)
tab_sup <- cbind.data.frame(dom_sup, mean_sup, class_sup) #junta vetores como colunas de uma data frame 1-Domínio, 2-Média valores suportados, 3-classificação
tab_sup$dom_sup <- as.factor(tab_sup$dom_sup)   #Transforma coluna de domínios em fator. Fica mais fácil de lidar nos gráficos e tabelas depois
tab_sup$class_sup <- as.factor(tab_sup$class_sup)  #Transforma classificação do autor em fator


hist_sup_domain <- ggplot(tab_sup, aes(x=mean_sup, fill=dom_sup)) +
    geom_histogram(binwidth=0.1, alpha=0.6, position="identity")+
    labs(x=expression(paste("[NaCl] mol.L"^-{1})),
    y="", fill="")+
    theme_bw()+
    scale_x_continuous(breaks = c(0,1,2,3,4,5))+
    theme(
    plot.title = element_text(family="", face="bold", size=16, hjust=.5), #hjust coloca título para o centro
    plot.subtitle = element_text(face="bold", size=12, hjust=0.5, colour="dark gray"),
    legend.text = element_text(size = 12),
    axis.text = element_text(colour = "black", size=12),
    axis.title = element_text(size=12),
    axis.title.x = element_text(margin = margin(t = 0, r = 0, b = 20, l = 0)))+
    scale_color_manual(values=c("#6B2639", "#70846B"))+
    scale_fill_manual(values=c("#6B2639", "#70846B"), labels=c("Archaea", "Bacteria"))

hist_sup_domain
ggsave('hist_sup_domain_mean.png', hist_sup_domain, device='png', unit='cm', width = 12, height = 6.5, dpi=300)




hist_sup_mean <- ggplot(tab_sup) +
    geom_histogram(aes(x=mean_sup), binwidth=0.1, fill="#112225")+  #quanto maior alpha, mais opacas as cores
    labs(x=expression(paste("[NaCl] mol.L"^-{1})), y="Número de espécies")+
    coord_cartesian(xlim = c(0,5))+
    scale_x_continuous(breaks = c(0,1,2,3,4,5))+
    theme_bw()+
    theme(
    plot.title = element_text(family="", face="bold", size=16, hjust=.5), #hjust coloca título para o centro
    plot.subtitle = element_text(face="bold", size=12, hjust=0.5, colour="dark gray"),
    legend.text = element_text(size = 12),
    axis.text = element_text(colour = "black", size=12),
    axis.title = element_text(size=12),
    axis.title.y = element_text(margin = margin(t = 0, r = 20, b = 0, l = 0)),
    axis.title.x = element_text(margin = margin(t = 0, r = 0, b = 20, l = 0)))
hist_sup_mean
ggsave('hist_sup_mean.png', hist_sup_mean,  device='png', unit='cm', width = 12, height = 6.5, dpi=300)



#********** ÓTIMAS **********

dom_opt <- as.vector(halo_opt[,2]) #Cria vetor com domínios
dom_opt <- as.character(dom_opt)

class_opt <- as.vector(halo_opt[,5]) #Cria vetor com classificação dada pelo autor
class_opt <- as.character(class_opt)

mean_opt <- c()  #Vetor onde serão colocadas as médias entre concentração ótimas máxima e mínima

for (n in seq_along(halo_opt[,4])){    #pega cada bicho
    media <- (halo_opt[n,3]+halo_opt[n,4])/2    #faz a média entre o mínimo e o máximo, tornando um único valor
    mean_opt <-append(mean_opt, media, after=length(mean_opt))  #coloca valor da média em um vetor na mesma ordem da tabela
    }

mean_opt <- as.numeric(mean_opt)
tab_opt <- cbind.data.frame(dom_opt, mean_opt,class_opt) #junta vetores como colunas de uma data frame 1-Domínio, 2-Média valores suportados, 3-classificação
tab_opt$dom_opt <- as.factor(tab_opt$dom_opt)   #Transforma coluna de domínios em fator. Fica mais fácil de lidar nos gráficos e tabelas depois
tab_opt$class_opt <- as.factor(tab_opt$class_opt)  #Transforma classificação do autor em fator


hist_opt_domain <- ggplot(tab_opt, aes(x=mean_opt, fill=dom_opt)) +
    geom_histogram(binwidth=0.1, alpha=0.6, position="identity")+
    labs(x=expression(paste("[NaCl] mol.L"^-{1})),
    y="", fill="")+
    scale_x_continuous(breaks = c(0,1,2,3,4,5))+
    theme_bw()+
    theme(
    plot.title = element_text(family="", face="bold", size=16, hjust=.5), #hjust coloca título para o centro
    plot.subtitle = element_text(face="bold", size=12, hjust=0.5, colour="dark gray"),
    legend.text = element_text(size = 12),
    axis.text = element_text(colour = "black", size=12),
    axis.title = element_text(size=12),
    axis.title.x = element_text(margin = margin(t = 0, r = 0, b = 20, l = 0)))+
    scale_fill_manual(values=c("#973E3A", "#EC9C61"), labels=c("Archaea", "Bacteria"))

hist_opt_domain
	ggsave('hist_opt_domain_mean.png', hist_opt_domain, device='png', unit='cm', width = 12, height = 6.5, dpi=300)




hist_opt_mean <- ggplot(tab_opt) +
      geom_histogram(aes(x=mean_opt), binwidth=0.1, fill="#3C1933")+  #quanto maior alpha, mais opacas as cores
      labs(x=expression(paste("[NaCl] mol.L"^-{1})), y="Número de espécies")+
      coord_cartesian(xlim = c(0,5))+
      scale_x_continuous(breaks = c(0,1,2,3,4,5))+
      theme_bw()+
      theme(
      plot.title = element_text(family="", face="bold", size=16, hjust=.5), #hjust coloca título para o centro
      plot.subtitle = element_text(face="bold", size=12, hjust=0.5, colour="dark gray"),
      legend.text = element_text(size = 12),
      axis.text = element_text(colour = "black", size=12),
      axis.title = element_text(size=12),
      axis.title.y = element_text(margin = margin(t = 0, r = 20, b = 0, l = 0)),
      axis.title.x = element_text(margin = margin(t = 0, r = 0, b = 20, l = 0)))
  hist_opt_mean
  ggsave('hist_opt_mean.png', hist_opt_mean,  device='png', unit='cm', width = 12, height = 6.5, dpi=300)



#********** TOTAL = OPT+SUP **********

mean_total <- c(mean_sup, mean_opt)
dom_total <- c(dom_sup, dom_opt)
class_total <- c(class_sup, class_opt)

s <- rep("Suportado", times=length(dom_sup)) #cria vetor com Suportado com número de vezes igual ao tamanho dos outros dois vetores
o <- rep("Otimo", times=length(dom_opt))

inter <- c(s, o) #O vetor inter diz se aquela linha corresponde a um ponto do suportado ou do ótimo

tab_total <- cbind.data.frame(dom_total, mean_total, inter, class_total) #tabela que junta as informações de Domínio, Concentração, ótimo/suportado


hist_total <- ggplot(tab_total, aes(x=mean_total, fill=inter)) +
    geom_histogram(binwidth=0.1, alpha=.55, position="identity")+
    labs(x=expression(paste("[NaCl] mol.L"^-{1})), y="Número de espécies", fill="")+
    theme_bw()+
    scale_x_continuous(breaks = c(0,1,2,3,4,5))+
    theme(
    plot.title = element_text(family="", face="bold", size=16, hjust=.5), #hjust coloca título para o centro
    plot.subtitle = element_text(face="bold", size=12, hjust=0.5, colour="dark gray"),
    legend.text = element_text(size = 12),
    axis.text = element_text(colour = "black", size=12),
    axis.title = element_text(size=12),
    axis.title.y = element_text(margin = margin(t = 0, r = 20, b = 0, l = 0)))+
    scale_fill_manual(values=c('#3C1933', '#112225'), labels=c("Ótimo", "Suportado"))
hist_total
ggsave('hist_total_mean.png', hist_total,  device='png', unit='cm', width = 12, height = 6.5, dpi=300)

#********** GRID **********
grid_mean <- grid.arrange(hist_sup_mean, hist_sup_domain,hist_opt_mean, hist_opt_domain,hist_total,
                      layout_matrix=matrix(c(1,1,1,2,2,2,2,3,3,3,4,4,4,4,NA,5,5,5,5,5,NA), byrow=TRUE, ncol=7))

ggsave('grid_mean.png', grid_mean,  device='png', unit='cm', width = 18, height = 27, dpi=300)
