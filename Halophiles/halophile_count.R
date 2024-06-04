#Código em R para analizar levantamento de dados sobre halófilos. Foram procurados papers com a termo "halophilic", coloquei na tabela nome da espécie, domínio (representado por A=Archea, B=Bacteria, E=Eukarya),
#mínimo da concentração ótima, máximo da concentração ótima, mínimo da concentração suportada, máximo da concentração suportada, como ela é classificada (extreme, slight,moderate)
#meio de cultura, referência do paper. A ideia é ver a distribuição: quantos bixos sobrevivem/tem ótimo em cera quantidade.
#Escrito em JULHO/2019 para doutorado


#**********PACOTES E BIBLIOTECAS*********
install.packages("ggplot2") #Linha necessária apenas se não tiver pacote instalado
library(ggplot2)
install.packages('plyr')
library(plyr)

#lê as duas tabelas, coloca em objetos no R
sup <- read.table("halo_sup.tsv", sep="", dec=".", header=TRUE, na.strings="NA") #lê tabela de concentração suportada
opt <- read.table("halo_opt.tsv", sep="", dec=".", header=TRUE, na.strings="NA") #lê tabela de concentração ótima

halo_sup <- na.omit(sup) #tira linhas com valores nulos
halo_opt <- na.omit(opt)


#*********** SUPORTADO **********

dom_sup <- as.vector(halo_sup[,2]) #cria um vetor de dominios correspondentes
dom_sup <- as.character(dom_sup) #transforam em vetor de character, pra depois na tabela ficar no formato certo

class_sup <- as.vector(halo_sup[,5])
class_sup <- as.character(class_sup)


valor_max_sup <- max(halo_sup[,4]) #dá o valor máximo de concentração


name_sup <- c(mode="character") #vetor emparelhado do graf_sup com os dominios de cada concentração no graf_sup
graf_sup <- c() #vetor que vai ser usado para fazer histograma
clas_sup <- c(mode="character")



concent_sup <- seq(0, valor_max_sup, 0.1) #posição dado por i. Vetor com o valor das concentrações, passo a passo. 0.01 porque é a precisão dos dados, então não perco nenhum
#contador <- vector(mode="integer", length=length(concent)) # Vetor com os contadores associados a cada concentração, então contador[i] tem a quantidade de bichos que sobrevivem na concentração concent[i]

for (n in seq_along(halo_sup[,4])){    #pega cada bicho
    for (i in seq_along(concent_sup)){  #testa cada concentração
        if (concent_sup[i] >= halo_sup[n,3] && concent_sup[i] <= halo_sup[n,4]){   #se a concentração testada está dentro dos limites do bicho, adiciona esse bicho no número de bichos que pode viver naquela condição
            #contador[i]=contador[i]+1
            graf_sup <- append(graf_sup, concent_sup[i], after=length(graf_sup))  #coloca valor de concentração no vetor do histograma. Coloca toda vez que valor está dentro do intervalo do bicho.
            name_sup <- append(name_sup, dom_sup[n], after=length(name_sup))  #para cada concentração no vetor do histograma, coloca o dominio daquele bicho no mesmo lugar no vetor name
            clas_sup <- append(clas_sup, class_sup[n], after=length(clas_sup)) #assim como name_sup, mas coloca a classificação dada pelo autor
            }
        }
    }

graf_sup <- as.numeric(graf_sup)
name_sup <- name_sup[-1] #name[1] por algum motivo bizonho era 'mode "character"', então essa linha retira esse valor, deixando name e graf com o mesmo tamanho
clas_sup <- clas_sup[-1] #
tab_sup <- cbind.data.frame(name_sup, graf_sup, clas_sup)
tab_sup$name_sup <- as.factor(tab_sup$name_sup)
tab_sup$clas_sup <- as.factor(tab_sup$clas_sup)

max(graf_sup)
hist_sup_domain <- ggplot(tab_sup, aes(x=graf_sup, fill=name_sup)) +
     geom_histogram(binwidth=0.1, alpha=0.6, position="identity")+
     labs(x=expression(paste("[NaCl] mol.L"^-{1})), y="", fill="")+
     coord_cartesian(xlim = c(0,6.5))+
     theme_bw()+
     scale_x_continuous(breaks = c(0,1,2,3,4,5,6,7))+
     theme(
     plot.title = element_text(family="", face="bold", size=16, hjust=.5), #hjust coloca título para o centro
     plot.subtitle = element_text(face="bold", size=12, hjust=0.5, colour="dark gray"),
     legend.text = element_text(size = 12),
     axis.text = element_text(colour = "black", size=12),
     axis.title = element_text(size=12),
     axis.title.x = element_text(margin = margin(t = 0, r = 0, b = 20, l = 0)))+
     scale_fill_manual(values=c("#6B2639", "#70846B"), labels=c("Archaea", "Bacteria"))

 hist_sup_domain
 ggsave('hist_sup_domain_count.png', hist_sup_domain, device='png', unit='cm', width = 12, height = 6.5, dpi=300)


 hist_sup_count <- ggplot(tab_sup) +
     geom_histogram(aes(x=graf_sup), binwidth=0.1, fill="#112225")+  #quanto maior alpha, mais opacas as cores
     labs(x=expression(paste("[NaCl] mol.L"^-{1})), y="Número de espécies")+
     coord_cartesian(xlim = c(0,6.5))+
     scale_x_continuous(breaks = c(0,1,2,3,4,5,6))+
     theme_bw()+
     theme(
     plot.title = element_text(family="", face="bold", size=16, hjust=.5), #hjust coloca título para o centro
     plot.subtitle = element_text(face="bold", size=12, hjust=0.5, colour="dark gray"),
     legend.text = element_text(size = 12),
     axis.text = element_text(colour = "black", size=12),
     axis.title = element_text(size=12),
     axis.title.y = element_text(margin = margin(t = 0, r = 20, b = 0, l = 0)),
     axis.title.x = element_text(margin = margin(t = 0, r = 0, b = 20, l = 0)))
 hist_sup_count
 ggsave('hist_sup_count.png', hist_sup_count,  device='png', unit='cm', width = 12, height = 6.5, dpi=300)


#********** ÓTIMAS **********

dom_opt <- as.vector(halo_opt[,2])
dom_opt <- as.character(dom_opt)

class_opt <- as.vector(halo_opt[,5])
class_opt <- as.character(class_opt)


valor_max_opt <- max(halo_opt[,4])


name_opt <- c(mode="character")
graf_opt <- c()
clas_opt <- c(mode="character")


concent_opt <- seq(0, valor_max_opt, 0.1) #posição dado por i. Vetor com o valor das concentrações, passo a passo. 0.01 porque é a precisão dos dados, então não perco nenhum
#contador <- vector(mode="integer", length=length(concent)) # Vetor com os contadores associados a cada concentração, então contador[i] tem a quantidade de bichos que sobrevivem na concentração concent[i]

for (n in seq_along(halo_opt[,4])){    #pega cada bicho
    for (i in seq_along(concent_opt)){  #testa cada concentração
        if (concent_opt[i] >= halo_opt[n,3] && concent_opt[i] <= halo_opt[n,4]){   #se a concentração testada está dentro dos limites do bicho, adiciona esse bicho no número de bichos que pode viver naquela condição
            #contador[i]=contador[i]+1
            graf_opt = append(graf_opt, concent_opt[i], after=length(graf_opt))  #coloca valor de concentração no vetor do histograma. Coloca toda vez que valor está dentro do intervalo do bicho.
            name_opt = append(name_opt, dom_opt[n], after=length(name_opt))  #para cada concentração no vetor do histograma, coloca o dominio daquele bicho no mesmo lugar no vetor name
            clas_opt = append(clas_opt, class_opt[n], after=length(clas_opt)) #assim como name_opt, mas coloca a classificação dada pelo autor
            }
        }
    }

graf_opt <- as.numeric(graf_opt)
name_opt <- name_opt[-1] #name[1] por algum motivo bizonho era 'mode "character"', então essa linha retira esse valor, deixando name e graf com o mesmo tamanho
clas_opt <- clas_opt[-1]
tab_opt <- cbind.data.frame(name_opt, graf_opt, clas_opt)
tab_opt$name_opt <- as.factor(tab_opt$name_opt)
tab_opt$clas_opt <- as.factor(tab_opt$clas_opt)
tab_opt$clas_opt <- factor(tab_opt$clas_opt, levels=c('Tolerant', 'Slight', 'Moderate', 'Extreme', 'No Class'))
tab_opt$clas_opt <- revalue(tab_opt$clas_opt, c('Tolerant'='Tolerante', 'Slight'='Leve', 'Moderate'='Moderado', 'Extreme'='Extremo', 'No Class'='Sem classificação'))


quartil <- quantile(graf_opt, probs = c(0,0.25,0.5,0.75,1))


 #geom_vline(xintercept=0.2, color="darkgray", size=1, linetype="dashed")+
 #geom_vline(xintercept=0.85, color="darkgray", size=1, linetype="dashed")+
 #geom_vline(xintercept=3.4, color="darkgray", size=1, linetype="dashed")+
 #geom_vline(xintercept=5.1, color="darkgray", size=1, linetype="dashed")+
max(graf_opt)

 hist_opt_domain <- ggplot(tab_opt, aes(x=graf_opt, fill=name_opt)) +
      geom_histogram(binwidth=0.1, alpha=0.6, position="identity")+
      labs(x=expression(paste("[NaCl] mol.L"^-{1})), y="", fill="")+
      coord_cartesian(xlim = c(0,5.2))+
      theme_bw()+
      scale_x_continuous(breaks = c(0,1,2,3,4,5))+
      theme(
      plot.title = element_text(family="", face="bold", size=16, hjust=.5), #hjust coloca título para o centro
      plot.subtitle = element_text(face="bold", size=12, hjust=0.5, colour="dark gray"),
      legend.text = element_text(size = 12),
      axis.text = element_text(colour = "black", size=12),
      axis.title = element_text(size=12),
      axis.title.x = element_text(margin = margin(t = 0, r = 0, b = 20, l = 0)))+
      scale_fill_manual(values=c("#973E3A", "#EC9C61"), labels=c("Archaea", "Bacteria"))

  hist_opt_domain
  ggsave('hist_opt_domain_count.png', hist_opt_domain, device='png', unit='cm', width = 12, height = 6.5, dpi=300)


  hist_opt_count <- ggplot(tab_opt) +
      geom_histogram(aes(x=graf_opt), binwidth=0.1, fill="#3C1933")+  #quanto maior alpha, mais opacas as cores
      geom_vline(xintercept=0.2, color="gray", size=.9, linetype="dotted")+
      geom_vline(xintercept=0.5, color="gray", size=.9, linetype="dotted")+
      geom_vline(xintercept=2.5, color="gray", size=.9, linetype="dotted")+
      geom_vline(xintercept=0.7, color="#70846B", size=.9, linetype="dashed")+
      geom_vline(xintercept=1.2, color="#70846B", size=.9, linetype="dashed")+
      geom_vline(xintercept=2.1, color="#70846B", size=.9, linetype="dashed")+
      labs(x=expression(paste("[NaCl] mol.L"^-{1})), y="Número de espécies")+
      coord_cartesian(xlim = c(0,5.2))+
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
  hist_opt_count
  ggsave('hist_opt_count_2.png', hist_opt_count,  device='png', unit='cm', width = 12, height = 6.5, dpi=300)

hist_opt_class <- ggplot(tab_opt, aes(x=graf_opt, fill=clas_opt)) + # Histograma concentração ótima separado por domínio da vida
    geom_histogram(binwidth=0.1, alpha=0.95, position="identity")+  #quanto maior alpha, mais opacas as cores
    geom_vline(xintercept=0.2, color="gray", size=.9, linetype="dashed")+
    geom_vline(xintercept=0.5, color="gray", size=.9, linetype="dashed")+
    geom_vline(xintercept=2.5, color="gray", size=.9, linetype="dashed")+
    facet_wrap(~clas_opt, dir='v', nrow=5)+
    labs(x=expression(paste("[NaCl] mol.L"^-{1})), y="Número de espécies", fill="Classificação pelos autores", title="", subtitle="")+
    coord_cartesian(xlim = c(0,5.2))+
    scale_x_continuous(breaks = c(0,1,2,3,4,5))+
    #scale_y_continuous(breaks = c(0, 30, 60))+
    theme_bw()+
    theme(legend.position = 'none',
      strip.text = element_text(size=14),
      axis.text = element_text(colour = "black", size=12),
      axis.title = element_text(size=14),
      panel.background=element_rect(fill='transparent', color=NA))+
      scale_fill_manual(values=c('#ec9c61', '#70846b', "#973e3a", "#6b2639", "#301631"))
hist_opt_class
ggsave('hist_opt_class.png', width=20, height=25, unit='cm', dpi=320, bg='transparent')

hist_opt_class_mine <- ggplot(tab_opt, aes(x=graf_opt, fill=clas_opt)) + # Histograma concentração ótima separado por domínio da vida
    geom_histogram(binwidth=0.1, alpha=0.95, position="identity")+  #quanto maior alpha, mais opacas as cores
    geom_vline(xintercept=0.7, color="gray", size=.9, linetype="dashed")+
    geom_vline(xintercept=1.2, color="gray", size=.9, linetype="dashed")+
    geom_vline(xintercept=2.1, color="gray", size=.9, linetype="dashed")+
    facet_wrap(~clas_opt, dir='v', nrow=5)+
    labs(x=expression(paste("[NaCl] mol.L"^-{1})), y="", fill="Classificação pelos autores", title="", subtitle="")+
    coord_cartesian(xlim = c(0,5.2))+
    scale_x_continuous(breaks = c(0,1,2,3,4,5))+
    #scale_y_continuous(breaks = c(0, 30, 60))+
    theme_bw()+
    theme(legend.position = 'none',
      strip.text = element_text(size=14),
      axis.text = element_text(colour = "black", size=12),
      axis.title = element_text(size=14),
      panel.background=element_rect(fill='transparent', color=NA))+
      scale_fill_manual(values=c('#ec9c61', '#70846b', "#973e3a", "#6b2639", "#301631"))
hist_opt_class_mine
ggsave('hist_opt_class_mine.png', width=20, height=25, unit='cm', dpi=320, bg='transparent')



#********** TOTAL = OPT+SUP **********

graf_total <- c(graf_sup, graf_opt)
name_total <- c(name_sup, name_opt)
class_total <- c(class_sup, class_opt)

s <- rep("Suportado", times=length(name_sup)) #cria vetor com Suportado com número de vezes igual ao tamanho dos outros dois vetores
o <- rep("Ótimo", times=length(name_opt))

inter <- c(s, o) #O vetor inter diz se aquela linha corresponde a um ponto do suportado ou do ótimo

tab_total <- cbind.data.frame(name_total, graf_total, inter) #tabela que junta as informações de Domínio, Concentração, ótimo/suportado



hist_total <- ggplot(tab_total, aes(x=graf_total, fill=inter)) +
    geom_histogram(binwidth=0.1, alpha=.55, position="identity")+
    labs(x=expression(paste("[NaCl] mol.L"^-{1})), y="Número de espécies", fill="")+
    coord_cartesian(xlim = c(0,6.5))+
    scale_x_continuous(breaks = c(0,1,2,3,4,5,6))+
    theme_bw()+
    theme(
    plot.title = element_text(family="", face="bold", size=16, hjust=.5), #hjust coloca título para o centro
    plot.subtitle = element_text(face="bold", size=12, hjust=0.5, colour="dark gray"),
    legend.text = element_text(size = 12),
    axis.text = element_text(colour = "black", size=12),
    axis.title = element_text(size=12),
    axis.title.y = element_text(margin = margin(t = 0, r = 20, b = 0, l = 0)))+
    scale_fill_manual(values=c('#3C1933', '#112225'))
hist_total
ggsave('hist_total_mean.png', hist_total,  device='png', unit='cm', width = 12, height = 6.5, dpi=300)

#********** GRID **********
grid_count <- grid.arrange(hist_sup_count, hist_sup_domain,hist_opt_count, hist_opt_domain,hist_total,
                      layout_matrix=matrix(c(1,1,1,2,2,2,2,3,3,3,4,4,4,4,NA,5,5,5,5,5,NA), byrow=TRUE, ncol=7))

ggsave('grid_count.png', grid_count,  device='png', unit='cm', width = 18, height = 27, dpi=300)

grid_class <- grid.arrange(hist_opt_class, hist_opt_class_mine, ncol=2)
