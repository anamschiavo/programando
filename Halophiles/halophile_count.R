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



concent_sup <- seq(0, valor_max_sup, 0.01) #posição dado por i. Vetor com o valor das concentrações, passo a passo. 0.01 porque é a precisão dos dados, então não perco nenhum
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


caixa_sup <-ggplot(tab_sup, aes(x=name_sup, y=graf_sup)) +
  geom_boxplot()+
  labs(x="Domínio", y="[NaCl] (mol/L)")+
  ggtitle("NaCl Survival Range")+
  theme_bw()




hist_sup_domain <- ggplot(tab_sup, aes(x=graf_sup, fill=name_sup)) +
    geom_histogram(binwidth=0.01, alpha=0.6, position="identity")+
    labs(x=expression(paste("[NaCl] mol.L"^-{1})),
    y="", fill="Domínio",
    title="Concentração suportada NaCl por halófilos")+
    coord_cartesian(xlim = c(0,6.5))+
    scale_x_continuous(breaks = c(0,1,2,3,4,5,6,7))+
    theme_bw()+
    theme(
    plot.title = element_text(family="", face="bold", size=16, hjust=.5), #hjust coloca título para o centro
    plot.subtitle = element_text(face="bold", size=12, hjust=0.5, colour="dark gray"),
    legend.text = element_text(size = 12),
    axis.text = element_text(colour = "black", size=12),
    axis.title = element_text(size=12))+
    scale_color_manual(values=c('#301631', '#973e3a'))+
    scale_fill_manual(values=c('#301631', '#973e3a'), labels=c("Archaea", "Bacteria"))
hist_sup_domain
	ggsave('hist_sup_domain.png', hist_sup_domain, device='png', unit='cm', width=16, height=14)


hist_sup <- ggplot(tab_sup, aes(x=graf_sup)) +
 #geom_vline(xintercept=0.2, color="darkgray", size=1, linetype="dashed")+
 #geom_vline(xintercept=0.85, color="darkgray", size=1, linetype="dashed")+
 #geom_vline(xintercept=3.4, color="darkgray", size=1, linetype="dashed")+
 #geom_vline(xintercept=5.1, color="darkgray", size=1, linetype="dashed")+
    geom_histogram(binwidth=0.01, color="#973e3a", fill="#EC9C61", position="identity")+  #quanto maior alpha, mais opacas as cores
    labs(x=expression(paste("[NaCl] mol.L"^-{1})), y="N° de espécies", title="Concentração suportada NaCl por halófilos", subtitle="Método contagem")+
    #geom_density(alpha=.45, fill="#f6e8a1ff", color="#5e5e5eff")+
    coord_cartesian(xlim = c(0,6.5))+
    scale_x_continuous(breaks = c(0,1,2,3,4,5,6,7))+
    theme_bw()+
    theme(plot.title=element_text(family="", face="bold", size=16, hjust=.5),
    plot.subtitle = element_text(face="bold", size=12, hjust=0.5, colour="dark gray"),
    legend.text = element_text(size = 12),
    legend.title = element_text(face = "bold", size=12),
    axis.text=element_text(size=10),
    axis.title=element_text(size=12))
hist_sup
	ggsave('hist_sup.png', hist_sup, device='png', unit='cm', width=16, height=14)


#********** ÓTIMAS **********

dom_opt <- as.vector(halo_opt[,2])
dom_opt <- as.character(dom_opt)

class_opt <- as.vector(halo_opt[,5])
class_opt <- as.character(class_opt)


valor_max_opt <- max(halo_opt[,4])


name_opt <- c(mode="character")
graf_opt <- c()
clas_opt <- c(mode="character")


concent_opt <- seq(0, valor_max_opt, 0.01) #posição dado por i. Vetor com o valor das concentrações, passo a passo. 0.01 porque é a precisão dos dados, então não perco nenhum
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



caixa_opt_domain <-ggplot(tab_opt, aes(x=name_opt, y=graf_opt))+ #boxplot das concentrações ótimas separada por domínio
  geom_boxplot()+
  labs(x="Domain", y="[NaCl] (mol/L)")+
  ggtitle("NaCl Concentration Optimal Growth")+
  theme_bw()

hist_opt_domain <- ggplot(tab_opt, aes(x=graf_opt, fill=name_opt)) + # Histograma concentração ótima separado por domínio da vida
 # geom_vline(xintercept=0.2, color="darkgray", size=1, linetype="dashed")+
 # geom_vline(xintercept=0.85, color="darkgray", size=1, linetype="dashed")+
 # geom_vline(xintercept=3.4, color="darkgray", size=1, linetype="dashed")+
 # geom_vline(xintercept=5.1, color="darkgray", size=1, linetype="dashed")+
    geom_histogram(binwidth=0.01, alpha=0.6, position="identity")+  #quanto maior alpha, mais opacas as cores
    labs(x=expression(paste("[NaCl] mol.L"^-{1})), y="", fill="Domínio", title="Concentração ótima de crescimento de halófilos")+
    coord_cartesian(xlim = c(0,5.5))+
    scale_x_continuous(breaks = c(0,1,2,3,4,5))+
    theme_bw()+
    theme(plot.title = element_text(family="", face="bold", size=16, hjust=.5), #hjust coloca título para o centro
    plot.subtitle = element_text(face="bold", size=12, hjust=0.5, colour="dark gray"),
    legend.text = element_text(size = 12),
    axis.text = element_text(colour = "black", size=12),
    axis.title = element_text(size=12),
    panel.background=element_rect(fill='transparent', color=NA))+
    scale_color_manual(values=c('#6b2639', '#70846b'))+
    scale_fill_manual(values=c('#6b2639', '#70846b'), labels=c("Archaea", "Bacteria"))
hist_opt_domain
	ggsave('hist_opt_domain.png', hist_opt_domain, device='png', unit='cm', width=16, height=14)



hist_opt_class <- ggplot(tab_opt, aes(x=graf_opt, fill=clas_opt)) + # Histograma concentração ótima separado por domínio da vida
    geom_vline(xintercept=0.2, color="gray", size=.9, linetype="dashed")+
    geom_vline(xintercept=0.85, color="gray", size=.9, linetype="dashed")+
    geom_vline(xintercept=3.4, color="gray", size=.9, linetype="dashed")+
    geom_vline(xintercept=5.1, color="gray", size=.9, linetype="dashed")+
    geom_histogram(binwidth=0.01, alpha=0.8, position="identity")+  #quanto maior alpha, mais opacas as cores
    facet_wrap(~clas_opt, dir='v', nrow=5)+
    labs(x=expression(paste("[NaCl] mol.L"^-{1})), y="", fill="Classificação pelos autores", title="", subtitle="")+
    coord_cartesian(xlim = c(0,5.5))+
    scale_x_continuous(breaks = c(0,1,2,3,4,5))+
    scale_y_continuous(breaks = c(0, 30, 60))+
    theme_bw()+
    theme(legend.position = 'none',
      strip.text = element_text(size=16, face='bold'),
      axis.text = element_text(colour = "black", size=12),
      axis.title = element_text(size=16),
      panel.background=element_rect(fill='transparent', color=NA))+
      scale_fill_manual(values=c('#ec9c61', '#70846b', "#973e3a", "#6b2639", "#301631"))
hist_opt_class
ggsave('hist_opt_class.png', width=20, height=25, unit='cm', dpi=320, bg='transparent')

hist_opt_nova_class <- ggplot(tab_opt, aes(x=graf_opt, fill=clas_opt)) + # Histograma concentração ótima separado por domínio da vida
    geom_vline(xintercept=0.65, color="gray", size=.9, linetype="dashed")+
    geom_vline(xintercept=1.21, color="gray", size=.9, linetype="dashed")+
    geom_vline(xintercept=1.99, color="gray", size=.9, linetype="dashed")+
    geom_vline(xintercept=5.2, color="gray", size=.9, linetype="dashed")+
    geom_histogram(binwidth=0.01, alpha=0.8, position="identity")+  #quanto maior alpha, mais opacas as cores
    facet_wrap(~clas_opt, dir='v', nrow=5)+
    labs(x=expression(paste("[NaCl] mol.L"^-{1})), y="", fill="Classificação pelos autores", title="", subtitle="")+
    coord_cartesian(xlim = c(0,5.5))+
    scale_x_continuous(breaks = c(0,1,2,3,4,5))+
    scale_y_continuous(breaks = c(0, 30, 60))+
    theme_bw()+
    theme(legend.position = 'none',
      strip.text = element_text(size=16, face='bold'),
      axis.text = element_text(colour = "black", size=12),
      axis.title = element_text(size=16),
      panel.background=element_rect(fill='transparent', color=NA))+
      scale_fill_manual(values=c('#ec9c61', '#70846b', "#973e3a", "#6b2639", "#301631"))
hist_opt_nova_class
ggsave('hist_opt_nova_class.png', width=20, height=25, unit='cm', dpi=320, bg='transparent')


hist_opt <- ggplot(tab_opt, aes(x=graf_opt)) +
    # geom_vline(xintercept=0.2, color="darkgray", size=1, linetype="dashed")+
    # geom_vline(xintercept=0.85, color="darkgray", size=1, linetype="dashed")+
    # geom_vline(xintercept=3.4, color="darkgray", size=1, linetype="dashed")+
    # geom_vline(xintercept=5.1, color="darkgray", size=1, linetype="dashed")+
    geom_histogram(binwidth=0.01, color="#6b2639", fill="#ec9c61", position="identity")+  #quanto maior alpha, mais opacas as cores
    labs(x=expression(paste("[NaCl] mol.L"^-{1})), y="Frequência", title="Concentração ótima NaCl de halófilos", subtitle="Método contagem")+
    #geom_density(alpha=.35, fill="#f6e8a1ff", color="#5e5e5eff")+
    coord_cartesian(xlim = c(0,5.5))+
    scale_x_continuous(breaks = c(0,1,2,3,4,5))+
    theme_bw()+
    theme(plot.title=element_text(family="", face="bold", size=16, hjust=.5),
    plot.subtitle = element_text(face="bold", size=12, hjust=0.5, colour="dark gray"),
    legend.text = element_text(size = 12),
    legend.title = element_text(face = "bold", size=12),
    axis.text=element_text(size=10),
    axis.title=element_text(size=12))
hist_opt

	ggsave('hist_opt.png', hist_opt, device='png', unit='cm', width=16, height=14)



#********** TOTAL = OPT+SUP **********

graf_total <- c(graf_sup, graf_opt)
name_total <- c(name_sup, name_opt)
class_total <- c(class_sup, class_opt)

s <- rep("Suportado", times=length(name_sup)) #cria vetor com Suportado com número de vezes igual ao tamanho dos outros dois vetores
o <- rep("Otimo", times=length(name_opt))

inter <- c(s, o) #O vetor inter diz se aquela linha corresponde a um ponto do suportado ou do ótimo

tab_total <- cbind.data.frame(name_total, graf_total, inter) #tabela que junta as informações de Domínio, Concentração, ótimo/suportado



hist_total <- ggplot(tab_total, aes(x=graf_total, fill=inter)) +
 geom_histogram(binwidth=0.01, alpha=.7, position="identity")+
 labs(x=expression(paste("[NaCl] moles.L"^-{1})), y="Count", fill="")+
 ggtitle("NaCl Concentration Optimal Growth Range")+
 theme_bw()+
 theme(plot.title=element_text(family="", face="bold", size=16), legend.text = element_text(size = 12),
legend.title = element_text(face = "bold", size=12), axis.text=element_text(size=10), axis.title=element_text(size=12))+
 scale_color_manual(values=c('#e07761ff', '#f9c78cff'))+
 scale_fill_manual(values=c('#e07761ff', '#f9c78cff'))

	ggsave('hist_total.png')


caixa_optsup <- ggplot(tab_total, aes(x=name_total, y= graf_total, fill=inter))+
    geom_boxplot(outlier.size=0.4, alpha=1)+
    labs(x="Domínio", y=expression(paste("[NaCl] mol.L"^-{1})), fill="", title="", subtitle="")+
    theme_bw()+
    scale_fill_manual(values=c('#35b779', '#31688e'), labels=c("Crescimento Ótimo", "Faixa de Sobrevivência"))+
    scale_x_discrete(labels=c('Archaea', 'Bacteria'))+
    theme(
    legend.position='top',
    plot.title = element_text(face="bold", size=16, hjust=.5), #hjust coloca título para o centro
    plot.subtitle = element_text(face="bold", size=12, hjust=0.5, colour="dark gray"),
    legend.text = element_text(size = 12),
    axis.text = element_text(size=10),
    axis.title = element_text(size=12))

	ggsave('caixa_optsup.png', width=12, height=14, unit='cm', dpi=320, bg='transparent')


  violin_optsup <- ggplot(tab_total, aes(x=name_total, y= graf_total, fill=inter))+
      geom_violin(outlier.size=0.4, alpha=1)+
      labs(x="Domínio", y=expression(paste("[NaCl] mol.L"^-{1})), fill="", title="", subtitle="")+
      theme_bw()+
      scale_fill_manual(values=c('#35b779', '#31688e'), labels=c("Crescimento Ótimo", "Faixa de Sobrevivência"))+
      scale_x_discrete(labels=c('Archaea', 'Bacteria'))+
      theme(
      legend.position='top',
      plot.title = element_text(face="bold", size=16, hjust=.5), #hjust coloca título para o centro
      plot.subtitle = element_text(face="bold", size=12, hjust=0.5, colour="dark gray"),
      legend.text = element_text(size = 12),
      axis.text = element_text(size=10),
      axis.title = element_text(size=12))

  	ggsave('violin_optsup.png', width=12, height=14, unit='cm', dpi=320, bg='transparent')
