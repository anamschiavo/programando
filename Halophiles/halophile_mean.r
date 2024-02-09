#Código em R para analizar levantamento de dados sobre halófilos. Foram procurados papers com a termo "halophilic", coloquei na tabela nome da espécie, domínio (representado por A=Archea, B=Bacteria, E=Eukarya),
#mínimo da concentração ótima, máximo da concentração ótima, mínimo da concentração suportada, máximo da concentração suportada, como ela é classificada (extreme, slight,moderate)
#meio de cultura, referência do paper. A ideia é ver a distribuição: quantos bixos sobrevivem/tem ótimo em cera quantidade.
#Escrito em JULHO/2019 para doutorado


#**********PACOTES E BIBLIOTECAS*********
install.packages("ggplot2") #Linha necessária apenas se não tiver pacote instalado
library(ggplot2)

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
tab_sup <- cbind.data.frame(dom_sup, mean_sup,class_sup) #junta vetores como colunas de uma data frame 1-Domínio, 2-Média valores suportados, 3-classificação
tab_sup$dom_sup <- as.factor(tab_sup$dom_sup)   #Transforma coluna de domínios em fator. Fica mais fácil de lidar nos gráficos e tabelas depois
tab_sup$class_sup <- as.factor(tab_sup$class_sup)  #Transforma classificação do autor em fator

caixa_sup <-ggplot(tab_sup, aes(x=dom_sup, y=mean_sup)) +  # box plot da concentração suportada média dividido por domínio
  geom_boxplot()+
  labs(x="Domínio", y="[NaCl] (mol/L)")+
  ggtitle("NaCl Survival Range")+
  theme_bw()




hist_sup_domain <- ggplot(tab_sup, aes(x=mean_sup, fill=dom_sup, stat(density))) +
    geom_histogram(binwidth=0.1, alpha=0.4, position="identity")+
    labs(x=expression(paste("[NaCl] mol.L"^-{1})),
    y="Frequência", fill="Domínio",
    title="Faixa de sobreviência de micro-organismos halófilos",
    subtitle="Por Domínio")+
    theme_bw()+
    theme(
    plot.title = element_text(family="", face="bold", size=16, hjust=.5), #hjust coloca título para o centro
    plot.subtitle = element_text(face="bold", size=12, hjust=0.5, colour="dark gray"),
    legend.text = element_text(size = 12),
    axis.text = element_text(colour = "black", size=12),
    axis.title = element_text(size=12))+
    scale_color_manual(values=c("#719154", "#963B36", "#278577"))+
    scale_fill_manual(values=c("#719154", "#963B36", "#278577"), labels=c("Archaea", "Bacteria", "Eukarya"))

	ggsave('hist_sup_domain_mean.png')




hist_sup <- ggplot(tab_sup, aes(x=mean_sup, stat(density))) +
    #geom_vline(xintercept=0.2, color="darkgray", size=1, linetype="dashed")+
    #geom_vline(xintercept=0.85, color="darkgray", size=1, linetype="dashed")+
    #geom_vline(xintercept=3.4, color="darkgray", size=1, linetype="dashed")+
    #geom_vline(xintercept=5.1, color="darkgray", size=1, linetype="dashed")+
    geom_histogram(binwidth=0.1, color="#e6d994ff", fill="#e6d994ff", position="identity")+  #quanto maior alpha, mais opacas as cores
    labs(x=expression(paste("[NaCl] mol.L"^-{1})), y="Frequência", title="Concentração suportada NaCl por halófilos", subtitle="Método média")+
    geom_density(alpha=.45, fill="#e6d994ff", color="#5e5e5eff")+
    coord_cartesian(xlim = c(0,6.5), ylim=c(0, 1))+
    scale_x_continuous(breaks = c(0,1,2,3,4,5,6,7))+ 
    theme_bw()+
    theme(plot.title=element_text(family="", face="bold", size=16, hjust=.5),
    plot.subtitle = element_text(face="bold", size=12, hjust=0.5, colour="#5e5e5eff"),
    legend.text = element_text(size = 12),
    legend.title = element_text(face = "bold", size=12),
    axis.text=element_text(size=10),
    axis.title=element_text(size=12))

	ggsave('hist_sup_mean.png')



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



caixa_opt_domain <-ggplot(tab_opt, aes(x=dom_opt, y=mean_opt))+ #boxplot das concentrações ótimas separada por domínio
  geom_boxplot()+
  labs(x="Domain", y="[NaCl] (mol/L)")+
  ggtitle("NaCl Concentration Optimal Growth")+
  theme_bw()



hist_opt_domain <- ggplot(tab_opt, aes(x=mean_opt, fill=dom_opt, stat(density))) + # Histograma concentração ótima separado por domínio da vida
 # geom_vline(xintercept=0.2, color="darkgray", size=1, linetype="dashed")+
 # geom_vline(xintercept=0.85, color="darkgray", size=1, linetype="dashed")+
 # geom_vline(xintercept=3.4, color="darkgray", size=1, linetype="dashed")+
 # geom_vline(xintercept=5.1, color="darkgray", size=1, linetype="dashed")+
 geom_histogram(binwidth=0.1, alpha=0.65, position="identity")+  #quanto maior alpha, mais opacas as cores
 labs(x=expression(paste("[NaCl] mol.L"^-{1})), y="Frequência", fill="Domínio", title="Concentração ótima de crescimento", subtitle="Por Domínio")+
 theme_bw()+
 theme(   plot.title = element_text(family="", face="bold", size=16, hjust=.5), #hjust coloca título para o centro
    plot.subtitle = element_text(face="bold", size=12, hjust=0.5, colour="dark gray"),
    legend.text = element_text(size = 12),
    axis.text = element_text(colour = "black", size=12),
    axis.title = element_text(size=12))+
    scale_color_manual(values=c("#719154", "#963B36", "#278577"))+
    scale_fill_manual(values=c("#719154", "#963B36", "#278577"), labels=c("Archaea", "Bacteria", "Eukarya"))

	ggsave('hist_opt_domain_mean.png')
 
 

 
hist_opt <- ggplot(tab_opt, aes(x=mean_opt, stat(density))) +
    # geom_vline(xintercept=0.2, color="darkgray", size=1, linetype="dashed")+
    # geom_vline(xintercept=0.85, color="darkgray", size=1, linetype="dashed")+  
    # geom_vline(xintercept=3.4, color="darkgray", size=1, linetype="dashed")+
    # geom_vline(xintercept=5.1, color="darkgray", size=1, linetype="dashed")+
    geom_histogram(binwidth=0.1, color="#e6d994ff", fill="#e6d994ff", position="identity")+  #quanto maior alpha, mais opacas as cores
    labs(x=expression(paste("[NaCl] mol.L"^-{1})), y="Frequência", title="Concentração ótima NaCl de halófilos", subtitle="Método média")+
    geom_density(alpha=.45, fill="#e6d994ff", color="#5e5e5eff")+
    coord_cartesian(xlim = c(0,6.5), ylim=c(0, 1))+
    scale_x_continuous(breaks = c(0,1,2,3,4,5,6,7))+ 
    theme_bw()+
    theme(plot.title=element_text(family="", face="bold", size=16, hjust=.5),
    plot.subtitle = element_text(face="bold", size=12, hjust=0.5, colour="#5e5e5eff"),
    legend.text = element_text(size = 12),
    legend.title = element_text(face = "bold", size=12),
    axis.text=element_text(size=10),
    axis.title=element_text(size=12))
 
	ggsave('hist_opt_mean.png')


 

#********** TOTAL = OPT+SUP **********
 
mean_total <- c(mean_sup, mean_opt)
dom_total <- c(dom_sup, dom_opt)
class_total <- c(class_sup, class_opt)

s <- rep("Suportado", times=length(dom_sup)) #cria vetor com Suportado com número de vezes igual ao tamanho dos outros dois vetores
o <- rep("Otimo", times=length(dom_opt))

inter <- c(s, o) #O vetor inter diz se aquela linha corresponde a um ponto do suportado ou do ótimo

tab_total <- cbind.data.frame(dom_total, mean_total, inter, class_total) #tabela que junta as informações de Domínio, Concentração, ótimo/suportado




hist_total <- ggplot(tab_total, aes(x=mean_total, fill=inter)) +
 geom_histogram(binwidth=0.1, alpha=.45, position="identity")+
 labs(x=expression(paste("[NaCl] moles.L"^-{1})), y="Count", fill="")+
 ggtitle("NaCl Concentration Optimal Growth Range")+
 theme_bw()+
 theme(plot.title=element_text(family="", face="bold", size=16), legend.text = element_text(size = 12),
legend.title = element_text(face = "bold", size=12), axis.text=element_text(size=10), axis.title=element_text(size=12))+
 scale_color_manual(values=c('#00332b', '#278577'))+
 scale_fill_manual(values=c('#00332b', '#278577'))
 
	ggsave('hist_total_mean.png')




caixa_optsup <- ggplot(tab_total, aes(x=dom_total, y= mean_total, fill=inter))+ 
    geom_boxplot(outlier.size=0.3, alpha=1)+
    labs(x="Domínio", y=expression(paste("[NaCl] mol.L"^-{1})), fill="", title="Crescimento de micro-organismos halófilos", subtitle="Por Domínio")+
    theme_bw()+
    scale_fill_manual(values=c("#278577", "#719154"), labels=c("Ótimo", "Suportado"))+
    theme(
    plot.title = element_text(face="bold", size=16, hjust=.5), #hjust coloca título para o centro
    plot.subtitle = element_text(face="bold", size=12, hjust=0.5, colour="dark gray"),
    legend.text = element_text(size = 12),
    axis.text = element_text(colour = "black", size=12),
    axis.title = element_text(size=12))

	ggsave('caixa_optsup_mean.png')
  
# hist_total_optsup <- ggplot(tab_total, aes(x=graf_total, color=inter, stat(density))) +    # histograma das concentrações todas separado por classificação
    # geom_histogram(fill="white", position="dodge", binwidth=0.2)+
    # theme_bw()  #dodge serve pra colocar os histogramas um do lado do outro (em vez de um em cima do outro)
    
    
# hist_total_optsup <- ggplot(tab_total, aes(x=graf_total, color=inter, stat(density))) +    # histograma das concentrações todas separado por classificação
    # geom_histogram(fill="white", position="identity", alpha=0.3, binwidth=0.2)+
    # theme_bw()  # alpha: quanto maior valor, mais opaco
    