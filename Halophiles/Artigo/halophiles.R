# Scrip de pré-processamento da tabela de halófilos da literatura

# REGRAS DE ENTRADA:
# ***********ORDEM DAS COLUNAS: Espécie|Domínio|Opt min|Opt max|Sup sup|Sup max|Classificação|É NaCl?|Unidade de concentração***************
# Concentrações não podem ter símbolo %, se for dado assim no paper, coloque "p" na coluna 7 (Unidade de concentração)
# Decimal é dado por vírgula
# Coluna 7 (Unidade de concentração) pode ser p = % massa/volume | g = g/L | m = mol/L
# Coluna 6 (É NaCl?) pode ser "s" para sim. Qualquer outra entrada diz que é outro sal ou combinação: ISSO FARÁ COM QUE SEJA RETIRADO!
# Coluna 2 (Dominío) pode ser B = bacteria | A = archaea | E = eukarya

#********** PACOTES E BIBLOTECAS **********
install.packages('plyr')
library(plyr)
install.packages("dplyr")   #Linha necessária apenas se não tiver pacote já instalado
library(dplyr)
install.packages("ggplot2") #Linha necessária apenas se não tiver pacote instalado
library(ggplot2)
install.packages('viridis')
library(viridis)
install.packages('gridExtra')
library(gridExtra)
install.packages('fmsb')
library(fmsb)
install.packages('janitor')
library(janitor)


#********** LEITURA **********
halo <- read.table(file='halophiles.tsv', sep='\t', dec=',', quote='', header=TRUE, fill=TRUE) #Lê tabela

# ********** PAPER POR ANO **********
paper_ano <- count(halo, Domain, Year)
tab_paper_ano <- cbind.data.frame(paper_ano$Year, paper_ano$Domain, paper_ano$n)
colnames(tab_paper_ano) <- c('Year', 'Domain', 'n')
graf_paper_ano <- ggplot(tab_paper_ano, aes(x=Year, y=n, fill=Domain))+
    geom_bar(position='stack', stat='identity')+
    theme_bw()+
    labs(y='Number of species')+
    theme(
    plot.title = element_text(family="", face="bold", size=16, hjust=.5), #hjust coloca título para o centro
    plot.subtitle = element_text(face="bold", size=12, hjust=0.5, colour="dark gray"),
    axis.title = element_text(size=12),
    axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0)),
    axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0)),
    panel.grid.minor=element_blank(),
    legend.position='top')+
    scale_fill_manual(values=c('#49C1ADFF', '#357BA2FF', '#3E356BFF'), labels=c('Archaea', 'Bacteria', 'Eukarya'))+
    scale_x_continuous(breaks = seq(min(tab_paper_ano$Year), max(tab_paper_ano$Year), by = 1))+
    coord_flip()
graf_paper_ano
ggsave('Species x Year.png', graf_paper_ano, device='png', unit='mm', width=119, height=200, dpi=600)
ggsave('Species x Year.tiff', graf_paper_ano, device='tiff', unit='mm', width=119, height=200, dpi=600)


# ********** QUANTIDADE POR CLADO (estilo artigo Chen et al 2025) **********
fam <- count(halo, Family)
x <- rep(c('1'), times=length(fam$Family))
fam <- cbind.data.frame(fam, x)
graf_fam <- ggplot(fam, aes(x=x, y=n, fill=Family))+
    geom_bar(position='stack', stat='identity', color='black')+
    theme_void()+
    labs(y='(%)', x='', fill='' )+
    theme(
    plot.title = element_text(family="", face="bold", size=16, hjust=.5), #hjust coloca título para o centro
    plot.subtitle = element_text(face="bold", size=12, hjust=0.5, colour="dark gray"),
    axis.title = element_text(size=12),
    axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0)),
    axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0)),
    panel.grid.minor=element_blank(),
    legend.position="none")+
    coord_flip()
graf_fam

cla <- count(halo, Class)
x <- rep(c('1'), times=length(cla$Class))
cla <- cbind.data.frame(cla, x)
graf_cla <- ggplot(cla, aes(x=x, y=n, fill=Class))+
    geom_bar(position='stack', stat='identity', color='black')+
    theme_void()+
    labs(y='(%)', x='', fill='' )+
    theme(
    plot.title = element_text(family="", face="bold", size=16, hjust=.5), #hjust coloca título para o centro
    plot.subtitle = element_text(face="bold", size=12, hjust=0.5, colour="dark gray"),
    axis.title = element_text(size=12),
    axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0)),
    axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0)),
    panel.grid.minor=element_blank(),
    legend.position="none")+
        coord_flip()
graf_cla

# ********** PRÉ-PROCESSAMENTO ***********
# Transformar em número as colunas de concentrações
halo$Min.opt <- as.numeric(sub(",", ".", halo$Min.opt, fixed=TRUE))
halo$Max.opt <- as.numeric(sub(",", ".", halo$Max.opt, fixed=TRUE))
halo$Min.sup <- as.numeric(sub(",", ".", halo$Min.sup, fixed=TRUE))
halo$Max.sup <- as.numeric(sub(",", ".", halo$Max.sup, fixed=TRUE))

bad_halo <- c()   # Cria vetor com linhas que sal não é NaCl
for(n in seq_along(halo$Min.opt)){    # Passa por cada espécie
    if(halo$NaCl.[n]!='s'){    # Testa de sal é NaCl
        bad_halo <- append(bad_halo, n, after=length(bad_halo))}    # Se não é NaCl, armazena linha num vetor

    else if(halo$Unit[n] == 'p'){    # p=porcentual m/v | loop faz conversão para mol/L
        halo$Min.opt[n] <- ((halo$Min.opt[n]*10)/58.44)
        halo$Max.opt[n] <- ((halo$Max.opt[n]*10)/58.44)
        halo$Min.sup[n] <- ((halo$Min.sup[n]*10)/58.44)
        halo$Max.sup[n] <- ((halo$Max.sup[n]*10)/58.44)}

    else if(halo$Unit[n] == 'g'){   # g=g/L | loop faz conversão para mol/L
      halo$Min.opt[n] <- (halo$Min.opt[n]/58.44)
      halo$Max.opt[n] <- (halo$Max.opt[n]/58.44)
      halo$Min.sup[n] <- (halo$Min.sup[n]/58.44)
      halo$Max.sup[n] <- (halo$Max.sup[n]/58.44)}
}

halotab <- halo[-bad_halo,]

sup <- halotab[,-c(3, 4, 10)] #Retira colunas de concentração ótima e outras informações desnecessárias
opt <- halotab[,-c(5, 6, 10)] #Retira colunas de concentração suportada e outras informações desnecessárias

halo_sup <- na.omit(sup) #Tira linhas sem informações de concentração
halo_opt <- na.omit(opt)

write.table(halo_opt, file="halo_opt.tsv") #Salva arquivo .tsv da tabela limpa de concentrações ótimas
table(halo_opt$Domain)

write.table(halo_sup, file="halo_sup.tsv") #Salva arquivo .tsv da tabela limpa de concentrações ótimas
table(halo_sup$Domain)

#********** CONTAGEM CONCENTRAÇÃO ÓTIMA DE CRESCIMENTO **********

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


quartil <- quantile(graf_opt, probs = c(0,0.25,0.5,0.75,1))


# *********** CLASSIFICAÇÕES **********
# *--------* POR MÉDIA *--------*
halotab <- na.omit(halotab)
mean.opt <- (1/2*halotab$Min.opt)+(1/2*halotab$Max.opt)
#----------- Classificação 1 (DasSarma and Arora 2001) ----------
class1 <- c()

for(n in seq_along(halotab$Species)){
    if(mean.opt[n]<0.2 & halotab$Max.sup[n]>=0.2){
        class1 <- append(class1, 'Tolerant', after=length(class1))}
    else if(0.2<=mean.opt[n] & mean.opt[n]<0.85){
        class1 <- append(class1, 'Slight', after=length(class1))}
    else if(0.85<=mean.opt[n] & mean.opt[n]<3.4){
        class1 <- append(class1, 'Moderate', after=length(class1))}
    else if(3.4<=mean.opt[n]){
        class1 <- append(class1, 'Extreme', after=length(class1))}
    else{
        class1 <- append(class1, 'Not', after=length(class1))}
}

#----------- Classificação 2 (Kushner and Kamemura 1988) ----------
class2 <- c()

for(n in seq_along(halotab$Species)){
    if(mean.opt[n]<0.2 & halotab$Max.sup[n]>=0.2){
        class2 <- append(class2, 'Tolerant', after=length(class1))}
    else if(0.2<=mean.opt[n] & mean.opt[n]<0.5){
        class2 <- append(class2, 'Slight', after=length(class1))}
    else if(0.5<=mean.opt[n] & mean.opt[n]<2.5){
        class2 <- append(class2, 'Moderate', after=length(class1))}
    else if(2.5<=mean.opt[n]){
        class2 <- append(class2, 'Extreme', after=length(class1))}
    else{
        class2 <- append(class2, 'Not', after=length(class1))}
}

#----------- Classificação 3 (The Prokaryotes 2006) ----------
class3 <- c()

for(n in seq_along(halotab$Species)){
    if(mean.opt[n]<20/58.44 & halotab$Max.sup[n]>=20/58.44){
        class3 <- append(class3, 'Tolerant', after=length(class1))}
    else if(20/58.44<=mean.opt[n] & mean.opt[n]<50/58.44){
        class3 <- append(class3, 'Slight', after=length(class1))}
    else if(50/58.44<=mean.opt[n] & mean.opt[n]<200/58.44){
        class3 <- append(class3, 'Moderate', after=length(class1))}
    else if(200/58.44<=mean.opt[n]){
        class3 <- append(class3, 'Extreme', after=length(class1))}
    else{
        class3 <- append(class3, 'Not', after=length(class1))}
}

#----------- Classificação 4 (Proposta) ----------
class4 <- c()

for(n in seq_along(halotab$Species)){
    if(mean.opt[n]<0.6 & halotab$Max.sup[n]>=0.6){
        class4 <- append(class4, 'Tolerant', after=length(class1))}
    else if(0.6<=mean.opt[n] & mean.opt[n]<1.3){
        class4 <- append(class4, 'Slight', after=length(class1))}
    else if(1.3<=mean.opt[n] & mean.opt[n]<2.3){
        class4 <- append(class4, 'Moderate', after=length(class1))}
    else if(2.3<=mean.opt[n]){
        class4 <- append(class4, 'Extreme', after=length(class1))}
    else{
        class4 <- append(class4, 'Not', after=length(class1))}
}

# ---------- Comparação de classificação por média -----------
clastab <- cbind.data.frame(halotab$Classif, class1, class2, class3, class4)
colnames(clastab) <- c('author', 'clas1', 'clas2', 'clas3', 'clas4')
clastab_clas <- filter(clastab, author!='No Class')

clas.compare1 <- clastab_clas$author==clastab_clas$clas1
clas.compare2 <- clastab_clas$author==clastab_clas$clas2
clas.compare3 <- clastab_clas$author==clastab_clas$clas3
clas.compare4 <- clastab_clas$author==clastab_clas$clas4

table(clas.compare1)
table(clas.compare2)
table(clas.compare3)
table(clas.compare4)

clas1.t <- (sum(clas.compare1, na.rm=TRUE)/length(clas.compare1))*100
clas1.f <- 100-clas1.t
clas2.t <- (sum(clas.compare2, na.rm=TRUE)/length(clas.compare2))*100
clas2.f <- 100-clas2.t
clas3.t <- (sum(clas.compare3, na.rm=TRUE)/length(clas.compare3))*100
clas3.f <- 100-clas3.t
clas4.t <- (sum(clas.compare4, na.rm=TRUE)/length(clas.compare4))*100
clas4.f <- 100-clas4.t

clas.compare <- c(clas1.t, clas2.t, clas3.t, clas4.t, clas1.f, clas2.f, clas3.f, clas4.f)
status <- rep(c('Equal', 'Different'), each=4)
classification <- rep(c('1', '2', '3','4'), times=2)
compare.tab <- cbind.data.frame(classification, status, clas.compare)

graf_compare <- ggplot(compare.tab, aes(x=classification, y=clas.compare, fill=status))+
    geom_bar(position='stack', stat='identity')+
    theme_bw()+
    labs(y='(%)', x='Classification System', fill='' )+
    theme(
    plot.title = element_text(family="", face="bold", size=16, hjust=.5), #hjust coloca título para o centro
    plot.subtitle = element_text(face="bold", size=12, hjust=0.5, colour="dark gray"),
    axis.title = element_text(size=12),
    axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0)),
    axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0)),
    panel.grid.minor=element_blank())
graf_compare
ggsave('Compare_mean.png', graf_compare, device='png', unit='cm', width=15, height=12, dpi=600)

# *--------* POR CONTAGEM *--------*
#----------- Classificação 1 (DasSarma and Arora 2001) ----------
Tolerant <- c()
Slight <- c()
Moderate <- c()
Extreme <- c()
Not <- c()

for(n in seq_along(halotab[,3])){
  opt_range <- seq(halotab$Min.opt[n], halotab$Max.opt[n], 0.1)
  sup_max <- halotab$Max.sup[n]
  t <- 0
  s <- 0
  m <- 0
  e <- 0
  nt <- 0
  for(i in seq_along(opt_range)){
    if(0.2<=opt_range[i] & opt_range[i]<0.85){
      s <- s+1}
    else if(0.85<=opt_range[i] & opt_range[i]<3.4){
      m <- m+1}
    else if(3.4<=opt_range[i]){
      e <- e+1}
    else if(opt_range[i]<0.2){
      if(sup_max>=0.2){
        t <- t+1}
      else{nt <- nt+1}
    }
  }
  Tolerant <- append(Tolerant, t, after=length(Tolerant))
  Slight <- append(Slight, s, after=length(Slight))
  Moderate <- append(Moderate, m, after=length(Moderate))
  Extreme <- append(Extreme, e, after=length(Extreme))
  Not <- append(Not, nt , after=length(Not))
}

clastab1 <- cbind.data.frame(Tolerant, Slight, Moderate, Extreme, Not)
max_col <- max.col(clastab1)
clas1_count <- colnames(clastab1)[max_col]
halo1 <- cbind.data.frame(halotab, clas1_count)

#----------- Classificação 2 (Kushner and Kamemura 1988) ----------
Tolerant <- c()
Slight <- c()
Moderate <- c()
Extreme <- c()
Not <- c()

for(n in seq_along(halotab[,3])){
  opt_range <- seq(halotab$Min.opt[n], halotab$Max.opt[n], 0.1)
  sup_max <- halotab$Max.sup[n]
  t <- 0
  s <- 0
  m <- 0
  e <- 0
  nt <- 0
  for(i in seq_along(opt_range)){
    if(0.2<=opt_range[i] & opt_range[i]<0.5){
      s <- s+1}
    else if(0.5<=opt_range[i] & opt_range[i]<2.5){
      m <- m+1}
    else if(2.5<=opt_range[i]){
      e <- e+1}
    else if(opt_range[i]<0.2){
      if(sup_max>=0.2){
        t <- t+1}
      else{nt <- nt+1}
    }
  }
  Tolerant <- append(Tolerant, t, after=length(Tolerant))
  Slight <- append(Slight, s, after=length(Slight))
  Moderate <- append(Moderate, m, after=length(Moderate))
  Extreme <- append(Extreme, e, after=length(Extreme))
  Not <- append(Not, nt , after=length(Not))
}

clastab2 <- cbind.data.frame(Tolerant, Slight, Moderate, Extreme, Not)
max_col <- max.col(clastab2)
clas2_count <- colnames(clastab2)[max_col]
halo2 <- cbind.data.frame(halotab, clas2_count)

#----------- Classificação 3 (The Prokaryotes 2006) ----------
Tolerant <- c()
Slight <- c()
Moderate <- c()
Extreme <- c()
Not <- c()

for(n in seq_along(halotab[,3])){
  opt_range <- seq(halotab$Min.opt[n], halotab$Max.opt[n], 0.1)
  sup_max <- halotab$Max.sup[n]
  t <- 0
  s <- 0
  m <- 0
  e <- 0
  nt <- 0
  for(i in seq_along(opt_range)){
    if(0.34<=opt_range[i] & opt_range[i]<0.86){
      s <- s+1}
    else if(0.86<=opt_range[i] & opt_range[i]<3.42){
      m <- m+1}
    else if(3.42<=opt_range[i]){
      e <- e+1}
    else if(opt_range[i]<0.34){
      if(sup_max>=0.34){
        t <- t+1}
      else{nt <- nt+1}
    }
  }
  Tolerant <- append(Tolerant, t, after=length(Tolerant))
  Slight <- append(Slight, s, after=length(Slight))
  Moderate <- append(Moderate, m, after=length(Moderate))
  Extreme <- append(Extreme, e, after=length(Extreme))
  Not <- append(Not, nt , after=length(Not))
}

clastab3 <- cbind.data.frame(Tolerant, Slight, Moderate, Extreme, Not)
max_col <- max.col(clastab3)
clas3_count <- colnames(clastab3)[max_col]
halo3 <- cbind.data.frame(halotab, clas3_count)

#----------- Classificação 4 (Proposta) ----------
Tolerant <- c()
Slight <- c()
Moderate <- c()
Extreme <- c()
Not <- c()

for(n in seq_along(halotab[,3])){
  opt_range <- seq(halotab$Min.opt[n], halotab$Max.opt[n], 0.1)
  sup_max <- halotab$Max.sup[n]
  t <- 0
  s <- 0
  m <- 0
  e <- 0
  nt <- 0
  for(i in seq_along(opt_range)){
    if(0.6<=opt_range[i] & opt_range[i]<1.3){
      s <- s+1}
    else if(1.3<=opt_range[i] & opt_range[i]<2.3){
      m <- m+1}
    else if(2.3<=opt_range[i]){
      e <- e+1}
    else if(opt_range[i]<0.6){
      if(sup_max>=0.6){
        t <- t+1}
      else{nt <- nt+1}
    }
  }
  Tolerant <- append(Tolerant, t, after=length(Tolerant))
  Slight <- append(Slight, s, after=length(Slight))
  Moderate <- append(Moderate, m, after=length(Moderate))
  Extreme <- append(Extreme, e, after=length(Extreme))
  Not <- append(Not, nt , after=length(Not))
}

clastab4 <- cbind.data.frame(Tolerant, Slight, Moderate, Extreme, Not)
max_col <- max.col(clastab4)
clas4_count <- colnames(clastab4)[max_col]
halo4 <- cbind.data.frame(halotab, clas4_count)

# ---------- Comparação de classificação por contagem -----------
clastab <- cbind.data.frame(halotab$Classif, clas1_count, clas2_count, clas3_count, clas4_count)
colnames(clastab) <- c('author', 'clas1', 'clas2', 'clas3', 'clas4')
clastab_clas <- filter(clastab, author!='No Class')

clas.compare1 <- clastab_clas$author==clastab_clas$clas1
clas.compare2 <- clastab_clas$author==clastab_clas$clas2
clas.compare3 <- clastab_clas$author==clastab_clas$clas3
clas.compare4 <- clastab_clas$author==clastab_clas$clas4

table(clas.compare1)
table(clas.compare2)
table(clas.compare3)
table(clas.compare4)

clas1.t <- (sum(clas.compare1, na.rm=TRUE)/length(clas.compare1))*100
clas1.f <- 100-clas1.t
clas2.t <- (sum(clas.compare2, na.rm=TRUE)/length(clas.compare2))*100
clas2.f <- 100-clas2.t
clas3.t <- (sum(clas.compare3, na.rm=TRUE)/length(clas.compare3))*100
clas3.f <- 100-clas3.t
clas4.t <- (sum(clas.compare4, na.rm=TRUE)/length(clas.compare4))*100
clas4.f <- 100-clas4.t

clas.compare <- c(clas1.t, clas2.t, clas3.t, clas4.t, clas1.f, clas2.f, clas3.f, clas4.f)
status <- rep(c('Equal', 'Different'), each=4)
classification <- rep(c('1', '2', '3', '4'), times=2)
compare.tab <- cbind.data.frame(classification, status, clas.compare)

graf_compare_count <- ggplot(compare.tab, aes(x=classification, y=clas.compare, fill=status))+
    geom_bar(position='stack', stat='identity')+
    theme_bw()+
    labs(y='(%)', x='Classification System', fill='' )+
    theme(
    plot.title = element_text(family="", face="bold", size=16, hjust=.5), #hjust coloca título para o centro
    plot.subtitle = element_text(face="bold", size=12, hjust=0.5, colour="dark gray"),
    axis.title = element_text(size=12),
    axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0)),
    axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0)),
    panel.grid.minor=element_blank())+
    scale_fill_manual(values=c('#3497A9FF', '#382A54FF'))
graf_compare_count
ggsave('Compare_count.png', graf_compare_count, device='png', unit='mm', width=100, height=90, dpi=600)
ggsave('Compare_count.tiff', graf_compare_count, device='tiff', unit='mm', width=100, height=90, dpi=600)

# RADAR CHART
table.author <- tabyl(clastab$author)
colnames(table.author) <- c('class', 'n', 'percent')
tableauthor <- filter(table.author, class!='No Class')
table1 <- tabyl(clas1_count)
table1 <- filter(table1, clas1_count!='Not')
table2 <- tabyl(clas2_count)
table2 <- filter(table2, clas2_count!='Not')
table3 <- tabyl(clas3_count)
table3 <- filter(table3, clas3_count!='Not')
table4 <- tabyl(clas4_count)
table4 <- filter(table4, clas4_count!='Not')

# Tudo junto
tabclas.count <- rbind.data.frame(tableauthor$n, table1$n, table2$n, table3$n, table4$n)
rownames(tabclas.count) <- c("Author's", 'DasSarma and Arora 2001', 'Kushner and Kamekura 1988', 'Larsen 1986', 'Proposed')
colnames(tabclas.count) <- c('Extreme', 'Moderate', 'Slight', 'Tolerant')


tabclas.count <- rbind.data.frame(tableauthor$n, table1$n, table2$n, table3$n, table4$n)
colnames(tabclas.count) <- c('Extreme', 'Moderate', 'Slight', 'Tolerant')
#max_min <- data.frame(Extreme=c(max(tabclas_count$Extreme), 0), Moderate=c(max(tabclas_count$Moderate), 0), Slight=c(max(tabclas_count$Slight), 0), Tolerant=c(max(tabclas_count$Tolerant), 0))
max_min <- data.frame(Extreme=c(600, 0), Moderate=c(600, 0), Slight=c(600, 0), Tolerant=c(600, 0))
tabclascount <- rbind.data.frame(max_min, tabclas.count)
rownames(tabclascount) <- c('Max', 'Min', "Author's description", 'DasSarma and Arora 2001', 'Kushner and Kamekura 1988', 'Larsen 1986', 'Proposed')

borda <- viridis_pal(option = "mako")(5)
dentro <- alpha(borda, 0.15)
radarchart(tabclascount, axistype=0, maxmin=T,
    pcol=borda, pfcol=dentro, plwd=4, plty=1,
    cglcol="grey", cglty=1, axislabcol="black", cglwd=0.8,
    vlcex=0.8)
legend(x=0.7, y=1, legend = rownames(tabclascount[-c(1,2),]), bty = "n", pch=20 , col=borda , text.col = "grey", cex=1.2, pt.cex=3)

# Proposed vs. authors
tabclas_count <- rbind.data.frame(tableauthor$n, table4$n)
colnames(tabclas_count) <- c('Extreme', 'Moderate', 'Slight', 'Tolerant')
#max_min <- data.frame(Extreme=c(max(tabclas_count$Extreme), 0), Moderate=c(max(tabclas_count$Moderate), 0), Slight=c(max(tabclas_count$Slight), 0), Tolerant=c(max(tabclas_count$Tolerant), 0))
max_min <- data.frame(Extreme=c(500, 0), Moderate=c(500, 0), Slight=c(500, 0), Tolerant=c(500, 0))
tabclascount <- rbind.data.frame(max_min,tabclas_count)
rownames(tabclascount) <- c('Max', 'Min', "Author's description", 'Proposed')

#viridis_pal(option = "mako")(5)
borda <- c('#49C1ADFF', '#3E356BFF')
dentro <- alpha(borda, 0.3)
radarchart(tabclascount, axistype=0, maxmin=T,
    pcol=borda, pfcol=dentro, plwd=4, plty=1,
    cglcol="grey", cglty=1, axislabcol="black", cglwd=0.8,
    vlcex=0.8)
legend(x=0.7, y=1, legend = rownames(tabclascount[-c(1,2),]), bty = "n", pch=20 , col=borda , text.col = "darkgrey", cex=1, pt.cex=3)


# Proposta vs. DasSarma and Arora 2001
tabclas_1 <- rbind.data.frame(table1$n, table4$n)
colnames(tabclas_1) <- c('Extreme', 'Moderate', 'Slight', 'Tolerant')
#max_min <- data.frame(Extreme=c(max(tabclas_count$Extreme), 0), Moderate=c(max(tabclas_count$Moderate), 0), Slight=c(max(tabclas_count$Slight), 0), Tolerant=c(max(tabclas_count$Tolerant), 0))
max_min <- data.frame(Extreme=c(500, 0), Moderate=c(500, 0), Slight=c(500, 0), Tolerant=c(500, 0))
tabclas1 <- rbind.data.frame(max_min,tabclas_1)
rownames(tabclas1) <- c('Max', 'Min', "DasSarma and Arora 2001", 'Proposed')

#viridis_pal(option = "mako")(5)
borda <- c('#49C1ADFF', '#3E356BFF')
dentro <- alpha(borda, 0.3)
radarchart(tabclas1, axistype=0, maxmin=T,
    pcol=borda, pfcol=dentro, plwd=4, plty=1,
    cglcol="grey", cglty=1, axislabcol="black", cglwd=0.8,
    vlcex=0.8)
legend(x=0.7, y=1, legend = rownames(tabclas1[-c(1,2),]), bty = "n", pch=20 , col=borda , text.col = "darkgrey", cex=1, pt.cex=3)

# Proposta vs. Kushner and Kamekura 1988
tabclas_2 <- rbind.data.frame(table2$n, table4$n)
colnames(tabclas_2) <- c('Extreme', 'Moderate', 'Slight', 'Tolerant')
#max_min <- data.frame(Extreme=c(max(tabclas_count$Extreme), 0), Moderate=c(max(tabclas_count$Moderate), 0), Slight=c(max(tabclas_count$Slight), 0), Tolerant=c(max(tabclas_count$Tolerant), 0))
max_min <- data.frame(Extreme=c(500, 0), Moderate=c(600, 0), Slight=c(500, 0), Tolerant=c(500, 0))
tabclas2 <- rbind.data.frame(max_min,tabclas_2)
rownames(tabclas2) <- c('Max', 'Min', "DasSarma and Arora 2001", 'Proposed')

#viridis_pal(option = "mako")(5)
borda <- c('#49C1ADFF', '#3E356BFF')
dentro <- alpha(borda, 0.3)
radarchart(tabclas2, axistype=0, maxmin=T,
    pcol=borda, pfcol=dentro, plwd=4, plty=1,
    cglcol="grey", cglty=1, axislabcol="black", cglwd=0.8,
    vlcex=0.8)
legend(x=0.7, y=1, legend = rownames(tabclas2[-c(1,2),]), bty = "n", pch=20 , col=borda , text.col = "darkgrey", cex=1, pt.cex=3)

# Proposta vs. Larsen 1986
tabclas_3 <- rbind.data.frame(table3$n, table4$n)
colnames(tabclas_3) <- c('Extreme', 'Moderate', 'Slight', 'Tolerant')
#max_min <- data.frame(Extreme=c(max(tabclas_count$Extreme), 0), Moderate=c(max(tabclas_count$Moderate), 0), Slight=c(max(tabclas_count$Slight), 0), Tolerant=c(max(tabclas_count$Tolerant), 0))
max_min <- data.frame(Extreme=c(500, 0), Moderate=c(500, 0), Slight=c(500, 0), Tolerant=c(500, 0))
tabclas3 <- rbind.data.frame(max_min,tabclas_3)
rownames(tabclas3) <- c('Max', 'Min', "DasSarma and Arora 2001", 'Proposed')

#viridis_pal(option = "mako")(5)
borda <- c('#49C1ADFF', '#3E356BFF')
dentro <- alpha(borda, 0.3)
radarchart(tabclas3, axistype=0, maxmin=T,
    pcol=borda, pfcol=dentro, plwd=4, plty=1,
    cglcol="grey", cglty=1, axislabcol="black", cglwd=0.8,
    vlcex=0.8)
legend(x=0.7, y=1, legend = rownames(tabclas3[-c(1,2),]), bty = "n", pch=20 , col=borda , text.col = "darkgrey", cex=1, pt.cex=3)


#Individuais
max_min <- data.frame(Extreme=c(700, 0), Moderate=c(700, 0), Slight=c(700, 0), Tolerant=c(700, 0))
tab.prop <- rbind.data.frame(max_min, table4$n)
colnames(tab.prop) <- c('Extreme', 'Moderate', 'Slight', 'Tolerant')
rownames(tab.prop) <- c('Max', 'Min', 'Proposed')
borda <- c('#0B0405FF')
dentro <- alpha(borda, 0.3)
png(filename = "radar_proposta.png", width = 13, height = 10, units = "cm", res = 600)
op <- par(mar = c(1, 0.5, 1, 1))
radarchart(tab.prop, axistype=0, maxmin=T,
    pcol=borda, pfcol=dentro, plwd=4, plty=1,
    cglcol="grey", cglty=1, axislabcol="black", cglwd=0.8,
    vlcex=0.7, title= 'Schiavo et al. 2025')
par(op)
dev.off()

tab1 <- rbind.data.frame(max_min, table1$n)
colnames(tab1) <- c('Extreme', 'Moderate', 'Slight', 'Tolerant')
rownames(tab1) <- c('Max', 'Min', 'DasSarma and Arora 2001')
borda <- c('#382A54FF')
dentro <- alpha(borda, 0.3)
png(filename = "radar_dassarma.png", width = 13, height = 10, units = "cm", res = 600)
op <- par(mar = c(1, 0.5, 1, 1))
radarchart(tab1, axistype=0, maxmin=T,
    pcol=borda, pfcol=dentro, plwd=4, plty=1,
    cglcol="grey", cglty=1, axislabcol="black", cglwd=0.8,
    vlcex=0.7, title= 'DasSarma and DasSarma 2012')
par(op)
dev.off()

tab2 <- rbind.data.frame(max_min, table2$n)
colnames(tab2) <- c('Extreme', 'Moderate', 'Slight', 'Tolerant')
rownames(tab2) <- c('Max', 'Min', 'Kushner and Kamekura 1988')
borda <- c('#395D9CFF')
dentro <- alpha(borda, 0.3)
png(filename = "radar_Kushner.png", width = 13, height = 10, units = "cm", res = 600)
op <- par(mar = c(1, 0.5, 1, 1))
radarchart(tab2, axistype=0, maxmin=T,
    pcol=borda, pfcol=dentro, plwd=4, plty=1,
    cglcol="grey", cglty=1, axislabcol="black", cglwd=0.8,
    vlcex=0.7, title= 'Kushner and Kamekura 1988')
par(op)
dev.off()

tab3 <- rbind.data.frame(max_min, table3$n)
colnames(tab3) <- c('Extreme', 'Moderate', 'Slight', 'Tolerant')
rownames(tab3) <- c('Max', 'Min', 'Larsen 1986')
borda <- c('#3497A9FF')
dentro <- alpha(borda, 0.3)
png(filename = "radar_larsen.png", width = 13, height = 10, units = "cm", res = 600)
op <- par(mar = c(1, 0.5, 1, 1))
rc3 <- radarchart(tab3, axistype=0, maxmin=T,
    pcol=borda, pfcol=dentro, plwd=4, plty=1,
    cglcol="grey", cglty=1, axislabcol="black", cglwd=0.8,
    vlcex=0.7, title= 'Larsen 1986')
par(op)
dev.off()

tab.aut <- rbind.data.frame(max_min, tableauthor$n)
colnames(tab.aut) <- c('Extreme', 'Moderate', 'Slight', 'Tolerant')
rownames(tab.aut) <- c('Max', 'Min', "Original description")
borda <- c('#60CEACFF')
dentro <- alpha(borda, 0.3)
png(filename = "radar_original.png", width = 13, height = 10, units = "cm", res = 600)
op <- par(mar = c(1, 0.5, 1, 1))
radarchart(tab.aut, axistype=0, maxmin=T,
    pcol=borda, pfcol=dentro, plwd=4, plty=1,
    cglcol="grey", cglty=1, axislabcol="black", cglwd=0.8,
    vlcex=0.7, title= 'Original description')
par(op)
dev.off()

# ********** HEATMAP CLADO POR CLASSIFICAÇÃO DADA PELO AUTOR **********
halotab <- na.omit(halotab)
hm_author <- count(halotab, Classif, Class)
hm_author$Classif <- factor(hm_author$Classif, levels=c('Tolerant', 'Slight', 'Moderate', 'Extreme', 'No Class'))
graf_hm_author_c <-ggplot(hm_author, aes(x=Classif, y=Class, fill=n))+
    geom_tile()+
    theme_bw()+
    labs(x='', title='a')+
    theme(
      panel.grid=element_blank())+
    scale_fill_viridis(option='mako', direction=-1)
graf_hm_author_c
ggsave('heatmap_author_classe.png', graf_hm_author_c, device='png', unit='cm', width=15, height=20, dpi=600)

hm_author <- count(halotab, Classif, Phylum)
hm_author$Classif <- factor(hm_author$Classif, levels=c('No Class', 'Tolerant', 'Slight', 'Moderate', 'Extreme'))
graf_hm_author_f <-ggplot(hm_author, aes(x=Classif, y=Phylum, fill=n))+
    geom_tile()+
    theme_bw()+
    labs(x='')+
    theme(
      panel.grid=element_blank())+
    scale_fill_viridis(option='mako', direction=-1)
graf_hm_author_f
ggsave('heatmap_author_filo.png', graf_hm_author_f, device='png', unit='cm', width=15, height=20, dpi=600)

hm_author <- count(halotab, Classif, Kingdom)
hm_author$Classif <- factor(hm_author$Classif, levels=c('No Class', 'Tolerant', 'Slight', 'Moderate', 'Extreme'))
graf_hm_author_r <-ggplot(hm_author, aes(x=Classif, y=Kingdom, fill=n))+
    geom_tile()+
    theme_bw()+
    labs(x='')+
    theme(
      panel.grid=element_blank())+
    scale_fill_viridis(option='mako', direction=-1)
graf_hm_author_r
ggsave('heatmap_author_reino.png', graf_hm_author_r, device='png', unit='cm', width=15, height=20, dpi=600)

hm_author <- count(halotab, Classif, Order)
hm_author$Classif <- factor(hm_author$Classif, levels=c('No Class','Tolerant', 'Slight', 'Moderate', 'Extreme'))
graf_hm_author_o <-ggplot(hm_author, aes(x=Classif, y=Order, fill=n))+
    geom_tile()+
    theme_bw()+
    labs(x='')+
    theme(
      panel.grid=element_blank())+
    scale_fill_viridis(option='mako', direction=-1)
graf_hm_author_o
ggsave('heatmap_author_ordem.png', graf_hm_author_o, device='png', unit='cm', width=15, height=22, dpi=600)


# ********** HEATMAP CLADO POR CLASSIFICAÇÃO PROPOSTA **********
haloclas <- cbind.data.frame(halotab, clas4_count)
hm_prop <- count(haloclas, clas4_count, Class)
hm_prop$clas4_count <- factor(hm_prop$clas4_count, levels=c('Tolerant', 'Slight', 'Moderate', 'Extreme', 'Not'))
graf_hm_prop_c <-ggplot(hm_prop, aes(x=clas4_count, y=Class, fill=n))+
    geom_tile()+
    theme_bw()+
    labs(x='', title='b')+
    theme(
      panel.grid=element_blank())+
    scale_fill_viridis(option='mako', direction=-1)
graf_hm_prop_c
ggsave('heatmap_prop_classe.png', graf_hm_prop_c, device='png', unit='cm', width=15, height=20, dpi=600)

hm_prop <- count(haloclas, clas4_count, Phylum)
hm_prop$clas4_count <- factor(hm_prop$clas4_count, levels=c('Not', 'Tolerant', 'Slight', 'Moderate', 'Extreme'))
graf_hm_prop_f <-ggplot(hm_prop, aes(x=clas4_count, y=Phylum, fill=n))+
    geom_tile()+
    theme_bw()+
    labs(x='')+
    theme(
      panel.grid=element_blank())+
    scale_fill_viridis(option='mako', direction=-1)
graf_hm_prop_f
ggsave('heatmap_prop_filo.png', graf_hm_prop_f, device='png', unit='cm', width=15, height=20, dpi=600)

hm_prop <- count(haloclas, clas4_count, Kingdom)
hm_prop$clas4_count <- factor(hm_prop$clas4_count, levels=c('Not', 'Tolerant', 'Slight', 'Moderate', 'Extreme'))
graf_hm_prop_r <-ggplot(hm_prop, aes(x=clas4_count, y=Kingdom, fill=n))+
    geom_tile()+
    theme_bw()+
    labs(x='')+
    theme(
      panel.grid=element_blank())+
    scale_fill_viridis(option='mako', direction=-1)
graf_hm_prop_r
ggsave('heatmap_prop_reino.png', graf_hm_prop_r, device='png', unit='cm', width=15, height=20, dpi=600)

hm_prop <- count(haloclas, clas4_count, Order)
hm_prop$clas4_count <- factor(hm_prop$clas4_count, levels=c('Not', 'Tolerant', 'Slight', 'Moderate', 'Extreme'))
graf_hm_prop_o <-ggplot(hm_prop, aes(x=clas4_count, y=Order, fill=n))+
    geom_tile()+
    theme_bw()+
    labs(x='')+
    theme(
      panel.grid=element_blank())+
    scale_fill_viridis(option='mako', direction=-1)
graf_hm_prop_o
ggsave('heatmap_prop_ordem.png', graf_hm_prop_o, device='png', unit='cm', width=15, height=22, dpi=600)

# ---------- Criação de grid -----------
grid_hm_class <- grid.arrange(graf_hm_author_c, graf_hm_prop_c, ncol=1)
ggsave('grid_hm_class.tiff', grid_hm_class,  device='tiff', unit='mm', width = 150, height = 300, dpi=600)
ggsave('grid_hm_class.png', grid_hm_class,  device='png', unit='mm', width = 133, height = 300, dpi=600)







manynot <- count(haloclas, clas4_count)
many.author <- count(halotab, Classif)
sp_domain <- count(halotab, Domain)
