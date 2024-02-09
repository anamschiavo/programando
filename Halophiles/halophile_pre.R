# Scrip de pré-processamento da tabela de halófilos da literatura

# REGRAS DE ENTRADA:
# ***********ORDEM DAS COLUNAS: Espécie|Domínio|Opt min|Opt max|Sup sup|Sup max|Classificação|É NaCl?|Unidade de concentração***************
# Concentrações não podem ter símbolo %, se for dado assim no paper, coloque "p" na coluna 7 (Unidade de concentração)
# Decimal é dado por vírgula
# Coluna 7 (Unidade de concentração) pode ser p = % massa/volume | g = g/L | m = mol/L
# Coluna 6 (É NaCl?) pode ser "s" para sim. Qualquer outra entrada diz que é outro sal ou combinação: ISSO FARÁ COM QUE SEJA RETIRADO!
# Coluna 2 (Dominío) pode ser B = bactérias | A = archaea | E = eukarya

#********** PACOTES E BIBLOTECAS **********
install.packages("dplyr")   #Linha necessária apenas se não tiver pacote já instalado
library(dplyr)



#********** LEITURA E PRIMEIRO PRÉ-PROCESSAMENTO **********
halo <- read.table(file='halophile_bruto.tsv', sep='\t', dec=',', quote='', header=TRUE, fill=TRUE) #Lê tabela

halotab <- filter(halo, journal=='Micro Society')


sup <- halotab[,-c(3, 4, 10, 11, 12, 13)] #Retira colunas de concentração ótima e outras informações desnecessárias
opt <- halotab[,-c(5, 6, 10, 11, 12 ,13)] #Retira colunas de concentração suportada e outras informações desnecessárias



sup[,3] <- as.numeric(sub(",", ".", sup[,3], fixed=TRUE))   #Transforma mínimo suportado em números, considerando "," como decimal e mantendo dígitos decimais
sup[,4] <- as.numeric(sub(",", ".", sup[,4], fixed=TRUE))   #Transforma máximo suportado em números, considerando "," como decimal e mantendo dígitos decimais

opt[,3] <- as.numeric(sub(",", ".", opt[,3], fixed=TRUE))   #Transforma mínimo ótimo em números, considerando "," como decimal e mantendo dígitos decimais
opt[,4] <- as.numeric(sub(",", ".", opt[,4], fixed=TRUE))   #Transforma máximo ótimo em números, considerando "," como decimal e mantendo dígitos decimais

halo_sup <- na.omit(sup) #Tira linhas sem informações de concentração
halo_opt <- na.omit(opt)



#********** CONCENTRAÇÕES ÓTIMAS **********

bad_opt <- c() #Cria vetor em que será colocado quais linhas devem ser retiradas por não serem NaCl

for (n in seq_along(halo_opt[,1])){ #passa por cada bicho
    if (halo_opt[n,6] != "s"){ #Vê quais linhas não correspondem a NaCl | s= sim, é NaCl. Qualquer coisa diferente quer dizer outros sais
        bad_opt = append(bad_opt, n, after=length(bad_opt))} #Armazena no vetor qual linha deve ser retirada

    else if (halo_opt[n,7] == "p"){ # p=porcentual m/v | loop faz conversão para mol/L
        halo_opt[n,3] <- ((halo_opt[n,3]*10)/58.44)
        halo_opt[n,4] <- ((halo_opt[n,4]*10)/58.44)}

    else if (halo_opt[n,7] == "g"){ # g=g/L | loop faz conversão para mol/L
        halo_opt[n,3] <- (halo_opt[n,3]/58.44)
        halo_opt[n,4] <- (halo_opt[n,4]/58.44)}
}

halo_opt <- halo_opt[-bad_opt,] #Retira linhas marcadas
#halo_opt <- na.omit(halo_opt) #Retira qualquer linha vazia (provavelmente redundante)
write.table(halo_opt, file="halo_opt.tsv") #Salva arquivo .tsv da tabela limpa de concentrações ótimas

table(halo_opt[,2])


#********** CONCENTRAÇÕES SUPORTADAS **********

bad_sup <- c() #Cria vetor em que será colocado quais linhas devem ser retiradas por não serem NaCl

for (n in seq_along(halo_sup[,1])){ #passa por cada bicho
    if (halo_sup[n,6] != "s"){  #Vê quais linhas não correspondem a NaCl | s= sim, é NaCl. Qualquer coisa diferente quer dizer outros sais. Vê se classificação está como "not" ou "Not", que significa que não é halófilo.
        bad_sup = append(bad_sup, n, after=length(bad_sup))}

    else if (halo_sup[n,7] == "p"){ # p=porcentual m/v | loop faz conversão para mol/L
        halo_sup[n,3] <- ((halo_sup[n,3]*10)/58.44)
        halo_sup[n,4] <- ((halo_sup[n,4]*10)/58.44)}

    else if (halo_sup[n,7] == "g"){ # g=g/L | loop faz conversão para mol/L
        halo_sup[n,3] <- (halo_sup[n,3]/58.44)
        halo_sup[n,4] <- (halo_sup[n,4]/58.44)}
}

halo_sup <- halo_sup[-bad_sup,] #Retira linhas marcadas
#halo_sup <- na.omit(halo_sup) #Retira qualquer linha vazia (provavelmente redundante)
write.table(halo_sup, file="halo_sup.tsv") #Salva arquivo .tsv da tabela limpa de concentrações suportadas

table(halo_sup[,2])
