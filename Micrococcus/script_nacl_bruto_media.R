#-------MANIPULAÇÃO DAS TABELAS DO PLATE READER E PLOTS-----------
#------------COLB perclorato de magnésio (Mg(ClO4)2-----------------------

# Carregar pacotes necessários
install.packages("readxl")
library(readxl)
library(dplyr)
install.packages("openxlsx")
library(openxlsx)

# Carregar os dados brutos
raw_data <- read_excel("/Users/isabe/Downloads/nepalensis/Curva_perc_DO/colB_22_7/colB_perc.xlsx", sheet = "Plate 2 - Raw Data")

# Remover as colunas cujos nomes terminam com 1 (nessa placa era água)
dados_filtrados <- raw_data %>% select(-matches("^[A-H]1$"))

# Remover todas as colunas que correspondem à linha A da placa
dados_filtrados <- dados_filtrados %>% select(-starts_with("A"))


# Converter todas as colunas (exceto Time) para numérico
dados_filtrados[,-1] <- lapply(dados_filtrados[,-1], as.numeric)

# Manter a coluna de tempo separada para que ela não seja alterada
time_col <- dados_filtrados$`Time`

# Subtrair o valor do blank correspondente para cada linha e concentração
dados_corrigidos <- dados_filtrados

# Iterar sobre cada coluna de blank (H2 a H12) e subtrair das colunas correspondentes
for (blank_col in grep("^H[2-9]|^H1[0-2]$", colnames(dados_filtrados), value = TRUE)) {
  
  # Identificar o padrão da concentração associada ao blank (ex: H3 -> colunas com sufixo 3)
  conc_pattern <- gsub("H", "", blank_col)
  
  # Encontrar todas as colunas que correspondem à mesma concentração
  target_cols <- grep(paste0("[B-G]", conc_pattern, "$"), colnames(dados_corrigidos), value = TRUE)
  
  # Subtrair o valor do blank correspondente de cada linha
  for (col in target_cols) {
    dados_corrigidos[[col]] <- dados_corrigidos[[col]] - dados_filtrados[[blank_col]]
  }
}

# Adicionar de volta a coluna de tempo sem alterar o formato
dados_corrigidos <- cbind(Time = time_col, dados_corrigidos[,-1])

# Verificar o resultado
head(dados_corrigidos)

write.xlsx(dados_corrigidos, file = "/Users/isabe/Downloads/nepalensis/Curva_perc_DO/colB_22_7/dados_corrigidos.xlsx", rowNames = FALSE)

#-----GRÁFICO FACETADO REPLICATAS---------------

library(ggplot2)
library(tidyr)
library(dplyr)

#-----------TABELA PARA PLOT-------------
#time como numérico
dados_corrigidos$Time <- as.numeric(dados_corrigidos$Time)

# Transformar os dados em formato longo, agora focando nas colunas B2 a G12
dados_long <- dados_corrigidos %>%
  select(Time, B2:G12) %>%  # Seleciona explicitamente as colunas B2 a G12 junto com o tempo
  pivot_longer(
    cols = B2:G12,  # Seleciona as colunas de interesse (B2 a G12)
    names_to = c("Replicata", "Concentracao"),
     names_pattern = "([B-G])(2|3|4|5|6|7|8|9|10|11|12)",  # Define explicitamente os números de 2 a 12
    values_to = "DO"
  ) %>%
  mutate(
    DO = abs(DO),  # Transformar valores de DO em módulo
    Concentracao = case_when(
      Concentracao == "2" ~ "0",
      Concentracao == "3" ~ "0.1",
      Concentracao == "4" ~ "0.2",
      Concentracao == "5" ~ "0.3",
      Concentracao == "6" ~ "0.4",
      Concentracao == "7" ~ "0.5",
      Concentracao == "8" ~ "0.6",
      Concentracao == "9" ~ "0.7",
      Concentracao == "10" ~ "0.8",
      Concentracao == "11" ~ "0.9",
      Concentracao == "12" ~ "1.0"
    )
  )# Verificar o resultado
head(dados_long)
write.xlsx(dados_long, file = "/Users/isabe/Downloads/nepalensis/Curva_perc_DO/colB_22_7/dados_long.xlsx", rowNames = FALSE)



# Criar o gráfico com eixo y em escala logarítmica
plot_replicatas_perc <- ggplot(dados_long, aes(x = Time, y = DO, color = Concentracao, group = Concentracao)) +
  geom_point(size=1) +
  facet_wrap(~ Replicata) +
  scale_y_log10() +  # Aplicar transformação logarítmica ao eixo y
  labs(
    title = "Curva de Crescimento por Replicata e Concentração de Perclorato",
    x = "Tempo",
    y = "Densidade Óptica (DO - Escala Log)"
  ) +
  theme_minimal()

plot_replicatas_perc

ggsave("/Users/isabe/Downloads/nepalensis/Curva_perc_DO/colB_22_7/plot_replicatas.png", plot = last_plot(), width = 8, height = 6, dpi = 300, bg = "white")


#--------GRÁFICO DAS MÉDIAS----------------------

#TABELA DAS MÉDIAS
library(dplyr)
library(readxl)
library(openxlsx)

# Carregar os dados corrigidos do arquivo Excel
dados_corrigidos <- read_excel("/Users/isabe/Downloads/nepalensis/Curva_nacl_DO/colB_22_7/dados_corrigidos.xlsx")


# Transformar os dados em formato longo, agora focando nas colunas B a G
dados_long <- dados_corrigidos %>%
  select(Time, B3:G10) %>%  # Seleciona explicitamente as colunas B3 a G10 junto com o tempo
  pivot_longer(
    cols = B3:G10, # Seleciona as colunas de interesse (B3 a G10)
    names_to = c("Replicata", "Concentracao"),
    names_pattern = "([B-G])([0-9]+)",
    values_to = "DO"
  ) %>%
  mutate(
    DO = abs(DO),  # Transformar valores de DO em módulo
    Concentracao = case_when(
      Concentracao == "3" ~ "0",
      Concentracao == "4" ~ "0.5",
      Concentracao == "5" ~ "1",
      Concentracao == "6" ~ "1.5",
      Concentracao == "7" ~ "2",
      Concentracao == "8" ~ "3",
      Concentracao == "9" ~ "3.5",
      Concentracao == "10" ~ "4"
    )
  )



dados_long

# Filtrar os dados para excluir a replicata E na concentração 4
dados_long_filtrado <- dados_long %>%
  filter(!(Replicata == "E" & Concentracao == "4") & !(Replicata == "F" & Concentracao == "4"))

 
# Verificar o resultado
print(tail(dados_long_filtrado, 20))


# Calcular a média e o desvio padrão das replicatas por tempo e por concentração
dados_medias <- dados_long_filtrado %>%
  group_by(Time, Concentracao) %>%
  summarise(
    Media_DO = mean(DO, na.rm = TRUE),
    DesvioPadrao_DO = sd(DO, na.rm = TRUE)
  )

# Verificar o resultado
print(head(dados_medias))

# Salvar os dados de média e desvio padrão em um novo arquivo Excel
write.xlsx(dados_medias, file = "/Users/isabe/Downloads/nepalensis/Curva_nacl_DO/colB_22_7/dados_media.xlsx", rowNames = FALSE)


#PLOT MEDIAS

library(ggplot2)
library(dplyr)
library(readxl)
library(openxlsx)

# Carregar os dados de média e desvio padrão do arquivo Excel
dados_medias <- read_excel("/Users/isabe/Downloads/nepalensis/Curva_nacl_DO/colB_22_7/dados_media.xlsx")

# Garantir que os valores de DO estejam em módulo
dados_medias <- dados_medias %>%
  mutate(
    Media_DO = abs(Media_DO), 
    DesvioPadrao_DO = abs(DesvioPadrao_DO)
  )



# Criar o gráfico com eixo y em escala logarítmica
plot_medias <- ggplot(dados_medias, aes(x = Time, y = Media_DO, color = Concentracao, group = Concentracao)) +
  geom_point() +
  geom_errorbar(aes(ymin = Media_DO - DesvioPadrao_DO, ymax = Media_DO + DesvioPadrao_DO), width = 0.2) +
  scale_y_log10() +  # Aplicar transformação logarítmica ao eixo y
  labs(
    title = "Curva de Crescimento Média por Concentração de NaCl",
    x = "Tempo",
    y = "Densidade Óptica (DO - Escala Log)"
  ) +
  theme_minimal()

plot_medias <- ggsave("/Users/isabe/Downloads/nepalensis/Curva_nacl_DO/colB_22_7/plot_medias.png", plot = last_plot(), width = 8, height = 6, dpi = 300, bg = "white")























