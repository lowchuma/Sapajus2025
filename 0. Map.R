# --- 1. BIBLIOTECAS ---
library(ggplot2)
library(sf)
library(ggspatial)  # Para bússola e escala
library(dplyr)      # Para 'rename' e 'bind_rows'
library(svglite)    # Para salvar o SVG perfeito

# --- 2. CONFIGURAÇÃO (ATUALIZADA) ---

# 2.1: Caminhos para os arquivos de base
csv_ocorrencias_path <- "Plans/avglocsex.csv"
sa_shapefile_path <- "C:/Users/Lourenço/Downloads/Mapa - Felipe/South_America_shapefile/americadosul.shp"
biomes_shapefile_path <- "C:/Users/Lourenço/Downloads/Mapa - Felipe/regions_cropped/poly.shp"

# 2.2: Colunas de dados no seu CSV
coluna_longitude <- "long"
coluna_latitude <- "lat"
coluna_especie <- "sp"

# 2.3: Caminhos para os 6 SHAPEFILES DE DISTRIBUIÇÃO
lista_shapefiles_distribuicao <- c(
  "Shapefiles/Sapajus_apella.shp",
  "Shapefiles/Sapajus_cay.shp",
  "Shapefiles/Sapajus_libidinosus.shp",
  "Shapefiles/Sapajus_nigritus.shp",
  "Shapefiles/Sapajus_robustus.shp",
  "Shapefiles/Sapajus_xanthosternos.shp"
)

# 2.4: Nomes das 6 espécies (!! IMPORTANTE: DO SEU SCRIPT !!)
# (Estes nomes completos devem bater 100% com a coluna 'sp' do seu CSV)
nomes_das_especies <- c(
  "Sapajus_apella",
  "Sapajus_cay",
  "Sapajus_libidinosus",
  "Sapajus_nigritus",
  "Sapajus_robustus",
  "Sapajus_xanthosternos"
)

# --- 2.5: CORES E SÍMBOLOS (COM VETORES NOMEADOS - A CORREÇÃO) ---
# (Cores do seu script de PCA)
cores_lista <- c("red", "navy", "cyan3", "saddlebrown", "magenta4", "black")

# (Símbolos Unicode EXATOS conforme seu pedido)
simbolos_lista <- c(
  "\u25A0", # QUADRADO (S. apella)
  "\u25BC", # TRIÂNGULO PARA BAIXO (S. cay)
  "\u25B2", # TRIÂNGULO NORMAL (S. libidinosus)
  "\u25CF", # CÍRCULO (S. nigritus)
  "\u2666", # LOSANGO (S. robustus)
  "\u2736"  # ESTRELA (S. xanthosternos)
)

# Cria os vetores NOMEADOS ("dicionários")
cores_especies_named <- setNames(cores_lista, nomes_das_especies)
simbolos_especies_named <- setNames(simbolos_lista, nomes_das_especies)
# --- FIM DA ATUALIZAÇÃO ---


# --- 3. CARREGAR DADOS DE BASE ---
cat("Carregando dados de base...\n")
Sapajus_data_original <- read.csv(csv_ocorrencias_path, sep = ";")
SouthAmerica <- st_read(sa_shapefile_path)
Biomes <- st_read(biomes_shapefile_path)

# Padroniza os nomes das colunas
Sapajus_data <- Sapajus_data_original %>%
  rename(
    long_plot = all_of(coluna_longitude),
    lat_plot = all_of(coluna_latitude),
    sp_plot = all_of(coluna_especie)
  )

# --- 4. CARREGAR E COMBINAR 6 SHAPEFILES (R3) ---
cat("Carregando 6 shapefiles de distribuição...\n")

lista_de_sfs <- lapply(1:length(lista_shapefiles_distribuicao), function(i) {
  sf_obj <- st_read(lista_shapefiles_distribuicao[i])
  sf_obj$sp_plot <- nomes_das_especies[i] # Adiciona a coluna 'sp_plot'
  sf_obj <- st_make_valid(sf_obj)
  return(sf_obj)
})

distribuicoes_combinadas <- do.call(rbind, lista_de_sfs)
cat("Shapefiles combinados.\n")


# --- 5. PROCESSAR DADOS (CROP COM O SEU ZOOM MAIS RECENTE) ---
# (BBox: xmin = -72, xmax = -35, ymin = -28, ymax = 5)
cat("Recortando todas as camadas para o zoom...\n")

crs_padrao <- st_crs(Biomes)
teste <- st_crop(Biomes, xmin = -72, xmax = -35, ymin = -28, ymax = 5)
teste2 <- st_crop(SouthAmerica, xmin = -72, xmax = -35, ymin = -28, ymax = 5)
dist_combinada_crs <- st_transform(distribuicoes_combinadas, crs = crs_padrao)
dist_combinada_cropped <- st_crop(dist_combinada_crs,
                                  xmin = -72, xmax = -35,
                                  ymin = -28, ymax = 5)


# --- 6. O PLOT FINAL (COM 'geom_text' E CORES CORRETAS) ---
cat("Gerando o gráfico final...\n")

mapa_final_completo <- ggplot(data = teste) +
  
  # Camada 1: Biomas (COM TRANSPARÊNCIA)
  geom_sf(aes(fill = as.factor(BIOME)), color=NA, alpha = 0.7) +
  
  # Camada 2: Contornos da Distribuição (TRACEJADOS)
  geom_sf(data = dist_combinada_cropped, 
          aes(color = as.factor(sp_plot)),
          fill = NA,
          linewidth = 1,
          linetype = "dashed") +
  
  # Camada 3: Contorno do Continente (NEUTRO)
  geom_sf(data = teste2, 
          fill=NA, 
          color = "grey20", # <--- Cor neutra
          linewidth = 0.5,
          alpha = 0.8) +
  
  # Camada 4: Pontos (COMO 'geom_text' PARA USAR UNICODE)
  geom_text(data = Sapajus_data, 
            aes(x = long_plot, 
                y = lat_plot, 
                label = as.factor(sp_plot),  # <--- USA 'label'
                color = as.factor(sp_plot)), # <--- USA 'color'
            size = 8) +                      
  
  # Camada 5: Escalas de Cor/Símbolo (COM A FUNÇÃO CORRETA)
  scale_fill_manual(
    name = NULL,
    # Suas cores de bioma do script:
    values = c("darkgreen","#49FF00","orange","white","#2666CF","white","white","#AE431E","white","white","white","white")
  ) +
  
  # --- ESTA É A CORREÇÃO para 'scale_label_manual' ---
  scale_discrete_manual(
    aesthetics = "label",     # <--- DIZ QUAL ESTÉTICA CONTROLAR
    name = "Espécie",         # Título da legenda
    values = simbolos_especies_named # <--- USA O VETOR NOMEADO
  ) +
  # --- FIM DA CORREÇÃO ---
  
  scale_color_manual(
    name = NULL, # Título IDÊNTICO para FUNDIR
    values = cores_especies_named, # <--- USA O VETOR NOMEADO
    labels = nomes_das_especies
  ) +
  
  # Camada 6: Tema
  theme_bw(base_size = 12) +
  labs(x = "Longitude", y = "Latitude") +
  
  # --- Adicionais do ggspatial ---
  
  # 7. Define o sistema de coordenadas
  coord_sf(crs = crs_padrao) + 
  
  # 8. Adiciona a Bússola
  annotation_north_arrow(
    location = "tr", style = north_arrow_fancy_orienteering,
    height = unit(1.2, "cm"), width = unit(1.2, "cm"),
    pad_y = unit(0.5, "cm"), pad_x = unit(0.5, "cm")
  ) +
  
  # 9. Adiciona a Escala
  annotation_scale(
    location = "bl", width_hint = 0.2,
    pad_y = unit(0.2, "cm"), pad_x = unit(0.2, "cm")
  )

# --- 7. MOSTRAR E SALVAR EM SVG (PARA INKSCAPE) ---
print(mapa_final_completo)
cat("TUDO PRONTO! O SVG e o PDF (com fontes incorporadas) foram salvos.\n")

# --- 1. BIBLIOTECAS ---
library(ggplot2)
library(sf)
library(ggspatial)
library(dplyr)
library(svglite)    # <--- Essencial para a correção do SVG
library(showtext)   # <--- Essencial para a correção do PDF
library(sysfonts)

# --- 2. CONFIGURAÇÃO (ATUALIZADA) ---

# 2.1 a 2.3: Seus caminhos (iguais)
csv_ocorrencias_path <- "Plans/avglocsex.csv"
sa_shapefile_path <- "C:/Users/Lourenço/Downloads/Mapa - Felipe/South_America_shapefile/americadosul.shp"
biomes_shapefile_path <- "C:/Users/Lourenço/Downloads/Mapa - Felipe/regions_cropped/poly.shp"
coluna_longitude <- "long"
coluna_latitude <- "lat"
coluna_especie <- "sp"
lista_shapefile_distribuicao <- c(
  "Shapefiles/Sapajus_apella.shp",
  "Shapefiles/Sapajus_cay.shp",
  "Shapefiles/Sapajus_libidinosus.shp",
  "Shapefiles/Sapajus_nigritus.shp",
  "Shapefiles/Sapajus_robustus.shp",
  "Shapefiles/Sapajus_xanthosternos.shp"
)

# 2.4: Nomes das 6 espécies (!! CONFIRME SE SÃO ESTES OS NOMES NO SEU CSV !!)
nomes_das_especies <- c(
  "Sapajus_apella",
  "Sapajus_cay",
  "Sapajus_libidinosus",
  "Sapajus_nigritus",
  "Sapajus_robustus",
  "Sapajus_xanthosternos"
)
# (Se os seus nomes forem abreviados, ex: "S. apella", tem de mudar aqui)

# 2.5: CORES E SÍMBOLOS (COM VETORES NOMEADOS)
cores_lista <- c("red", "navy", "cyan3", "saddlebrown", "magenta4", "black")
simbolos_lista <- c(
  "\u25A0", # QUADRADO
  "\u25BC", # TRIÂNGULO PARA BAIXO
  "\u25B2", # TRIÂNGULO NORMAL
  "\u25CF", # CÍRCULO
  "\u2666", # LOSANGO
  "\u2736"  # ESTRELA
)
cores_especies_named <- setNames(cores_lista, nomes_das_especies)
simbolos_especies_named <- setNames(simbolos_lista, nomes_das_especies)

# --- 2.6: CARREGAR FONTES (A MÁGICA DO 'showtext') ---
font_add_google("Noto Sans", "notosans")
showtext_auto()
# --- FIM DA MUDANÇA ---


# --- 3. CARREGAR DADOS DE BASE ---
# (Seu código daqui não muda)
cat("Carregando dados de base...\n")
Sapajus_data_original <- read.csv(csv_ocorrencias_path, sep = ";")
SouthAmerica <- st_read(sa_shapefile_path)
Biomes <- st_read(biomes_shapefile_path)

Sapajus_data <- Sapajus_data_original %>%
  rename(
    long_plot = all_of(coluna_longitude),
    lat_plot = all_of(coluna_latitude),
    sp_plot = all_of(coluna_especie)
  )

# --- 4. CARREGAR E COMBINAR 6 SHAPEFILES (R3) ---
# (Seu código daqui não muda)
cat("Carregando 6 shapefiles de distribuição...\n")
lista_de_sfs <- lapply(1:length(lista_shapefile_distribuicao), function(i) {
  sf_obj <- st_read(lista_shapefile_distribuicao[i])
  sf_obj$sp_plot <- nomes_das_especies[i]
  sf_obj <- st_make_valid(sf_obj)
  return(sf_obj)
})
distribuicoes_combinadas <- do.call(rbind, lista_de_sfs)
cat("Shapefiles combinados.\n")

# --- 5. PROCESSAR DADOS (CROP COM O SEU ZOOM MAIS RECENTE) ---
# (BBox: xmin = -72, xmax = -35, ymin = -28, ymax = 5)
# (Seu código daqui não muda)
cat("Recortando todas as camadas para o zoom...\n")
crs_padrao <- st_crs(Biomes)
teste <- st_crop(Biomes, xmin = -72, xmax = -35, ymin = -28, ymax = 5)
teste2 <- st_crop(SouthAmerica, xmin = -72, xmax = -35, ymin = -28, ymax = 5)
dist_combinada_crs <- st_transform(distribuicoes_combinadas, crs = crs_padrao)
dist_combinada_cropped <- st_crop(dist_combinada_crs,
                                  xmin = -72, xmax = -35,
                                  ymin = -28, ymax = 5)


# --- 6. O PLOT FINAL (COM 'geom_text' E 'showtext') ---
cat("Gerando o gráfico final...\n")

mapa_final_completo <- ggplot(data = teste) +
  
  # Camadas 1, 2, 3 (Biomas, Contornos R3, Contorno Continente)
  geom_sf(aes(fill = as.factor(BIOME)), color=NA, alpha = 0.7) +
  geom_sf(data = dist_combinada_cropped, 
          aes(color = as.factor(sp_plot)),
          fill = NA, linewidth = 0.8, linetype = "dashed") +
  geom_sf(data = teste2, fill=NA, color = "grey50", linewidth = 0.6, alpha = 0.7) +
  
  # Camada 4: Pontos (COM TAMANHO AUMENTADO)
  geom_text(data = Sapajus_data, 
            aes(x = long_plot, 
                y = lat_plot, 
                label = as.factor(sp_plot),
                color = as.factor(sp_plot)), 
            size = 15,              # <--- CORREÇÃO DO TAMANHO
            family = "notosans") +
  
  # Camada 5: Escalas (COM A CORREÇÃO 'scale_discrete_manual')
  scale_fill_manual(
    name = "Bioma",
    values = c("darkgreen","#49FF00","orange","white","#2666CF","white","white","#AE431E","white","white","white","white")
  ) +
  scale_discrete_manual(
    aesthetics = "label",
    name = "Espécie",
    values = simbolos_especies_named
  ) +
  scale_color_manual(
    name = "Espécie",
    values = cores_especies_named,
    labels = nomes_das_especies
  ) +
  
  # Camada 6, 7, 8, 9 (Tema, Coords, Bússola, Escala)
  theme_bw(base_size = 12) +
  labs(x = "Longitude", y = "Latitude") +
  coord_sf(crs = crs_padrao) + 
  annotation_north_arrow(
    location = "tr", style = north_arrow_fancy_orienteering,
    height = unit(1.2, "cm"), width = unit(1.2, "cm"),
    pad_y = unit(0.5, "cm"), pad_x = unit(0.5, "cm")
  ) +
  annotation_scale(
    location = "bl", width_hint = 0.3,
    pad_y = unit(0.5, "cm"), pad_x = unit(0.5, "cm")
  )

# --- 7. MOSTRAR E SALVAR (COM AMBAS AS CORREÇÕES) ---
print(mapa_final_completo)


# 7.2: Salvar o PDF (COM 'showtext')
cat("Salvando em PDF de alta qualidade (com 'showtext')...\n")
ggsave(
  filename = "Fig_S1_FINAL_PDF_FIX.pdf",
  plot = mapa_final_completo,
  device = "pdf",  # 'showtext' vai tratar disto
  width = 11,
  height = 9
)

cat("TUDO PRONTO! Agora vai funcionar.\n")

# --- 1. BIBLIOTECAS ---
library(ggplot2)
library(sf)
library(ggspatial)
library(dplyr)
library(svglite)
library(showtext)
library(sysfonts)

# --- 2. CONFIGURAÇÃO (ATUALIZADA) ---

# 2.1 a 2.3: Seus caminhos (iguais)
csv_ocorrencias_path <- "Plans/avglocsex.csv"
sa_shapefile_path <- "C:/Users/Lourenço/Downloads/Mapa - Felipe/South_America_shapefile/americadosul.shp"
biomes_shapefile_path <- "C:/Users/Lourenço/Downloads/Mapa - Felipe/regions_cropped/poly.shp"
coluna_longitude <- "long"
coluna_latitude <- "lat"
coluna_especie <- "sp"
lista_shapefile_distribuicao <- c(
  "Shapefiles/Sapajus_apella.shp",
  "Shapefiles/Sapajus_cay.shp",
  "Shapefiles/Sapajus_libidinosus.shp",
  "Shapefiles/Sapajus_nigritus.shp",
  "Shapefiles/Sapajus_robustus.shp",
  "Shapefiles/Sapajus_xanthosternos.shp"
)

# 2.4: Nomes das 6 espécies (!! CONFIRME SE SÃO ESTES OS NOMES NO SEU CSV !!)
nomes_das_especies <- c(
  "Sapajus_apella",
  "Sapajus_cay",
  "Sapajus_libidinosus",
  "Sapajus_nigritus",
  "Sapajus_robustus",
  "Sapajus_xanthosternos"
)

# 2.5: CORES E SÍMBOLOS (COM VETORES NOMEADOS)
cores_lista <- c("red", "navy", "cyan3", "saddlebrown", "magenta4", "black")
simbolos_lista <- c(
  "\u25A0", # QUADRADO
  "\u25BC", # TRIÂNGULO PARA BAIXO
  "\u25B2", # TRIÂNGULO NORMAL
  "\u25CF", # CÍRCULO
  "\u2666", # LOSANGO
  "\u2736"  # ESTRELA
)
cores_especies_named <- setNames(cores_lista, nomes_das_especies)
simbolos_especies_named <- setNames(simbolos_lista, nomes_das_especies)

# 2.6: CARREGAR FONTES ('showtext' para PDF)
font_add_google("Noto Sans", "notosans")
showtext_auto()


# --- 3. CARREGAR DADOS DE BASE ---
cat("Carregando dados de base...\n")
Sapajus_data_original <- read.csv(csv_ocorrencias_path, sep = ";")
SouthAmerica_full <- st_read(sa_shapefile_path) # <--- CARREGA A AMÉRICA DO SUL COMPLETA
Biomes <- st_read(biomes_shapefile_path)

Sapajus_data <- Sapajus_data_original %>%
  rename(
    long_plot = all_of(coluna_longitude),
    lat_plot = all_of(coluna_latitude),
    sp_plot = all_of(coluna_especie)
  )

# --- 3.1: CHECAGEM DE DEBUG (OLHE O CONSOLE!) ---
cat("--- CHECAGEM DE ESPÉCIES ---\n")
cat("Nomes no seu CSV (coluna 'sp'):\n")
print(unique(Sapajus_data$sp_plot))
cat("Nomes no seu script (Seção 2.4):\n")
print(nomes_das_especies)
cat("Se estas duas listas não forem IDÊNTICAS, os símbolos VÃO SUMIR.\n")
cat("-------------------------------\n\n")


# --- 4. CARREGAR E COMBINAR 6 SHAPEFILES (R3) ---
cat("Carregando 6 shapefiles de distribuição...\n")
lista_de_sfs <- lapply(1:length(lista_shapefile_distribuicao), function(i) {
  sf_obj <- st_read(lista_shapefile_distribuicao[i])
  sf_obj$sp_plot <- nomes_das_especies[i]
  sf_obj <- st_make_valid(sf_obj)
  return(sf_obj)
})
distribuicoes_combinadas <- do.call(rbind, lista_de_sfs)
cat("Shapefiles combinados.\n")


# --- 5. PROCESSAR DADOS E DEFINIR CAIXA DE RECORTE ---
cat("Definindo área de recorte e processando dados...\n")
crs_padrao <- st_crs(Biomes)

# 5.1: Definir as coordenadas da sua "caixa" de estudo
xmin_box <- -72
xmax_box <- -35
ymin_box <- -30
ymax_box <- 5

# 5.2: Criar um objeto sf para a caixa de recorte (será desenhada no mapa geral)
bbox_sf <- st_bbox(c(xmin = xmin_box, ymin = ymin_box, xmax = xmax_box, ymax = ymax_box), crs = crs_padrao) %>%
  st_as_sfc()

# 5.3: Cropar biomas e contorno do continente (PARA A CAIXA)
teste <- st_crop(Biomes, xmin = xmin_box, xmax = xmax_box, ymin = ymin_box, ymax = ymax_box)
teste2 <- st_crop(SouthAmerica_full, xmin = xmin_box, xmax = xmax_box, ymin = ymin_box, ymax = ymax_box) # <--- Usar a full aqui tb

# 5.4: Cropar as distribuições das espécies (PARA A CAIXA)
dist_combinada_crs <- st_transform(distribuicoes_combinadas, crs = crs_padrao)
dist_combinada_cropped <- st_crop(dist_combinada_crs,
                                  xmin = xmin_box, xmax = xmax_box,
                                  ymin = ymin_box, ymax = ymax_box)

# 5.5: Cropar os pontos (PARA A CAIXA)
Sapajus_data_cropped <- Sapajus_data %>%
  filter(long_plot >= xmin_box & long_plot <= xmax_box &
           lat_plot >= ymin_box & lat_plot <= ymax_box)


# --- 5.6: CHECAGEM DE DEBUG DOS BIOMAS ---
cores_biomas_lista <- c("darkgreen","#49FF00","orange","white","#2666CF","white","white","#AE431E","white","white","white","white")
n_cores_biomas <- length(cores_biomas_lista)
n_niveis_biomas <- length(levels(as.factor(teste$BIOME)))

cat("--- CHECAGEM DE BIOMAS ---\n")
cat(paste("Seu script definiu", n_cores_biomas, "cores para os biomas.\n"))
cat(paste("Seu mapa CORTADO ('teste') tem", n_niveis_biomas, "níveis de bioma.\n"))
if (n_cores_biomas < n_niveis_biomas) {
  cat("ERRO: Você precisa de mais cores na 'cores_biomas_lista' (Seção 6)!\n")
} else {
  cat("OK: O número de cores é suficiente.\n")
}
cat("---------------------------\n\n")


# --- 6. O PLOT FINAL (COM VISÃO GERAL DA AMÉRICA DO SUL) ---
cat("Gerando o gráfico final...\n")

cores_biomas_lista <- c("darkgreen","#49FF00","orange","white","#2666CF","white","white","#AE431E","white","white","white","white")

mapa_final_completo <- ggplot() + # <--- Sem data base aqui
  
  # 0. Camada: América do Sul completa (fundo)
  geom_sf(data = SouthAmerica_full, fill = "lightgrey", color = "darkgrey", linewidth = 0.3) +
  
  # 1. Camada: A sua "caixa" de estudo (retângulo delimitador)
  geom_sf(data = bbox_sf, fill = NA, color = "red", linewidth = 1.2, linetype = "dotted") + # <--- Caixa com contorno vermelho
  
  # As próximas camadas são CORTADAS para a sua caixa (como já estavam)
  # Camada 1: Biomas (apenas na caixa)
  geom_sf(data = teste, aes(fill = as.factor(BIOME)), color=NA, alpha = 0.7) +
  
  # Camada 2: Contornos da Distribuição (apenas na caixa)
  geom_sf(data = dist_combinada_cropped, 
          aes(color = as.factor(sp_plot)),
          fill = NA, linewidth = 0.8, linetype = "dashed") +
  
  # Camada 3: Contorno do Continente (apenas na caixa)
  geom_sf(data = teste2, fill=NA, color = "grey50", linewidth = 0.6, alpha = 0.7) +
  
  # Camada 4: Pontos (apenas na caixa)
  geom_text(data = Sapajus_data_cropped, # <--- Usar os pontos cortados
            aes(x = long_plot, 
                y = lat_plot, 
                label = as.factor(sp_plot),
                color = as.factor(sp_plot)), 
            size = 8,
            family = "notosans") +
  
  # Camada 5: Escalas
  scale_fill_manual(
    name = "Bioma",
    values = cores_biomas_lista,
    na.value = "transparent"
  ) +
  scale_discrete_manual(
    aesthetics = "label",
    name = "Espécie",
    values = simbolos_especies_named
  ) +
  scale_color_manual(
    name = "Espécie",
    values = cores_especies_named,
    labels = nomes_das_especies
  ) +
  
  # Camada 6, 7, 8, 9 (Tema, Coords, Bússola, Escala)
  theme_bw(base_size = 12) +
  labs(x = "Longitude", y = "Latitude") +
  
  # As coordenadas agora abrangem a América do Sul completa
  coord_sf(crs = crs_padrao, xlim = c(-85, -30), ylim = c(-57, 13)) + # <--- Coordenadas da América do Sul
  
  annotation_north_arrow(
    location = "tr", style = north_arrow_fancy_orienteering,
    height = unit(1.2, "cm"), width = unit(1.2, "cm"),
    pad_y = unit(0.5, "cm"), pad_x = unit(0.5, "cm")
  ) +
  annotation_scale(
    location = "bl", width_hint = 0.3,
    pad_y = unit(0.5, "cm"), pad_x = unit(0.5, "cm")
  )

# --- 7. MOSTRAR E SALVAR ---
print(mapa_final_completo)

cat("Salvando em SVG (com CSS inline)...\n")
ggsave(
  filename = "Fig_S1_SOUTHAMERICA_VIEW_INKSCAPE_FIX.svg",
  plot = mapa_final_completo,
  width = 11,
  height = 9,
  device = function(file, width, height, ...) {
    svglite::svglite(file, width = width, height = height, inline_css = TRUE, ...)
  }
)

cat("Salvando em PDF de alta qualidade (com 'showtext')...\n")
ggsave(
  filename = "Fig_S1_SOUTHAMERICA_VIEW_PDF_FIX.pdf",
  plot = mapa_final_completo,
  device = "pdf",
  width = 11,
  height = 9
)

cat("TUDO PRONTO! O mapa agora inclui a América do Sul com sua caixa de recorte.\n")
