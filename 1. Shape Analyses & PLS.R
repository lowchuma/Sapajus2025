### Parte básica do script ###

### SetWD - carrega as pastas onde estão os dados###
setwd("C:/Users/Lourenço/OneDrive/UFSM/LabMastozoo/SAPAJUS/Data - From forests to savannas")

### Carregar pacotes ###

# Pacotes não diretamente instaláveis pelo R (via github ou CRAN)

#install_github("geomorphR/geomorph", ref = "Stable", build_vignettes = TRUE)
#install_github("mlcollyer/RRPP")
#install.packages("geiger")
#install.packages("RColorBrewer")
#install.packages("vegan")
#install.packages("tidyverse")
#install.packages ("rgl")

# Desenvolvimento de Pacotes
library(devtools)   # Ferramentas para facilitar o desenvolvimento de pacotes R

# Manipulação e Análise de Dados
library(tidyverse)  # Conjunto de pacotes para ciência de dados, incluindo ggplot2, dplyr, tidyr, readr, purrr e outros
library(reshape2)   # Funções para reorganizar dados entre formatos amplos e longos

# Estatística e Modelagem
library(MASS)       # Funções e conjuntos de dados para estatísticas aplicadas
library(psych)      # Ferramentas para análise psicológica e psicométrica
library(klaR)       # Funções para classificação e visualização, incluindo análise discriminante regularizada

# Visualização de Dados
library(ggplot2)       # Sistema para criar gráficos baseados em camadas
library(RColorBrewer)  # Paletas de cores para visualização
library(ggord)         # Criação de gráficos de ordenação com ggplot2
library(dotwhisker)    # Gráficos de coeficientes para modelos estatísticos

# Ecologia e Biologia Evolutiva
library(vegan)     # Métodos de ordenação e análise de diversidade para ecologia comunitária
library(ape)       # Análises de filogenética e evolução
library(phytools)  # Ferramentas para biologia comparativa filogenética
library(picante)   # Integração de filogenias e ecologia
library(geiger)    # Métodos estatísticos para análise de dados filogenéticos

# Morfometria Geométrica
library(geomorph)  # Análises de morfometria geométrica de dados de pontos de referência 2D e 3D
library(Morpho)    # Ferramentas para morfometria geométrica e processamento de malhas
library(shapes)    # Rotinas para análise estatística de formas baseadas em pontos de referência
library(Rvcg)      # Processamento e análise de malhas 3D

# Modelagem de Distribuição de Espécies
library(usdm)      # Análise de incerteza para modelos de distribuição de espécies

# Estatísticas Base
library(stats)     # Funções estatísticas básicas do R

### Ler os dados brutos e criar os fatores de agrupamento ###

## Amostras, nesse caso, são dados de cada localidade amostrada 
# (onde haviam mais de um espécime foi feito a média)

dados <- readland.tps("Datasets/tps/avglocsex.tps", specID =  "ID")

fatores <- read.csv("Planilhas/avglocsex.csv", sep=";")

names(fatores)

species <- as.factor(fatores$sp)
summary(species)

biome <- as.factor(fatores$biome)
summary(biome)

sex <- as.factor(fatores$sex)
summary (sex)

latitude <-as.numeric(fatores$lat)

longitude <- as.numeric(fatores$long)

size <- as.numeric(fatores$LogCS)

### Execute a GPA - Alinha tudo e tira o impacto da dimensionalidade dos dados crus ###

plot(dados) # dados crus
gpa <- gpagen(dados)

link<-read.table("Datasets/tps/link.txt") ## Links entre os landmarks
plot(gpa,link=link)

shp <- two.d.array(gpa$coords) # criando uma matriz para a forma
cov_matrix <- cov(shp)  
print(cov_matrix)

### Quando necessário, busque o espécime médio para desenhar o outline

tps_0 <- readland.tps("Datasets/tps/jaw_18LM.tps") ## dados de todos os espécimes para achar a foto
Y.gpa <- gpagen(tps_0)
findMeanSpec(Y.gpa$coords) # 90

plotOutliers(Y.gpa$coords, inspect.outliers = TRUE) ## verificação

### Dados separados por sexo para análises separadas (existem outros meios)

tps_M <- readland.tps("Datasets/tps/avgloc_M.tps")
M_gpa <- gpagen(tps_M)
Males <- read.csv("Planilhas/avglocsex_M.csv", sep =";")
tps_F <- readland.tps("Datasets/tps/avgloc_F.tps")
F_gpa <- gpagen(tps_F)
Females <- read.csv("Planilhas/avglocsex_F.csv", sep =";")

### Carregar o outline
drawinglandmark<-readland.tps("Datasets/outline2/outline.tps")
outline<-read.table("Datasets/outline2/outline.txt", header=FALSE)
summary(drawinglandmark)

mshape<-mshape(Y.gpa$coords)
Sapajusoutline<-warpRefOutline(file = "Datasets/outline2/outline.txt", drawinglandmark[,,1],mshape)#run dev.off() in case of error message ##set outline configuration

grid.pars <- gridPar(grid.col = "white") ## se quiser só os grids "invisíveis"


GP <- gridPar(n.col.cell = 100, pt.bg = "gray", pt.size = 0.8, tar.pt.bg = "cyan", 
              tar.pt.size = 0.8,tar.out.col = "gray10", tar.out.cex = 0.5, 
              grid.col = "white", grid.lwd = 0.5,txt.pos = 1, txt.col = "steelblue") ## grids personalizados

plotRefToTarget(mshape,Y.gpa$coords[,,90],outline = Sapajusoutline$outline,method="points", gridPars=GP)

dev.off()

############ ORDENAÇÃO & EXPLORAÇÃO ###############

################## PCA #####################

# Exploratória

PCA <- gm.prcomp(gpa$coords)
summary(PCA)


### Basicamente os parâmetros para personalizar os plots

symbols <- c(
  Sapajus_apella = "\u25A0",        # Quadrado 
  Sapajus_cay = "\u25BC",           # Triângulo apontando para baixo
  Sapajus_libidinosus = "\u25B2",   # Triângulo apontando para cima
  Sapajus_nigritus = "\u25CF",      # Círculo 
  Sapajus_robustus = "\u2726",      # Estrela com 6 pontas
  Sapajus_xanthosternos = "\u2736"  # Estrela com 8 pontas
)
symbols <- symbols[as.character(species)]

Fac <- as.factor(fatores$fac) ## Fator de ambientes

levels_fac <- levels(Fac)

cores_fac <- c("green3", "darkgreen", "goldenrod1")
names(cores_fac) <- levels_fac

indv <- 1:nrow(PCA$x)

pca_df <- data.frame(
  PC1 = PCA$x[, "Comp1"],
  PC2 = PCA$x[, "Comp2"],
  PC3 = PCA$x[, "Comp3"],
  Species = species, 
  Individuo = 1:nrow(PCA$x), Fac = Fac, Sex = sex)

species_colors <- c("red","navy","cyan3","saddlebrown","magenta4","black")
names(species_colors) <- levels(pca_df$Species)

pch <- as.numeric(fatores$pch)
uPCH <- unique(pch)
names(uPCH) <- unique(pca_df$Species)
cores_fac <- c(AF = "green3", AM = "darkgreen", SV = "goldenrod1")[as.numeric(as.factor(Fac))]

# Plot PCA
mat<-matrix(c(4,5,0,1,1,2,1,1,3),3) # Dividir a janela do plot
layout(mat, widths=c(1,1,1), heights=c(1,1,0.6))
par(mar=c(4, 4, 1, 1))
xlab <- paste("PC1 (", round(100 * PCA$sdev[1]^2 / sum(PCA$sdev^2), 1), "%)", sep = "")
ylab <- paste("PC2 (", round(100 * PCA$sdev[2]^2 / sum(PCA$sdev^2), 1), "%)", sep = "")

plot(PCA$x[,1],PCA$x[,2],cex=4,pch = 21, bg=cores_fac, xlab=xlab,ylab=ylab, 
     col = cores_fac)
legend("topright",legend= unique(Fac),pch=19,col = unique(cores_fac), cex = 2, pt.cex = 3)

plotRefToTarget(mshape,PCA$shapes$shapes.comp1$min, mag = 2,outline = Sapajusoutline$outline,method="points",gridPars=GP)
plotRefToTarget(mshape,PCA$shapes$shapes.comp1$max, mag = 2,outline = Sapajusoutline$outline,method="points",gridPars=GP)
plotRefToTarget(mshape,PCA$shapes$shapes.comp2$max, mag = 2,outline = Sapajusoutline$outline,method="points",gridPars=GP)
plotRefToTarget(mshape,PCA$shapes$shapes.comp2$min, mag = 2,outline = Sapajusoutline$outline,method="points",gridPars=GP)
par(mfrow=c(1,1)) # janela default  

# Gera elipses de 95%

plot(PCA$x[,1],PCA$x[,2],cex=4,pch = 19,col = cores_fac,xlab=xlab,ylab=ylab)

ordiellipse(PCA$x,group = Fac, kind="sd", conf=0.95, col = unique(cores_fac))

ordihull(PCA$x,group = Fac, col = unique(Fac))

########### Plots complexos com "ggplot2" ########

compute_hull <- function(pca_df) pca_df[chull(pca_df$PC1, pca_df$PC2), ]
hulls <- pca_df %>%
  group_by(Species) %>%
  do(compute_hull(.))

ggplot(pca_df, aes(x = PC1, y = PC2)) +
  geom_point(aes(color = Species, shape = Species), size = 5) +  
  geom_text(aes(label = Individuo), size = 3, vjust = 1.5) +  
  geom_polygon(data = hulls, aes(group = Species, fill = Species, color = Species), alpha = 0.15) +  
  labs(title = "Shape's Principal Components", x = paste("PC1 (", round(100 * PCA$sdev[1]^2 / sum(PCA$sdev^2), 1), "%)", sep = ""), y = paste("PC2 (", round(100 * PCA$sdev[2]^2 / sum(PCA$sdev^2), 1), "%)", sep = ""), color = "Species", shape = "Biomes") +
  theme_minimal() +
  scale_color_manual(values = species_colors, name = "Species") +  
  scale_fill_manual(values = species_colors, name = "Species") +  
  scale_shape_manual(values = uPCH, name = "Species")

compute_hull <- function(pca_df) pca_df[chull(pca_df$PC1, pca_df$PC2), ]
hulls <- pca_df %>%
  group_by(Fac) %>%
  do(compute_hull(.))
ggplot(pca_df, aes(x = PC1, y = PC2)) +
  geom_point(aes(shape = Species, color = Fac), size = 6) +  
  geom_text(aes(label = Individuo), size = 2, vjust = 2) +  
  geom_polygon(data = hulls, aes(x = PC1, y = PC2, group = Fac, fill = Fac, color = Fac), alpha = 0.2) +  
  labs(title = "Shape's Principal Components", 
       x = paste("PC1 (", round(100 * PCA$sdev[1]^2 / sum(PCA$sdev^2), 1), "%)", sep = ""), 
       y = paste("PC2 (", round(100 * PCA$sdev[2]^2 / sum(PCA$sdev^2), 1), "%)", sep = ""), 
       shape = "Species", color = "Fac") +
  theme_minimal() +
  scale_shape_manual(values = uPCH, labels = levels(pca_df$Species), name = "Species") +  
  scale_color_manual(values = cores_fac, name = "Environment") +  
  scale_fill_manual(values = cores_fac, name = "Environment")

compute_hull <- function(pca_df) pca_df[chull(pca_df$PC2, pca_df$PC3), ]
hulls <- pca_df %>%
  group_by(Fac) %>%
  do(compute_hull(.))
ggplot(pca_df, aes(x = PC2, y = PC3)) +
  geom_point(aes(shape = Species, color = Fac), size = 8) +  
  geom_text(aes(label = Individuo), size = 3, vjust = 2) +  
  geom_polygon(data = hulls, aes(x = PC2, y = PC3, group = Fac, fill = Fac, color = Fac), alpha = 0.2) +  
  labs(title = NULL, 
       x = paste("PC2 (", round(100 * PCA$sdev[2]^2 / sum(PCA$sdev^2), 1), "%)", sep = ""), 
       y = paste("PC3 (", round(100 * PCA$sdev[3]^2 / sum(PCA$sdev^2), 1), "%)", sep = ""), 
       shape = "Species", color = "Fac") +
  theme_minimal() +
  scale_shape_manual(values = uPCH, labels = levels(pca_df$Species), name = "Species") +  
  scale_color_manual(values = cores_fac, name = "Environment") +  
  scale_fill_manual(values = cores_fac, name = "Environment")


levels_sex <- levels(sex)

cores_sex <- c("salmon","steelblue")
names(cores_sex) <- levels_sex

compute_hull <- function(pca_df) pca_df[chull(pca_df$PC1, pca_df$PC2), ]
hulls <- pca_df %>%
  group_by(Sex) %>%
  do(compute_hull(.))

ggplot(pca_df, aes(x = PC1, y = PC2)) +
  geom_point(aes(shape = Species, color = Sex), size = 10) +  
  geom_text(aes(label = Individuo), size = 2, vjust = 2) +  
  geom_polygon(data = hulls, aes(x = PC1, y = PC2, group = Sex, fill = Sex, color = Sex), alpha = 0.2) +  
  labs(
    title = "Shape's Principal Components", 
    x = paste("PC1 (", round(100 * PCA$sdev[1]^2 / sum(PCA$sdev^2), 1), "%)", sep = ""), 
    y = paste("PC2 (", round(100 * PCA$sdev[2]^2 / sum(PCA$sdev^2), 1), "%)", sep = ""), 
    shape = "Species", color = "Sex"
  ) +
  theme_minimal() +
  scale_shape_manual(values = uPCH, labels = levels(pca_df$Species), name = "Species") +  
  scale_color_manual(values = cores_sex, name = "Sex") +  # Defina `cores_sexo` com as cores correspondentes aos sexos
  scale_fill_manual(values = cores_sex, name = "Sex")

PC1 <- PCA$x[, 1] # Scores para PC1
PC2 <- PCA$x[, 2] # Scores para PC2
PC3 <- PCA$x[, 3]

preds_comb <- shape.predictor(
  gpa$coords,
  x = cbind(PC1, PC2), # Usar PC1 e PC2 juntos
  Intercept = FALSE,
  pred1 = c(min(PC1), min(PC2)), # Min PC1 + Min PC2
  pred2 = c(max(PC1), min(PC2)), # Max PC1 + Min PC2
  pred3 = c(min(PC1), max(PC2)), # Min PC1 + Max PC2
  pred4 = c(max(PC1), max(PC2)),  # Max PC1 + Max PC2
  pred5 = c(min(PC2), min(PC3)),
  pred6 = c(min(PC2), max(PC3)),
  pred7 = c(max(PC2), min(PC3)),
  pred8 = c(max(PC2), max(PC3))
)

M <- mshape(gpa$coords)

# Visualizar as formas associadas a cada combinação dos extremos das componentes

par(mfrow = c(2, 2)) 

plotRefToTarget(M, preds_comb$pred3, main = "Min PC1 + Max PC2", mag = 2, outline = Sapajusoutline$outline,gridPars = GP,  method = "points")
plotRefToTarget(M, preds_comb$pred4, main = "Max PC1 + Max PC2", mag = 2, outline = Sapajusoutline$outline,gridPars = GP,  method = "points")
plotRefToTarget(M, preds_comb$pred1, main = "Min PC1 + Min PC2", mag = 2, outline = Sapajusoutline$outline,gridPars = GP,  method = "points")
plotRefToTarget(M, preds_comb$pred2, main = "Max PC1 + Min PC2", mag = 2, outline = Sapajusoutline$outline,gridPars = GP,  method = "points")

### Para PC3 
plotRefToTarget(M, preds_comb$pred6, main = "Min PC2 + Max PC3", mag = 2, outline = Sapajusoutline$outline,gridPars = GP,  method = "points")
plotRefToTarget(M, preds_comb$pred8, main = "Max PC2 + Max PC3", mag = 2, outline = Sapajusoutline$outline,gridPars = GP,  method = "points")
plotRefToTarget(M, preds_comb$pred5, main = "Min PC2 + Min PC3", mag = 2, outline = Sapajusoutline$outline,gridPars = GP,  method = "points")
plotRefToTarget(M, preds_comb$pred7, main = "Max PC2 + Min PC3", mag = 2, outline = Sapajusoutline$outline,gridPars = GP,  method = "points")

par(mfrow = c(1, 1))

####### PCA original separada por sexo ###
compute_hull <- function(pca_df) pca_df[chull(pca_df$PC1, pca_df$PC2), ]

hulls <- pca_df %>%
  group_by(Sex, Fac) %>%  # Separar por `Sex` e `Fac`
  do(compute_hull(.)) %>%
  ungroup()

hulls <- hulls %>%
  mutate(Sex = factor(Sex), Fac = factor(Fac))

ggplot(pca_df, aes(x = PC1, y = PC2)) +
  geom_point(aes(shape = Species, color = Fac), size = 6) +  
  geom_text(aes(label = Individuo), size = 2, vjust = 2) +  
  geom_polygon(
    data = hulls,
    aes(x = PC1, y = PC2, group = interaction(Sex, Fac), fill = Fac, color = Fac),
    alpha = 0.2,   # Opacidade para o preenchimento
    linewidth = 1  # Espessura do contorno para destaque
  ) +  
  labs(
    title = "Shape's Principal Components (Separated by Sex)",
    x = paste("PC1 (", round(100 * PCA$sdev[1]^2 / sum(PCA$sdev^2), 1), "%)", sep = ""), 
    y = paste("PC2 (", round(100 * PCA$sdev[2]^2 / sum(PCA$sdev^2), 1), "%)", sep = ""), 
    shape = "Species", color = "Environment", fill = "Environment"
  ) +
  theme_minimal() +
  scale_shape_manual(values = uPCH, labels = levels(pca_df$Species), name = "Species") +  
  scale_color_manual(values = cores_fac, name = "Environment") +  
  scale_fill_manual(values = cores_fac, name = "Environment") +
  facet_wrap(~Sex)  # Separar por sexo


par(mfrow=c(1,2))

scorespc1<-PCA$x[,1] #scores for pc1
scorespc2<-PCA$x[,2] #scores for pc2
scorespc3<-PCA$x[,3] #scores for pc2

### Visualizar as formas relacionadas aos valores extremos de cada PC

predsPC1 <- shape.predictor(gpa$coords, scorespc1, Intercept = TRUE, predmin = min(scorespc1), predmax = max(scorespc1)) #estimating shape configurations based on pc scores
plotRefToTarget(mshape, predsPC1$predmin, mag=2, outline = Sapajusoutline$outline,gridPars = GP,  method = "points") #shape related to small pc scores
plotRefToTarget(mshape, predsPC1$predmax, mag=2, outline = Sapajusoutline$outline,gridPars = GP,  method = "points")#shape related to large pc scores

predsPC2 <- shape.predictor(gpa$coords, scorespc2, Intercept = TRUE, predmin = min(scorespc2), predmax = max(scorespc2)) #estimating shape configurations based on pc scores
plotRefToTarget(mshape, predsPC2$predmin, mag=2, outline = Sapajusoutline$outline,gridPars = GP,  method = "points") #shape related to small pc scores
plotRefToTarget(mshape, predsPC2$predmax, mag=2, outline = Sapajusoutline$outline,gridPars = GP,  method = "points")#shape related to large pc scores

predsPC3<- shape.predictor(gpa$coords, scorespc3, Intercept = TRUE, predmin = min(scorespc3), predmax = max(scorespc3)) #estimating shape configurations based on pc scores
plotRefToTarget(mshape, predsPC3$predmin, mag=2, outline = Sapajusoutline$outline,gridPars = GP,  method = "points") #shape related to small pc scores
plotRefToTarget(mshape, predsPC3$predmax, mag=2, outline = Sapajusoutline$outline,gridPars = GP,  method = "points")#shape related to large pc scores

par(mfrow=c(1,1))

dev.off()


####### CORRELAÇÃO DA PCA ##########

### Pegar as 3 pcs e correlacionar ao tamanho ###

PC1 <- pca_df$PC1
PC2 <- pca_df$PC2
PC3 <- pca_df$PC3

correlacao <- cor(cbind(PC1, PC2, PC3, size))
print(correlacao)

correlacao <- cor(cbind(PC1, PC2, PC3, biome))
print(correlacao)

summary(PCA)
main_components <- PCA$x

correlacao_bio <- cor(cbind(main_components),as.numeric(biome))
print(correlacao)

correlacao <- cor(cbind(main_components),size)
print(correlacao)


pca_data <- data.frame(pca_scores = main_components, Size = fatores$LogCS)

size_correlations <- correlacao
bio_cor <- correlacao_bio

names(bio_cor) <- c("Comp1", "Comp2", "Comp3", "Comp4", "Comp5", "Comp6", 
                    "Comp7", "Comp8", "Comp9", "Comp10", "Comp11", "Comp12",
                    "Comp13", "Comp14", "Comp15", "Comp16", "Comp17", "Comp18", 
                    "Comp19", "Comp20", "Comp21", "Comp22", "Comp23", "Comp24", 
                    "Comp25", "Comp26", "Comp27", "Comp28", "Comp29", "Comp30", 
                    "Comp31", "Comp32")
names(size_correlations) <- c("Comp1", "Comp2", "Comp3", "Comp4", "Comp5", "Comp6", 
                              "Comp7", "Comp8", "Comp9", "Comp10", "Comp11", "Comp12",
                              "Comp13", "Comp14", "Comp15", "Comp16", "Comp17", "Comp18", 
                              "Comp19", "Comp20", "Comp21", "Comp22", "Comp23", "Comp24", 
                              "Comp25", "Comp26", "Comp27", "Comp28", "Comp29", "Comp30", 
                              "Comp31", "Comp32")


correlation_df <- data.frame(PC = names(bio_cor), Bio = bio_cor, Size = size_correlations)

correlation_df$PC <- factor(correlation_df$PC, 
                            levels = c("Comp1", "Comp2", "Comp3", "Comp4", "Comp5", 
                                       "Comp6", "Comp7", "Comp8", "Comp9", "Comp10", 
                                       "Comp11", "Comp12", "Comp13", "Comp14", 
                                       "Comp15", "Comp16", "Comp17", "Comp18", 
                                       "Comp19", "Comp20", "Comp21", "Comp22", 
                                       "Comp23", "Comp24", "Comp25", "Comp26", 
                                       "Comp27", "Comp28", "Comp29", "Comp30", 
                                       "Comp31", "Comp32"))

ggplot(correlation_df, aes(x = PC, y = Bio, fill = Bio)) +
  geom_bar(stat = "identity") +
  scale_fill_gradientn(colors = c("red", "gray80", "steelblue"), 
                       values = scales::rescale(c(-0.5, 0, 0.5)),
                       name = "Cor") +
  theme_classic() +
  labs(title = "Correlogram: Main Components vs. Biome",
       x = "Main Components",
       y = "Correlation") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

correlation_df$PC <- factor(correlation_df$PC, 
                            levels = c("Comp1", "Comp2", "Comp3", "Comp4", "Comp5", 
                                       "Comp6", "Comp7", "Comp8", "Comp9", "Comp10", 
                                       "Comp11", "Comp12", "Comp13", "Comp14", 
                                       "Comp15", "Comp16", "Comp17", "Comp18", 
                                       "Comp19", "Comp20", "Comp21", "Comp22", 
                                       "Comp23", "Comp24", "Comp25", "Comp26", 
                                       "Comp27", "Comp28", "Comp29", "Comp30", 
                                       "Comp31", "Comp32"))

ggplot(correlation_df, aes(x = PC, y = Size, fill = Size)) +
  geom_bar(stat = "identity") +
  scale_fill_gradientn(colors = c("red", "gray80", "steelblue"), 
                       values = scales::rescale(c(-0.5, 0, 0.5)),
                       name = "Cor") +
  theme_classic() +
  labs(title = "Correlogram: Main Components vs. Size",
       x = "Main Components",
       y = "Correlation") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


### isso nos diz se é a correlação é verossímil

resultado_teste <- cor.test(PC1, size, method = "pearson")
print(resultado_teste)

par(mfrow=c(3,1))

plot(PC1, size, main = "Relação entre PC1 e LogCS", xlab = "PC1", ylab = "LogCS",pch = 21, bg = cores_fac, col = cores_fac, cex = 2,)
abline(lm(size ~ PC1), col = "red")

resultado_teste <- cor.test(PC2, size, method = "pearson")
print(resultado_teste)

plot(PC2, size, main = "Relação entre PC2 e LogCS", xlab = "PC2", ylab = "LogCS",pch = 21, bg = cores_fac, col = cores_fac, cex = 2,)
abline(lm(size ~ PC2), col = "red")

resultado_teste <- cor.test(PC3, size, method = "pearson")
print(resultado_teste)

plot(PC3, size, main = "Relação entre PC3 e LogCS", xlab = "PC3", ylab = "LogCS",pch = 21, bg = cores_fac, col = cores_fac, cex = 2,)
abline(lm(size ~ PC3), col = "blue")
par(mfrow=c(1,1))


####################################### CVA ############################################


if (require(shapes)) {
  alldat<-gpa$coords
  # create factors
  groups<-as.factor(fatores$fac)
  # perform CVA and test Mahalanobis distance
  # between groups with permutation test by 100 rounds)            
  cvall<-CVA(alldat,groups,rounds=10000)     
  ## visualize a shape change from score -5 to 5:
  cvvis5 <- 5*matrix(cvall$CVvis[,1],nrow(cvall$Grandm),ncol(cvall$Grandm))+cvall$Grandm
  cvvisNeg5 <- -5*matrix(cvall$CVvis[,1],nrow(cvall$Grandm),ncol(cvall$Grandm))+cvall$Grandm
  plot(cvvis5,asp=1,pch = 21, bg = "black")
  points(cvvisNeg5,col=2, pch = 25, cex = 2, bg = "red")
  for (i in 1:nrow(cvvisNeg5))
    lines(rbind(cvvis5[i,],cvvisNeg5[i,]))
}

cva.1=CVA(gpa$coords, groups=Fac)

## get the typicality probabilities and resulting classifications - tagging
## all specimens with a probability of < 0.01 as outliers (assigned to no class)

typprobs <- typprobClass(cva.1$CVscores,groups=Fac)
print(typprobs)

## visualize the CV scores by their groups estimated from (cross-validated)
## typicality probabilities:
ffinCV <- as.factor(typprobs$groupaffinCV)

levels_ffinCV <- unique(typprobs$groupaffinCV)

cores_CVA <- c(AF = "green3", AM = "darkgreen", SV = "goldenrod1")

if (require(car)) {
  scatterplot(cva.1$CVscores[,1],cva.1$CVscores[,2],groups=typprobs$groupaffinCV, pch = c(21,22,25) ,col = unique(cores_CVA), cex = 4, bg = unique(cores_CVA))
}
text(cva.1$CVscores, as.character(ffinCV), col = "black", cex = 0.5)

cva.1$CVscores[,1]

# plot the CVA
plot(cva.1$CVscores, pch = 21, bg = cores_fac, col= cores_fac, asp=1, cex = 3,
     xlab=paste("1st canonical axis", paste(round(cva.1$Var[1,2],1),"%")),
     ylab=paste("2nd canonical axis", paste(round(cva.1$Var[2,2],1),"%")))

text(cva.1$CVscores, as.character(Fac), col = "black", cex=.7)

# add chull (merge groups)
for(jj in 1:length(levels(Fac))){
  ii=levels(Fac)[jj]
  kk=chull(cva.1$CVscores[Fac==ii,1:2])
  lines(cva.1$CVscores[Fac==ii,1][c(kk, kk[1])],
        cva.1$CVscores[Fac==ii,2][c(kk, kk[1])], col = (1:3))
}

# add 80% ellipses
if (require(car)) {
  for(ii in 1:length(levels(Fac))){
    dataEllipse(cva.1$CVscores[Fac==levels(Fac)[ii],1],
                cva.1$CVscores[Fac==levels(Fac)[ii],2], 
                add=TRUE,levels=.95, col=unique(cores_CVA)[ii])}
}
# histogram per group
if (require(lattice)) {
  histogram(~cva.1$CVscores[,1]|Fac,
            layout=c(1,length(levels(Fac))),
            xlab=paste("1st canonical axis", paste(round(cva.1$Var[1,2],1),"%")))
  histogram(~cva.1$CVscores[,2]|Fac, layout=c(1,length(levels(Fac))),
            xlab=paste("2nd canonical axis", paste(round(cva.1$Var[2,2],1),"%")))
} 
# plot Mahalahobis
dendroS=hclust(cva.1$Dist$GroupdistMaha)
dendroS$labels=levels(Fac)
par(mar=c(4,4.5,1,1))
dendroS=as.dendrogram(dendroS)
plot(dendroS, main='',sub='', xlab="Ecoregions",
     ylab='Mahalahobis distance')

fac <- Fac
col.group<-rainbow(length(levels(Fac))) # criar vetor de cores para grupos
names(col.group)<-levels(Fac)
col.group<-col.group[match(Fac,names(col.group))]

PCA

cva<-CVA(PCA$x[,1:32],groups = Fac,cv=TRUE) # Jackknife Cross-validation
cva
plot(cva$CVscores, col = cores_fac, cex = 3)

cva<-CVA(gpa$coords,groups = Fac,cv=FALSE) # CVA validation 
cva
plot(cva$CVscores, col = cores_fac, pch = 21, bg = cores_fac, cex = 3)

# LDA (Linear Discriminant Analysis) #2 alternativa com apenas dois grupos

cva2<-lda(PCA$x[,1:32],biome)
plot(cva2)

cva2<-lda(PCA$x[,1:32],biome,CV=T) #LDA com Jackknife cross validation
plot(cva2$posterior, col = unique(cores_fac), pch = 21, bg = unique(cores_fac), cex = 1)

tab<-table(biome,cva2$class) 
lda.p<-diag(tab)/summary(biome)*100
lda.p # proporção de classificação correta para cada grupo

# CVA (Canonical Variate Analysis)
cva$Var # Porcentagem de explica??o dos eixos da CVA
plot(cva$CVscores,col = col.group,bg=col.group,pch=21,cex=2,xlab=paste("CV1  (",paste(round(cva$Var[1,2],1),"%)"),sep=""),ylab=paste("CV2  (",paste(round(cva$Var[2,2],1),"%)"),sep=""))
legend("topleft",legend=unique(fac),pch=19,col=unique(col.group)) 

# gera elipses de confiançaa 95%
ordiellipse(cva$CVscores,group=fac, kind="sd", conf=0.95) 
ordihull(cva$CVscores,group=fac)

# Visualização da forma nos eixos da CVA
cv1max<-max(cva$CVscores[,1])*matrix(cva$CVvis[,1],nrow(cva$Grandm),ncol(cva$Grandm))+cva$Grandm
cv1min<-min(cva$CVscores[,1])*matrix(cva$CVvis[,1],nrow(cva$Grandm),ncol(cva$Grandm))+cva$Grandm
cv2max<-max(cva$CVscores[,2])*matrix(cva$CVvis[,2],nrow(cva$Grandm),ncol(cva$Grandm))+cva$Grandm
cv2min<-min(cva$CVscores[,2])*matrix(cva$CVvis[,2],nrow(cva$Grandm),ncol(cva$Grandm))+cva$Grandm

plotRefToTarget(cv1min,cv1max,outline = Sapajusoutline$outline,gridPars = GP,  method = "points", mag = 1.5) 
plotRefToTarget(mshape,cv1min,outline = Sapajusoutline$outline,gridPars = GP,  method = "points", mag = 1.5)  
plotRefToTarget(mshape,cv2max,outline = Sapajusoutline$outline,gridPars = GP,  method = "points", mag = 1.5) 
plotRefToTarget(mshape,cv2min,outline = Sapajusoutline$outline,gridPars = GP,  method = "points", mag = 1.5) 

# CVA Estilizada
mat<-matrix(c(4,5,0,1,1,2,1,1,3),3) # Dividir a janela do plot
layout(mat, widths=c(1,1,1), heights=c(1,1,0.6))
par(mar=c(4, 4, 1, 1))
plot(cva$CVscores,bg=col.group,pch=21,cex=2,xlab=paste("CV1  (",paste(round(cva$Var[1,2],1),"%)"),sep=""),ylab=paste("CV2  (",paste(round(cva$Var[2,2],1),"%)"),sep=""))
legend(+3,+3,legend=unique(Fac),pch=19,col=unique(col.group))
ordiellipse(cva$CVscores,group=fac, kind="sd", conf=0.95) 

plotRefToTarget(mshape,cv1min,outline = Sapajusoutline$outline,gridPars = GP,  method = "points", mag = 2) 
plotRefToTarget(mshape,cv1max,outline = Sapajusoutline$outline,gridPars = GP,  method = "points", mag = 2)  
plotRefToTarget(mshape,cv2max,outline = Sapajusoutline$outline,gridPars = GP,  method = "points", mag = 2) 
plotRefToTarget(mshape,cv2min,outline = Sapajusoutline$outline,gridPars = GP,  method = "points", mag = 2) 
par(mfrow=c(1,1))

### CVA Corrigida ###
layout(matrix(c(4,5,0,1,1,2,1,1,3), 3), widths=c(1,1,1), heights=c(1,1,0.6))
par(mar=c(4, 4, 1, 1))

plot(
  cva$CVscores,
  bg = cores_fac[as.character(Fac)],
  pch = 21,
  cex = 4,
  xlab = paste("CV1  (", round(cva$Var[1,2], 1), "%)", sep = ""),
  ylab = paste("CV2  (", round(cva$Var[2,2], 1), "%)", sep = ""),
  col = cores_fac[as.character(Fac)],
  lwd = 2
)

legend(+2.5,+3,legend=unique(Fac),pch=19,col=unique(cores_fac),title = "Ecoregions",
       cex = 1.2,
       pt.cex = 3,)

# add 95% ellipses
if (require(car)) {
  for(ii in 1:length(levels(Fac))){
    dataEllipse(cva.1$CVscores[Fac==levels(Fac)[ii],1],
                cva.1$CVscores[Fac==levels(Fac)[ii],2], 
                add=TRUE,levels=.95, col=unique(cores_CVA)[ii])}
}

#ordiellipse(cva$CVscores,group = fac,kind = "sd",conf = 0.95,col = legend_colors,lwd = 2)

text(
  cva$CVscores,
  labels = rownames(cva$CVscores),
  cex = 0.7,
  col = "black"
)

plotRefToTarget(mshape,cv1min,outline = Sapajusoutline$outline,gridPars = GP,  method = "points", mag = 2) 
plotRefToTarget(mshape,cv1max,outline = Sapajusoutline$outline,gridPars = GP,  method = "points", mag = 2)  
plotRefToTarget(mshape,cv2max,outline = Sapajusoutline$outline,gridPars = GP,  method = "points", mag = 2) 
plotRefToTarget(mshape,cv2min,outline = Sapajusoutline$outline,gridPars = GP,  method = "points", mag = 2) 
par(mfrow=c(1,1))

# Comparação das formas médias dos grupos entre si # 1 = cinza, 2 = preto # Exemplos:
par(mfrow=c(2,2))

plotRefToTarget(cva$groupmeans[,,"SV"],cva$groupmeans[,,"AM"],outline = Sapajusoutline$outline,gridPars = GP,  method = "points", mag = 3)
plotRefToTarget(cva$groupmeans[,,"AM"],cva$groupmeans[,,"SV"],outline = Sapajusoutline$outline,gridPars = GP,  method = "points", mag = 3)
plotRefToTarget(cva$groupmeans[,,"SV"],cva$groupmeans[,,"AF"],outline = Sapajusoutline$outline,gridPars = GP,  method = "points", mag = 3)
plotRefToTarget(cva$groupmeans[,,"AF"],cva$groupmeans[,,"SV"],outline = Sapajusoutline$outline,gridPars = GP,  method = "points", mag = 3)


# Plotar o Scree Plot (gráfico que mostra a dispersão da variância acumulada pela PCA)
par(mfrow=c(1,2))

scree_data <- data.frame(Principal_Component = 1:length(PCA$sdev), Variance_Explained = PCA$sdev^2 / sum(PCA$sdev^2))

plot(scree_data$Principal_Component, scree_data$Variance_Explained, type = "b", pch = 19, col = "blue", xlab = "Componente Principal", ylab = "Porcentagem de Variância Explicada", main = "Original")

scree_data <- data.frame(Principal_Component = 1:length(pca_res$sdev), Variance_Explained = pca_res$sdev^2 / sum(pca_res$sdev^2))

plot(scree_data$Principal_Component, scree_data$Variance_Explained, type = "b", pch = 19, col = "red", xlab = "Componente Principal", ylab = "Porcentagem de Variância Explicada", main = "Resíduos Alométricos")

par(mfrow=c(1,1))

dev.off()


############################### ANÁLISES DE VARIÂNCIA ####################################


####################### Analises de Variância Univariadas Qualitativas ###################
df <- data.frame(Species = species, Latitude = latitude, Biome = biome, Envi = Fac, Size = size, Sex = sex)

aov_res <- aov(Size~Species, data = df)
summary(aov_res)
TukeyHSD(aov_res)

aov_res <- aov(Size~Sex, data = df)
summary(aov_res)

aov_res <- aov(Size~Biome, data = df)
summary(aov_res)

aov_res <- aov(Size~Envi, data = df)
summary(aov_res)
TukeyHSD(aov_res)

### Multivariadas (para forma)
## Criar um geomorph dataframe ##

gdf <- geomorph.data.frame(Shape = gpa$coords, Size = gpa$Csize,Species = species, Latitude = latitude, Biome = biome, Fac = fatores$fac, Sex = sex)

fit.size <- procD.lm(gpa$coords~log(gpa$Csize), data = gdf) 
summary(fit.size)

fitl <- procD.lm(Shape~Size*Latitude, data = gdf, iter = 999) ### A forma não está associada à variações latitudinais
summary(fitl)

fitl1 <- procD.lm(Shape~log(Size)*Biome, data = gdf, iter = 999) 
fitl2 <- procD.lm(Shape~log(Size)*Species, data = gdf, iter = 999) ### Dimorfismo sexual 
fitl3 <- procD.lm(Shape~log(Size)*Fac, data = gdf, iter = 999)
fitl4 <- procD.lm(Shape~log(Size)*Sex, data = gdf, iter = 999)
fitl5 <- procD.lm(Shape~log(Size)*Fac*Sex, data = gdf, iter = 999) ### Interação CS : Biome (?)

summary(fitl1)
summary(fitl2)
summary(fitl3)
summary(fitl4)
summary(fitl5)

## Não tendo interações com o fator Sexo, as análises podem ser feitas com ambos juntos

########################## Phylo Analisys #########################

# Métodos Filogenéticos Comparativos
# Carregar arquivo tps da mandíbula
tps<-readland.tps("Datasets/tps/avglocsex.tps",specID = "ID", readcurves = FALSE)
gpa<-gpagen(tps)
shape<-gpa$coords
size<-gpa$Csize
ref.mand<-mshape(shape)

shape.2d<-two.d.array(shape)
shape.2d <- as.matrix(shape.2d,ncol(32))

# Carregar classificadores a partir de lista externa
plan<-read.csv("Planilhas/avglocsex.csv",sep=";")
fac <- as.factor(plan$fac)
species <- as.factor(plan$sp)
names <- plan$sp_ives
sexis <- as.factor(plan$sex)

names(fac) <- names
names(species) <- names
names(sexis) <- names

###### Minha análise filogenética pra vários spécimes #####
# Carregar árvore filogen?tica no formato nexus
tree<-read.tree("Datasets/trees/AVGLOCSEX_Lima.tre")
plot(tree)
tree<-compute.brlen(tree,1) # definir comprimento dos ramos = 1

# Plotar Filogenia no espaço de forma (phylomorphoespace)

tree$tip.label

dimnames(shape.2d)[[1]] <- paste(names)

pms <- gm.prcomp(shape.2d,phy=tree)
pms
plot(pms,phylo=TRUE)

# PCA filogenética
pPCA_BM<-phyl.pca(tree,pms$x[,1:16],method="BM")
pPCA_BM
phylomorphospace(tree,pPCA_BM$S[,1:2],label="horizontal")

pPCA_L<-phyl.pca(tree,pms$x[,1:16],method="lambda")
pPCA_L
phylomorphospace(tree,pPCA_L$S[,1:2],label="horizontal")

GP <- gridPar(n.col.cell = 100, pt.bg = "gray", pt.size = 0.8, tar.pt.bg = "cyan", 
              tar.pt.size = 0.8,tar.out.col = "gray10", tar.out.cex = 0.5, 
              grid.col = "white", grid.lwd = 0.5,txt.pos = 1, txt.col = "steelblue") ## grids personalizados

plot(tree)
nodelabels()

sinal.k<-physignal(shape.2d,tree,iter=999)
sinal.k
physignal.z(shape.2d, tree)

# Sinal Filogenético para tamanho do crânio
names(size)<-names
sinal.k.size<-physignal(size,tree,iter=999)
sinal.k.size

# Valor observado de K
K_obs <- sinal.k.size$phy.signal

# Permutações (nula)
K_null <- sinal.k.size$random.K

# Effect size como Z-score
effect_size <- (K_obs - mean(K_null)) / sd(K_null)
effect_size
library(phytools)

# lnCS é um vetor nomeado com os nomes das espécies
lambda_result <- phylosig(tree, size, method = "lambda", test = TRUE)
print(lambda_result)

### com esses códigos da pra fazer análises entre espécies (médias), mas como o meu foco era a análise entre ambientes considerando toda a variação, não foquei nisso

fac <- as.factor(plan$fac)
names(fac)<-names

fit<-procD.lm(shape ~ fac, iter = 999)
summary(fit)

# Execute o modelo (forma como resposta, grupo como preditor)

fit <- procD.pgls(shape.2d ~ Size, phy = tree, data=gdf,iter = 999) 
summary(fit)

# as politomias entre os grupos são muito próximas, o controle filogenético pode estar "apagando" qualquer diferença na forma entre eles.
fit <- procD.pgls(shape.2d ~ Size*Fac, phy = tree, data=gdf, iter = 999) 
summary(fit)

dimnames(gdf$Shape)[[3]] <- paste (names)

fitl1 <- procD.pgls(Shape~log(Size)*Biome, data = gdf, , phy = tree, iter = 999) 
fitl2 <- procD.pgls(Shape~log(Size)*Species, data = gdf, , phy = tree, iter = 999) 
fitl3 <- procD.pgls(Shape~log(Size)*Fac, data = gdf, , phy = tree, iter = 999)
fitl4 <- procD.pgls(Shape~log(Size)*Sex, data = gdf, , phy = tree, iter = 999)
fitl5 <- procD.pgls(Shape~log(Size)*Sex*Biome, data = gdf, , phy = tree, iter = 999)

summary(fitl1)
summary(fitl2)
summary(fitl3)
summary(fitl4)
summary(fitl5)

### ANOVA Filogenética - Essa parte do script não está muito correta
dim(shape)  # Verifique o número de linhas e colunas de shape
length(names)  # Verifique o comprimento de names
formula <- (size~Fac)
names(Fac) <- names

library(nlme)
library(ape)

# Matriz de covariância filogenética
C <- vcv(tree, corr = TRUE)

# Ajustar modelo GLS multivariado com essa matriz
model <- gls(cbind(PC) ~ Fac, correlation = corSymm(value = C[lower.tri(C)], fixed = TRUE))
summary(model)


###### PGLS de acordo com o R. Maestri
########################################################################
# Calcular forma média por grupos
shape.2d<-two.d.array(gpa$coords)
shape.2d.means<-rowsum(shape.2d,species)/as.vector(table(species))
shape.means<-arrayspecs(shape.2d.means,dim(gpa$coords)[1],dim(gpa$coords)[2]) #arrayspecs transforma matriz em tps
shape.means

size.means<-rowsum(size,species)/as.vector(table(species))
size.means <- log(size.means)

# Carregar árvore filogen?tica no formato nexus
tree<-read.tree("Datasets/trees/MCC_Lima_calibrated.tre")
plot(tree)
tree<- drop.tip(tree,c("Sapajus_macrocephalus","Sapajus_flavius" ,"Cebus_albifrons", "Cebus_capucinus"))
tree<-compute.brlen(tree,1) # definir comprimento dos ramos = 1

# Plotar Filogenia no espaço de forma (phylomorphoespace)
# pms<-plotGMPhyloMorphoSpace(tree,shape.means,plot.param=list(t.bg="green",txt.col="black",n.bg="black",n.cex=1,lwd=2,l.col="blue")) #não funciona mais
pms <- gm.prcomp(shape.means,phy=tree)
pms
plot(pms,phylo=TRUE)

# 3d
pca<-gm.prcomp(shape.means)
phylomorphospace(tree,pms$x[,1:2])
phylomorphospace3d(tree,pms$x[,1:3])

# PCA filogenética
pPCA_BM<-phyl.pca(tree,pms$x[,1:5],method="BM")
pPCA_BM
phylomorphospace(tree,pPCA_BM$S[,1:2],label="horizontal")

pPCA_L<-phyl.pca(tree,pca$x[,1:5],method="lambda")
pPCA_L
phylomorphospace(tree,pPCA_L$S[,1:2],label="horizontal")

GP <- gridPar(n.col.cell = 100, pt.bg = "gray", pt.size = 0.8, tar.pt.bg = "cyan", 
              tar.pt.size = 0.8,tar.out.col = "gray10", tar.out.cex = 0.5, 
              grid.col = "white", grid.lwd = 0.5,txt.pos = 1, txt.col = "steelblue") ## grids personalizados

plot(tree)
nodelabels()

# Visualização das formas ancestrais # n1 = raiz # 1st forma em cinza, 2nd em preto
ancestral.shapes<-arrayspecs(pms$ancestors,18,2)
ancestral.shapes
mshp <- mshape(shape.means)
par(mfrow = c(1, 1))
plotRefToTarget(mshp,ancestral.shapes[,,1],method="points", outline = Sapajusoutline$outline, gridPars=GP, mag = 4)
plotRefToTarget(ancestral.shapes[,,1],ancestral.shapes[,,2],method="points", outline = Sapajusoutline$outline, gridPars=GP, mag = 4)
plotRefToTarget(ancestral.shapes[,,1],ancestral.shapes[,,3],method="points", outline = Sapajusoutline$outline, gridPars=GP, mag = 4)
plotRefToTarget(ancestral.shapes[,,1],ancestral.shapes[,,4],method="points", outline = Sapajusoutline$outline, gridPars=GP, mag = 4)
plotRefToTarget(ancestral.shapes[,,1],ancestral.shapes[,,5],method="points", outline = Sapajusoutline$outline, gridPars=GP, mag = 4)

par(mfrow = c(2, 3))
plotRefToTarget(ancestral.shapes[,,1],shape.means[,,"Sapajus_xanthosternos"],method="points", outline = Sapajusoutline$outline, gridPars=GP, mag = 3)
plotRefToTarget(ancestral.shapes[,,1],shape.means[,,"Sapajus_robustus"],method="points", outline = Sapajusoutline$outline, gridPars=GP, mag = 3)
plotRefToTarget(ancestral.shapes[,,1],shape.means[,,"Sapajus_nigritus"],method="points", outline = Sapajusoutline$outline, gridPars=GP, mag = 3)
plotRefToTarget(ancestral.shapes[,,1],shape.means[,,"Sapajus_libidinosus"],method="points", outline = Sapajusoutline$outline, gridPars=GP, mag = 3)
plotRefToTarget(ancestral.shapes[,,1],shape.means[,,"Sapajus_cay"],method="points", outline = Sapajusoutline$outline, gridPars=GP, mag = 3)
plotRefToTarget(ancestral.shapes[,,1],shape.means[,,"Sapajus_apella"],method="points", outline = Sapajusoutline$outline, gridPars=GP, mag = 3)

# Sinal Filogenético para forma do crânio
sinal.k<-physignal(shape.means,tree,iter=999)
sinal.k
physignal.z(shape.means,tree)

# Sinal Filogenético para tamanho do crânio
size.means<-rowsum(size,species)/as.vector(table(species)) # m?dia por esp?cie
sinal.k.size<-physignal(size.means,tree,iter=999)
sinal.k.size

# Visualizar tamanho na filogenia
size.means1<-as.vector(size.means)
names(size.means1)=rownames(size.means)
tree1<-compute.brlen(tree,method="Grafen") # make ultrametric
# cores
contMap(tree1,size.means1)

# Phylogenetic Generalized Least Squares (PGLS) 

# Variáveis hipotéticas para usar como exemplo

nomes<-levels(as.factor(plan$sp))
temp<-as.numeric(plan$bio1)
prec<-as.numeric(plan$bio12)
bio.means <- c("Sapajus_xanthosternos"="Forest","Sapajus_robustus"="Forest","Sapajus_nigritus"="Forest",
               "Sapajus_libidinosus"="Savanna","Sapajus_cay"="Savanna","Sapajus_apella"="Forest")
temp.means<-scale(rowsum(temp,species)/as.vector(table(species)))
prec.means<-scale(rowsum(prec,species)/as.vector(table(species)))
spec <- levels(as.factor(plan$sp))
names(temp.means)<-nomes
names(prec.means)<-nomes
names(bio.means)<-nomes
names(spec) <- nomes

# PGLS
gdf.m<-geomorph.data.frame(Shape=shape.means,Size=size.means,Species = spec,Temp=temp.means,Prec=prec.means,Biome=bio.means,tree=tree)
dimnames(gdf.m$Shape)[[3]] <- paste (nomes)

fitl1 <- procD.lm(Shape~Size*Biome, data = gdf.m, iter = 999) 
fitl2 <- procD.lm(Shape~Species, data = gdf.m, iter = 999) 

summary(fitl1)
summary(fitl2)

fitl1 <- procD.pgls(Shape~log(Size)*Biome, data = gdf.m, phy = tree, iter = 999) 
fitl2 <- procD.pgls(Shape~log(Size)*Species, data = gdf.m, phy = tree, iter = 999) 

summary(fitl1)
summary(fitl2)

p<-procD.lm(shape.means~size.means*prec,data=gdf.m)
summary(p)

pgls.shape<-procD.pgls(shape.means~size.means*prec,tree,data=gdf.m,iter=999) #inclui a árvore filogenética #procrustes
summary(pgls.shape)

# Plot
p.plot<-
  plot(pgls.shape,type="regression",predictor=as.numeric(prec.means),reg.type="RegScore",xlab="Precipitation")
reg.score<-as.vector(p.plot$RegScore)
names(reg.score)=rownames(p.plot$RegScore)
phylomorphospace(tree,cbind(prec.means,reg.score))

# Variação de forma associada
preds <- shape.predictor(shape.means, x = prec.means, Intercept = TRUE,
                         predmin = min(prec.means),
                         predmax = max(prec.means))

plotRefToTarget(preds$predmin,preds$predmax,method = "points", outline = Sapajusoutline$outline, gridPars=GP, mag = 2)

# PGLS size
pgls.size<-procD.pgls(size.means~temp.means,tree,data=gdf,iter=999)
summary(pgls.size)
plot(size.means~temp.means)

pgls.biome <- procD.pgls(shape.means~bio.means,tree,data=gdf,iter=999)
summary(pgls.biome)

### ANOVA Filogenética - Essa parte do script não está muito correta, não funciona
sz <- as.matrix(size.means, ncol(1))
grp <- as.factor(bio.means)
names(grp)=nomes
names(sz)=nomes

x = aov.phylo(sz~grp,tree, nsim=999, test="Pillai")


################# Modelos de regressão linear simples automatizado ###########################

variables <- c("lat", "long", "bio1","bio2", "bio3","bio4","bio5","bio6","bio7","bio8", "bio9","bio10","bio11", "bio12","bio13","bio14", "bio15","bio16","bio17", "bio18", "bio19", 
               "npp", "humid", "soilmoist")
# Lista para armazenar os resultados
model_results <- list()

names(fatores)
fatores[,c(4,5,11:32)] <- scale(fatores[,c(4,5,11:32)])

# Loop para criar um modelo linear para cada variável
for (var in variables) {
  # Montando a fórmula para o modelo linear
  formula <- as.formula(paste("CS ~", var))
  
  # Ajustando o modelo linear
  model <- lm(formula, data = fatores)  # 'locations' é o dataframe com seus dados
  
  # Armazenando os resultados na lista
  model_results[[var]] <- summary(model)
}

# Exibindo os resultados
model_results
r_squared <- sapply(model_results, function(model) model$r.squared)
p_values <- sapply(model_results, function(model) {
  coef_table <- model$coefficients
  coef_table[2, "Pr(>|t|)"]
})
beta_coefficients <- sapply(model_results, function(model) {
  coef_table <- model$coefficients
  coef_table[2, "Estimate"]
})

plot_data <- data.frame(
  Variable = names(beta_coefficients),
  Beta = beta_coefficients,
  P_Value = p_values,
  R_Squared = r_squared
)

plot_data$Significant <- plot_data$P_Value < 0.05

plot_data <- plot_data %>%
  mutate(Significant = P_Value < 0.05)

ggplot(plot_data, aes(x = Beta, y = reorder(Variable, Beta), color = Significant)) +
  geom_point(size = 4) +  # Adiciona os pontos para os coeficientes Beta
  geom_errorbarh(aes(xmin = Beta - 1.96 * R_Squared, xmax = Beta + 1.96 * R_Squared), height = 0.2) +  # Intervalos de confiança
  labs(
    title = "Comparação de Coeficientes Beta com Intervalos de Confiança",
    subtitle = "Intervalos de confiança com coeficientes beta das variáveis preditoras",
    x = "Coeficiente Beta",
    y = "Variáveis",
    caption = "Fonte: Autor"
  ) +
  theme_minimal() +  # Usando um tema minimalista
  scale_color_manual(values = c("TRUE" = "deepskyblue3", "FALSE" = "rosybrown")) +  # Cores para modelos significativos e não significativos
  theme(
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),  # Centralizando e estilizando o título
    plot.subtitle = element_text(hjust = 0.5, size = 12, face = "italic"),  # Subtítulo em itálico
    axis.title = element_text(size = 12),  # Tamanho do título dos eixos
    axis.text = element_text(size = 10),   # Tamanho do texto dos eixos
    plot.caption = element_text(size = 10, face = "italic", hjust = 1)  # Estilizando a legenda
  ) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black") +  # Linha vertical no zero
  coord_cartesian(xlim = c(-0.555, 0.55)) +  # Definindo um limite fixo para o eixo X, centrando o zero
  theme(legend.position = "none")  # Remover a legenda se não for necessária

### Shapiro - Wilk - Aplicar onde necessário conferir a homogeneidade/normalidade das variáveis (geralmente feito nos resíduos de modelos)

Size <- log(size)
head (Size)

shapiro.test(Size)
summary(Size)

hist(Size)
qqnorm(Size) 
qqline(Size)

dev.off()

# Violin plot

ggplot(df, aes(x = Envi, y = Size, fill = Envi)) +
  geom_violin(trim = FALSE) +
  geom_boxplot(width = 0.1, fill = "white", color = "black") +
  scale_fill_manual(values = cores_fac, labels = c("Amazon","Atlantic Forest", "Savanna")) +  # Usando as cores definidas acima
  labs(x = "Biome", y = "LogCS", title = "Violin Plot of Size by Biome") +
  theme_minimal()

ggplot(df, aes(x = Biome, y = Size, fill = Biome)) +
  geom_violin(trim = FALSE) +
  geom_boxplot(width = 0.1, fill = "white", color = "black") +
  scale_fill_manual(values = c("darkgreen", "goldenrod1"), labels = c("Forest", "Savanna")) +  # Usando as cores definidas acima
  labs(x = "Biome", y = "LogCS", title = "Violin Plot of Size by Biome") +
  theme_minimal()

# Ajuste para Regressão Linear

modelobio <- lm(Size ~ Biome, data = df) # Adjusted R-squared:  0.1412
summary(modelobio)

modelo <- lm(Size ~ Species, data = df) #S. cay e S. lib aparecem com maiores coeficientes de redução
summary(modelo)

modelosex <- lm(Size ~ sex, data = df) #Machos maiores
summary(modelo)

model <- lm (log(CS) ~ lat, data = fatores)
summary(model)

main_components <- PCA$x
biome_pca <- cbind(PCs = main_components, Biome = biome, Sex = sex, Species = as.factor(fatores$sp), Size = log(gpa$Csize))

head(biome_pca)

biome_pca <- as.data.frame(biome_pca)

## Gráfico de Efeito

#install.packages("effects")
library(effects)

plot(allEffects(modelobio)["Biome"])
plot(predictorEffects(modelobio,residuals=T))


########################################################################
# MANOVA (Multivariate Analysis of Variance) ## Método R. Maestri

manova<-procD.lm(shape~species, iter=999, RRPP= TRUE) # P/ espécies
summary(manova)

manova.sex<-procD.lm(shape~sex, iter=999, RRPP= TRUE) # P/ sexo
summary(manova.sex)

manova.bio<-procD.lm(shape~biome, iter=999, RRPP= TRUE) # P/ bioma
summary(manova.bio)

# MANOVA Wilks's lambda #Mais prox de 0 mais diferença entre grupos #Teste parecido com o de p. Para comparar entre os grupos

PCA
manova.w<-manova(PCA$x[,1:16]~biome)
summary(manova.w,test="Wilks")

manova.w<-manova(PCA$x[,1:16]~Fac)
summary(manova.w,test="Wilks")

manova.w.sex<-manova(PCA$x[,1:32]~sex)
summary(manova.w.sex,test="Wilks")

manova.w.spec<-manova(PCA$x[,1:16]~log(gpa$Csize)*sex)
summary(manova.w.spec,test="Wilks")

manova.w<-manova(PCA$x[,1:16]~ log(gpa$Csize)*biome)
summary(manova.w,test="Wilks")

manova.w<-manova(PCA$x[,1:16]~log(gpa$Csize)*Fac)
summary(manova.w,test="Wilks")

# Pairwise comparisons
require(RRPP)
manova.pairwise<-pairwise(manova.w.spec,groups=species)
summary(manova.pairwise)

# Fenograma (árvore de distância morfológica), distância de Procrustes, agrupamento
#de Neighbor-Joining
par(mfrow=c(1,1))

obj<-summary(manova.pairwise)
obj$pairwise.tables$D
means.dist<-obj$pairwise.tables$D
fen<-nj(means.dist)
plot(fen)
plot(fen,type="unrooted")


df <- data.frame(Shape = PCA$x, Species = species, Latitude = latitude, Biome = biome, Envi = Fac, Size = size)

# Maior complexidade na morfologia: As PCs add Maior complexidade na morfologia

# As PCs adicionais (PC4 a PC8) podem capturar aspectos mais sutis e complexos da morfologia que são relevantes para diferenciar entre os biomas. Isso pode indicar que a morfologia dos indivíduos varia de forma mais sutil do que o que é capturado pelas três primeiras PCs.

# Maior poder de discriminação: A inclusão de mais PCs aumenta o poder estatístico do modelo para detectar diferenças significativas entre os biomas. Isso pode ser especialmente verdadeiro se houver uma grande quantidade de variação morfológica entre os biomas, que só pode ser capturada por uma quantidade maior de PCs.


# Boxplot

# Gráfico de dispersão #

ggplot(df, aes(x = biome, y = Size, color = sex)) +
  geom_point(position = position_jitter(width = 0.12), size = 3) +  
  geom_smooth(method = "lm", se = TRUE, aes(group = sex), color = "black", linetype = "dashed") +
  labs(title = "Centroid Size relation with Biome by sex",
       x = "Biome",
       y = "Centroid Size",
       color = "Sex") +
  theme_minimal()

# Direfença estatísticamente significativa entre as espécies de Floresta e Savana (Cerrado + Caatinga)

# Mais Boxplots...

ggplot(df, aes(x = paste(Fac,sex), y = Size, fill = Fac)) +  # Use 'biome' para mapear as cores
  geom_boxplot() +
  scale_fill_manual(values = c("green3","darkgreen","goldenrod"), labels = c("Atlantic Forest","Amazon Forest", "Savanna")) +  # Usando as cores definidas acima
  labs(x = "Enviroment", y = "Centroid Size", fill = "Biome") +
  ggtitle("Size by Ecoregions")


##################################### Alometria ##########################################

gdf <- geomorph.data.frame(Shape = gpa$coords, Size = gpa$Csize,Species = species, Latitude = latitude, Biome = biome, Fac = fatores$fac, Sex = sex)

fit.size <- procD.lm(gpa$coords~log(gpa$Csize), data = gdf, iter = 999)
summary(fit.size) ### Resultado geral

fit.size <- procD.lm(gpa$coords~log(gpa$Csize)*Fac, data = gdf, iter = 999)
summary(fit.size) ### Resultado considerando os diferentes slopes

plotAllometry(fit.size, 
              size = gpa$Csize, 
              logsz = TRUE, 
              method = "PredLine", 
              pch = as.numeric(pch), bg = cores_fac, 
              cex = 4, 
              col = cores_fac)


legend("bottomright", legend = unique(biome), pch = unique(as.numeric(biome) + 22),title = "Biomes", cex = 2, col = "black")

# Size-Shape PCA

pc.plot <- plotAllometry(fit.size, size = (gpa$Csize), logsz = TRUE, 
                         method = "size.shape", 
                         pch = as.numeric(pch), bg = cores_fac, 
                         cex = 5, 
                         col = cores_fac)

AlloPCA <- pc.plot$size.shape.PCA
pc.plot$plot_args

scorespc1<-AlloPCA$x[,1] #scores for pc1
scorespc2<-AlloPCA$x[,2] #scores for pc2
cor.test(scorespc2, as.numeric(sex))

AlloPCs <- AlloPCA$x

mshape (gpa$coords)

predmin <- min((scorespc1))
predmax <- max((scorespc1))
preds <- shape.predictor(gpa$coords, scorespc1, Intercept = TRUE, predmin = predmin, predmax = predmax)

par(mfrow=c(1,2))
plotRefToTarget(mshape, preds$predmin, mag = 2, outline = Sapajusoutline$outline, gridPars = GP,  method = "points")  # Forma associada ao menor score
plotRefToTarget(mshape, preds$predmax, mag = 2, outline = Sapajusoutline$outline, gridPars = GP,  method = "points")  # Forma associada ao maior score
par(mfrow=c(1,1))

# Comparação em pares entre espécies - Ordinary Least Squares

gdf <- geomorph.data.frame(Shape = gpa$coords, Size = as.numeric(log(gpa$Csize)), Species = as.factor(fatores$sp), Sex = as.factor(fatores$sex), Biome = as.factor(fatores$biome), Fac = as.factor(fatores$fac))

fit.common <- procD.lm(gpa$coords ~ log(gpa$Csize)*species, data = gdf, print.progress = FALSE)
summary (fit.common)

fit.null <- procD.lm(gpa$coords ~ log(gpa$Csize)+ Sex, data = gdf, print.progress = FALSE)
summary(fit.null)

PW <- pairwise(fit.common,fit.null, group = species, print.progress = FALSE)
summary(PW, confidence = 0.95, test.type = "dist",stat.table=FALSE) #distance between vectors - no differences
summary(PW, confidence = 0.95, test.type = "DL",stat.table=FALSE) # no differences in vector lengths
summary(PW, confidence = 0.95, test.type = "VC", angle.type = "deg",stat.table=FALSE)# =vectors have similar orientation
summary(PW, confidence = 0.95, test.type = "var",stat.table=FALSE) # vectors differ in variance (amount of dispersion around estimated values for groups)

# Primeiro, monte a matriz de distâncias Procrustes:
sp <- c(
  "Sapajus_apella", "Sapajus_cay", "Sapajus_libidinosus",
  "Sapajus_nigritus", "Sapajus_robustus", "Sapajus_xanthosternos"
)

distances <- matrix(
  c(
    0.00000000, 0.02053338, 0.02600128, 0.03254954, 0.03394679, 0.04158019,
    0.02053338, 0.00000000, 0.02467312, 0.04457783, 0.03912553, 0.05200746,
    0.02600128, 0.02467312, 0.00000000, 0.03265863, 0.03809999, 0.04999526,
    0.03254954, 0.04457783, 0.03265863, 0.00000000, 0.04101024, 0.03830677,
    0.03394679, 0.03912553, 0.03809999, 0.04101024, 0.00000000, 0.03271396,
    0.04158019, 0.05200746, 0.04999526, 0.03830677, 0.03271396, 0.00000000
  ),
  nrow = 6, byrow = TRUE
)
rownames(distances) <- sp
colnames(distances) <- sp

# Em seguida, use ggplot2 para desenhar o heatmap:
if (!requireNamespace("ggplot2", quietly=TRUE)) install.packages("ggplot2")
if (!requireNamespace("reshape2", quietly=TRUE)) install.packages("reshape2")

library(reshape2)
library(ggplot2)

# "Derreta" a matriz para um data.frame longo
df_heat <- melt(distances)
colnames(df_heat) <- c("Species1", "Species2", "ProcrustesDist")

# Plot
ggplot(df_heat, aes(x = Species1, y = Species2, fill = ProcrustesDist)) +
  geom_tile(color = "white") +
  scale_fill_viridis_c(option = "plasma", name = "Distância\nProcrustes") +
  coord_fixed() +
  theme_minimal(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.title = element_blank()
  ) +
  labs(title = "Heatmap das Distâncias Procrustes entre Espécies")

# Disparidade morfológica comparada

morphol.disparity(Shape ~ 1, groups = ~ Species, data = gdf, iter = 999, print.progress = FALSE)
morphol.disparity(Shape ~ Size, groups = ~ Species, data = gdf, iter = 999, print.progress = FALSE)
morphol.disparity(Shape ~ Size, groups = ~ Species + Sex, data = gdf, iter = 999, print.progress = FALSE)
morphol.disparity(Shape ~ Size, groups = ~ Biome, data = gdf, iter = 999, print.progress = FALSE)
morphol.disparity(Shape ~ Size, groups = ~ Fac + Sex, data = gdf, iter = 999, print.progress = FALSE)
morphol.disparity(Shape ~ Size + Sex, groups = ~ Species, data = gdf, iter = 999, print.progress = FALSE)

### RESÍDUOS DA ALOMETRIA

# Obtendo resíduos formaXtamanho (allometry-free shapes)
alometria<-procD.lm(gpa$coords~log(gpa$Csize),iter=999)
dim(gpa$coords)
shape.resid<-
  arrayspecs(alometria$residuals,p=dim(gpa$coords)[1],k=dim(gpa$coords)[2]) ## extrai dos resultados os residuos da associação de forma e tamanho #size adjusted residuals
adj.shape<-shape.resid+array(gpa$consensus, dim(shape.resid)) # allometry-free ## variável de forma livre de efeito de alometria
dim(adj.shape)
#shapes--somar a média retorna os números para o espaço de forma original e permite
#visualização de mudanças de forma

gdf <- geomorph.data.frame(Shape = adj.shape, Species = spec, Biome = biome, 
                           Size = log(gpa$Csize), Sex = sex)


gdf$residuals <- adj.shape

pca_res <- gm.prcomp(adj.shape) ## forma dos resíduos alométricos
summary(pca_res)

residual_components <- pca_res$x

colnames(pca_res$x)
pca_res_df <- data.frame(
  PC1 = pca_res$x[, "Comp1"],
  PC2 = pca_res$x[, "Comp2"],
  PC3 = pca_res$x[, "Comp3"],
  pch = as.numeric(fatores$pch),
  Species = as.factor(species),  # Supondo que "species" seja o vetor com os fatores das espécies
  Individuo = 1:nrow(PCA$x), Fac = Fac, residual_comps = pca_res$x)

PC1 <- pca_res_df$PC1
PC2 <- pca_res_df$PC2
PC3 <- pca_res_df$PC3

correlacao <- cor(cbind(PC1, PC2, PC3, size))
print(correlacao)
correlacao <- cor(cbind(PC1, PC2, PC3, biome))
print(correlacao)

summary(pca_res)
main_components <- pca_res$x
residual_components <- pca_res$x

correlacao_bio <- cor(cbind(main_components),as.numeric(biome))
print(correlacao)

correlacao <- cor(cbind(main_components),size)
print(correlacao)

### PCA dos Resíduos alométricos

pca_res_df$Sex <- as.factor(fatores$sex)

compute_hull <- function(pca_res_df) pca_res_df[chull(pca_res_df$PC1, pca_res_df$PC2), ]
hulls <- pca_res_df %>%
  group_by(Sex) %>%
  do(compute_hull(.))

ggplot(pca_res_df, aes(x = PC1, y = PC2)) +
  geom_point(aes(shape = Species, color = Sex), size = 6) +  
  geom_text(aes(label = Individuo), size = 2, vjust = 2) +  
  geom_polygon(data = hulls, aes(x = PC1, y = PC2, group = Sex, fill = Sex, color = Sex), alpha = 0.2) +  
  labs(
    title = "Residual Principal Components", 
    x = paste("PC1 (", round(100 * pca_res$sdev[1]^2 / sum(pca_res$sdev^2), 1), "%)", sep = ""), 
    y = paste("PC2 (", round(100 * pca_res$sdev[2]^2 / sum(pca_res$sdev^2), 1), "%)", sep = ""), 
    shape = "Species", color = "Sex"
  ) +
  theme_minimal() +
  scale_shape_manual(values = uPCH, labels = levels(pca_res_df$Species), name = "Species") +  
  scale_color_manual(values = cores_sex, name = "Sex") +  # Defina `cores_sexo` com as cores correspondentes aos sexos
  scale_fill_manual(values = cores_sex, name = "Sex")


compute_hull <- function(pca_res_df) pca_res_df[chull(pca_res_df$PC1, pca_res_df$PC2), ]
hulls <- pca_res_df %>%
  group_by(Fac) %>%
  do(compute_hull(.))

ggplot(pca_res_df, aes(x = PC1, y = PC2)) +
  geom_point(aes(shape = Species, color = Fac), size = 6) +  
  geom_text(aes(label = Individuo), size = 2, vjust = 2) +  
  geom_polygon(data = hulls, aes(x = PC1, y = PC2, group = Fac, fill = Fac, color = Fac), alpha = 0.2) +  
  labs(title = "Residuals Principal Components", 
       x = paste("PC1 (", round(100 * pca_res$sdev[1]^2 / sum(pca_res$sdev^2), 1), "%)", sep = ""), 
       y = paste("PC2 (", round(100 * pca_res$sdev[2]^2 / sum(pca_res$sdev^2), 1), "%)", sep = ""), 
       shape = "Species", color = "Fac") +
  theme_minimal() +
  scale_shape_manual(values = uPCH, labels = levels(pca_res_df$Species), name = "Species") +  
  scale_color_manual(values = cores_fac, name = "Environment") +  
  scale_fill_manual(values = cores_fac, name = "Environment")

compute_hull <- function(pca_res_df) pca_res_df[chull(pca_res_df$PC2, pca_res_df$PC3), ]
hulls <- pca_res_df %>%
  group_by(Fac) %>%
  do(compute_hull(.))
ggplot(pca_res_df, aes(x = PC2, y = PC3)) +
  geom_point(aes(shape = Species, color = Fac), size = 6) +  
  geom_text(aes(label = Individuo), size = 2, vjust = 2) +  
  geom_polygon(data = hulls, aes(x = PC2, y = PC3, group = Fac, fill = Fac, color = Fac), alpha = 0.2) +  
  labs(title = "Shape's Principal Components", 
       x = paste("PC2 (", round(100 * pca_res$sdev[2]^2 / sum(pca_res$sdev^2), 1), "%)", sep = ""), 
       y = paste("PC3 (", round(100 * pca_res$sdev[3]^2 / sum(pca_res$sdev^2), 1), "%)", sep = ""), 
       shape = "Species", color = "Fac") +
  theme_minimal() +
  scale_shape_manual(values = uPCH, labels = levels(pca_res_df$Species), name = "Species") +  
  scale_color_manual(values = cores_fac, name = "Environment") +  
  scale_fill_manual(values = cores_fac, name = "Environment")



compute_hull <- function(pca_res_df) pca_res_df[chull(pca_res_df$PC1, pca_res_df$PC2), ]

hulls <- pca_res_df %>%
  group_by(Sex, Fac) %>%  # Separar por `Sex` e `Fac`
  do(compute_hull(.)) %>%
  ungroup()

hulls <- hulls %>%
  mutate(Sex = factor(Sex), Fac = factor(Fac))

ggplot(pca_res_df, aes(x = PC1, y = PC2)) +
  geom_point(aes(shape = Species, color = Fac), size = 6) +  
  geom_text(aes(label = Individuo), size = 2, vjust = 2) +  
  geom_polygon(
    data = hulls,
    aes(x = PC1, y = PC2, group = interaction(Sex, Fac), fill = Fac, color = Fac),
    alpha = 0.2,   # Opacidade para o preenchimento
    linewidth = 1  # Espessura do contorno para destaque
  ) +  
  labs(
    title = "Resitual's Principal Components (Separated by Sex)",
    x = paste("PC1 (", round(100 * pca_res$sdev[1]^2 / sum(pca_res$sdev^2), 1), "%)", sep = ""), 
    y = paste("PC2 (", round(100 * pca_res$sdev[2]^2 / sum(pca_res$sdev^2), 1), "%)", sep = ""), 
    shape = "Species", color = "Environment", fill = "Environment"
  ) +
  theme_minimal() +
  scale_shape_manual(values = uPCH, labels = levels(pca_res_df$Species), name = "Species") +  
  scale_color_manual(values = cores_fac, name = "Environment") +  
  scale_fill_manual(values = cores_fac, name = "Environment") +
  facet_wrap(~Sex)  # Separar por sexo


par(mfrow=c(1,2))
scorespc1<-pca_res$x[,1] #scores for pc1
scorespc2<-pca_res$x[,2] #scores for pc2

mshape <- mshape(adj.shape)
predsPC1 <- shape.predictor(adj.shape, scorespc1, Intercept = TRUE, predmin = min(scorespc1), predmax = max(scorespc1)) #estimating shape configurations based on pc scores
plotRefToTarget(mshape, predsPC1$predmin, mag=2, outline = Sapajusoutline$outline,gridPars = GP,  method = "points") #shape related to small pc scores
plotRefToTarget(mshape, predsPC1$predmax, mag=2, outline = Sapajusoutline$outline,gridPars = GP,  method = "points")#shape related to large pc scores

predsPC2 <- shape.predictor(adj.shape, scorespc2, Intercept = TRUE, predmin = min(scorespc2), predmax = max(scorespc2)) #estimating shape configurations based on pc scores
plotRefToTarget(mshape, predsPC2$predmin, mag=2, outline = Sapajusoutline$outline,gridPars = GP,  method = "points") #shape related to small pc scores
plotRefToTarget(mshape, predsPC2$predmax, mag=2, outline = Sapajusoutline$outline,gridPars = GP,  method = "points")#shape related to large pc scores

scorespc3<-pca_res$x[,3] #scores for pc3
par(mfrow=c(1,2))

predsPC3 <- shape.predictor(gpa$coords, scorespc3, Intercept = TRUE, predmin = min(scorespc3), predmax = max(scorespc3)) #estimating shape configurations based on pc scores
plotRefToTarget(mshape, predsPC3$predmin, mag = 2, outline = Sapajusoutline$outline,gridPars = GP,  method = "points") #shape related to small pc scores
plotRefToTarget(mshape, predsPC3$predmax, mag = 2, outline = Sapajusoutline$outline,gridPars = GP,  method = "points")#shape related to large pc scores

par(mfrow=c(1,1))

dev.off()


############################## Análises separadas por sexo ####################################

### Sexed M

Species_M <- as.factor(Males$sp)
Biome_M <- as.factor(Males$biome)
symbols_M <- symbols[as.character(Males$sp)]
M_mshape <- mshape(M_gpa$coords)

Fac_M <- as.factor(Males$fac)

levels_fac_M <- unique(Males$fac)

cores_fac_M <- c(AF = "green3", AM = "darkgreen", SV ="goldenrod1")[as.factor(Fac)]

M_pca <- gm.prcomp(M_gpa$coords)
summary(M_pca)
indv <- 1:nrow(M_pca$x)

M_pca_df <- data.frame(
  PC1 = M_pca$x[, "Comp1"],
  PC2 = M_pca$x[, "Comp2"],
  PC3 = M_pca$x[, "Comp3"],
  pch = as.numeric(Males$pch),
  Species = as.factor(Males$sp),  # Supondo que "species" seja o vetor com os fatores das espécies
  Individuo = 1:nrow(M_pca$x), Fac = Fac_M)

species_colors_M <- c("red","navy","cyan3","saddlebrown","magenta4","black")
names(species_colors_M) <- levels(M_pca_df$Species)

M_fit.size <- procD.lm(M_gpa$coords~log(M_gpa$Csize), iter = 9999)
summary(M_fit.size)

# Obtendo resíduos formaXtamanho (allometry-free shapes)
alometria<-procD.lm(M_gpa$coords~log(M_gpa$Csize),iter=999)
shape.resid<-
  arrayspecs(alometria$residuals,p=dim(M_gpa$coords)[1],k=dim(M_gpa$coords)[2]) ## extrai dos resultados os residuos da associação de forma e tamanho #size adjusted residuals
M_adj.shape<-shape.resid+array(M_gpa$consensus, dim(shape.resid)) # allometry-free ## variável de forma livre de efeito de alometria
#shapes--somar a média retorna os números para o espaço de forma original e permite
#visualização de mudanças de forma

M_residuals <- M_adj.shape

M_fit_bio <- procD.lm(M_gpa$coords~log(M_gpa$Csize)*Biome_M, iter = 9999)
summary(M_fit_bio)

M_fit_bio <- procD.lm(M_residuals~Biome_M, iter = 9999)
summary(M_fit_bio)

M_pca_res <- gm.prcomp(M_residuals)
summary(M_pca_res)

M_pch <- as.numeric(Males$pch)
M_uPCH <- unique(M_pch)
names(M_uPCH) <- unique(M_pca_df$Species)


plotAllometry(M_fit_bio, 
              size = log(M_gpa$Csize), 
              logsz = FALSE, 
              method = "PredLine", 
              pch = c(15,18)[as.numeric(as.factor(Biome_M))], 
              cex = 2, 
              col = c("green3", "goldenrod")[as.numeric(as.factor(Biome_M))])

legend("bottomright",  legend = unique(Biome_M), pch = 16, title = "Biome", cex = 1.2, col = c("green3", "goldenrod")[unique(as.factor(Biome_M))])

dev.off()

### Sexed F 
Species_F<- as.factor(Females$sp)
Biome_F<- as.factor(Females$biome)
symbols_F <- symbols[as.character(Females$sp)]

Fac_F <- as.factor(Females$fac)

levels_fac_F <- unique(Females$fac)

cores_fac_F <-c(AF = "green3", AM = "darkgreen", SV ="goldenrod1")[as.factor(Fac)]

F_pca <- gm.prcomp(F_gpa$coords)
summary(F_pca)
indv <- 1:nrow(F_pca$x)

F_pca_df <- data.frame(
  PC1 = F_pca$x[, "Comp1"],
  PC2 = F_pca$x[, "Comp2"],
  PC3 = F_pca$x[, "Comp3"], PC5 = F_pca$x[, "Comp5"],
  pch = as.numeric(Females$pch),
  Species = as.factor(Females$sp),  # Supondo que "species" seja o vetor com os fatores das espécies
  Individuo = 1:nrow(F_pca$x), Fac = Fac_F)

species_colors_F <- c("red","navy","cyan3","saddlebrown","magenta4","black")
names(species_colors_F) <- levels(F_pca_df$Species)

# Obtendo resíduos formaXtamanho (allometry-free shapes)
alometria<-procD.lm(F_gpa$coords~log(F_gpa$Csize),iter=999)
shape.resid<-
  arrayspecs(alometria$residuals,p=dim(F_gpa$coords)[1],k=dim(F_gpa$coords)[2]) ## extrai dos resultados os residuos da associação de forma e tamanho #size adjusted residuals
F_adj.shape<-shape.resid+array(F_gpa$consensus, dim(shape.resid)) # allometry-free ## variável de forma livre de efeito de alometria
#shapes--somar a média retorna os números para o espaço de forma original e permite
#visualização de mudanças de forma

F_residuals <- F_adj.shape

F_fit_bio <- procD.lm(F_gpa$coords~log(F_gpa$Csize)*Biome_F, iter = 9999)
summary(F_fit_bio)

F_fit_bio <- procD.lm(F_residuals~Biome_F, iter = 9999)
summary(F_fit_bio)

F_pca_res <- gm.prcomp(F_residuals)
summary(F_pca_res)

F_pch <- as.numeric(Females$pch)
F_uPCH <- unique(F_pch)
names(F_uPCH) <- unique(F_pca_df$Species)
F_pca_res_df <- data.frame(
  PC1 = F_pca_res$x[, "Comp1"],
  PC2 = F_pca_res$x[, "Comp2"],
  PC3 = F_pca_res$x[, "Comp3"],
  pch = as.factor(Females$pch),
  Species = as.factor(Females$sp),  # Supondo que "species" seja o vetor com os fatores das espécies
  Individuo = 1:nrow(F_pca$x), Fac = Fac_F, residual_comps = F_pca_res$x)


plotAllometry(F_fit_bio, 
              size = log(F_gpa$Csize), 
              logsz = FALSE, 
              method = "PredLine", 
              pch = c(15,18)[as.numeric(as.factor(Biome_F))], 
              cex = 2, 
              col = c("green3", "goldenrod")[as.numeric(as.factor(Biome_F))])

legend("bottomright",  legend = unique(Biome_F), pch = 16, title = "Biome", cex = 1.2, col = c("green3", "goldenrod")[unique(as.factor(Biome_F))])

# Forma média geral dos machos
mshapeM <- mshape(M_gpa$coords)  

## Separando dados por bioma ##
biomes <- unique(Males$biome)  

# Separar coordenadas de machos por bioma
sppsplitpanthM_biome <- coords.subset(tps_M, Males$biome)

# Configurar layout para múltiplos gráficos (1 linha por número de biomas)
par(mfrow = c(1, length(biomes)), oma = c(0, 0, 0, 0), mar = c(2, 2, 2, 2), mai = c(0.3, 0.3, 0.3, 0.3))

for (i in biomes) {
  print(i)  # Para acompanhar o progresso
  
  # GPA para machos por bioma
  GPA_M <- gpagen(sppsplitpanthM_biome[[i]])  
  
  # Forma média dos machos no bioma específico
  mshapeM_biome <- mshape(GPA_M$coords)
  
  # Plotar comparação entre a forma média geral dos machos e a do bioma
  plotRefToTarget(mshape, mshapeM_biome, mag=3, outline = Sapajusoutline$outline, gridPars = GP,  method = "points")
  title(paste("Machos -", i))
}

# Resetar layout para padrão
par(mfrow = c(1,1))


# Forma média geral das fêmeas
mshapeF <- mshape(F_gpa$coords)  

## Separando dados por bioma ##
biomes <- unique(Females$biome)  

# Separar coordenadas de fêmeas por bioma
sppsplitpanthF_biome <- coords.subset(tps_F, Females$biome)

# Configurar layout para múltiplos gráficos (1 linha por número de biomas)
par(mfrow = c(1, length(biomes)), oma = c(0, 0, 0, 0), mar = c(2, 2, 2, 2), mai = c(0.3, 0.3, 0.3, 0.3))

for (i in biomes) {
  print(i)  # Para acompanhar o progresso
  
  # GPA para fêmeas por bioma
  GPA_F <- gpagen(sppsplitpanthF_biome[[i]])  
  
  # Forma média das fêmeas no bioma específico
  mshapeF_biome <- mshape(GPA_F$coords)
  
  # Plotar comparação entre a forma média geral das fêmeas e a do bioma
  plotRefToTarget(mshape, mshapeF_biome, mag=3, outline = Sapajusoutline$outline, gridPars = GP,  method = "points")
  title(paste("Fêmeas -", i))
}

par(mfrow = c(1,1))


################################### Climatic Variables - AVGLOCSEX ########################

############ Partial Least Squares ###################

names(fatores)

bioclim <- as.matrix(fatores[,c(11:31)])
#write.csv(bioclim, "bioclim_avg_std.csv")
bioclim <- scale(bioclim,center=T,scale=T)

vifstep_result <- vifstep(bioclim,keep = c("bio12")) ### Vif é usado pra verificar se as variáveis não causam multicolinearidade
print(vifstep_result)

selected_vars <- vifstep_result@results$Variables
#bioclim <- as.matrix(fatores[,selected_vars])
bioclim <- scale(bioclim,center=T,scale=T)

rownames(bioclim) <- paste(fatores$sp, fatores$lat,fatores$long, fatores$sex)
dimnames(gpa$coords)[[3]] <- paste(fatores$sp, fatores$lat,fatores$long, fatores$sex)

#rownames(bioclim) <- paste(1:106)
#dimnames(gpa$coords)[[3]] <- paste(1:106)

fatores$fac <- as.factor(fatores$fac)
fatores$biome <- as.factor(fatores$biome)
levels_fac <- levels(fatores$fac)
levels_biome <- levels(fatores$biome)
cores_fac <- c("green3", "darkgreen", "goldenrod3")
names(cores_fac) <- levels_fac
bioma_shapes <- c("Atlantic Forest" = 16, "Amazon" = 15, "SV" = 17)
#?pch
names(bioma_shapes) <- levels_fac

PLS <- two.b.pls(A1 = bioclim, A2 = gpa$coords,iter = 9999, print.progress = TRUE)
summary (PLS)

PLS$A1
PLS$r.pls
PLS$P.value
PLS$Z

#PLS1 = 0.4765
#PLS2 = 0.2652
#PLS3 = 0.1367
#PLS4 = 0.0587

data_plot <- na.omit(data.frame(
  PLS_Block_1 = PLS$XScores[, 1],
  PLS_Block_2 = PLS$YScores[, 1],
  fac = fatores$fac,
  biome = fatores$biome, fac = fatores$fac))

P <-ggplot(data_plot, aes(x = PLS_Block_1, y = PLS_Block_2)) +
  geom_point(aes(color = fac, shape = fac), size = 10) +
  scale_color_manual(values = cores_fac) +
  scale_shape_manual(values = bioma_shapes) +
  theme_minimal() +
  labs(title = "PLS1 Plot", x = "Block 1 (Climatic variables)", y = "Block 2 (Shape)", color = "Ecoregions", shape = NULL) +
  theme(
    legend.position = "bottom",
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 8),
    plot.title = element_text(hjust = 0.5, size = 14)
  ) +
  guides(color = guide_legend(override.aes = list(size = 4)),
         shape = guide_legend(override.aes = list(size = 4)))
print(P)
par(mfrow=c(1,1))

mshape <- mshape(gpa$coords)
symbols <- c(
  Sapajus_apella = "\u25A0",        # Quadrado 
  Sapajus_cay = "\u25BC",           # Triângulo apontando para baixo
  Sapajus_libidinosus = "\u25B2",   # Triângulo apontando para cima
  Sapajus_nigritus = "\u25CF",      # Círculo 
  Sapajus_robustus = "\u2726",      # Estrela com 6 pontas
  Sapajus_xanthosternos = "\u2736"  # Estrela com 8 pontas
)
symbols <- symbols[as.character(species)]
fac <- as.factor(fatores$fac)
cex <- ifelse(species %in% c("Sapajus_cay", "Sapajus_libidinosus", "Sapajus_apella", "Sapajus_nigritus"), 1, 1.2)
P <- plot (PLS, col = fac, pch = as.numeric(fatores$pch), cex =cex + 1.5)
legend("bottomright", legend = unique(fac), col = unique(fac), pch = 19, title = "Enviroments", cex = 0.8)
legend("topleft", legend = unique(fatores$sp), col = "skyblue", pch = unique(as.numeric(fatores$pch)), title = "Species", cex = 0.8)


miny <- min(P$plot_args$y)
maxy <- max(P$plot_args$y)
mean <- mean(P$plot_args$y)

preds <- shape.predictor(P$A2, x = P$plot.args$y, min = miny, max = maxy, mean = mean)
grid.pars <- gridPar(grid.col = "white") 
GP <- gridPar(n.col.cell = 100, pt.bg = "gray", pt.size = 0.8, tar.pt.bg = "cyan", 
               tar.pt.size = 1.5,tar.out.col = "gray10", tar.out.cex = 0.8, grid.col = "gray", grid.lwd = 0.5,txt.pos = 1, txt.col = "steelblue")


par(mfrow=c(1,2))
plotRefToTarget(mshape(P$A2), preds$min, outline = Sapajusoutline$outline,gridPars = GP, mag = 2, method = "points")
plotRefToTarget(mshape(P$A2), preds$max, outline = Sapajusoutline$outline,gridPars = GP, mag = 2, method = "points")
par(mfrow=c(1,1))
dev.off()

## PLS alometria ##

names(gpa$Csize)<- paste(fatores$sp, fatores$lat,fatores$long, fatores$sex)
PLS <- two.b.pls(A1 = log(gpa$Csize), A2 = gpa$coords, print.progress = TRUE)
summary (PLS)

data_plot <- na.omit(data.frame(
  PLS_Block_1 = PLS$XScores[, 1],
  PLS_Block_2 = PLS$YScores[, 1],
  fac = fatores$fac,
  biome = fatores$biome, fac = fatores$fac))

P <-ggplot(data_plot, aes(x = PLS_Block_1, y = PLS_Block_2)) +
  geom_point(aes(color = fac, shape = fac), size = 4) +
  scale_color_manual(values = cores_fac) +
  scale_shape_manual(values = bioma_shapes) +
  theme_minimal() +
  labs(title = "PLS1 Plot",x = "Block 1 (Ln Centroid Size)", y = "Block 2 (Shape)", color = "Ecoregions", shape = NULL) +
  theme(
    legend.position = "bottom",
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 8),
    plot.title = element_text(hjust = 0.5, size = 14)
  ) +
  guides(color = guide_legend(override.aes = list(size = 4)),
         shape = guide_legend(override.aes = list(size = 4)))
print(P)

mshape (gpa$coords)

fac <- as.factor(fatores$fac)
P <- plot (PLS, col = fac, pch = symbols, cex = cex + 0.)
legend("bottomright", legend = unique(fac), col = unique(fac), title = "Enviroments", pch = as.numeric(fatores$pch),cex = 0.8)
legend("topleft", legend = unique(fatores$sp), col = c("red", "lightgreen", "lightgreen","black","black","black"), pch = as.numeric(fatores$pch), title = "Species", cex = 0.8)

PLS$right.pls.vectors

minx <- min(P$plot_args$y)
maxx <- max(P$plot_args$y)
mean <- mean(P$plot_args$y)

preds <- shape.predictor(P$A2, x = P$plot.args$y, min = minx, max = maxx, mean = mean)

par(mfrow=c(1,3))
plotRefToTarget(mshape(P$A2), preds$min, outline = Sapajusoutline$outline,gridPars = GP,  method = "points", mag = 2)
plotRefToTarget(mshape(P$A2), preds$mean, outline = Sapajusoutline$outline,gridPars = GP,  method = "points", mag = 2)
plotRefToTarget(mshape(P$A2), preds$max, outline = Sapajusoutline$outline,gridPars = GP,  method = "points", mag = 2)
par(mfrow=c(1,1))
dev.off()

### Jamile's PLS ##
names(fatores)
env <- bioclim
env <- as.matrix(scale(fatores[,c(11:31)]))

fit.size <- procD.lm(gpa$coords ~ log(gpa$Csize), iter = 999)

rownames(env) <- paste(fatores$sp, fatores$lat,fatores$long,fatores$sex)
dimnames(gpa$coords)[[3]] <- paste(fatores$sp, fatores$lat,fatores$long, fatores$sex)

pls<-two.b.pls(env, gpa$coords,iter=9999)
summary(pls)
plot (pls, col = fac, pch = as.numeric(fatores$pch), cex = cex + 0.)
pls$left.pls.vectors
predpls<-pls$YScores[,1]

preds <- shape.predictor(gpa$coords, predpls, Intercept = TRUE, predmin = min(predpls), predmax = max(predpls)) #estimating shape configurations based on pc scores

envvector<-data.frame(pls$left.pls.vectors[,1])
names(envvector)[names(envvector) == "pls.left.pls.vectors...1."] <- "envvar"
envvector$id <- row.names(envvector)
envvector$colour <- ifelse(envvector$envvar < 0, "Negative","Positive")
envvector$hjust <- ifelse(envvector$envvar > 0, 1.3, -0.3)

ggplot(data=envvector,aes(id,envvar,label="",hjust=hjust))+
  geom_text(aes(y=0,colour=colour))+
  geom_bar(stat="identity",position="identity",aes(fill = colour))+theme_classic()

require(ape)

#non-allometric component

alometria<-procD.lm(gpa$coords~log(gpa$Csize),iter=999)
shape.resid<-
  arrayspecs(alometria$residuals,p=dim(gpa$coords)[1],k=dim(gpa$coords)[2]) ## extrai dos resultados os residuos da associação de forma e tamanho #size adjusted residuals
adj.shape<-shape.resid+array(gpa$consensus, dim(shape.resid)) # allometry-free ## variável de forma livre de efeito de alometria

#shapes--somar a média retorna os números para o espaço de forma original e permite
#visualização de mudanças de forma

plsAlo<-two.b.pls(env,adj.shape,iter=9999)
summary(plsAlo)
plot (plsAlo, col = as.factor(fatores$fac), pch = symbols, cex = cex + 0.)
print(plsAlo$svd)
d <- plsAlo$svd$d

plot(plsAlo$svd$u[,1], plsAlo$svd$u[,2], main = "PLS1 vs PLS2 (Bloco 1)", xlab = "PLS1", ylab = "PLS2")
barplot(plsAlo$svd$vt[ ,1], las=2, main="Pesos das variáveis no PLS1", horiz=TRUE)

round(d^2 / sum(d^2), 3)

# Scores dos blocos
XScores <- plsAlo$XScores
YScores <- plsAlo$YScores
n_axes <- ncol(XScores)
obs.cor <- sapply(1:n_axes, function(i) cor(XScores[, i], YScores[, i]))
# Tamanho do vetor de permutações
n_permutations <- length(plsAlo$random.r)

# Inicializar os p-valores
pvals <- sapply(1:n_axes, function(j) {
  # Calcula o p-valor para o eixo j
  mean(abs(plsAlo$random.r) >= abs(obs.cor[j]))
})

# Resultados em tabela
pls.results <- data.frame(
  Eixo = paste0("PLS", 1:n_axes),
  Cor_Observada = round(obs.cor, 3),
  P_Valor = round(pvals, 4)
)

# Exibir resultados
print(pls.results)

envvector<-data.frame(plsAlo$left.pls.vectors[,1])
names(envvector)[names(envvector) == "plsAlo.left.pls.vectors...1."] <- "envvar"
envvector$id <- row.names(envvector)
envvector$colour <- ifelse(envvector$envvar < 0, "Negative", "Positive")
envvector$hjust <- ifelse(envvector$envvar > 0, 1.3, -0.3)

ggplot(data=envvector,aes(id,envvar,label="",hjust=hjust))+
  geom_text(aes(y=0,colour=colour))+
  geom_bar(stat="identity",position="identity",aes(fill = colour))+theme_classic()


data_plot <- na.omit(data.frame(
  PLS_Block_1 = plsAlo$XScores[, 1],
  PLS_Block_2 = plsAlo$YScores[, 1],
  fac = fatores$fac,
  biome = fatores$biome, fac = fatores$fac))

P <-ggplot(data_plot, aes(x = PLS_Block_1, y = PLS_Block_2)) +
  geom_point(aes(color = fac, shape = fac), size = 10) +
  scale_color_manual(values = cores_fac) +
  scale_shape_manual(values = bioma_shapes) +
  theme_minimal() +
  labs(title = "Non Allometric PLS", x = "Block 1 (Climatic variables)", y = "Block 2 (Shape)", color = "Ecoregions", shape = NULL) +
  theme(
    legend.position = "bottom",
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 8),
    plot.title = element_text(hjust = 0.5, size = 14)
  ) +
  guides(color = guide_legend(override.aes = list(size = 4)),
         shape = guide_legend(override.aes = list(size = 4)))
print(P)
par(mfrow=c(1,1))

P <- plot (plsAlo, col = fac, pch = symbols, cex = cex + 0.)
legend("bottomright", legend = unique(fac), col = unique(fac), title = "Enviroments", pch = as.numeric(fatores$pch),cex = 0.8)
legend("topleft", legend = unique(fatores$sp), col = c("red", "lightgreen", "lightgreen","black","black","black"), pch = as.numeric(fatores$pch), title = "Species", cex = 0.8)

plsAlo$right.pls.vectors

minx <- min(P$plot_args$y)
maxx <- max(P$plot_args$y)
mean <- mean(P$plot_args$y)
preds <- shape.predictor(P$A2, x = P$plot.args$y, min = minx, max = maxx, mean = mean)

par(mfrow=c(1,2))
plotRefToTarget(mshape(P$A2), preds$min, outline = Sapajusoutline$outline,gridPars = GP,  method = "points", mag = 2)
plotRefToTarget(mshape(P$A2), preds$max, outline = Sapajusoutline$outline,gridPars = GP,  method = "points", mag = 2)
par(mfrow=c(1,1))
dev.off()


### Aqui estão códigos para salvar os dados originais de diferentes formas

##### Average specimens by species+locality+sex

tps_0 <- readland.tps("Datasets/tps/jaw_18LM.tps")
Y.gpa <- gpagen(tps_0)
fatores <- read.csv("global.csv", sep = ";")

data<-fatores
sapajus<-tps_0

avgterm<-as.factor(paste(data$sp,data$lat,data$long,data$sex)) #creating averaging factor
x <- two.d.array(sapajus)#convert to 2d array

means <- rowsum(x, avgterm)/as.vector(table(avgterm)) # rowsum() is a simple base function to get the sum of the rows, while table() is a great function for getting group sizes (group n).

Y <- arrayspecs(means,dim(sapajus)[1],dim(sapajus)[2], sep=NULL)# then all you have to do is put the averaged data back into an 3D array.

writeland.tps(Y,file="Datasets/tps/avglocsex.tps")

names(data)
datatoaverage<-data[,c(6,7)]

meanenv <- rowsum(datatoaverage, avgterm)/as.vector(table(avgterm)) #average env data

write.table(meanenv,file="meanenv.txt")

### Mean Spec

avgterm<-as.factor(paste(data$sp)) #creating averaging factor
x <- two.d.array(sapajus)#convert to 2d array

means <- rowsum(x, avgterm)/as.vector(table(avgterm)) # rowsum() is a simple base function to get the sum of the rows, while table() is a great function for getting group sizes (group n).

Y <- arrayspecs(means,dim(sapajus)[1],dim(sapajus)[2], sep=NULL)# then all you have to do is put the averaged data back into an 3D array.

writeland.tps(Y,file="Datasets/tps/avgspec.tps")

names(data)
datatoaverage<-data[,c(6,7)]

meanenv <- rowsum(datatoaverage, avgterm)/as.vector(table(avgterm)) #average env data

write.table(meanenv,file="meanspec.txt")


##### Média por espécie + localidade + sexo separadamente #####


##### Average specimens by species+locality+sex

# Carregar os dados
tps_0 <- readland.tps("Datasets/tps/avglocsex.tps")
Y.gpa <- gpagen(tps_0)  # Processamento de GPA

# Carregar fatores de diferenciação
fatores <- read.csv("Planilhas/avglocsex.csv", sep = ";")

# Definir dados
data <- fatores
sapajus <- tps_0

# Criar fator de média para diferentes combinações de localização e sexo
avgterm <- as.factor(paste(data$sp, data$lat, data$long, data$sex))

# Converter dados para matriz 2D
x <- two.d.array(sapajus)

# Calcular as médias para cada grupo
means <- rowsum(x, avgterm) / as.vector(table(avgterm))

# Converter os dados de volta para uma matriz 3D
Y <- arrayspecs(means, dim(sapajus)[1], dim(sapajus)[2], sep = NULL)

# Filtrar os dados por sexo (Feminino e Masculino) e salvar os resultados
data_F <- subset(data, sex == "F")
data_M <- subset(data, sex == "M")

# Separando o array de "sapajus" para cada sexo, garantindo que estamos usando índices válidos
# A correspondência entre os dados e o array 'sapajus' é feita de acordo com a ordem das espécies.
# Vamos assumir que a ordem de 'data_F$sp' e 'data_M$sp' corresponde à ordem das linhas em 'sapajus'.

# Filtrando por sexo feminino
indices_F <- which(data$sex == "F")
tps_F <- sapajus[,,indices_F]

# Filtrando por sexo masculino
indices_M <- which(data$sex == "M")
tps_M <- sapajus[,,indices_M]

# Calcular as médias para cada sexo
avgterm_F <- as.factor(paste(data_F$sp, data_F$lat, data_F$long, data_F$sex))
avgterm_M <- as.factor(paste(data_M$sp, data_M$lat, data_M$long, data_M$sex))

# Para as mulheres
x_F <- two.d.array(tps_F)
means_F <- rowsum(x_F, avgterm_F) / as.vector(table(avgterm_F))
Y_F <- arrayspecs(means_F, dim(tps_F)[1], dim(tps_F)[2], sep = NULL)

# Para os homens
x_M <- two.d.array(tps_M)
means_M <- rowsum(x_M, avgterm_M) / as.vector(table(avgterm_M))
Y_M <- arrayspecs(means_M, dim(tps_M)[1], dim(tps_M)[2], sep = NULL)

# Salvar os arquivos TPS para cada grupo
writeland.tps(Y_F, file = "Datasets/tps/avgloc_F.tps")
writeland.tps(Y_M, file = "Datasets/tps/avgloc_M.tps")

names(data_F)

# Salvar as tabelas de médias de ambiente para cada sexo
meanenv_F <- rowsum(data_F[, c(7:8,11:32)], avgterm_F) / as.vector(table(avgterm_F))
meanenv_M <- rowsum(data_M[, c(7:8,11:32)], avgterm_M) / as.vector(table(avgterm_M))

write.table(meanenv_F, file = "meanenv_F.txt")
write.table(meanenv_M, file = "meanenv_M.txt")
