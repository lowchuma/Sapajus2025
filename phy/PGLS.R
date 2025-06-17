setwd("C:/Users/Lourenço/OneDrive/UFSM/LabMastozoo/Sapajus/Data")

library(ape)
library(phytools)
library(phangorn)
library(tidyverse)
library(geiger)
library(nlme)


# Criar uma árvore filogenética aleatória
set.seed(123)  # Para garantir a reprodutibilidade

# Visualizar a árvore filogenética

arvore <- read.tree("trees/AVGLOCSEX_Upham_2024.tre")
arvore

mcc <- mcc(arvore, tree = TRUE, part = NULL, rooted = TRUE)
plot(mcc)
plot.phylo (mcc, type = "phylogram", show.tip.label = TRUE, 
            show.node.label = TRUE, edge.color = "black", edge.width = 1.8, 
            tip.color = "black", cex = 0.6, label.offset = 0.08)

dados <- read.csv("avglocsex.csv", sep = ";")
head(dados)

# Criar a matriz de design com a variável preditora lncsize
x <- model.matrix(LogCS ~ bio8, data = dados)

# Criar a matriz de covariância filogenética a partir da árvore "mcc"
V <- vcv(mcc)

# Definir a variável de resposta como uma matriz de coluna
y <- as.matrix(dados$LogCS)

rownames(dados) <- dados$sp_ives
#dados <- dados[,c(2,3,4,5,6,7)]

mcc<- drop.tip(mcc,c("Cebus_albifrons", "Cebus_capucinus", "Sapajus_macrocephalus"))
obj<-name.check(mcc,dados)
obj

library(ggplot2)
modelo <- lm(LogCS ~ bio12, data = dados)
summary(modelo)

ggplot(dados, aes(x = bio12, y = LogCS)) +
  geom_point(color = "blue", alpha = 0.6) +  # Adiciona pontos com transparência
  geom_smooth(method = "lm", se = FALSE, color = "red") +  # Adiciona linha de tendência linear
  labs(x = "Mean Diurnal Range", y = "Size",  # Adiciona rótulos aos eixos
       title = "Mean Diurnal Range vs. Body Size") +  # Adiciona título
  theme_minimal() +  # Define um tema minimalista
  theme(plot.title = element_text(hjust = 0.5)) 

# Extract columns
size <- dados[, "LogCS"]
MDR <- dados[, "bio12"]

# Give them names
names(size) <- names(MDR) <- rownames(dados)

# Calculate PICs
hPic <- pic(size, mcc)
aPic <- pic(MDR, mcc)

# Make a model
picModel <- lm(hPic ~ aPic - 1)

# Yes, significant
summary(picModel)
