library(ggplot2)
library(geomorph)
library(phytools)
library(vegan)
library(geiger)
library(phangorn)
library(tidyverse)
library(usdm)

phy1 <-read.tree("Datasets/trees/AVGLOCSEX_Lima.tre")
phy2 <-read.tree("Datasets/trees/AVGLOCSEX_Wright.tre")
phy3 <-read.tree("Datasets/trees/AVGLOCSEX_Upham.tre")

library(caper)  # Pacote para rodar o PGLS
library(devtools)
#install_github("cran/difconet")
library(difconet)
library(car)


# Supondo que suas 19 variáveis estejam no data frame 'bio_data' e que as filogenias sejam 'phy1', 'phy2', e 'phy3'
data <- read.csv("Planilhas/avglocsex.csv", sep = ";")
response_variable <- as.data.frame(data$CS)
rownames(response_variable) <-data$sp_ives
names(data)
bio_data <- as.data.frame(scale(data[,11:32]))
CS <- as.numeric(data$CS)
rownames(bio_data) <-data$sp_ives

## verificar as correlações

M <- cor(bio_data)  # Remove a variável dependente "CS"
library(corrplot)
corrplot(M,method = 'circle', col = COL2(n=200), order = 'hclust', hclust.method = 'ward.D2', addrect = 4,tl.srt = 45)

### PCA para reduzir dimensionalidade das variáveis

pca <- prcomp(bio_data, scale. = TRUE)
summary(pca)
pca_data <- as.data.frame(pca$x)  # Extrai os scores dos componentes principais
pca_model <- lm(CS ~., data = cbind(pca_data))
summary(pca_model)
step_model <- step(pca_model, direction = "both")
summary(step_model)
vif(step_model)

plot(pca_model, pch = as.numeric(data$pch), col = as.factor(data$fac))
plot(step_model, pch = as.numeric(data$pch), col = as.factor(data$fac))

shapiro.test(pca_model$residuals)
shapiro.test(step_model$residuals)

plot(
  pca$x, 
  xlab = "PC1", 
  ylab = "PC2", 
  main = "PCA Environment Plot", 
  pch = as.numeric(data$pch), 
  col = c("green3","forestgreen", "goldenrod")[as.factor(data$fac)], cex = 2, bg = "goldenrod"
)
text(pca$x, labels = rownames(pca$x), pos = 4, cex = 0.3)

## Regressão linear "lm"

model <- lm(log(CS) ~ ., data=bio_data)
summary (model)
step_model <- step(model, direction = "both")
summary(step_model)

model_reduc <- lm(log(CS) ~ bio4+bio8+bio11+bio12+bio16+bio17, data=bio_data)
summary(model_reduc)
vif(model_reduc)

plot(model_reduc, pch = as.numeric(data$pch), col = as.factor(data$fac))
plot(model_reduc$fitted.values~model_reduc$residuals, pch = as.numeric(data$pch), col = as.factor(data$fac))
abline(model_reduc)
shapiro.test(model_reduc$residuals)

# Carregar pacotes necessários
#install.packages("mpm")
library(mpmcorrelogram)
require(car)
require(AICcmodavg)
require(reshape2)
require(ape)
require(vegan)

# Ajuste do modelo linear múltiplo com as variáveis selecionadas
fit <- lm(log(CS) ~ bio4+bio8+bio11+bio12+bio16+bio17, data = data)
vif(fit) # Verificação de multicolinearidade

# Modelos lineares simples para as variáveis selecionadas
fit1 <- lm(log(CS) ~ bio12, data = data)
summary(fit1)

fit2 <- lm(log(CS) ~ bio8, data = data)
summary(fit2)

fit3 <- lm(log(CS) ~ bio17, data = data)
summary(fit3)

fit4 <- lm(log(CS) ~ bio16, data = data)
summary(fit4)


# Comparação de modelos com AIC
models <- list(fit1, fit2, fit3, fit4)
model.names <- c('bio12', 'bio8', 'bio17','bio16')
aictab(cand.set = models, modnames = model.names)

# Aplicar a função vifstep() para remover automaticamente variáveis com VIF alto
names(bio_data)
bioclim <- as.matrix(bio_data)

vifstep_result <- vifstep(bioclim,th = 4, keep = "bio12")
vifstep_result

# "bio2","bio3","bio8","bio12","bio15","bio18","npp", "soilmoist"

# Ajustar um modelo apenas com as variáveis selecionadas
selected_vars <- vifstep_result@results$Variables
model_reduced <- lm(log(CS) ~ ., data = bio_data[, selected_vars])
summary(model_reduced)
step_model <- step(model_reduced, direction = "both")
summary(step_model) #bio8 e bio12
vif(step_model)

### Após testar com modelos lineares convencionais, é hora de selecionar a PGLS -> Best phylo ################

models_phy1 <- pgls.SEy(CS ~ bio2+bio3+bio8+bio12+bio15+bio18+npp+soilmoist ,data=bio_data,
                        tree=phy1)
summary(models_phy1)

models_phy2 <- pgls.SEy(CS ~  bio2+bio3+bio8+bio12+bio15+bio18+npp+soilmoist ,data=bio_data,
                        tree=phy2)
summary(models_phy2)

models_phy3 <- pgls.SEy(CS ~  bio2+bio3+bio8+bio12+bio15+bio18+npp+soilmoist ,data=bio_data,
                        tree=phy3)
summary(models_phy3)

# Compbio_data# Comparar os valores de AIC
aic_phy1 <- AIC(models_phy1)
aic_phy2 <- AIC(models_phy2)
aic_phy3 <- AIC(models_phy3)

# Imprimir os AICs
print(c("AIC - Lima et al (2018)" = aic_phy1, "AIC - Wright et al (2015)" = aic_phy2, "AIC - Upham et al (2019)" = aic_phy3))

# Comparar os valores de log-likelihood (opcional)
logLik_phy1 <- logLik(models_phy1)
logLik_phy2 <- logLik(models_phy2)
logLik_phy3 <- logLik(models_phy3)

print(c("Log-likelihood - Lima et al (2018)" = logLik_phy1, "Log-likelihood - Wright et al (2015)" = logLik_phy2, "Log-likelihood - Upham et al (2019)" = logLik_phy3))

# Modelos univariados
model_bio2 <- pgls.SEy(CS ~ bio2, data = bio_data, tree = phy1)
model_bio3 <- pgls.SEy(CS ~ bio3, data = bio_data, tree = phy1)
model_bio8 <- pgls.SEy(CS ~ bio8, data = bio_data, tree = phy1)
model_bio12 <- pgls.SEy(CS ~ bio12, data = bio_data, tree = phy1)
model_bio18 <- pgls.SEy(CS ~ bio18, data = bio_data, tree = phy1)
model_soil <- pgls.SEy(CS ~ soilmoist, data = bio_data, tree = phy1)

# Resumos dos modelos
summary(model_bio2)
summary(model_bio3)#
summary(model_bio8)
summary(model_bio12)#
summary(model_bio18)#
summary(model_soil)#

model_bio2 <- pgls.SEy(CS ~ bio2 + pcnm1 + pcnm5, data = bio_data, tree = phy1)
model_bio3 <- pgls.SEy(CS ~ bio3 + pcnm1 + pcnm5, data = bio_data, tree = phy1)
model_bio10 <- pgls.SEy(CS ~ bio10 + pcnm1 + pcnm5, data = bio_data, tree = phy1)
model_bio12 <- pgls.SEy(CS ~ bio12 + pcnm1 + pcnm5, data = bio_data, tree = phy1)
model_bio18 <- pgls.SEy(CS ~ bio18 + pcnm1 + pcnm5, data = bio_data, tree = phy1)
model_soil <- pgls.SEy(CS ~ soilmoist + pcnm1 + pcnm5, data = bio_data, tree = phy1)

# Resumos dos modelos
summary(model_bio2)
summary(model_bio3)
summary(model_bio10)#
summary(model_bio12)#
summary(model_bio18)#
summary(model_soil)

# Carregar pacotes necessários
library(ggplot2)
#install.packages("caret")
library(caret)

# 1. Comparação Gráfica dos Modelos
# Criar um dataframe com os dados observados e preditos para cada modelo
predictions_phy1 <- predict(models_phy1)
predictions_phy2 <- predict(models_phy2)
predictions_phy3 <- predict(models_phy3)

# Plotar os resíduos dos modelos
par(mfrow = c(3, 1)) # Organizar os gráficos em uma coluna
plot(residuals(models_phy1), main = "Resíduos - Lima (2018)", ylab = "Resíduos", xlab = "Observações")
abline(h = 0, col = "red")
plot(residuals(models_phy2), main = "Resíduos - Wright (2015)", ylab = "Resíduos", xlab = "Observações")
abline(h = 0, col = "red")
plot(residuals(models_phy3), main = "Resíduos - Upham (2019)", ylab = "Resíduos", xlab = "Observações")
abline(h = 0, col = "red")


# 2

# Analisar a normalidade dos resíduos com um histograma e QQ-plot
par(mfrow = c(3, 2)) # Organizar os gráficos
hist(residuals(models_phy1), main = "Histograma dos Resíduos - Lima (2018)", xlab = "Resíduos")
qqnorm(residuals(models_phy1))
qqline(residuals(models_phy1), col = "red")

hist(residuals(models_phy2), main = "Histograma dos Resíduos - Wright (2015)", xlab = "Resíduos")
qqnorm(residuals(models_phy2))
qqline(residuals(models_phy2), col = "red")

hist(residuals(models_phy3), main = "Histograma dos Resíduos - Upham (2019)", xlab = "Resíduos")
qqnorm(residuals(models_phy3))
qqline(residuals(models_phy3), col = "red")

# Modelo Browniano
model_brownian <- pgls.SEy(CS ~ bio3 + bio10 + bio12 + bio18 + soilmoist, data = bio_data, tree = phy1, corClass=corBrownian)
library(ape)

bio_data$sp <- as.factor(data$sp_ives)

# Modelo Ornstein-Uhlenbeck
model_ou <- pgls.SEy(CS ~ bio3 + bio10 + bio12 + bio18 + soilmoist, data = bio_data, tree = phy1, corClass=corMartins)

model_pagel <- pgls.SEy(CS ~ bio3 + bio10 + bio12 + bio18 + soilmoist, data = bio_data, tree = phy1, corClass=corPagel)

# AIC dos modelos
aic_brownian <- AIC(model_brownian)
aic_ou <- AIC(model_ou)
aic_pagel <- AIC(model_pagel)

# Log-likelihood dos modelos
logLik_brownian <- logLik(model_brownian)
logLik_ou <- logLik(model_ou)
logLik_pagel <- logLik(model_pagel)

# Imprimir os resultados
print(c("AIC - Modelo Browniano" = aic_brownian, "AIC - Modelo OU" = aic_ou, "AIC - Modelo Pagel" = aic_pagel))
print(c("Log-likelihood - Modelo Browniano" = logLik_brownian, "Log-likelihood - Modelo OU" = logLik_ou, "Log-likelihood - Modelo Pagel" = logLik_pagel))

# Pagel's !!!
c("bio4", "bio10", "bio12", "bio13", "bio14", "bio16", "bio17", "bio18", "humid", "soilmoist")

