### Testes Geomorph (Adams & Otárola, 2013)
#devtools::install_github("akiopteryx/lambda")
library(LaMBDA)
library(geomorph)
library(phytools)
library(picante)


setwd("C:/Users/Lourenço/OneDrive/UFSM/LabMastozoo/Sapajus/Data")

data <- read.table("meanspec.txt", sep = ",",header = TRUE)
head (data)
land <- readland.tps ("tps/avgspec.tps")

Y.gpa <- gpagen(land)
links = read.table("tps/link.txt")
plotAllSpecimens(Y.gpa$coords, links = links)
y<-two.d.array(Y.gpa$coords)
dev.off()
#resultados_lasec <- lasec(coord.data = y, n.dim = 2, iter = 1000, show.progress = TRUE)
#resultados_lasec

fit <-procD.lm(Y.gpa$coords ~ Y.gpa$Csize * as.factor(data$biome), iter = 9999)
summary (fit)


phy <-read.tree("trees/MCC_Lima_calibrated.tre")
phy <- drop.tip(phy,c("Sapajus_macrocephalus", "Cebus_albifrons", "Cebus_capucinus", "Sapajus_flavius"))
class(phy)

dimnames(Y.gpa$coords)[[3]] <- data$spec
names(Y.gpa$Csize) <- data$spec
CS <- as.data.frame(data$LogCS)
rownames(CS) <- data$spec
CS <- as.matrix(CS)


phy.sig <- physignal(Y.gpa$Csize,phy,iter=9999)
phy.sig
phy.sig$phy.signal

phy.sig <- physignal(Y.gpa$coords,phy,iter=9999)
phy.sig
phy.sig$phy.signal

effect_size <- physignal.z(Y.gpa$coords, phy)
print(effect_size)

hist(phy.sig$random.K)
phy.sig$pvalue

rownames(y) <- data$spec
lambda_result <- phylosig(phy, Y.gpa$Csize, method = "lambda")
print(lambda_result)

library(geiger)

lambda_result <- phylosig(phy, y, method = "lambda")
print(lambda_result)

K_result <- phylosig(phy, Y.gpa$Csize, method = "K", nsim =9999)
print(K_result)

K_result <- phylosig(phy, y, method = "K", nsim =9999)
print(K_result)

citation("geomorph")
loadedNamespaces()

plotOutliers(Y.gpa$coords, groups = NULL, inspect.outliers = TRUE)

globalIntegration (Y.gpa$coords, ShowPlot = TRUE)

land.gp <- factor(c(rep("Corpo da Mandíbula", 4), rep("Coronoide", 6), rep("Côndilo", 3),rep("Goníaco", 3), rep("Corpo da Mandíbula", 2)))

EMR <- compare.multi.evol.rates(A = Y.gpa$coords, gp = land.gp, Subset = TRUE, phy = phy)
summary(EMR)
plot(EMR)

library(ggplot2)

set.seed(123) 
taxas_simuladas <- data.frame(
  Grupo = rep(c("Côndilo", "Coronoide", "Corpo da Mandíbula", "Goníaco"), each = 100),
  Taxa = c(rnorm(100, mean = 9.34e-05, sd = 1e-06),
           rnorm(100, mean = 0.0001016, sd = 1e-06),
           rnorm(100, mean = 8.86e-05, sd = 1e-06),
           rnorm(100, mean = 0.0002231, sd = 1e-06))
)

ggplot(taxas_simuladas, aes(x = Grupo, y = Taxa, fill = Grupo)) +
  geom_boxplot() +
  theme_minimal() +
  labs(title = "Distribuição das Taxas de Evolução por Grupo", y = "Taxa de Evolução", x = "Grupo") +
  theme(legend.position = "none")


IT <- integration.test(Y.gpa$coords, partition.gp = land.gp)
summary(IT) 

MT <- modularity.test(Y.gpa$coords, land.gp, CI = FALSE)
summary(MT)
plot(MT) 

IT <- phylo.integration(Y.gpa$coords, partition.gp = land.gp, phy = phy)
summary(IT)

MT <- phylo.modularity(Y.gpa$coords, partition.gp = land.gp, phy = phy)
summary(MT)

PCA.w.phylo <- gm.prcomp(Y.gpa$coords, phy = phy)
summary(PCA.w.phylo)
plot(PCA.w.phylo, phylo = TRUE, main = "PCA.w.phylo")

phylo.PCA <- gm.prcomp(Y.gpa$coords, phy = phy, GLS = TRUE)
summary(phylo.PCA)
plot(phylo.PCA, phylo = TRUE, main = "phylo PCA")

phylo.tPCA <- gm.prcomp(Y.gpa$coords, phy = phy, 
                        GLS = TRUE, transform = TRUE)
summary(phylo.tPCA)
plot(phylo.tPCA, phylo = TRUE, main = "phylo PCA")

# OLS method (rotation of PCA)
PaCA.ols <- gm.prcomp(Y.gpa$coords, phy = phy, 
                      align.to.phy = TRUE)
summary(PaCA.ols)
plot(PaCA.ols, phylo = TRUE, main = "PaCA using OLS")

# GLS method (rotation of Phylogenetic PCA)
PaCA.gls <- gm.prcomp(Y.gpa$coords, phy = phy, 
                      align.to.phy = TRUE, GLS = TRUE)
summary(PaCA.gls)
plot(PaCA.gls, phylo = TRUE, main = "PaCA using GLS")

# GLS method (rotation of Phylogenetic PCA with transformed data)
PaCA.gls <- gm.prcomp(Y.gpa$coords, phy = phy, 
                      align.to.phy = TRUE, GLS = TRUE, transform = TRUE)
summary(PaCA.gls)
plot(PaCA.gls, phylo = TRUE, 
     main = "PaCA using GLS and transformed projection")


### Advanced Plotting
gps <- as.factor(data$biome) # Two random groups
par(mar=c(2, 2, 2, 2))
plot(PaCA.ols, pch=22, cex = 1.5, bg = gps, phylo = TRUE) 

# Modify options as desired
#  Add things as desired using standard R plotting
text(par()$usr[1], 0.1*par()$usr[3], labels = "PC1 - 45.64%", 
     pos = 4, font = 2)
text(0, 0.95*par()$usr[4], labels = "PC2 - 18.80%", pos = 4, font = 2)
legend("topleft", pch=22, pt.bg = unique(gps), legend = levels(gps))

## 3D plot with a phylogeny and time on the z-axis
plot(PCA.w.phylo, time.plot = TRUE)
plot(PCA.w.phylo, time.plot = TRUE, bg = "red", 
     phylo.par = list(tip.labels = TRUE, 
                      tip.txt.cex = 0.5, edge.color = "blue", edge.width = 2))

##### With averaged specimens

data <- read.csv("avglocsex.csv", sep = ";")
head (data)
land <- readland.tps ("tps/avglocsex.tps")

phy <- read.tree("trees/AVGLOCSEX_Lima.tre")

Y.gpa <- gpagen(land)
y<-two.d.array(Y.gpa$coords)
               
coords <- Y.gpa$coords
size <- log(Y.gpa$Csize)
sex <- as.factor(data$sex)
biome <- as.factor(data$biome)

dimnames(Y.gpa$coords)[[3]] <- data$sp_ives
names(Y.gpa$Csize) <- data$sp_ives
CS <- as.data.frame(data$LogCS)
rownames(CS) <- data$sp_ives

phy.sig <- physignal(Y.gpa$Csize,phy,iter=9999)
phy.sig
phy.sig$phy.signal

phy.sig <- physignal(Y.gpa$coords,phy,iter=9999)
phy.sig
phy.sig$phy.signal

effect_size <- physignal.z(Y.gpa$coords, phy)
print(effect_size)

hist(phy.sig$random.K)
phy.sig$pvalue

rownames(y) <- data$sp_ives
lambda_result <- phylosig(phy, Y.gpa$Csize, method = "lambda")
print(lambda_result)

library(geiger)

lambda_result <- phylosig(phy, y, method = "lambda")
print(lambda_result)

K_result <- phylosig(phy, Y.gpa$Csize, method = "K", nsim =9999)
print(K_result)

K_result <- phylosig(phy, y, method = "K", nsim =9999)
print(K_result)


globalIntegration (Y.gpa$coords, ShowPlot = TRUE)

land.gp <- factor(c(rep("Corpo da Mandíbula", 4), rep("Coronoide", 6), rep("Côndilo", 3),rep("Goníaco", 3), rep("Corpo da Mandíbula", 2)))

EMR <- compare.multi.evol.rates(A = Y.gpa$coords, gp = land.gp, Subset = TRUE, phy = phy)
summary(EMR)
plot(EMR)

library(ggplot2)

set.seed(123) 
taxas_simuladas <- data.frame(
  Grupo = rep(c("Côndilo", "Coronoide", "Corpo da Mandíbula", "Goníaco"), each = 100),
  Taxa = c(rnorm(100, mean = 9.34e-05, sd = 1e-06),
           rnorm(100, mean = 0.0001016, sd = 1e-06),
           rnorm(100, mean = 8.86e-05, sd = 1e-06),
           rnorm(100, mean = 0.0002231, sd = 1e-06))
)

ggplot(taxas_simuladas, aes(x = Grupo, y = Taxa, fill = Grupo)) +
  geom_boxplot() +
  theme_minimal() +
  labs(title = "Distribuição das Taxas de Evolução por Grupo", y = "Taxa de Evolução", x = "Grupo") +
  theme(legend.position = "none")


IT <- integration.test(Y.gpa$coords, partition.gp = land.gp)
summary(IT) 

MT <- modularity.test(Y.gpa$coords, land.gp, CI = FALSE)
summary(MT)
plot(MT) 

IT <- phylo.integration(Y.gpa$coords, partition.gp = land.gp, phy = phy)
summary(IT)

MT <- phylo.modularity(Y.gpa$coords, partition.gp = land.gp, phy = phy)
summary(MT)

### T.A, groups = Biomes

fit <- lm.rrpp(y ~ biome * size * sex, data = data, iter =999)
summary(fit)

TA <- trajectory.analysis(fit, groups = data$biome, 
                          traj.pts = data$sex, print.progress = FALSE)

# Magnitude difference (absolute difference between path distances)
summary(TA, attribute = "MD", show.trajectories = TRUE) 

# Correlations (angles) between trajectories
summary(TA, attribute = "TC", angle.type = "deg", show.trajectories = TRUE) 

# No shape differences between vectors
summary(TA, attribute = "SD") 

# Retain results
TA.summary <- summary(TA, attribute = "MD")
TA.summary$summary.table

# Plot results
dev.off()
biome_levels <- as.numeric(as.factor(data$biome))  # Converte biomas para números
sex_levels <- as.numeric(as.factor(data$sex))

TP <- plot(TA, pch = c(21,25)[biome_levels], bg = c("palegreen3", "khaki1")[biome_levels], cex = 3, col = "darkgrey")

add.trajectories(
  TP,
  traj.pch = c(21,25),         
  traj.col = 1,          
  traj.lty = 1,          
  traj.lwd = 2,          
  traj.bg = 1,           
  start.bg = "lightblue",
  end.bg = "lightcoral", 
  traj.cex = 2.5)

legend("topright", legend = unique(data$biome), 
       pch = c(21,25), pt.bg = c("palegreen3","khaki1"), title = "Biome", cex = 1.2)

legend("bottomright", legend = unique(data$sex), 
       pt.bg = c("lightblue", "lightcoral"), pch = c(23), title = "Sex", cex= 1.5)

### T.A, groups = Sexes

fit <- lm.rrpp(y ~ sex * size * biome, data = data, iter =999)
summary(fit)

TA <- trajectory.analysis(fit, groups = data$sex, 
                          traj.pts = data$biome, print.progress = FALSE)

# Magnitude difference (absolute difference between path distances)
summary(TA, attribute = "MD", show.trajectories = TRUE) 

# Correlations (angles) between trajectories
summary(TA, attribute = "TC", angle.type = "deg", show.trajectories = TRUE) 

# No shape differences between vectors
summary(TA, attribute = "SD") 

# Retain results
TA.summary <- summary(TA, attribute = "MD")
TA.summary$summary.table

# Plot results

biome_levels <- as.numeric(as.factor(data$biome))  # Forest = 1; Savanna = 2
sex_levels <- as.numeric(as.factor(data$sex)) # 1 = F; 2 = M

TP <- plot(TA, pch = c(23,24)[sex_levels], bg = c("lightcoral", "lightblue")[sex_levels], cex = 3, col = "darkgrey")

add.trajectories(
  TP,
  traj.pch = c(23,24),         
  traj.col = 1,          
  traj.lty = 1,          
  traj.lwd = 2,          
  traj.bg = 1,           
  start.bg = "palegreen3",
  end.bg = "khaki1", 
  traj.cex = 2.5)

legend("topright", legend = unique(data$biome), 
       pch = c(21), pt.bg = c("palegreen3","khaki1"), title = "Biome", cex = 1.2)

legend("bottomright", legend = unique(data$sex), 
       pt.bg = c("lightblue", "lightcoral"), pch = unique(c(23,24)[sex_levels]), title = "Sex", cex= 1.5)



gdf <- geomorph.data.frame(shape = Y.gpa$coords, size = log(Y.gpa$Csize), phy = phy, biome = as.factor(data$biome))

pgls <- procD.pgls(shape ~ biome * size, phy = phy, lambda = 0.8,seed = NULL, SS.type = "I",data = gdf)
predict(pgls)
summary(pgls)

colors <- c("forestgreen", "goldenrod")[as.numeric(as.factor(gdf$biome))]
bio_pch <- c(15, 17)[as.numeric(as.factor(gdf$biome))]

plot(pgls, 
     type = "PC", 
     reg.type = "RegScore", 
     predictor = as.numeric(gdf$size), 
     pch = bio_pch, 
     col = colors)

legend("topright",legend = levels(as.factor(gdf$biome)), 
       col = c("forestgreen", "goldenrod"), 
       pch = unique(bio_pch))

attributes(pgls) # Note the PGLS object
attributes(pgls$PGLS) # PGLS details embedded within PGLS object
pgls$LM$Pcov # the projection matrix derived from the 

# phylogenetic covariance matrix
pgls$pgls.fitted # the PGLS fitted values 
pgls$GM$pgls.fitted # The same fitted values, in a 3D array

# Changing lambda value

pgls2 <- procD.pgls(shape ~ size, phy = phy, lambda = 0.5, data = gdf)
anova(pgls2)
summary(pgls2)

PCA <- gm.prcomp(Y.gpa$coords)
PCS <- PCA$x

anova_PC1 <- aov(PCS [,1] ~ biome * sex, data = data)
summary(anova_PC1)

anova_PC2 <- aov(PCS [,2] ~ biome * sex, data = data)
summary(anova_PC2)

anova_PC3 <- aov(PCS[,3] ~ biome * sex ,data= data)
summary(anova_PC3)

anova_PC4 <- aov(PCS[,4] ~  biome * sex, data = data)
summary(anova_PC4)

anova_PCS <- aov(PCS[,5] ~  biome * sex, data = data)
summary(anova_PCS)

anova_PCS <- aov(PCS[,6] ~  biome * sex, data = data)
summary(anova_PCS)

anova_PCS <- aov(PCS[,7] ~  biome * sex, data = data)
summary(anova_PCS)

anova_PCS <- aov(PCS[,8] ~  biome * sex, data = data)
summary(anova_PCS)

anova_PCS <- aov(PCS[,9] ~  biome * sex, data = data)
summary(anova_PCS)

