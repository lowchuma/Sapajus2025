#devtools::install_github("geomorphR/geomorph", ref = "Develop", build_vignettes = TRUE)
library(ggplot2)
library(geomorph)
library(phytools)
library(vegan)
library(geiger)
library(phangorn)
library(tidyverse)
library (RColorBrewer)
library(usdm)
library(spaMM)
library(foreach)
library(doParallel)
setwd("C:/Users/Lourenço/OneDrive/UFSM/LabMastozoo/SAPAJUS/Data - From forests to savannas")

#############gm analises
sapajus<-readland.tps("Datasets/tps/avglocsex.tps", specID="ID")
gpa<-gpagen(sapajus)
dim(sapajus)

data<-read.csv("Planilhas/avglocsex.csv", sep=";")
dim(data)
names(data)

#link<-define.links(gpa$consensus, ptsize = 2, links = NULL)
#write.table(link,file="link.txt")

link <- read.table("Datasets/tps/link.txt")

gdf<-geomorph.data.frame(shape=gpa$coords, size = log(gpa$Csize), biome = as.factor(data$biome), facgp = as.factor(data$fac), spp = as.factor(data$sp))

drawinglandmark<-readland.tps("Datasets/outline2/outline.tps")
outline<-read.table("Datasets/outline2/outline.txt", header=FALSE)
summary(drawinglandmark)
mshape<-mshape(gpa$coords)
summary(mshape)
summary(outline)
Sapajusoutline<-warpRefOutline(file = "Datasets/outline2/outline.txt", drawinglandmark[,,1],mshape)#run dev.off() in case of error message ##set outline configuration

dev.off()

########################## PARVAR ##############################################
names(data)
coords<-data[,4:5]

geographicaldistance <- as.matrix(dist(coords))
pcnmcoord<-prcomp(geographicaldistance)
summary(pcnmcoord)
pcnm<-pcnmcoord$x
rownames(pcnm) <-data$sp_ives

#ou

pcnmcoord<- pcnm(dist(coords))
pcnm<-as.matrix(pcnmcoord$vectors)
rownames(pcnm) <-data$sp_ives

#write.table(pcnm, "pcnm.xls", row.names= TRUE, col.names = TRUE, sep = " ")

### Environment
names(data)
bioclim <- data[,c(11:32)]
bioclim <- scale(bioclim,center=T,scale=T)
head(bioclim)

vifstep_result <-vifstep (bioclim, keep = c("bio12"))
print(vifstep_result)
vifstep_result@results$Variables # "bio2", "bio3", "bio8", "bio12","bio14","bio15","bio18","bio19","npp", "humid","soilmoist"

data[,c(11:32)] <- scale(data[,c(11:32)])
rownames(bioclim) <-data$sp_ives

envpca<-prcomp(bioclim)
summary(envpca)
env<-as.matrix(envpca$x[,1:6])
rownames(env) <-data$sp_ives

### Shape (Y)
# Alometria simples ##Regressão multivariada
alometria<-procD.lm(gpa$coords~log(gpa$Csize),iter=999, RRPP = TRUE)
summary(alometria)
plot(alometria,type="regression",predictor=log(gpa$Csize),
     reg.type = "RegScore", pch=19, xlab="log(Size)")

shape.resid<-
  arrayspecs(alometria$residuals,p=dim(gpa$coords)[1],k=dim(gpa$coords)[2]) ## extrai dos resultados os residuos da associação de forma e tamanho #size adjusted residuals
adj.shape<-shape.resid+array(gpa$consensus, dim(shape.resid)) # allometry-free ## variável de forma livre de efeito de alometria
shape.2d <- two.d.array(adj.shape)
dimnames(shape.2d)[[1]] <- paste(as.factor(data$sp_ives))

pcamorph<-gm.prcomp(gpa$coords)
summary(pcamorph)
plot(pcamorph, pch = as.numeric(data$pch), col = as.factor(data$fac))
pcsgm<-as.matrix(pcamorph$x[,1:16]) #95%
rownames(pcsgm) <-data$sp_ives

# Size
shapiro.test(data$CS)
shapiro.test(log(data$CS))

size<-as.numeric(log(data$CS))
size <-as.matrix(size)
names(size)<-data$sp_ives
shapiro.test((size))

### Selection

#install.packages("packfor", repos="http://R-Forge.R-project.org")
library(packfor)
forward.sel(pcsgm,pcnm) ## 31, 5, 2, 7 /// 51,80,4,54,79
forward.sel(pcsgm,env) ## 3,2, 16
forward.sel(size,pcnm) ## 1,5,13 /// 69,27,24
forward.sel(size,env) ## 2,15

pcnm_axes <- pcnm[,c(2,31,1)]
env_axes <- env[,c(3,2)]

### Phylogeny
tree<-read.tree("Datasets/trees/AVGLOCSEX_Lima.tre")
rownames(data) <- data$sp_ives
name.check(tree,setNames(size,rownames(data)))
plot(tree, cex = 0.5)
axisPhylo()

cov<-vcv.phylo(tree)
phylo<-prcomp(cov)
summary(phylo)
phy<-as.matrix(phylo$x[,1:3])

## Elements

#size
pcnm1<- pcnm[,1]
pcnm5<-pcnm[,5]

#shape
pcnm31 <- pcnm[,31]
pcnm7 <- pcnm[,7]
pcnm2 <- pcnm[,2]

### Plotting

# Verificando nomes de linha
all.equal(rownames(pcnm), rownames(data))       # Deve retornar TRUE
all.equal(rownames(bioclim), rownames(data))   # Deve retornar TRUE
all.equal(rownames(pcsgm), rownames(data))     # Deve retornar TRUE
all.equal(rownames(env), rownames(data))       # Deve retornar TRUE

parvar<-varpart(pcsgm,~phy,~size,~env,~pcnm)

plot(parvar,cex=1.2,Xnames=c("Phylogeny","Spatial Structure","Size","Environment"),bg=c("yellow", "deeppink","navy","green3"))

y<-pcsgm
sf<-pcnm
env <-env

### for all parvar's variables #####

p_values <- c()

rda.result <- rda(y ~ phy)
anova_result <- anova(rda.result, step = 1000, perm.max = 2000)
p_values <- c(p_values, anova_result$`Pr(>F)`[1])

rda.result_alo <- rda(y ~ size)
anova_result <- anova(rda.result_alo, step = 1000, perm.max = 2000)
p_values <- c(p_values, anova_result$`Pr(>F)`[1])

rda.result <- rda(y ~ env)
anova_result <- anova(rda.result, step = 1000, perm.max = 2000)
p_values <- c(p_values, anova_result$`Pr(>F)`[1])

rda.result <- rda(y ~ sf)
anova_result <- anova(rda.result, step = 1000, perm.max = 2000)
p_values <- c(p_values, anova_result$`Pr(>F)`[1])

rda.result <- rda(y ~ phy + size)
anova_result <- anova(rda.result, step = 1000, perm.max = 2000)
p_values <- c(p_values, anova_result$`Pr(>F)`[1])

rda.result <- rda(y ~ phy + env)
anova_result <- anova(rda.result, step = 1000, perm.max = 2000)
p_values <- c(p_values, anova_result$`Pr(>F)`[1])

rda.result <- rda(y ~ phy + sf)
anova_result <- anova(rda.result, step = 1000, perm.max = 2000)
p_values <- c(p_values, anova_result$`Pr(>F)`[1])

rda.result <- rda(y ~ size + env)
anova_result <- anova(rda.result, step = 1000, perm.max = 2000)
p_values <- c(p_values, anova_result$`Pr(>F)`[1])

rda.result <- rda(y ~ size + sf)
anova_result <- anova(rda.result, step = 1000, perm.max = 2000)
p_values <- c(p_values, anova_result$`Pr(>F)`[1])

rda.result <- rda(y ~ env + sf)
anova_result <- anova(rda.result, step = 1000, perm.max = 2000)
p_values <- c(p_values, anova_result$`Pr(>F)`[1])

rda.result <- rda(y ~ phy + size + env)
anova_result <- anova(rda.result, step = 1000, perm.max = 2000)
p_values <- c(p_values, anova_result$`Pr(>F)`[1])

rda.result <- rda(y ~ phy + size + sf)
anova_result <- anova(rda.result, step = 1000, perm.max = 2000)
p_values <- c(p_values, anova_result$`Pr(>F)`[1])

rda.result <- rda(y ~ phy + env + sf)
anova_result <- anova(rda.result, step = 1000, perm.max = 2000)
p_values <- c(p_values, anova_result$`Pr(>F)`[1])

rda.result <- rda(y ~ size + env + sf)
anova_result <- anova(rda.result, step = 1000, perm.max = 2000)
p_values <- c(p_values, anova_result$`Pr(>F)`[1])

rda.result <- rda(y ~ env + size + phy + sf)
anova_result <- anova(rda.result, step = 1000, perm.max = 2000)
p_values <- c(p_values, anova_result$`Pr(>F)`[1])


### RDA condicionada

rda.result <- rda(y ~ phy + Condition(size) + Condition(env) + Condition(sf))
anova_result <- anova(rda.result, step = 1000, perm.max = 2000)
p_values <- c(p_values, anova_result$`Pr(>F)`[1])

rda.result <- rda(y ~ size + Condition(phy) + Condition(env) + Condition(sf))
anova_result <- anova(rda.result, step = 1000, perm.max = 2000)
p_values <- c(p_values, anova_result$`Pr(>F)`[1])

rda.result <- rda(y ~ env + Condition(phy) + Condition(size) + Condition(sf))
anova_result <- anova(rda.result, step = 1000, perm.max = 2000)
p_values <- c(p_values, anova_result$`Pr(>F)`[1])
plot(rda.result)

rda.result <- rda(y ~ sf + Condition(phy) + Condition(size) + Condition(env))
anova_result <- anova(rda.result, step = 1000, perm.max = 2000)
p_values <- c(p_values, anova_result$`Pr(>F)`[1])

results_df <- data.frame(
  p_value = p_values
)

print(results_df)

parvar$p_values <- results_df

fract <- parvar$part$fract
indfract <- parvar$part$indfract

results_fract <- data.frame(
  Partition = rownames(fract),
  Df = fract$Df,
  R_square = fract$R.square,
  Adj_R_square = fract$Adj.R.square,
  Testable = fract$Testable
)

results_indfract <- data.frame(
  Partition = rownames(indfract),
  Df = indfract$Df,
  R_square = indfract$R.square,
  Adj_R_square = indfract$Adj.R.square,
  Testable = indfract$Testable
)

final_results <- rbind(results_fract, results_indfract)

mapping <- c("X1" = "Phylogeny", 
             "X2" = "Size", 
             "X3" = "Environment", 
             "X4" = "Spatial Structure")

final_results$Partition <- gsub("X1", mapping["X1"], final_results$Partition)
final_results$Partition <- gsub("X2", mapping["X2"], final_results$Partition)
final_results$Partition <- gsub("X3", mapping["X3"], final_results$Partition)
final_results$Partition <- gsub("X4", mapping["X4"], final_results$Partition)
print(final_results)

p_values_complete <- c(p_values, rep(NA, nrow(final_results) - length(p_values)))
final_results$p_value <- p_values_complete
print(final_results)

write.csv(final_results, "parvar.csv", row.names = FALSE)

#### PGLS R #####

# Para o tamanho:

### PGLS.SEy
ecor <- as.factor(data$fac)
names(ecor) <- data$sp_ives

fit.allo <- pgls.SEy(y ~ size*ecor, data=data,tree=tree, corClass = corPagel)
summary(fit.allo)

fit<-pgls.SEy(size~bio1,data=data,tree=tree)
summary(fit)

fit<-pgls.SEy(size~bio2,data=data,tree=tree)
summary(fit) 

fit<-pgls.SEy(size~bio3,data=data,tree=tree)
summary(fit)#

fit<-pgls.SEy(size~bio4,data=data,tree=tree)
summary(fit)#

fit<-pgls.SEy(size~bio5,data=data,tree=tree)
summary(fit)

fit<-pgls.SEy(size~bio6,data=data,tree=tree)
summary(fit)

fit<-pgls.SEy(size~bio7,data=data,tree=tree)
summary(fit)

fit<-pgls.SEy(size~bio8,data=data,tree=tree)
summary(fit) #

fit<-pgls.SEy(size~bio9,data=data,tree=tree)
summary(fit)

fit<-pgls.SEy(size~bio10,data=data,tree=tree)
summary(fit)#####

fit<-pgls.SEy(size~bio11,data=data,tree=tree)
summary(fit)

fit<-pgls.SEy(size~bio12,data=data,tree=tree)
summary(fit) ##

fit<-pgls.SEy(size~bio13,data=data,tree=tree)
summary(fit) #

fit<-pgls.SEy(size~bio14,data=data,tree=tree)
summary(fit)

fit<-pgls.SEy(size~bio15,data=data,tree=tree)
summary(fit)

fit<-pgls.SEy(size~bio16,data=data,tree=tree)
summary(fit) #

fit<-pgls.SEy(size~bio17,data=data,tree=tree)
summary(fit)

fit<-pgls.SEy(size~bio18,data=data,tree=tree)
summary(fit) #

fit<-pgls.SEy(size~bio19,data=data,tree=tree)
summary(fit) 

fit<-pgls.SEy(size~humid,data=data,tree=tree)
summary(fit)

fit<-pgls.SEy(size~npp,data=data,tree=tree)
summary(fit) 

fit<-pgls.SEy(size~soilmoist,data=data,tree=tree)
summary(fit) #

##repeat model with filters

fit<-pgls.SEy(size~bio1+pcnm1+pcnm5,data=data,tree=tree)
summary(fit)

fit<-pgls.SEy(size~bio2+pcnm1+pcnm5,data=data,tree=tree)
summary(fit) 

fit<-pgls.SEy(size~bio3+pcnm1+pcnm5,data=data,tree=tree)
summary(fit) 

fit<-pgls.SEy(size~bio4+pcnm1+pcnm5,data=data,tree=tree)
summary(fit) 

fit<-pgls.SEy(size~bio5+pcnm1+pcnm5,data=data,tree=tree)
summary(fit)

fit<-pgls.SEy(size~bio6+pcnm1+pcnm5,data=data,tree=tree)
summary(fit)

fit<-pgls.SEy(size~bio7+pcnm1+pcnm5,data=data,tree=tree)
summary(fit)

fit<-pgls.SEy(size~bio8+pcnm1+pcnm5,data=data,tree=tree)
summary(fit)

fit<-pgls.SEy(size~bio9+pcnm1+pcnm5,data=data,tree=tree)
summary(fit)

fit<-pgls.SEy(size~bio10+pcnm1+pcnm5,data=data,tree=tree)
summary(fit)#

fit<-pgls.SEy(size~bio11+pcnm1+pcnm5,data=data,tree=tree)
summary(fit)

fit<-pgls.SEy(size~bio12+pcnm1+pcnm5,data=data,tree=tree)
summary(fit) ##

fit<-pgls.SEy(size~bio13+pcnm1+pcnm5,data=data,tree=tree)
summary(fit) 

fit<-pgls.SEy(size~bio14+pcnm1+pcnm5,data=data,tree=tree)
summary(fit)

fit<-pgls.SEy(size~bio15+pcnm1+pcnm5,data=data,tree=tree)
summary(fit)

fit<-pgls.SEy(size~bio16+pcnm1+pcnm5,data=data,tree=tree)
summary(fit) 

fit<-pgls.SEy(size~bio17+pcnm1+pcnm5,data=data,tree=tree)
summary(fit)

fit<-pgls.SEy(size~bio18+pcnm1+pcnm5,data=data,tree=tree)
summary(fit)#

fit<-pgls.SEy(size~bio19+pcnm1+pcnm5,data=data,tree=tree)
summary(fit)

fit<-pgls.SEy(size~npp+pcnm1+pcnm5,data=data,tree=tree)
summary(fit)

fit<-pgls.SEy(size~humid+pcnm1+pcnm5,data=data,tree=tree)
summary(fit)

fit<-pgls.SEy(size~soilmoist+pcnm1+pcnm5,data=data,tree=tree)
summary(fit)

#c("bio10", "bio12", "bio18") - #Lima (2018)

#####pagel

fit<-pgls.SEy(size~bio1,data=data,tree=tree, corClass=corPagel)
summary(fit)
fit<-pgls.SEy(size~bio2,data=data,tree=tree, corClass=corPagel)
summary(fit)
fit<-pgls.SEy(size~bio3,data=data,tree=tree, corClass=corPagel)
summary(fit) 
fit<-pgls.SEy(size~bio4,data=data,tree=tree, corClass=corPagel)
summary(fit) #
fit<-pgls.SEy(size~bio5,data=data,tree=tree, corClass=corPagel)
summary(fit) 
fit<-pgls.SEy(size~bio6,data=data,tree=tree, corClass=corPagel)
summary(fit)
fit<-pgls.SEy(size~bio7,data=data,tree=tree, corClass=corPagel) 
summary(fit) 
fit<-pgls.SEy(size~bio8,data=data,tree=tree, corClass=corPagel)
summary(fit) 
fit<-pgls.SEy(size~bio9,data=data,tree=tree, corClass=corPagel)
summary(fit) 
fit<-pgls.SEy(size~bio10,data=data,tree=tree, corClass=corPagel)
summary(fit) ##
fit<-pgls.SEy(size~bio11,data=data,tree=tree, corClass=corPagel)
summary(fit)
fit<-pgls.SEy(size~bio12,data=data,tree=tree, corClass=corPagel)
summary(fit) ##
fit<-pgls.SEy(size~bio13,data=data,tree=tree, corClass=corPagel)
summary(fit) #
fit<-pgls.SEy(size~bio14,data=data,tree=tree, corClass=corPagel)
summary(fit) #
fit<-pgls.SEy(size~bio15,data=data,tree=tree, corClass=corPagel)
summary(fit)
fit<-pgls.SEy(size~bio16,data=data,tree=tree, corClass=corPagel)
summary(fit) ##
fit<-pgls.SEy(size~bio17,data=data,tree=tree, corClass=corPagel)
summary(fit) ##
fit<-pgls.SEy(size~bio18,data=data,tree=tree, corClass=corPagel)
summary(fit) ##
fit<-pgls.SEy(size~bio19,data=data,tree=tree, corClass=corPagel)
summary(fit) 
fit<-pgls.SEy(size~npp,data=data,tree=tree, corClass=corPagel)
summary(fit)
fit<-pgls.SEy(size~humid,data=data,tree=tree, corClass=corPagel)
summary(fit) #
fit<-pgls.SEy(size~soilmoist,data=data,tree=tree, corClass=corPagel)
summary(fit) #.

### 4 mais importantes (decrescente): bio12 > bio16 > bio17 > bio10 > bio18 > bio4

fit<-pgls.SEy(size~bio1+pcnm1+pcnm5,data=data,tree=tree, corClass=corPagel)
summary(fit)

fit<-pgls.SEy(size~bio2+pcnm1+pcnm5,data=data,tree=tree, corClass=corPagel)
summary(fit) 

fit<-pgls.SEy(size~bio3+pcnm1+pcnm5,data=data,tree=tree, corClass=corPagel)
summary(fit) 

fit<-pgls.SEy(size~bio4+pcnm1+pcnm5,data=data,tree=tree, corClass=corPagel)
summary(fit) 

fit<-pgls.SEy(size~bio5+pcnm1+pcnm5,data=data,tree=tree, corClass=corPagel)
summary(fit)

fit<-pgls.SEy(size~bio6+pcnm1+pcnm5,data=data,tree=tree, corClass=corPagel)
summary(fit)

fit<-pgls.SEy(size~bio7+pcnm1+pcnm5,data=data,tree=tree, corClass=corPagel)
summary(fit)

fit<-pgls.SEy(size~bio8+pcnm1+pcnm5,data=data,tree=tree, corClass=corPagel)
summary(fit)

fit<-pgls.SEy(size~bio9+pcnm1+pcnm5,data=data,tree=tree, corClass=corPagel)
summary(fit)#######

fit<-pgls.SEy(size~bio10+pcnm1+pcnm5,data=data,tree=tree, corClass=corPagel)
summary(fit)

fit<-pgls.SEy(size~bio11+pcnm1+pcnm5,data=data,tree=tree, corClass=corPagel)
summary(fit)

fit<-pgls.SEy(size~bio12+pcnm1+pcnm5,data=data,tree=tree, corClass=corPagel)
summary(fit) #########

fit<-pgls.SEy(size~bio13+pcnm1+pcnm5,data=data,tree=tree, corClass=corPagel)
summary(fit) 

fit<-pgls.SEy(size~bio14+pcnm1+pcnm5,data=data,tree=tree, corClass=corPagel)
summary(fit) 

fit<-pgls.SEy(size~bio15+pcnm1+pcnm5,data=data,tree=tree, corClass=corPagel)
summary(fit)##### error

fit<-pgls.SEy(size~bio16+pcnm1+pcnm5,data=data,tree=tree, corClass=corPagel)
summary(fit) 

fit<-pgls.SEy(size~bio17+pcnm1+pcnm5,data=data,tree=tree, corClass=corPagel)
summary(fit)#

fit<-pgls.SEy(size~bio18+pcnm1+pcnm5,data=data,tree=tree, corClass=corPagel)
summary(fit)#

fit<-pgls.SEy(size~bio19+pcnm1+pcnm5,data=data,tree=tree, corClass=corPagel)
summary(fit)

fit<-pgls.SEy(size~npp+pcnm1+pcnm5,data=data,tree=tree, corClass=corPagel)
summary(fit) ############## FALSE CONVERGENCE

fit<-pgls.SEy(size~humid+pcnm1+pcnm5,data=data,tree=tree, corClass=corPagel)
summary(fit)

fit<-pgls.SEy(size~soilmoist+pcnm1+pcnm5,data=data,tree=tree, corClass=corPagel)
summary(fit)

# 4 mais importantes (decrescente): bio12 > bio17 bio18 


clim <- as.data.frame(bioclim[,c("bio12", "bio17","bio18")])
vif(clim)

vifstep_result <- vifstep(bioclim,th = 2, keep = c("bio12", "bio17","bio18"))
print(vifstep_result)

vifstep_result@results$Variables 
names(data)
bio_data <- as.data.frame(scale(data[,11:32]))
size <- as.numeric(data$size)

rownames(bio_data) <-data$sp_ives
names(size) <-data$sp_ives

selected_vars <- vifstep_result@results$Variables
bio12 <- data$bio12
names(bio12) <-data$sp_ives

gdf$bio18 <- data$bio18
dimnames(gdf$shape)[[3]] <- data$sp_ives

model_reduced <- lm(size~bio12+bio17+bio18+pcnm1+pcnm5,data=data)
summary(model_reduced)

AIC(model_reduced)
step_model <- step(model_reduced, direction = "both")
summary(step_model) 
#plot(step_model)
AIC(step_model) #

data$pcnm1 <- pcnm1

data$pcnm5 <- pcnm5

library(caper)
library(MuMIn)
print(names(data))

# Definir a variável de resposta e as variáveis explicativas
y <- "size"
x <- setdiff(selected_vars, "npp")


DREDGE_PARALLEL_LIMITED <- function(y, x, fn = pgls.SEy, data, tree, corClass = corClass, method = "ML", max_combinations = 100) {
  n <- length(x)
  ii <- list()
  
  # Limitar o número de combinações de variáveis para evitar sobrecarga
  for (i in 1:min(n, max_combinations)) {
    ii <- c(ii, combn(x, i, simplify = FALSE))
  }
  
  # Filtrar combinações que contenham pcnm1 e pcnm5
  ii <- Filter(function(vars) all(c("pcnm1", "pcnm5") %in% vars), ii)
  
  # Ajustar modelo nulo (apenas intercepto)
  fits <- list()
  fits[[1]] <- fn(as.formula(paste(y, "~ 1")), data = data, tree = tree, method = method, corClass = corClass)
  
  # Função para ajustar o modelo e capturar falhas de convergência
  fit_model <- function(vars) {
    model_formula <- as.formula(paste(y, "~", paste(vars, collapse = " + ")))
    tryCatch({
      fn(model_formula, data = data, tree = tree, method = method, corClass = corClass)
    }, error = function(e) {
      message("Erro no modelo com as variáveis: ", paste(vars, collapse = ", "), " - ", e$message)
      return(NULL)  # Retorna NULL se houver erro de convergência
    })
  }
  
  # Paralelizar o ajuste dos modelos, com controle de erros
  fits_combinados <- foreach(i = 1:length(ii), .packages = c("caper", "MASS", "MuMIn")) %dopar% {
    fit_model(ii[[i]])
  }
  
  # Adiciona os modelos combinados à lista que já tem o modelo nulo
  fits <- c(fits, fits_combinados)
  
  # Retorna a lista com todos os modelos ajustados
  return(fits)
}

# Executando com restrição e controle de combinações
fits.bm_parallel <- DREDGE_PARALLEL_LIMITED(y, c(x, "pcnm1", "pcnm5"), data = data, tree = tree, corClass = corBrownian, max_combinations = 1000)
fits.lambda <- DREDGE_PARALLEL_LIMITED(y, c(x, "pcnm1", "pcnm5"), data = data, tree = tree, method = "REML", corClass = corPagel, max_combinations = 1000)

# Combina os resultados dos dois métodos
fits <- c(fits.bm_parallel, fits.lambda)

# Seleção do melhor modelo com base no AIC
valid_fits <- Filter(Negate(is.null), fits)

# Calcular o AIC apenas para modelos válidos
aic_values <- sapply(valid_fits, function(model) AIC(model))
best_model_index <- which.min(aic_values)
best_model <- fits[[best_model_index]]
AIC(best_model)

# Verificando os coeficientes do melhor modelo
summary(best_model)
aic_values
write.table(aic_values, file = "Lima_AIC.txt", row.names = FALSE)

best_model <-fits[aic_values==min(aic_values)][[3]]
summary(best_model)
names (best_model)
best_model$varBeta

# Extrair as informações do modelo
coefficients <- best_model$coefficients
sigma <- best_model$sigma
logLik <- best_model$logLik
fitted_values <- best_model$fitted
residuals <- best_model$residuals

# Criar um data frame com as informações
model_summary <- data.frame(
  Coefficients = coefficients,
  Sigma = rep(sigma, length(coefficients)),   # Repete o valor de sigma para cada coeficiente
  LogLik = rep(logLik, length(coefficients)),  # Repete o valor de logLik para cada coeficiente
  Fitted = fitted_values[1:length(coefficients)],  # Ajuste para pegar os primeiros valores ajustados
  Residuals = residuals[1:length(coefficients)]  # Ajuste para pegar os primeiros resíduos
)

# Exportar para um arquivo de texto
write.table(model_summary, file = "Lima_best_model.txt", row.names = TRUE, sep = "\t")

tree <- read.tree("Datasets/trees/AVGLOCSEX_Lima.tre")

fitBM <- pgls.SEy(size ~ bio12 + pcnm1 + pcnm5, data = data, tree = tree, corClass = corBrownian)
summary(fitBM)

fitPG <- pgls.SEy(size ~ bio12 + pcnm1 + pcnm5, data = data, tree = tree, corClass = corPagel)
summary(fitPG)

AIC(fitPG, fitBM)
plot(fitPG, 
     pch = as.numeric(data$pch), 
     col = as.factor(data$fac), 
     main = "Ajuste PGLS", 
     xlab = "Valores Preditos", 
     ylab = "Resíduos", 
     cex = 1.2) 
# Transformar 'fac' em fator se não for
data$fac <- as.factor(data$fac)

# Verificar os níveis de 'fac'
levels_fac <- levels(data$fac)

# Definir cores para 'fac'
cores_fac <- c("AF" ="green3","AM" = "darkgreen", "SV" = "#FFDB58")
names(cores_fac) <- levels_fac
bioma_shapes <- c("Atlantic Forest" = 16,"Amazon" = 15, "SVcay" = 17, "SVlib" = 25)  # Substitua com os shapes corretos e correspondentes aos biomas


########### Plot modelos de menor AIC #############

fit<-pgls.SEy(size ~ bio18 + pcnm1 + pcnm5, data = data, tree = tree, corClass = corPagel)
coeff <- coef(fit)
summary(fit) 
res <- resid(fit, type="n") 
qqnorm(res) 
qqline(res)
shapiro.test(res)
AIC(fit)
anova(fit)

fitted_vals <- fitted(fit)
res <- resid(fit, type = "n") 

plot(fitted_vals, res, pch = as.numeric(data$pch), col = as.factor(data$fac),
     main = "Residual vs. Fitted (PGLS Model)",
     xlab = "Fitted Values", ylab = "Residuals")
grid()

legend("bottomleft", legend = unique(data$fac), col = unique(as.factor(data$fac)), 
       pch = 15, title = "Environment", cex = 1)

modelo_bio <- fit

# Extraia os valores ajustados e resíduos
plot_data <- data.frame(
  Fitted = fitted(fitPG),
  Residuals = residuals(fitPG),
  Group = as.factor(data$fac),  # Variável categórica para cores
  Symbol = as.numeric(data$pch) # Variável para os símbolos
)

# Criando o gráfico no ggplot
ggplot(plot_data, aes(x = Fitted, y = Residuals, color = Group, shape = Group)) +
  geom_point(size = 3) +   geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +  # Linha horizontal em zero
  geom_smooth(method = "lm", se = FALSE, linetype = "solid", color = "blue") +  # Linha de suavização LOESS
  # Pontos com símbolos e cores
  labs(
    title = "Ajuste PGLS para Lima et al. (2018)", subtitle = "Modelo: BIO18 + PCNM",
    x = "Valores Preditos",
    y = "Resíduos"
  ) + scale_color_manual(values = cores_fac) +
  theme_minimal() +  # Estilo limpo
  theme(
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 10)
  )

ggplot(plot_data, aes(x = Fitted, y = Residuals, color = Group, shape = Group)) +
  geom_point(size = 3, alpha = 0.8) +  # Reduzir opacidade para visualização melhor
  geom_hline(yintercept = 0, linetype = "dashed", color = "black", size = 0.8) +  # Linha horizontal mais evidente
  geom_smooth(method = "lm", se = FALSE, linetype = "dotted", color = "red", size = 1,aes(group = 1)) +  # Linha de tendência linear
  labs(
    title = "Ajuste PGLS para Lima et al. (2018)",  subtitle = "BIO18 + PCNMs",
    x = "Valores Preditos",
    y = "Resíduos"
  ) + scale_color_manual(values = cores_fac) +  theme_minimal() +
  theme( plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
         plot.subtitle = element_text(size = 12, face = "italic", hjust = 0.5),
    axis.title = element_text(size = 14),
    legend.title = element_blank(),
    legend.position = "bottom"
  )

############################################################################################

# Plotar os pontos e a linha ajustada
p_bio <- ggplot(data, aes(x = bio18 + pcnm1 + pcnm5, y = size, color = fac)) +
  geom_point(size = 7) +
  geom_smooth(method = "lm", formula = y ~ x, se = FALSE, color = "blue", linetype = "dashed") +
  geom_abline(intercept = coeff[1], slope = coeff[2], color = "red", linetype = "solid") +
  labs(x = "Precipitation in the warmest quarter (BIO18) + PCNM", y = "Mandibular size (size)",
       title = "Relationship between size and Bioclimatic Variables with PGLS",
       color = "Fac") + scale_color_manual(values = cores_fac) +  theme_minimal() +
  theme( plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
         plot.subtitle = element_text(size = 12, face = "italic", hjust = 0.5),
         axis.title = element_text(size = 14),
         legend.title = element_blank(),
         legend.position = "bottom"
  )

print(p_bio)

dev.off()


################### Para a Forma:


############################# Shape and Size ############################################

fit<-pgls.SEy(shape.2d~size+bio1,data=data,tree=tree, corClass=corPagel)
summary(fit) ###
fit<-pgls.SEy(shape.2d~size+bio2,data=data,tree=tree, corClass=corPagel)
summary(fit)
fit<-pgls.SEy(shape.2d~size+bio3,data=data,tree=tree, corClass=corPagel)
summary(fit) 
fit<-pgls.SEy(shape.2d~size+bio4,data=data,tree=tree, corClass=corPagel)
summary(fit) 
fit<-pgls.SEy(shape.2d~size+bio5,data=data,tree=tree)
summary(fit) ##################
fit<-pgls.SEy(shape.2d~size+bio6,data=data,tree=tree, corClass=corPagel)
summary(fit)
fit<-pgls.SEy(shape.2d~size+bio7,data=data,tree=tree, corClass=corPagel) 
summary(fit) 
fit<-pgls.SEy(shape.2d~size+bio8,data=data,tree=tree, corClass=corPagel)
summary(fit) ########################
fit<-pgls.SEy(shape.2d~size+bio9,data=data,tree=tree, corClass=corPagel)
summary(fit) 
fit<-pgls.SEy(shape.2d~size+bio10,data=data,tree=tree, corClass=corPagel)
summary(fit) #########################
fit<-pgls.SEy(shape.2d~size+bio11,data=data,tree=tree, corClass=corPagel)
summary(fit) 
fit<-pgls.SEy(shape.2d~size+bio12,data=data,tree=tree, corClass=corPagel)
summary(fit)##############
fit<-pgls.SEy(shape.2d~size+bio13,data=data,tree=tree, corClass=corPagel)
summary(fit) 
fit<-pgls.SEy(shape.2d~size+bio14,data=data,tree=tree, corClass=corPagel)
summary(fit) ##############
fit<-pgls.SEy(shape.2d~size+bio15,data=data,tree=tree, corClass=corPagel)
summary(fit) #############
fit<-pgls.SEy(shape.2d~size+bio16,data=data,tree=tree, corClass=corPagel)
summary(fit)
fit<-pgls.SEy(shape.2d~size+bio17,data=data,tree=tree, corClass=corPagel)
summary(fit) #############
fit<-pgls.SEy(shape.2d~size+bio18,data=data,tree=tree, corClass=corPagel)
summary(fit)
fit<-pgls.SEy(shape.2d~size+bio19,data=data,tree=tree, corClass=corPagel)
summary(fit) 
fit<-pgls.SEy(shape.2d~size+npp,data=data,tree=tree, corClass=corPagel)
summary(fit)
fit<-pgls.SEy(shape.2d~size+humid,data=data,tree=tree, corClass=corPagel)
summary(fit)
fit<-pgls.SEy(shape.2d~size+soilmoist,data=data,tree=tree, corClass=corPagel)
summary(fit) ############

##repeat model with filters

fit<-pgls.SEy(shape.2d~size+bio1+pcnm31+pcnm5+pcnm2+pcnm7,data=data,tree=tree, corClass=corPagel)
summary(fit)

fit<-pgls.SEy(shape.2d~size+bio2+pcnm31+pcnm5+pcnm2+pcnm7,data=data,tree=tree, corClass=corPagel)
summary(fit) ###################

fit<-pgls.SEy(shape.2d~size+bio3+pcnm31+pcnm5+pcnm2+pcnm7,data=data,tree=tree, corClass=corPagel)
summary(fit) 

fit<-pgls.SEy(shape.2d~size+bio4+pcnm31+pcnm5+pcnm2+pcnm7,data=data,tree=tree, corClass=corPagel)
summary(fit) 

fit<-pgls.SEy(shape.2d~size+bio5+pcnm31+pcnm5+pcnm2+pcnm7,data=data,tree=tree, corClass=corPagel)
summary(fit) ############################

fit<-pgls.SEy(shape.2d~size+bio6+pcnm31+pcnm5+pcnm2+pcnm7,data=data,tree=tree, corClass=corPagel)
summary(fit)

fit<-pgls.SEy(shape.2d~size+bio7+pcnm31+pcnm5+pcnm2+pcnm7,data=data,tree=tree, corClass=corPagel)
summary(fit) ########################

fit<-pgls.SEy(shape.2d~size+bio8+pcnm31+pcnm5+pcnm2+pcnm7,data=data,tree=tree, corClass=corPagel)
summary(fit) ########################

fit<-pgls.SEy(shape.2d~size+bio9+pcnm31+pcnm5+pcnm2+pcnm7,data=data,tree=tree, corClass=corPagel)
summary(fit)

fit<-pgls.SEy(shape.2d~size+bio10+pcnm31+pcnm5+pcnm2+pcnm7,data=data,tree=tree)
summary(fit) ##

fit<-pgls.SEy(shape.2d~size+bio11+pcnm31+pcnm5+pcnm2+pcnm7,data=data,tree=tree, corClass=corPagel)
summary(fit)

fit<-pgls.SEy(shape.2d~size+bio12+pcnm31+pcnm5+pcnm2+pcnm7,data=data,tree=tree, corClass=corPagel)
summary(fit) #######

fit<-pgls.SEy(shape.2d~size+bio13+pcnm31+pcnm5+pcnm2+pcnm7,data=data,tree=tree, corClass=corPagel)
summary(fit) 

fit<-pgls.SEy(shape.2d~size+bio14+pcnm31+pcnm5+pcnm2+pcnm7,data=data,tree=tree, corClass=corPagel)
summary(fit) ###############

fit<-pgls.SEy(shape.2d~size+bio15+pcnm31+pcnm5+pcnm2+pcnm7,data=data,tree=tree, corClass=corPagel)
summary(fit) #################

fit<-pgls.SEy(shape.2d~size+bio16+pcnm31+pcnm5+pcnm2+pcnm7,data=data,tree=tree, corClass=corPagel)
summary(fit) 

fit<-pgls.SEy(shape.2d~size+bio17+pcnm31+pcnm5+pcnm2+pcnm7,data=data,tree=tree, corClass=corPagel)
summary(fit) ###################

fit<-pgls.SEy(shape.2d~size+bio18+pcnm31+pcnm5+pcnm2+pcnm7,data=data,tree=tree, corClass=corPagel)
summary(fit)

fit<-pgls.SEy(shape.2d~size+bio19+pcnm31+pcnm5+pcnm2+pcnm7,data=data,tree=tree, corClass=corPagel)
summary(fit)

fit<-pgls.SEy(shape.2d~size+npp+pcnm31+pcnm5+pcnm2+pcnm7,data=data,tree=tree, corClass=corPagel)
summary(fit)

fit<-pgls.SEy(shape.2d~size+humid+pcnm31+pcnm5+pcnm2+pcnm7,data=data,tree=tree, corClass=corPagel)
summary(fit)

fit<-pgls.SEy(shape.2d~size+soilmoist+pcnm31+pcnm5+pcnm2+pcnm7,data=data,tree=tree, corClass=corPagel)
summary(fit)

####### BIO2,BIO5,BIO7,BIO8,BIO10,BIO12,BIO14,BIO15,BIO17


dimnames(gpa$coords)[[3]] <- paste(as.factor(data$sp_ives))
names(gpa$Csize) <- paste(as.factor(data$sp_ives))

fit_bio1 <- procD.pgls(gpa$coords ~ log(gpa$Csize) + bio1 + pcnm2 + pcnm31 + pcnm1 + pcnm7, phy = tree, data = data)
summary(fit_bio1)

fit_bio2 <- procD.pgls(gpa$coords ~ log(gpa$Csize) + bio2 + pcnm2 + pcnm31 + pcnm1 + pcnm7, phy = tree, data = data)
summary(fit_bio2)

fit_bio3 <- procD.pgls(gpa$coords ~ log(gpa$Csize) + bio3 + pcnm2 + pcnm31 + pcnm1 + pcnm7, phy = tree, data = data)
summary(fit_bio3)

fit_bio4 <- procD.pgls(gpa$coords ~ log(gpa$Csize) + bio4 + pcnm2 + pcnm31 + pcnm1 + pcnm7, phy = tree, data = data)
summary(fit_bio4)

fit_bio5 <- procD.pgls(gpa$coords ~ log(gpa$Csize) + bio5 + pcnm2 + pcnm31 + pcnm1 + pcnm7, phy = tree, data = data)
summary(fit_bio5)

fit_bio6 <- procD.pgls(gpa$coords ~ log(gpa$Csize) + bio6 + pcnm2 + pcnm31 + pcnm1 + pcnm7, phy = tree, data = data)
summary(fit_bio6)

fit_bio7 <- procD.pgls(gpa$coords ~ log(gpa$Csize) + bio7 + pcnm2 + pcnm31 + pcnm1 + pcnm7, phy = tree, data = data)
summary(fit_bio7)

fit_bio8 <- procD.pgls(gpa$coords ~ log(gpa$Csize) + bio8 + pcnm2 + pcnm31 + pcnm1 + pcnm7, phy = tree, data = data)
summary(fit_bio8)

fit_bio9 <- procD.pgls(gpa$coords ~ log(gpa$Csize) + bio9 + pcnm2 + pcnm31 + pcnm1 + pcnm7, phy = tree, data = data)
summary(fit_bio9)

fit_bio10 <- procD.pgls(gpa$coords ~ log(gpa$Csize) + bio10 + pcnm2 + pcnm31 + pcnm1 + pcnm7, phy = tree, data = data)
summary(fit_bio10)

fit_bio11 <- procD.pgls(gpa$coords ~ log(gpa$Csize) + bio11 + pcnm2 + pcnm31 + pcnm1 + pcnm7, phy = tree, data = data)
summary(fit_bio11)

fit_bio12 <- procD.pgls(gpa$coords ~ log(gpa$Csize) + bio12 + pcnm2 + pcnm31 + pcnm1 + pcnm7, phy = tree, data = data)
summary(fit_bio12)

fit_bio13 <- procD.pgls(gpa$coords ~ log(gpa$Csize) + bio13 + pcnm2 + pcnm31 + pcnm1 + pcnm7, phy = tree, data = data)
summary(fit_bio13)

fit_bio14 <- procD.pgls(gpa$coords ~ log(gpa$Csize) + bio14 + pcnm2 + pcnm31 + pcnm1 + pcnm7, phy = tree, data = data)
summary(fit_bio14)

fit_bio15 <- procD.pgls(gpa$coords ~ log(gpa$Csize) + bio15 + pcnm2 + pcnm31 + pcnm1 + pcnm7, phy = tree, data = data)
summary(fit_bio15)

fit_bio16 <- procD.pgls(gpa$coords ~ log(gpa$Csize) + bio16 + pcnm2 + pcnm31 + pcnm1 + pcnm7, phy = tree, data = data)
summary(fit_bio16)

fit_bio17 <- procD.pgls(gpa$coords ~ log(gpa$Csize) + bio17 + pcnm2 + pcnm31 + pcnm1 + pcnm7, phy = tree, data = data)
summary(fit_bio17)

fit_bio18 <- procD.pgls(gpa$coords ~ log(gpa$Csize) + bio18 + pcnm2 + pcnm31 + pcnm1 + pcnm7, phy = tree, data = data)
summary(fit_bio18)

fit_bio19 <- procD.pgls(gpa$coords ~ log(gpa$Csize) + bio19 + pcnm2 + pcnm31 + pcnm1 + pcnm7, phy = tree, data = data)
summary(fit_bio19)

fit_npp <- procD.pgls(gpa$coords ~ log(gpa$Csize) + npp + pcnm2 + pcnm31 + pcnm1 + pcnm7, phy = tree, data = data)
summary(fit_npp)

fit_humid <- procD.pgls(gpa$coords ~ log(gpa$Csize) + humid + pcnm2 + pcnm31 + pcnm1 + pcnm7, phy = tree, data = data)
summary(fit_humid)

fit_soilmoist <- procD.pgls(gpa$coords ~ log(gpa$Csize) + soilmoist + pcnm2 + pcnm31 + pcnm1 + pcnm7, phy = tree, data = data)
summary(fit_soilmoist)

library(tibble)
library(dplyr)
library(readr)

# Lista com os nomes dos objetos
nomes_modelos <- c(
  paste0("fit_bio", 1:19), "fit_npp", "fit_humid", "fit_soilmoist"
)

# Função para extrair informações de cada modelo
extrai_info <- function(nome_objeto) {
  fit <- get(nome_objeto)
  resumo <- summary(fit)
  
  # Nome da variável ambiental (assumindo padrão do nome do modelo)
  variavel <- gsub("fit_", "", nome_objeto)
  
  # Tenta pegar valor-p da variável ambiental (segunda posição no ANOVA)
  p_valor <- tryCatch({
    resumo$aov.table$`Pr(>F)`[2]
  }, error = function(e) NA)
  
  tibble(
    modelo = nome_objeto,
    variavel = variavel,
    Rsq = resumo$Rsq,
    AIC = AIC(fit),
    p_valor = p_valor
  )
}

# Aplica a função a todos os modelos
tabela_resultados <- purrr::map_dfr(nomes_modelos, extrai_info)

# Salva em CSV
write_csv(tabela_resultados, "resultados_modelos_procD.csv")


# Definir a variável de resposta e as variáveis explicativas

data$shape.2d <- shape.2d
y <- "shape.2d"

clim <- as.data.frame(bioclim[,c("bio7", "bio19")])
vif(clim)

vifstep_result <- vifstep(bioclim,th = 1.5, keep = c("bio7", "bio19"))
print(vifstep_result)

vifstep_result@results$Variables 
names(data)
bio_data <- as.data.frame(scale(data[,11:32]))
size <- as.numeric(data$size)

rownames(bio_data) <-data$sp_ives
names(size) <-data$sp_ives
names(bio12) <-data$sp_ives

selected_vars <- vifstep_result@results$Variables

data$pcnm31 <- pcnm31
data$pcnm7 <- pcnm7
data$pcnm2 <- pcnm2

x <- c("bio7", "bio19") #Lima (2018) - BM



DREDGE_PARALLEL_LIMITED <- function(y, x, fn = pgls.SEy, data, tree, corClass = corClass, method = "ML") {
  n <- length(x)
  ii <- list()
  
  # Gerar combinações de variáveis
  for (i in 1:n) ii <- c(ii, combn(x, i, simplify = FALSE))
  
  # Filtrar combinações que contenham pcnm1 e pcnm5
  ii <- Filter(function(vars) all(c("pcnm31", "pcnm7", "pcnm2") %in% vars), ii)
  
  # Ajustar modelo vazio
  fits <- list()
  fits[[1]] <- fn(as.formula(paste(y, "~ 1")), data = data, tree = tree, method = method)
  
  # Paralelizar o ajuste dos modelos
  fits <- foreach(i = 1:length(ii), .packages = c("caper", "MASS", "MuMIn")) %dopar% {
    model <- as.formula(paste(y, "~", paste(ii[[i]], collapse = "+")))
    fn(model, data = data, tree = tree, method = method)
  }
  
  # Retornar lista de modelos ajustados
  return(fits)
}

# Executando com restrição
fits.bm_parallel <- DREDGE_PARALLEL_LIMITED(y, c(x, "pcnm31", "pcnm7", "pcnm2"), data = data, tree = tree, method = "ML")
fits.lambda <- DREDGE_PARALLEL_LIMITED(y, c(x, "pcnm31", "pcnm7", "pcnm2"), data = data, corClass = corPagel, tree = tree, method = "ML")

fits <- c(fits.bm_parallel, fits.lambda)

# Seleção do melhor modelo com base no AIC
aic_values <- sapply(fits, function(model) AIC(model))
best_model_index <- which.min(aic_values)
best_model <- fits[[best_model_index]]

# Verificando os coeficientes do melhor modelo
summary(best_model)

## AIC values for each model
aic<-sapply(fits,AIC)
aic <- setdiff(aic,aic[1])
#write.table(aic, file = "Lima_AIC.txt", row.names = FALSE)

best_model <-fits[aic==min(aic)][[2]]
summary(best_model)
names (best_model)

# Extrair as informações do modelo
coefficients <- best_model$coefficients
sigma <- best_model$sigma
logLik <- best_model$logLik
fitted_values <- best_model$fitted
residuals <- best_model$residuals

# Criar um data frame com as informações
model_summary <- data.frame(
  Coefficients = coefficients,
  Sigma = rep(sigma, length(coefficients)),   # Repete o valor de sigma para cada coeficiente
  LogLik = rep(logLik, length(coefficients)),  # Repete o valor de logLik para cada coeficiente
  Fitted = fitted_values[1:length(coefficients)],  # Ajuste para pegar os primeiros valores ajustados
  Residuals = residuals[1:length(coefficients)]  # Ajuste para pegar os primeiros resíduos
)

# Exportar para um arquivo de texto
#write.table(model_summary, file = "Lima_best_model.txt", row.names = TRUE, sep = "\t")

tree <- read.tree("Datasets/trees/AVGLOCSEX_Lima.tre")

fitBM <- pgls.SEy(size ~ bio7 + bio19 + pcnm1 + pcnm5, data = data, tree = tree, corClass = corBrownian)
summary(fitBM)

fitPG <- pgls.SEy(size ~ bio7 + bio19 + pcnm1 + pcnm5, data = data, tree = tree, corClass = corPagel)
summary(fitPG)

# Aplicar AIC nos modelos listados
aic_bm <- print(AIC(fitBM))
aic_pg <- print(AIC(fitPG))

plot(fitPG, 
     pch = as.numeric(data$pch), 
     col = as.factor(data$fac), 
     main = "Ajuste PGLS", 
     xlab = "Valores Preditos", 
     ylab = "Resíduos", 
     cex = 1.2) 
# Transformar 'fac' em fator se não for
data$fac <- as.factor(data$fac)

# Verificar os níveis de 'fac'
levels_fac <- levels(data$fac)

# Definir cores para 'fac'
cores_fac <- c("AF" ="green3","AM" = "darkgreen", "SV" = "#FFDB58")
names(cores_fac) <- levels_fac
bioma_shapes <- c("Atlantic Forest" = 16,"Amazon" = 15, "SVcay" = 17, "SVlib" = 25)  # Substitua com os shapes corretos e correspondentes aos biomas


########### Plot modelos de menor AIC #############

fit<-pgls.SEy(size ~ bio7 + bio19 + pcnm1 + pcnm5, data = data, tree = tree, corClass = corPagel)
coef(fit)
summary(fit) 
res <- resid(fit, type="n") 
qqnorm(res) 
qqline(res)
shapiro.test(res)
AIC(fit)

fitted_vals <- fitted(fit)
res <- resid(fit, type = "n") 

plot(fitted_vals, res, pch = as.numeric(data$pch), col = as.factor(data$fac),
     main = "Residual vs. Fitted (PGLS Model)",
     xlab = "Fitted Values", ylab = "Residuals")
grid()

legend("bottomleft", legend = unique(data$fac), col = unique(as.factor(data$fac)), 
       pch = 15, title = "Environment", cex = 1)

modelo_bio <- fit

# Extraia os valores ajustados e resíduos
plot_data <- data.frame(
  Fitted = fitted(fitPG),
  Residuals = residuals(fitPG),
  Group = as.factor(data$fac),  # Variável categórica para cores
  Symbol = as.numeric(data$pch) # Variável para os símbolos
)

# Criando o gráfico no ggplot
ggplot(plot_data, aes(x = Fitted, y = Residuals, color = Group, shape = Group)) +
  geom_point(size = 3) +   geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +  # Linha horizontal em zero
  geom_smooth(method = "loess", se = FALSE, linetype = "solid", color = "blue") +  # Linha de suavização LOESS
  # Pontos com símbolos e cores
  labs(
    title = "Ajuste PGLS para Lima et al. (2018)", subtitle = "Modelo: BIO7 + BIO19 + PCNM",
    x = "Valores Preditos",
    y = "Resíduos"
  ) + scale_color_manual(values = cores_fac) +
  theme_minimal() +  # Estilo limpo
  theme(
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 10)
  )

ggplot(plot_data, aes(x = Fitted, y = Residuals, color = Group, shape = Group)) +
  geom_point(size = 3, alpha = 0.8) +  # Reduzir opacidade para visualização melhor
  geom_hline(yintercept = 0, linetype = "dashed", color = "black", size = 0.8) +  # Linha horizontal mais evidente
  geom_smooth(method = "lm", se = FALSE, linetype = "dotted", color = "red", size = 1,aes(group = 1)) +  # Linha de tendência linear
  labs(
    title = "Ajuste PGLS para Lima et al. (2018)",  subtitle = "BIO3 + BIO12 + BIO19 + PCNMs",
    x = "Valores Preditos",
    y = "Resíduos"
  ) + scale_color_manual(values = cores_fac) +  theme_minimal() +
  theme( plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
         plot.subtitle = element_text(size = 12, face = "italic", hjust = 0.5),
         axis.title = element_text(size = 14),
         legend.title = element_blank(),
         legend.position = "bottom"
  )

############################################################################################
fit
coef_bm <- coef(fit)

names(data)
data[,c(11:37)] <- scale(data[,c(11:30)],shape.2d)

y <- data_plot$shape.2d

# Plotar os pontos e a linha ajustada
p_bio12 <- ggplot(data, aes(x = bio7 + bio19, y = y[1], color = Fac)) +
  geom_point(size = 3) +
  geom_smooth(method = "lm", formula = y ~ x, se = FALSE, color = "blue", linetype = "dashed") +
  geom_abline(intercept = coef_bm[1], slope = coef_bm[2], color = "red", linetype = "solid") +
  labs(x = "BIO17", y = "Shape",
       title = "Relationship between size and Bioclimatic Variables with PGLS",
       color = "Fac") +
  scale_color_manual(values = cores_fac) +
  theme_minimal()

print(p_bio12)

dev.off()


