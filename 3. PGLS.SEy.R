#================================================== PARVAR + RDA and PGLS.Sey

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
setwd("~")

#=====================#=====================#=====================#gm analyses
sapajus<-readland.tps("tps/avglocsex.tps", specID="ID")
gpa<-gpagen(sapajus)
dim(sapajus)

data<-read.csv("Plans/avglocsex.csv", sep=";")
dim(data)
names(data)

#link<-define.links(gpa$consensus, ptsize = 2, links = NULL)
#write.table(link,file="link.txt")

link <- read.table("tps/link.txt")

gdf<-geomorph.data.frame(shape=gpa$coords, size = log(gpa$Csize), biome = as.factor(data$biome), facgp = as.factor(data$fac), spp = as.factor(data$sp))

drawinglandmark<-readland.tps("outline2/outline.tps")
outline<-read.table("outline2/outline.txt", header=FALSE)
summary(drawinglandmark)
mshape<-mshape(gpa$coords)
summary(mshape)
summary(outline)
Sapajusoutline<-warpRefOutline(file = "outline2/outline.txt", drawinglandmark[,,1],mshape)#run dev.off() in case of error message ##set outline configuration

dev.off()

#=====================#=====================#=====================#=====================#=====================#=====================## PARVAR #=====================#=====================#=====================#=====================#=====================#=====================#=====================#=====================#=====================#=====================#=====================##
names(data)
coords<-data[,4:5]

geographicaldistance <- as.matrix(dist(coords))
pcnmcoord<-prcomp(geographicaldistance)
summary(pcnmcoord)
pcnm<-pcnmcoord$x
rownames(pcnm) <-data$sp_ives

#or

pcnmcoord<- pcnm(dist(coords))
pcnm<-as.matrix(pcnmcoord$vectors)
rownames(pcnm) <-data$sp_ives

#write.table(pcnm, "pcnm.xls", row.names= TRUE, col.names = TRUE, sep = " ")

### Environment
names(data)
bioclim <- data[,c(11:32)]
bioclim <- scale(bioclim,center=T,scale=T)
head(bioclim)

data[,c(11:32)] <- scale(data[,c(11:32)])
rownames(bioclim) <-data$sp_ives

envpca<-prcomp(bioclim)
summary(envpca)
env<-as.matrix(envpca$x[,1:6]) # 95%
rownames(env) <-data$sp_ives

### Shape (Y)

allo<-procD.lm(gpa$coords~log(gpa$Csize),iter=999, RRPP = TRUE)
summary(allo)
plot(allo,type="regression",predictor=log(gpa$Csize),
     reg.type = "RegScore", pch=19, xlab="log(Size)")

shape.resid<-
  arrayspecs(allo$residuals,p=dim(gpa$coords)[1],k=dim(gpa$coords)[2]) 
adj.shape<-shape.resid+array(gpa$consensus, dim(shape.resid)) 
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
tree<-read.tree("trees/AVGLOCSEX_Lima.tre")
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
all.equal(rownames(pcnm), rownames(data))       # Deve retornar TRUE
all.equal(rownames(bioclim), rownames(data))   # Deve retornar TRUE
all.equal(rownames(pcsgm), rownames(data))     # Deve retornar TRUE
all.equal(rownames(env), rownames(data))       # Deve retornar TRUE

parvar<-varpart(pcsgm,~phy,~size,~env,~pcnm)

plot(parvar,cex=1.2,Xnames=c("Phylogeny","Spatial Structure","Size","Environment"),bg=c("yellow", "deeppink","navy","green3"))

y<-pcsgm
sf<-pcnm
env <-env

### for all parvar's variables #=====================#

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


### RDA condition

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

#===================== PGLS.SEy (Ives, 2007) R #=====================#

# This analyses set was the first used in this paper, but the scripts 4a and 4b proved to be more efficient

## a For can be done here

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
summary(fit)#=====================#

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

#===================== Pagel

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

# bio12 > bio16 > bio17 > bio10 > bio18 > bio4

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
summary(fit)#=====================###

fit<-pgls.SEy(size~bio10+pcnm1+pcnm5,data=data,tree=tree, corClass=corPagel)
summary(fit)

fit<-pgls.SEy(size~bio11+pcnm1+pcnm5,data=data,tree=tree, corClass=corPagel)
summary(fit)

fit<-pgls.SEy(size~bio12+pcnm1+pcnm5,data=data,tree=tree, corClass=corPagel)
summary(fit) #=====================#=====================#

fit<-pgls.SEy(size~bio13+pcnm1+pcnm5,data=data,tree=tree, corClass=corPagel)
summary(fit) 

fit<-pgls.SEy(size~bio14+pcnm1+pcnm5,data=data,tree=tree, corClass=corPagel)
summary(fit) 

fit<-pgls.SEy(size~bio15+pcnm1+pcnm5,data=data,tree=tree, corClass=corPagel)
summary(fit)#=====================# error

fit<-pgls.SEy(size~bio16+pcnm1+pcnm5,data=data,tree=tree, corClass=corPagel)
summary(fit) 

fit<-pgls.SEy(size~bio17+pcnm1+pcnm5,data=data,tree=tree, corClass=corPagel)
summary(fit)#

fit<-pgls.SEy(size~bio18+pcnm1+pcnm5,data=data,tree=tree, corClass=corPagel)
summary(fit)#

fit<-pgls.SEy(size~bio19+pcnm1+pcnm5,data=data,tree=tree, corClass=corPagel)
summary(fit)

fit<-pgls.SEy(size~npp+pcnm1+pcnm5,data=data,tree=tree, corClass=corPagel)
summary(fit) #=====================#=====================#=====================## FALSE CONVERGENCE

fit<-pgls.SEy(size~humid+pcnm1+pcnm5,data=data,tree=tree, corClass=corPagel)
summary(fit)

fit<-pgls.SEy(size~soilmoist+pcnm1+pcnm5,data=data,tree=tree, corClass=corPagel)
summary(fit)

#  bio12 > bio17 bio18 

#============= Here can be done a DREDGE, but this may not work 
y <- "size"
x <- setdiff(selected_vars, "npp")


DREDGE_PARALLEL_LIMITED <- function(y, x, fn = pgls.SEy, data, tree, corClass = corClass, method = "ML", max_combinations = 100) {
  n <- length(x)
  ii <- list()
  
  for (i in 1:min(n, max_combinations)) {
    ii <- c(ii, combn(x, i, simplify = FALSE))
  }
    ii <- Filter(function(vars) all(c("pcnm1", "pcnm5") %in% vars), ii)
    fits <- list()
  fits[[1]] <- fn(as.formula(paste(y, "~ 1")), data = data, tree = tree, method = method, corClass = corClass)
    fit_model <- function(vars) {
    model_formula <- as.formula(paste(y, "~", paste(vars, collapse = " + ")))
    tryCatch({
      fn(model_formula, data = data, tree = tree, method = method, corClass = corClass)
    }, error = function(e) {
      message("Error ", paste(vars, collapse = ", "), " - ", e$message)
      return(NULL)  
    })
  }
    fits_combinados <- foreach(i = 1:length(ii), .packages = c("caper", "MASS", "MuMIn")) %dopar% {
    fit_model(ii[[i]])
  }
    fits <- c(fits, fits_combinados)
    return(fits)
}
fits.bm_parallel <- DREDGE_PARALLEL_LIMITED(y, c(x, "pcnm1", "pcnm5"), data = data, tree = tree, corClass = corBrownian, max_combinations = 1000)
fits.lambda <- DREDGE_PARALLEL_LIMITED(y, c(x, "pcnm1", "pcnm5"), data = data, tree = tree, method = "REML", corClass = corPagel, max_combinations = 1000)

fits <- c(fits.bm_parallel, fits.lambda)
valid_fits <- Filter(Negate(is.null), fits)

aic_values <- sapply(valid_fits, function(model) AIC(model))
best_model_index <- which.min(aic_values)
best_model <- fits[[best_model_index]]
AIC(best_model)
summary(best_model)
aic_values
write.table(aic_values, file = "Lima_AIC.txt", row.names = FALSE)

best_model <-fits[aic_values==min(aic_values)][[3]]
summary(best_model)
names (best_model)
best_model$varBeta

coefficients <- best_model$coefficients
sigma <- best_model$sigma
logLik <- best_model$logLik
fitted_values <- best_model$fitted
residuals <- best_model$residuals

model_summary <- data.frame(
  Coefficients = coefficients,
  Sigma = rep(sigma, length(coefficients)),   
  LogLik = rep(logLik, length(coefficients)), 
  Fitted = fitted_values[1:length(coefficients)],
  Residuals = residuals[1:length(coefficients)]  
)

write.table(model_summary, file = "Lima_best_model.txt", row.names = TRUE, sep = "\t")
tree <- read.tree("trees/AVGLOCSEX_Lima.tre")
                     
fitBM <- pgls.SEy(size ~ bio12 + pcnm1 + pcnm5, data = data, tree = tree, corClass = corBrownian)
summary(fitBM)

fitPG <- pgls.SEy(size ~ bio12 + pcnm1 + pcnm5, data = data, tree = tree, corClass = corPagel)
summary(fitPG)

AIC(fitPG, fitBM)
plot(fitPG, 
     pch = as.numeric(data$pch), 
     col = as.factor(data$fac), 
     main = NULL, 
     xlab = "Fitted", 
     ylab = "Ressiduals", 
     cex = 1.2) 

data$fac <- as.factor(data$fac)
levels_fac <- levels(data$fac)

cores_fac <- c("AF" ="green3","AM" = "darkgreen", "SV" = "#FFDB58")
names(cores_fac) <- levels_fac
bioma_shapes <- c("Atlantic Forest" = 16,"Amazon" = 15, "SVcay" = 17, "SVlib" = 25) 

#=====================#=====================### Plot best AIC #=====================#

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
plot_data <- data.frame(
  Fitted = fitted(fitPG),
  Residuals = residuals(fitPG),
  Group = as.factor(data$fac),  
  Symbol = as.numeric(data$pch) 
)

ggplot(plot_data, aes(x = Fitted, y = Residuals, color = Group, shape = Group)) +
  geom_point(size = 3) +   geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") + 
  geom_smooth(method = "lm", se = FALSE, linetype = "solid", color = "blue") +  

  labs(
    title = "",
    x = "Fitted",
    y = "Residuals"
  ) + scale_color_manual(values = cores_fac) +
  theme_minimal() + 
  theme(
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 10)
  )

ggplot(plot_data, aes(x = Fitted, y = Residuals, color = Group, shape = Group)) +
  geom_point(size = 3, alpha = 0.8) + 
  geom_hline(yintercept = 0, linetype = "dashed", color = "black", size = 0.8) +  
  geom_smooth(method = "lm", se = FALSE, linetype = "dotted", color = "red", size = 1,aes(group = 1)) +  
  labs(
    title = "",  subtitle = "BIO18 + PCNMs",
    x = "Fitted",
    y = "Residuals"
  ) + scale_color_manual(values = cores_fac) +  theme_minimal() +
  theme( plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
         plot.subtitle = element_text(size = 12, face = "italic", hjust = 0.5),
    axis.title = element_text(size = 14),
    legend.title = element_blank(),
    legend.position = "bottom"
  )

#=====================#=====================#=====================#=====================#=====================#=====================#=====================#=====================#=====================#=====================#=====================#=====================#=====================#=====================#=====================#=====================#=====================#=====================#=====================#=====================#=====================#=====================#=====================

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


