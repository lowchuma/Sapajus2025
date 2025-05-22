#################### LDA & CVA #########################
require(devtools)
setwd("C:/Users/Lourenço/OneDrive/UFSM/LabMastozoo/SAPAJUS/Data")

#devtools::install_github("fawda123/ggord", force = TRUE)

require(klaR)
require(psych)
require(MASS)
require(ggord)
require(geomorph)
require(shapes)
require(ggplot2)
require(vegan)
require(devtools)
require(RColorBrewer)
require(tidyverse)
require(ape)
require(phytools)
require(picante)
require(Morpho)
require(Rvcg)
require(stats)
require(reshape2)
require(dotwhisker)
require(dplyr)
library(geiger)
library(phangorn)
library(usdm)


tps <- readland.tps("tps/avglocsex.tps", specID =  "ID")
classifiers <- read.csv("avglocsex.csv", sep=";")

gpa.obj <- gpagen(tps)
shape <- gpa.obj$coords
ref <- mshape(shape)

spec <- as.factor(classifiers$sp)
sex <- as.factor(classifiers$sex)
env <- as.factor(classifiers$fac)
bio <- as.factor(classifiers$biome)
l.env <- levels(env)
cores <- c("green3", "darkgreen", "goldenrod3")
names(cores) <- l.env
sym <- c("Atlantic Forest" = 16, "Amazon" = 15, "SV" = 17)
names(sym) <- l.env

### Carregar o outline
drawinglandmark<-readland.tps("outline2/outline.tps")
outline<-read.table("outline2/outline.txt", header=FALSE)
summary(drawinglandmark)
Sapajusoutline<-warpRefOutline(file = "outline2/outline.txt", drawinglandmark[,,1],ref)#run dev.off() in case of error message ##set outline configuration

### PCA ###

pca <- gm.prcomp(shape)

GP <- gridPar(n.col.cell = 100, pt.bg = "gray", pt.size = 1, tar.pt.bg = "cyan", 
              tar.pt.size = 1.5,tar.out.col = "gray10", tar.out.cex = 1, 
              grid.col = "white",grid.lwd = 0.5,txt.pos = 1, txt.col = "steelblue")

plotRefToTarget(pca$shapes$shapes.comp1$min,pca$shapes$shapes.comp1$max,outline = Sapajusoutline$outline,gridPars=GP,method="TPS") 
plotRefToTarget(pca$shapes$shapes.comp1$min,pca$shapes$shapes.comp1$max,outline = Sapajusoutline$outline,gridPars=GP,method="vector") 
plotRefToTarget(pca$shapes$shapes.comp1$min,pca$shapes$shapes.comp1$max,outline = Sapajusoutline$outline,gridPars=GP,method="points")

bioclim <- scale(classifiers[,11:30])
vs <-vifstep(bioclim, keep = "bio12")
select <- vs@results$Variables

clim <- as.matrix(bioclim[,select])

env.pca <- prcomp(bioclim)

??xlims

### LDA ###

linear <- lda(clim,env)
linear
ggord(linear,env, axes = c("1", "2"), grp_in = NULL, cols = cores, facet = FALSE, nfac = NULL,
  addpts = NULL, obslab = FALSE, ptslab = FALSE,  ellipse = TRUE,  ellipse_pro = 0.95,  poly = TRUE,
  polylntyp = "solid", hull = FALSE, arrow = 0.3, labcol = "black",veccol = "gray30", vectyp = "solid",
  veclsz = 1, ext = 1.3, repel = FALSE, vec_ext = 2, vec_lab = NULL, size = 6, sizelab = NULL, 
  addsize = size/2, addcol = "blue", addpch = 19, txt = 4,  alpha = 1,  alpha_el = 0.4,  xlims = c(-5,6),
  ylims = c(-5,4), var_sub = NULL, coord_fix = TRUE, parse = TRUE, grp_title = "Ecoregions", force = 1,
  max.overlaps = 10, exp = c(0, 0))

### CVA ###

cva<-CVA(shape,groups=env,cv=TRUE) # Jackknife Cross-validation
cva

plot(cva$CVscores,col = env, bg = env,pch= 21,cex=2,xlab=paste("CV1  (",paste(round(cva$Var[1,2],1),"%)"),sep=""),ylab=paste("CV2  (",paste(round(cva$Var[2,2],1),"%)"),sep=""))
legend("top",legend=unique(env),pch=19,col=unique(env)) 

# Comparação das formas médias dos grupos entre si # 1 = cinza, 2 = preto # Exemplos:
plotRefToTarget(cva$groupmeans[,,"Savanna"],cva$groupmeans[,,"Forest"],links=links.tuco,method="points", mag = 4)
plotRefToTarget(cva$groupmeans[,,"Forest"],cva$groupmeans[,,"Savanna"],links=links.tuco,method="points", mag = 4)

