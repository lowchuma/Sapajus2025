## Reconstruçãoo da evolução do atributo ao longo da filogenia ##

# Load required libraries
library(devtools)
library(gmodels)
library(ggplot2)
library(tidyverse)
library(phylosignal)
#install_github("uyedaj/treeplyr")
library(treeplyr)
library(ape)
library(caper)
library(picante)
library(phytools)
library(phylobase)
library(nlme)
library(geiger)
library(phangorn)
#windows(12,8)
## Reconstructing the evolution of eumelanic coloration throughout evolutionary history ##
citation ("phangorn")
## Importing color datasets and climate predictors
setwd("C:/Users/Lourenço/OneDrive/UFSM/LabMastozoo/Sapajus/Data")

data_color_pred<-read.csv("phy/sapajus_sp.csv", sep=";")
dim(data_color_pred)
head(data_color_pred)

##------- Importing phylogeny
#sapajus.tree <- read.tree("trees/MCC_Rtrees_calibrated.tre")
sapajus.tree <- read.tree("trees/MCC_Lima_calibrated.tre")
#sapajus.tree <- read.tree("trees/MCC_Wright.tre")

pgls_data <- data_color_pred
rownames(pgls_data) <- pgls_data$Species
# Importing only one consensus tree (MCC function)

plot(sapajus.tree)

tree <- sapajus.tree

# Maximum clade credibility tree

MCC<-maxCladeCred(sapajus.tree,rooted=TRUE)
random.tree<-sample(MCC,size=1)

##-------- Plotting phylogeny
plotTree(tree,ftype="i",fsize=0.8,lwd=2,type = "phylogram", cex = 0.8)

# Extract the color column (RGB/Lightness) 

color <- data_color_pred$LogCS

# Assign each species with the extracted variable 

color <- setNames(color,rownames(pgls_data))
color

# Distribution of the attribute throughout evolution

plotTree.barplot(tree,color)

tree<- drop.tip(tree,c("Sapajus_macrocephalus", "Cebus_albifrons", "Cebus_capucinus", "Sapajus_flavius"))
tree<- drop.tip(tree,c("Sapajus_macrocephalus"))
plot(tree)
# Ancestral reconstruction of continuous data 

# For each node attribute value

rec<-fastAnc(tree,color,vars=TRUE,CI=TRUE)
rec
print(rec,printlen=10)
dev.off()

plotTree(tree,fsize=1,ftype="i",ftype="i",lwd=1)
mean(rec$ace)

# Visualizing the reconstruction

rec.map<-contMap(tree,color, outline=FALSE)
rec.map

# customizing
library(viridis)
custom_palette <- plasma(1000)
rec.map <- setMap(rec.map, custom_palette)
layout(matrix(c(1, 2), ncol = 1), heights = c(6, 1))

# Margens mais controladas
par(mar = c(1, 4, 4, 2)) 
plot(rec.map,fsize=c(1,0.8),outline=FALSE,lwd=c(8,6),leg.txt="LogCS", type = "phylogram")
axisPhylo()
op <- par(no.readonly = TRUE)
par(op)
dev.off()

layout(matrix(c(1, 2), ncol = 1), heights = c(6, 1))
rec.map <- setMap(rec.map, custom_palette, invert = TRUE)
plot(rec.map,fsize=c(1,0.8),outline=FALSE,lwd=c(8,6),leg.txt="LogCS", type = "phylogram")
axisPhylo()
divergence_times <- branching.times(tree)
print(divergence_times)
nodelabels(round(divergence_times, 2), cex=0.8)

# Pacote RColorBrewer para paletas mais distintas
library(RColorBrewer)

# Definindo a paleta RdYlBu com mais contraste
custom_palette <- colorRampPalette(brewer.pal(n =11, name = "RdYlGn"))(400)
rec.map <- setMap(rec.map, custom_palette)

# Layout e margens ajustadas
layout(matrix(c(1, 2), ncol = 1), heights = c(6, 1))
par(mar = c(1, 4, 4, 2)) 
plot(rec.map, fsize = c(1, 0.8), outline = FALSE, lwd = c(8, 6), leg.txt = "LogCS", type = "phylogram")
axisPhylo()

# Reiniciando configurações para plot invertido, se desejado
par(op)
layout(matrix(c(1, 2), ncol = 1), heights = c(6, 1))
rec.map <- setMap(rec.map, custom_palette, invert = TRUE)
plot(rec.map, fsize = c(1, 0.8), outline = FALSE, lwd = c(8, 6), leg.txt = "LogCS", type = "phylogram")
axisPhylo()

# Exibir tempos de divergência
divergence_times <- branching.times(tree)
print(divergence_times)
nodelabels(round(divergence_times, 2), cex = 0.8)

# Carregar RColorBrewer e ajustar paleta sem o branco
library(RColorBrewer)

# Ajustar 'Spectral' com 9 cores e expandir sem branco central
custom_palette <- colorRampPalette(brewer.pal(n = 10, name = "Spectral")[c(-5)])(1000)
#custom_palette <- colorRampPalette(c("darkgreen", "#FFDB58"))(1000)  # Do verde ao amarelo

rec.map <- setMap(rec.map, custom_palette)

# Layout e ajustes de plot
layout(matrix(c(1, 2), ncol = 1), heights = c(6, 1))
par(mar = c(1, 4, 4, 2)) 
plot(rec.map, fsize = c(1.5, 1), outline = FALSE, lwd = c(15, 5), leg.txt = "LogCS", type = "phylogram")
axisPhylo()

# Reiniciando configurações para plot invertido, se desejado
par(op)
layout(matrix(c(1, 2), ncol = 1), heights = c(6, 1))
rec.map <- setMap(rec.map, custom_palette, invert = TRUE)
plot(rec.map, fsize = c(1, 0.8), outline = FALSE, lwd = c(8, 6), leg.txt = "LogCS", type = "phylogram")
axisPhylo()

# Exibir tempos de divergência
divergence_times <- branching.times(tree)
print(divergence_times)
nodelabels(round(divergence_times, 2), cex = 0.8)

# Adding a color error bar

rec.map<-contMap(tree,color,fsize=1.2,ftype="i", invert=TRUE, type = "phylogram")

obj<-contMap(tree,color,plot=FALSE)

plot(obj,fsize=c(0.26,0.6),outline=FALSE)
plot(obj,fsize=c(1,0.8),outline=FALSE,lwd=c(8,6),leg.txt="LogCS")
errorbar.contMap(rec.map)


#add.scale.bar()
#nodelabels()
#tiplabels()

## Reconstructing the evolution of feomelanic coloration throughout evolutionary history ##

## Importing color datasets (redness) and climate predictors

data_color_pred<-read.csv("Color_Redness_All_Species_final.csv")
dim(data_color_pred)
View(data_color_pred)

# Removing species with 1 phot

data_color_pred<-data_color_pred[data_color_pred$pic_num>1,]
dim(data_color_pred)


# Removing species that shed their coats in winter

table(data_color_pred$Troca.de.pelagem)
data_color_pred<-data_color_pred[!data_color_pred$Troca.de.pelagem=="S",]
table(data_color_pred$Troca.de.pelagem)
head(data_color_pred)
rownames(data_color_pred)<- data_color_pred$Species
pgls_data<- data_color_pred[,c(2,4)]
head(pgls_data)


##------- Importing phylogeny
Lagomorpha.tree<-read.nexus("Lagomorpha_tree_1000.nex") 
Lagomorpha.tree


# Importing only one consensus tree (MCC function)

plot(Lagomorpha.tree$tree_7555)


# Maximum clade credibility tree

MCC<-maxCladeCred(Accipitriformes.tree,rooted=FALSE)

##-------- Phylogeny parameters
str(MCC)
MCC$edge
MCC$tip.label
MCC$edge.length
MCC$Nnode

##-------- Plotting phylogeny
plotTree(MCC,ftype="i",fsize=0.44,lwd=1,type = "phylogram")


## Check correspondence between species in the phylogeny and in the data (name.check function)

name.check(Accipitriformes.tree[[1]],set_names(pgls_data$RGB_mean,rownames(pgls_data)))
names(pgls_data)

# Extract the coloring column (RGB/Lightness) 

color<-pgls_data$RGB_mean

# Assign each species with the extracted variable 

color<-setNames(color,rownames(pgls_data))
color

# Distribution of the attribute throughout evolution

plotTree.barplot(MCC,log(color))

# Ancestral reconstruction of continuous data 

# For each node attribute value

rec<-fastAnc(MCC,log(color),vars=TRUE,CI=TRUE)
rec
print(rec,printlen=10)
plotTree(MCC,fsize=0.5,ftype="i",ftype="i",lwd=1)

# Visualizing the reconstruction

rec.map<-contMap(MCC,log(color), outline=FALSE)
rec.map


# customizing
windows(12,8)
rec.map<-setMap(rec.map,invert=TRUE) # inverte cores

plot(rec.map,ftype="i",fsize=c(0.5,0.6),lwd=c(5,7), type = "phylogram", leg.txt = "Size")

# Adding a color error bar

rec.map<-contMap(MCC,log(color),fsize=0.4,ftype="i")

obj<-contMap(MCC,color,plot=FALSE)
obj<-setMap(obj,invert=TRUE)
plot(obj,fsize=c(0.26,0.6),outline=FALSE)
plot(obj,fsize=c(0.8,0.5),outline=FALSE,lwd=c(7,10),leg.txt="Centroid Size")
