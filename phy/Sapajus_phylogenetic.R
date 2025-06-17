####### Analises filogeneticas de Sapajus ########
#install_github("YuLab-SMU/ggtree")


### Carregar pacotes
library(phangorn)
library(phylosignal)
library(phytools)
library(picante)
library(vegan)
library(ggplot2)
library(devtools)
library(ggtree)
library(tidyverse)
library(tidytree)
library(geiger)
library(tibble)


### Carregar árvore filogenética

multiphy <- read.nexus("trees/tree_completed/output.nex")
MCC <- mcc(multiphy,tree = TRUE, part = NULL, rooted = TRUE)
plot.phylo (MCC, type = "phylogram", show.tip.label = TRUE, 
            show.node.label = TRUE, edge.color = "black", edge.width = 0.5, 
            tip.color = "black", cex = 0.8, label.offset = 0.1, main = "Maximum Credibility Tree")

axisPhylo()

### Carregar dados

pgls_data <- read.csv("phy/sapajus_sp.csv",sep=";")
head(pgls_data)

species <- as.factor(pgls_data$sp)

rownames(pgls_data) <- pgls_data$Species
head(pgls_data)
species
MCC<- drop.tip(MCC,c("Cebus_albifrons", "Cebus_capucinus","Saimiri_boliviensis","Sapajus_flavius", "Sapajus_macrocephalus"))

name.check(MCC,set_names(pgls_data$LogCS,rownames(pgls_data)))

spp_size <- pgls_data$LogCS
spp_size <- setNames(spp_size,rownames(pgls_data))

view(spp_size)
### Reconstrução de Atributo (traço)

rec<-fastAnc(MCC,spp_size,vars=TRUE,CI=TRUE)
print(rec,printlen=10)
mean(rec$ace)

obj<-contMap(MCC,spp_size,plot=FALSE)
plot(obj,fsize=c(0.8,0.6),outline=FALSE,lwd=c(6,10),leg.txt="Size")
axisPhylo()
dev.off()
