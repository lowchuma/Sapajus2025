## Análises filogenéticas ##
library(ape)
library(phytools)
library(devtools)
library(ggtree)
library(tidyverse)
library(phangorn)
library(vegan)
library(ggplot2)
library(GGally)
library(ggpubr)
library(picante)
library(geiger)
library(phyloregion)
library(pez)
library(reshape2)
library(betapart)
#install_github("daijiang/rtrees")
library(rtrees)

setwd("C:/Users/Lourenço/OneDrive/UFSM/LabMastozoo/Sapajus/Data")
test_tree = get_tree(sp_list = c("Sapajus_apella","Sapajus_cay","Sapajus_libidinosus","Sapajus_nigritus","Sapajus_robustus","Sapajus_xanthosternos"),
                     taxon = "mammal",
                     scenario = "at_basal_node",
                     show_grafted = TRUE)

credibility <- mcc(test_tree, tree = TRUE, part = NULL, rooted = TRUE)
credibility <- ladderize(credibility, right = FALSE)
plot.phylo (credibility, type = "phylogram", show.tip.label = TRUE, 
            show.node.label = TRUE, edge.color = "black", edge.width = 1.5, 
            tip.color = "gray15", cex = 1, label.offset = 0.15)
axisPhylo()

arvore <- read.nexus("trees/Upham_2024.nex")

mcc <- mcc(arvore, tree = TRUE, part = NULL, rooted = TRUE)
mcc<- drop.tip(mcc,c("Sapajus_macrocephalus", "Cebus_albifrons", "Cebus_capucinus"))
name.check (mcc, dados)
plot(mcc, show.tip.label = TRUE, 
     show.node.label = TRUE)

plot.phylo (mcc, type = "phylogram", show.tip.label = TRUE, 
            show.node.label = TRUE, edge.color = "black", edge.width = 1.5, 
            tip.color = "gray15", cex = 1, label.offset = 0.15)
axisPhylo()

names(mcc)
mcc$tip.label
mcc$edge.length
mcc$Nnode

strict_consensus <- consensus(test_tree)
majority_consensus <- consensus(test_tree, p=.5)
all_compat <- allCompat(test_tree)
max_clade_cred <- maxCladeCred(test_tree)

strict_consensus <- consensus(arvore)
majority_consensus <- consensus(arvore, p=.5)
all_compat <- allCompat(arvore)
max_clade_cred <- maxCladeCred(arvore)

old.par <- par(no.readonly = TRUE)
par(mfrow = c(2,2), mar = c(1,4,1,1))

plot(strict_consensus, main="Strict consensus tree")
plot(majority_consensus, main="Majority consensus tree")
plot(all_compat, main="Majority consensus tree with compatible splits")
plot(max_clade_cred, main="Maximum clade credibility tree")
par(old.par)

# compute clade credibility for trees given a prop.part object
pp <- prop.part(arvore)
tree <- rNNI(arvore[[1]], 20)
maxCladeCred(c(tree, arvore), tree=FALSE, part = pp)

# first value likely be -Inf


# Calcular mcc# Calcular a matriz de covariância da estrutura filogenética
phy_cov <- vcv(mcc)
dist_matrix <- cophenetic.phylo(mcc)

phy_cov <- vcv(credibility)
dist_matrix <- cophenetic.phylo(credibility)

# Executar a análise de clustering hierárquico
hc <- hclust(as.dist(dist_matrix))

# Plotar o dendrograma
plot(hc, main = "Clustering Hierárquico", cex = 0.6)

# Plotar a árvore para identificar os nós
plot(mcc)
plot(credibility)
nodelabels()

# Verificar se todos os nós têm valores numéricos
all_nodes <- c(tree$edge[,1], tree$edge[,2])
if (!all(is.numeric(all_nodes))) {
  stop("Todos os nós devem ser numéricos.")
}

# Definir calibração (supondo que o nó 7 é o MRCA de taxon1 e taxon2)
calibration <- makeChronosCalib(mcc, node=c(7,8,10, 11), 
                                age.min=c(2.65,2.36,0.84, 0.18), 
                                age.max=c(2.91,2.66,1.2, 0.22))

calibration <- makeChronosCalib(credibility, node=c(7,8,10, 11), 
                                age.min=c(2.65, 1.21, 0.18, 2.36), 
                                age.max=c(2.91,1.5,0.22, 2.66))

calibration <- makeChronosCalib(mcc, node=c(10,11,13,16, 17), 
                                age.min=c(6.6,2.44, 0.84, 0.18, 2.46), 
                                age.max=c(6.8,2.91,1.2,0.22,2.65))
# Aplicar a calibração à árvore
chronos_tree <- chronos(mcc, calibration=calibration, lambda = 0.1)
chronos_tree <- chronos(credibility, calibration=calibration, lambda = 0.1)

# Obter e visualizar os tempos de divergência
divergence_times <- branching.times(chronos_tree)
print(divergence_times)

# Plotar a árvore com tempos de divergência
plot(chronos_tree, show.tip.label=TRUE, cex = 0.8)
nodelabels(round(divergence_times, 2), cex=0.8)

axisPhylo()
dev.off()

write.tree(chronos_tree,"trees/MCC_Upham24_calibrated.tre")
write.tree(chronos_tree,"trees/MCC_Rtrees_calibrated.tre")
