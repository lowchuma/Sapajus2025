#### Consensus methods 
setwd("C:/Users/Lourenço/OneDrive/UFSM/LabMastozoo/Sapajus/Data")

library(phytools)
library(ape)
library(phangorn)

trees <- read.nexus("trees/Upham_2024.nex")
trees <- read.tree("trees/MCC_Lima_calibrated.tre")


mcctree <-maxCladeCred(trees, rooted = TRUE)
mcctree <- read.tree("trees/MCC_Upham24_calibrated.tre")
mcctree <- read.tree("trees/Sapajus_Lima.tre")

plot(mcctree, cex = 0.6)

tree<- drop.tip(mcctree,c("Sapajus_macrocephalus", "Cebus_albifrons", "Cebus_capucinus", "Sapajus_flavius"))

plot.phylo (tree, type = "phylogram", show.tip.label = TRUE, 
            show.node.label = TRUE, edge.color = "black", edge.width = 1.5, 
            tip.color = "gray15", cex = 1, label.offset = 0.25)

divergence_times <- branching.times(tree)
nodelabels(round(divergence_times, 2), cex=0.8)
axisPhylo()
dev.off()

write.tree(tree, file="trees/Sapajus_Upham.tre")

Tree <- read.tree("trees/Sapajus_Lima.tre")

TreeR <- read.tree("trees/MCC_Wright.tre")

plot(TreeR, type = "phylogram", cex = 0.8)
#plot(TreeR, type = "unrooted")   
axisPhylo()
### Sapajus

TreeU <- read.tree("trees/Sapajus_Upham.tre")
#plot(UTree,type = "unrooted")
plot(TreeU,type = "phylogram")
axisPhylo()

#write.tree(tree, file="MCC_Completed_Upham.tre")

### force.ultrametric method for ultrametric phylogenies that fail ape::is.ultrametric due to numerical precision (http://blog.phytools.org/2017/03/forceultrametric-method-for-ultrametric.html) 

is.ultrametric(TreeU)
is.ultrametric(TreeR)
is.ultrametric(Tree)

tree <- Tree

force.ultrametric<-function(tree,method=c("nnls","extend")){
  method<-method[1]
  if(method=="nnls") tree<-nnls.tree(cophenetic(tree),tree,
                                     rooted=TRUE,trace=0)
  else if(method=="extend"){
    h<-diag(vcv(tree))
    d<-max(h)-h
    ii<-sapply(1:Ntip(tree),function(x,y) which(y==x),
               y=tree$edge[,2])
    tree$edge.length[ii]<-tree$edge.length[ii]+d
  } else 
    cat("method not recognized: returning input tree\n\n")
  tree
}

ult.nnls<-force.ultrametric(tree) ## default method
is.ultrametric(ult.nnls)

ult.extend<-force.ultrametric(tree,method="extend")
is.ultrametric(ult.extend)

# Now let's compare the edge lengths of each of these two trees, as a function of edge height, to the edge lengths of our input tree:

ult.nnls<-reorder(ult.nnls)
h.nnls<-rowMeans(nodeHeights(ult.nnls))
h.extend<-rowMeans(nodeHeights(ult.extend))

plot(h.extend,tree$edge.length-ult.extend$edge.length,pch=21,
     ylim=c(-1e-6,1e-6),bg="grey",cex=1.5,xlab="edge height",
     ylab="difference between input & output edge lengths",
     main="force.ultrametric(...,method=\"extend\")")

## "nnls" method looks better

TreeU<-force.ultrametric(TreeU) ## default method
is.ultrametric(TreeU)

TreeR<-force.ultrametric(TreeR) ## default method
is.ultrametric(TreeR)

Tree<-force.ultrametric(tree) ## default method
is.ultrametric(Tree)

### Adding tips as a polytomy 
library(devtools)
#install_github("vanderleidebastiani/daee")
library(daee)

# Visualizar a nova árvore
Tree<- drop.tip(Tree,c("Sapajus_macrocephalus", "Cebus_albifrons", "Cebus_capucinus", "Sapajus_flavius"))
TreeR<- drop.tip(TreeR,c("Sapajus_macrocephalus"))

plot(TreeR, show.node.label=TRUE)

AddNode<-node.tree(Tree,m = 0,prefix = "N")$tree
plot(AddNode,show.node=T)
axisPhylo()
nodelabels()

#Tree<- drop.tip(Tree,c("Sapajus_macrocephalus", "Cebus_albifrons", "Cebus_capucinus"))
#write.tree(Tree,"Ultrametric_Upham_2024.tre")

AddNode<-node.tree(TreeU,m = 0,prefix = "N")$tree
plot(AddNode,show.node=T)
axisPhylo()

AddNode<-node.tree(TreeR,m = 0,prefix = "N")$tree
plot(AddNode,show.node=T)
axisPhylo()

### AvgLocSex

AddNode$tip.label
taxa<-matrix(data = NA,nrow = 77,ncol = 3)

taxa[,1]<- rep(c("Sapajus_nigritus","Sapajus_xanthosternos","Sapajus_cay","Sapajus_robustus","Sapajus_apella","Sapajus_libidinosus"),c(10,3,20,7,21,16))

taxa[1:10,2]<-paste("Sapajus_nigritus",2:11,sep = "_")
taxa[11:13,2]<-paste("Sapajus_xanthosternos",2:4,sep = "_")
taxa[14:33,2]<-paste("Sapajus_cay",2:21,sep = "_")
taxa[34:40,2]<-paste("Sapajus_robustus",2:8,sep = "_")
taxa[41:61,2]<-paste("Sapajus_apella",2:22,sep = "_")
taxa[62:77,2]<-paste("Sapajus_libidinosus",2:17,sep = "_")

plot(add.taxa.phylo(AddNode, taxa)$tree,cex=0.6)
axisPhylo()

write.tree(add.taxa.phylo(AddNode, taxa)$tree,"trees/AVGLOCSEX_Wright.tre")
write.tree(add.taxa.phylo(AddNode, taxa)$tree,"trees/AVGLOCSEX_Upham.tre")
write.tree(add.taxa.phylo(AddNode, taxa)$tree,"trees/AVGLOCSEX_Lima.tre")

### AVGLOCSEX

taxa <- matrix(data = NA, nrow = 55, ncol = 3)

taxa[,1]<- rep(c("Sapajus_nigritus","Sapajus_xanthosternos","Sapajus_cay","Sapajus_robustus","Sapajus_apella","Sapajus_libidinosus"),c(6,2,17,4,14,12))

taxa[1:6,2]<-paste("Sapajus_nigritus",2:7,sep = "_")
taxa[7:8,2]<-paste("Sapajus_xanthosternos",2:3,sep = "_")
taxa[9:25,2]<-paste("Sapajus_cay",2:18,sep = "_")
taxa[26:29,2]<-paste("Sapajus_robustus",2:5,sep = "_")
taxa[30:43,2]<-paste("Sapajus_apella",2:15,sep = "_")
taxa[44:55,2]<-paste("Sapajus_libidinosus",2:13,sep = "_")

plot(add.taxa.phylo(AddNode, taxa)$tree,cex=0.6)
axisPhylo()

write.tree(add.taxa.phylo(AddNode, taxa)$tree,"trees/avgloc_wright.tre")

### GLOBAL

taxa <- matrix(data = NA, nrow = 100, ncol = 3)

taxa[,1]<- rep(c("Sapajus_nigritus","Sapajus_xanthosternos","Sapajus_cay","Sapajus_robustus","Sapajus_apella","Sapajus_libidinosus"),c(14,5,27,11,23,20))

taxa[1:14,2]<-paste("Sapajus_nigritus",2:15,sep = "_")
taxa[15:19,2]<-paste("Sapajus_xanthosternos",2:6,sep = "_")
taxa[20:46,2]<-paste("Sapajus_cay",2:28,sep = "_")
taxa[47:57,2]<-paste("Sapajus_robustus",2:12,sep = "_")
taxa[58:80,2]<-paste("Sapajus_apella",2:24,sep = "_")
taxa[81:100,2]<-paste("Sapajus_libidinosus",2:21,sep = "_")

plot(add.taxa.phylo(AddNode, taxa)$tree,cex=0.6)
axisPhylo()

write.tree(add.taxa.phylo(AddNode, taxa)$tree,"trees/global_lima.tre")

### Sexed

factor <- read.csv("Global_M.csv", sep = ";")
summary (as.factor(factor$sp))
taxa <- matrix(data = NA, nrow = 51, ncol = 3)

taxa[,1]<- rep(c("Sapajus_nigritus","Sapajus_xanthosternos","Sapajus_cay","Sapajus_robustus","Sapajus_apella","Sapajus_libidinosus"),c(6,3,14,5,13,10))

taxa[1:6,2]<-paste("Sapajus_nigritus",2:7,sep = "_")
taxa[7:9,2]<-paste("Sapajus_xanthosternos",2:4,sep = "_")
taxa[10:23,2]<-paste("Sapajus_cay",2:15,sep = "_")
taxa[24:28,2]<-paste("Sapajus_robustus",2:6,sep = "_")
taxa[29:41,2]<-paste("Sapajus_apella",2:14,sep = "_")
taxa[42:51,2]<-paste("Sapajus_libidinosus",2:11,sep = "_")

plot(add.taxa.phylo(AddNode, taxa)$tree,cex=0.6)
axisPhylo()

write.tree(add.taxa.phylo(AddNode, taxa)$tree,"trees/Global_Lima_M.tre")
