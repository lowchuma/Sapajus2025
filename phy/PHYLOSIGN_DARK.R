
###---- SINAL FILOGENÉTICO P/ 1000 ÁRVORES ----###

library(vegan)
library(ggplot2)
library(GGally)
library(ggpubr)
library(picante)
library(phytools)
library(ape)
library(geiger)
library(phyloregion)
library(pez)
library(reshape2)
library(betapart)
library(phangorn)
library(ggtree)
library(treeio)
library(tidytree)
library(phylosignal)
library(tidyverse)
#install_github("YuLab-SMU/treeio")
library(rtrees)

#windows(12,8)

## Confere a filogenia pra não dar brete
setwd("C:/Users/Lourenço/OneDrive/UFSM/LabMastozoo/Sapajus/Data")

sapajus.tree <- read.nexus("trees/Upham_2024.nex")
pgls_data <- read.csv("phy/sapajus_sp.csv", sep=";")
names(pgls_data)
bioclim_data<-pgls_data[,c(1,4:22)]
bioclim_data[,2:20] <- scale(bioclim_data[,2:20])

rownames(pgls_data) <- pgls_data$Species
rownames(bioclim_data)<-pgls_data$Species

mcc <- mcc(sapajus.tree, tree = TRUE, part = NULL, rooted = TRUE)
mcc<- drop.tip(mcc,c("Sapajus_macrocephalus", "Cebus_albifrons", "Cebus_capucinus"))
name.check(mcc,set_names(bioclim_data$bio2,rownames(pgls_data)))


# Criando a matriz

matriz<- data.frame(matrix(NA,nrow = 1000,ncol = 3))

matriz[,1]<-seq(1:1000)
colnames(matriz)<-c("TreeNumber","K","P")
head(matriz)


library(phylosignal)


###---- Rodando o for pra lambda das 1000 árvores

i=1

for (i in 1:length(sapajus.tree)) {
  K_value<-phylosig(tree = sapajus.tree[[i]],x = setNames(pgls_data$LogCS,row.names(pgls_data)),method = "K",test = T)
  
  #matriz[i] <- lambda_value[1]
  
  matriz[i,2]<- K_value$K # valor de lambda
  matriz[i,3]<- K_value$P # valor do p 
  print(i)
  
  
}

print(K_value)
mean(matriz$K)
mean(matriz$P)
hist(matriz$P)
hist(matriz$K)


hist(matriz$K, xlab="Blomberg's K", ylab="Frequency", main= "K mean ",
     col = c("grey80"), border = FALSE)

abline(v=mean(matriz$K), col="blue",lwd=2.5,lty="dashed")

legend(x="topright", # posição da legenda
       c("Mean"), # nomes da legenda
       col=c("blue","red"), #cores
       lty=c(2,2), # estilo da linha
       lwd=c(2,2,4))# grossura das linhas



###---- FAZENDO O HISTOGRAMA PRA O LAMBDA

gmodels::ci(matriz$K)

ggplot(matriz, aes(P))+
  geom_histogram(fill= 'grey80', border = TRUE, bins = 8)+
  theme_classic()+
  labs(x="Blomberg's K", title = "Size", y = "Frequency")+
  geom_vline(xintercept =0.465382379, col="black",lwd=0.4, lty= "dashed")+
  geom_vline(xintercept =0.468819076 , col="blue",lwd=0.4, lty= "solid")+
  geom_vline(xintercept =0.472255772, col="black",lwd=0.4, lty="dashed")

###---- FAZENDO O HISTOGRAMA PRO VALOR DE P


gmodels::ci(matriz$P)


ggplot(matriz, aes(P))+
  geom_histogram(fill= 'grey80', border =TRUE, bins = 8)+
  theme_classic()+
  labs(x="P-value", y = "Frequency", title = "Size")+
  geom_vline(xintercept =  0.73447168, col="black",lwd=0.4, lty= "dashed")+
  geom_vline(xintercept = 0.73829280, col="blue",lwd=0.4, lty= "solid")+
  geom_vline(xintercept = 0.74211392, col="black",lwd=0.4, lty="dashed")

################################################################################

    ####------ sinal filogenético !!!!!!ÁRVORE CONSENSO!!!!! ----------####
MCC <- maxCladeCred(sapajus.tree,rooted=TRUE)


#------- Aqui vão ser os dados (SPP, COR, BIO's)-------##

pgls_data <- read.csv("phy/sapajus_sp.csv", sep=";")
dim(pgls_data)
head(pgls_data)
View(pgls_data)

##-- Mema fita pra arrumar os dados

rownames(pgls_data) <- pgls_data$Species


##-- Dá aquele check na filogenia de novo pra não dar brete

name.check(mcc,setNames(pgls_data$Species,rownames(pgls_data)))
pgls_data <- pgls_data[,c(2,3,4,5,6,7)]
head(pgls_data)
View(pgls_data)
dim(pgls_data)

##--- Faz o valor de lambda pra árvere consenso

lambda_value <- phylosig(tree = mcc,x = setNames(pgls_data$LogCS,row.names(pgls_data)),method = "K",test = T)

print(lambda_value)

