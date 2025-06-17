##  PGLS muultipla (MCC)

# Load packages
library(ape)
library(picante)
library(phytools)
library(phylobase)
library(nlme)
library(geiger)
library(phangorn)

##------- Importing phylogeny
Lagomorpha<-read.tree("Phylogenetics/Sapajus_10000_all_spp/Yves_Sapajus_Wright_ultrametric.tre")
Lagomorpha<-read.tree("Phylogenetics/Sapajus_10000_all_spp/MCC_Sapajus_NodeDated_10000_PGLSIves.tre")Lagomorpha
Lagomorpha<-read.tree("trees/MCC_Ives_Upham.tre")
plot(Lagomorpha, cex = 0.4)

# Maximum Clade Credibility tree (MCC)

MCC<-maxCladeCred(Lagomorpha,rooted=TRUE)
MCC<-Lagomorpha
##-------- Plotting phylogeny MCC
plot(MCC,edge.width = 2,label.offset = 0.05,type = "phylogram", cex = 0.3)
nodelabels()
tiplabels()
axisPhylo()

##------- Importing data-------##

Color_data<-read.csv("avglocmatrix.csv", sep=";")
dim(Color_data)
head(Color_data)
View(Color_data)
rownames(Color_data)<- Color_data$sp_ives
species <- as.factor(Color_data$sp)
pgls_data<- Color_data[,c(13,14:34)]
rownames(pgls_data)<- Color_data$sp_ives
head(pgls_data)
View(pgls_data)
dim(pgls_data)
names(pgls_data)

##-------- Checando se os dados da filogenia e coloracao batem ------- ###

name.check(MCC,setNames(pgls_data$LogCS,rownames(pgls_data)))
keep_tips <- function(tree, tip) drop.tip(tree, setdiff(tree$tip.label, tip))
Clean_Tree_temp<- keep.tip(phy = Lagomorpha, tip = rownames(pgls_data))
name.check(Clean_Tree_temp,setNames(pgls_data$LogCS,rownames(pgls_data)))

## size vs. Temperatura ## 

pglsModel_temp_BM<- gls(LogCS~bio4+bio12,correlation = corBrownian(phy=Clean_Tree_temp, form = ~species),data=data.frame(scale(pgls_data)),method = "ML")
#pglsModel_temp_BM<- gls(LogCS~BIO1+BIO12,correlation = corBrownian(phy=Lagomorpha, form = ~species),data=pgls_data,method = "ML")
anova(pglsModel_temp_BM)
coef(pglsModel_temp_BM)
summary(pglsModel_temp_BM)

## Pegamos a vari?vel precipita??o (BIO12) deixando o efeito da vari?vel temperatura constante.

pglsModel_prec_BM<- gls(LogCS~bio12,correlation = corBrownian(phy=Clean_Tree_temp, form = ~species),data=data.frame(scale(pgls_data)),method = "ML")
anova(pglsModel_prec_BM)

# Calculamos o valor do residuo, ou seja, o que n?o foi explicado pela vari?vel precipita??o 
pgls_data$residuals<-pglsModel_temp_BM$residuals
# Geramos um gr?fico associando o residuo da temperatura + efeito da vari?vel temperatura
plot(pglsModel_temp_BM$residuals~pgls_data$bio12)
abline (gls(residuals~bio12,correlation = corBrownian(phy=Lagomorpha, form = ~species),data=pgls_data,method = "ML"))

plot(LogCS~bio4, pch = as.numeric(Color_data$pch), col = as.factor(Color_data$biome),data = pgls_data)
abline (gls(LogCS~bio4,correlation = corBrownian(phy=Lagomorpha, form = ~species), data = pgls_data,method = "ML"))



## size vs. precipitaçao ## 

pglsModel_prec_BM<- gls(LogCS~bio12+bio4,correlation = corBrownian(phy=Clean_Tree_temp, form = ~species),data=data.frame(scale(pgls_data)),method = "ML")
pglsModel_prec_BM<- gls(LogCS~BIO12+BIO1,correlation = corBrownian(phy=Clean_Tree_temp, form = ~species),data=pgls_data,method = "ML")

anova(pglsModel_prec_BM)
coef(pglsModel_prec_BM)
summary(pglsModel_prec_BM)


## Pegamos a vari?vel temperatura (BIO1) deixando o efeito da vari?vel precipita??o constante.


pglsModel_prec_BM<- gls(LogCS~bio4,correlation = corBrownian(phy=Clean_Tree_temp, form = ~species),data=data.frame(scale(pgls_data)),method = "ML")

# Calculamos o valor do residuo, ou seja, o que n?o foi explicado pela vari?vel temperatura 
pgls_data$residualsPREC<-pglsModel_prec_BM$residuals
# Geramos um gr?fico associando o residuo da precipita??o + efeito da vari?vel precipita??o
plot(pglsModel_prec_BM$residuals~pgls_data$bio4)
abline(gls(residualsPREC~bio4,correlation = corBrownian(phy=Clean_Tree_temp, form = ~species),data=pgls_data,method = "ML"))

plot(LogCS~bio12, pch = as.numeric(Color_data$pch), col = as.factor(Color_data$biome),data = pgls_data)
abline (gls(LogCS~bio12,correlation = corBrownian(phy=Lagomorpha, form = ~species), data = pgls_data,method = "ML"))


