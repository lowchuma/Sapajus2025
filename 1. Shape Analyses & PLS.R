#================ Common basic part

#install_github("geomorphR/geomorph", ref = "Stable", build_vignettes = TRUE)

#install_github("mlcollyer/RRPP")
#install_github("fawda123/ggord")
#install.packages("geiger")
#install.packages("RColorBrewer")
#install.packages("vegan")
#install.packages("tidyverse")
#install.packages("rgl")

library(devtools)
library(tidyverse) # Set of packages for data science, including ggplot2, dplyr, tidyr, readr, purrr, and others
library(reshape2) # Functions for reshaping data between wide and long formats

# Statistics and Modeling
library(MASS) # Functions and datasets for applied statistics
library(psych) # Tools for psychological and psychometric analysis
library(klaR) # Functions for classification and visualization, including regularized discriminant analysis

# Data Visualization
library(ggplot2) # System for creating layer-based graphs
library(RColorBrewer) # Color palettes for visualization
library(ggord) # Creating ordination graphs with ggplot2
library(dotwhisker) # Coefficient graphs for statistical models

# Ecology and Evolutionary Biology
library(vegan) # Ordination methods and diversity analysis for community ecology
library(ape) # Phylogenetic and evolutionary analyses
library(phytools) # Tools for phylogenetic comparative biology
library(picante) # Integration of phylogenies and ecology
library(geiger) # Statistical methods for analyzing phylogenetic data

# Morphometrics Geometric
library(geomorph) # Geometric morphometric analysis of 2D and 3D landmark data
library(Morpho) # Tools for geometric morphometrics and mesh processing
library(shapes) # Routines for statistical analysis of landmark-based shapes
library(Rvcg) # 3D mesh processing and analysis

# Species Distribution Modeling
library(usdm) # Uncertainty analysis for species distribution models

### Read the raw data and create the clustering factors ###

## Samples, in this case, are data from each sampled location
# (where there was more than one specimen, the average was calculated)

dados <- readland.tps("tps/avglocsex.tps", specID =  "ID")

factors <- read.csv("Plans/avglocsex.csv", sep=";")

names(factors)

species <- as.factor(factors$sp)
summary(species)

biome <- as.factor(factors$biome)
summary(biome)

sex <- as.factor(factors$sex)
summary (sex)

latitude <-as.numeric(factors$lat)

longitude <- as.numeric(factors$long)

### Run GPA - Aligns everything and takes the impact of raw data dimensionality

plot(dados)
gpa <- gpagen(dados)
link<-read.table("tps/link.txt") 
plot(gpa,link=link)

shp <- two.d.array(gpa$coords) 
cov_matrix <- cov(shp)  
print(cov_matrix)

### When necessary, find the mean specimen to draw the outline.
tps_0 <- readland.tps("tps/jaw_18LM.tps")
Y.gpa <- gpagen(tps_0)
findMeanSpec(Y.gpa$coords) # 90

plotOutliers(Y.gpa$coords, inspect.outliers = TRUE) ## verification

### Separate data by sex for separate analyses (other methods exist)

tps_M <- readland.tps("tps/avgloc_M.tps")
M_gpa <- gpagen(tps_M)
Males <- read.csv("Plans/avglocsex_M.csv", sep =";")
tps_F <- readland.tps("tps/avgloc_F.tps")
F_gpa <- gpagen(tps_F)
Females <- read.csv("Plans/avglocsex_F.csv", sep =";")

### Load the outline
drawinglandmark<-readland.tps("outline2/outline.tps")
outline<-read.table("outline2/outline.txt", header=FALSE)
summary(drawinglandmark)

mshape<-mshape(Y.gpa$coords)
Sapajusoutline<-warpRefOutline(file = "outline2/outline.txt", drawinglandmark[,,1],mshape)#run dev.off() in case of error message ##set outline configuration

grid.pars <- gridPar(grid.col = "white") ## if you want only the "invisible" grids


GP <- gridPar(n.col.cell = 100, pt.bg = "gray", pt.size = 0.8, tar.pt.bg = "cyan",
tar.pt.size = 0.8, tar.out.col = "gray10", tar.out.cex = 0.5,
grid.col = "white", grid.lwd = 0.5, txt.pos = 1, txt.col = "steelblue") ## Custom grids

plotRefToTarget(mshape,Y.gpa$coords[,,90],outline = Sapajusoutline$outline,method="points", gridPars=GP)

dev.off()

#=============== ORDERING & EXPLORATION #================#

#================ PCA #================####

# Exploratory

PCA <- gm.prcomp(gpa$coords)
summary(PCA)

### Basically, the parameters to customize the plots

symbols <- c(
Sapajus_apella = "\u25A0", # Square
Sapajus_cay = "\u25BC", # Downward-pointing triangle
Sapajus_libidinosus = "\u25B2", # Upward-pointing triangle
Sapajus_nigritus = "\u25CF", # Circle
Sapajus_robustus = "\u2726", # 6-pointed star
Sapajus_xanthosternos = "\u2736" # 8-pointed star
)
symbols <- symbols[as.character(species)]

Fac <- as.factor(factors$fac) ## Environment factor

levels_fac <- levels(Fac)

colors_fac <- c("green3", "darkgreen", "goldenrod1")
names(colors_fac) <- levels_fac

indv <- 1:nrow(PCA$x)

pca_df <- data.frame( 
PC1 = PCA$x[, "Comp1"], 
PC2 = PCA$x[, "Comp2"], 
PC3 = PCA$x[, "Comp3"], 
Species = species, 
Individual = 1:nrow(PCA$x), Fac = Fac, Sex = sex)

species_colors <- c("red","navy","cyan3","saddlebrown","magenta4","black")
names(species_colors) <- levels(pca_df$Species)

pch <- as.numeric(factors$pch)
uPCH <- unique(pch)
names(uPCH) <- unique(pca_df$Species)
colors_fac <- c(AF = "green3", AM = "darkgreen", SV = "goldenrod1")[as.numeric(as.factor(Fac))]

# PCA Plot
mat<-matrix(c(4,5,0,1,1,2,1,1,3),3) # Split the plot window
layout(mat, widths=c(1,1,1), heights=c(1,1,0.6))
par(mar=c(4, 4, 1, 1))
xlab <- paste("PC1 (", round(100 * PCA$sdev[1]^2 / sum(PCA$sdev^2), 1), "%)", sep = "")
ylab <- paste("PC2 (", round(100 * PCA$sdev[2]^2 / sum(PCA$sdev^2), 1), "%)", sep = "")

plot(PCA$x[,1],PCA$x[,2],cex=4,pch = 21, bg=colors_fac, xlab=xlab,ylab=ylab, 
col = colors_fac)
legend("topright",legend= unique(Fac),pch=19,col = unique(colors_fac), cex = 2, pt.cex = 3)

plotRefToTarget(mshape,PCA$shapes$shapes.comp1$min, mag = 2,outline = Sapajusoutline$outline,method="points",gridPars=GP)
plotRefToTarget(mshape,PCA$shapes$shapes.comp1$max, mag = 2,outline = Sapajusoutline$outline,method="points",gridPars=GP)
plotRefToTarget(mshape,PCA$shapes$shapes.comp2$max, mag = 2,outline = Sapajusoutline$outline,method="points",gridPars=GP)
plotRefToTarget(mshape,PCA$shapes$shapes.comp2$min, mag = 2,outline = Sapajusoutline$outline,method="points",gridPars=GP)
par(mfrow=c(1,1)) # default window

# Generates 95% ellipses

plot(PCA$x[,1],PCA$x[,2],cex=4,pch = 19,col = colors_fac,xlab=xlab,ylab=ylab)

ordiellipse(PCA$x,group = Fac, kind="sd", conf=0.95, col = unique(colors_fac))

ordihull(PCA$x,group = Fac, col = unique(Fac))

########### Complex plots with "ggplot2" ########

compute_hull <- function(pca_df) pca_df[chull(pca_df$PC1, pca_df$PC2), ]
hulls <- pca_df %>%
  group_by(Species) %>%
  do(compute_hull(.))

ggplot(pca_df, aes(x = PC1, y = PC2)) +
  geom_point(aes(color = Species, shape = Species), size = 5) +  
  geom_text(aes(label = Individual), size = 3, vjust = 1.5) +  
  geom_polygon(data = hulls, aes(group = Species, fill = Species, color = Species), alpha = 0.15) +  
  labs(title = "Shape's Principal Components", x = paste("PC1 (", round(100 * PCA$sdev[1]^2 / sum(PCA$sdev^2), 1), "%)", sep = ""), y = paste("PC2 (", round(100 * PCA$sdev[2]^2 / sum(PCA$sdev^2), 1), "%)", sep = ""), color = "Species", shape = "Biomes") +
  theme_minimal() +
  scale_color_manual(values = species_colors, name = "Species") +  
  scale_fill_manual(values = species_colors, name = "Species") +  
  scale_shape_manual(values = uPCH, name = "Species")

compute_hull <- function(pca_df) pca_df[chull(pca_df$PC1, pca_df$PC2), ]
hulls <- pca_df %>%
  group_by(Fac) %>%
  do(compute_hull(.))
ggplot(pca_df, aes(x = PC1, y = PC2)) +
  geom_point(aes(shape = Species, color = Fac), size = 6) +  
  geom_text(aes(label = Individual), size = 2, vjust = 2) +  
  geom_polygon(data = hulls, aes(x = PC1, y = PC2, group = Fac, fill = Fac, color = Fac), alpha = 0.2) +  
  labs(title = "Shape's Principal Components", 
       x = paste("PC1 (", round(100 * PCA$sdev[1]^2 / sum(PCA$sdev^2), 1), "%)", sep = ""), 
       y = paste("PC2 (", round(100 * PCA$sdev[2]^2 / sum(PCA$sdev^2), 1), "%)", sep = ""), 
       shape = "Species", color = "Fac") +
  theme_minimal() +
  scale_shape_manual(values = uPCH, labels = levels(pca_df$Species), name = "Species") +  
  scale_color_manual(values = colors_fac, name = "Environment") +  
  scale_fill_manual(values = colors_fac, name = "Environment")

compute_hull <- function(pca_df) pca_df[chull(pca_df$PC2, pca_df$PC3), ]
hulls <- pca_df %>%
  group_by(Fac) %>%
  do(compute_hull(.))
ggplot(pca_df, aes(x = PC2, y = PC3)) +
  geom_point(aes(shape = Species, color = Fac), size = 8) +  
  geom_text(aes(label = Individual), size = 3, vjust = 2) +  
  geom_polygon(data = hulls, aes(x = PC2, y = PC3, group = Fac, fill = Fac, color = Fac), alpha = 0.2) +  
  labs(title = NULL, 
       x = paste("PC2 (", round(100 * PCA$sdev[2]^2 / sum(PCA$sdev^2), 1), "%)", sep = ""), 
       y = paste("PC3 (", round(100 * PCA$sdev[3]^2 / sum(PCA$sdev^2), 1), "%)", sep = ""), 
       shape = "Species", color = "Fac") +
  theme_minimal() +
  scale_shape_manual(values = uPCH, labels = levels(pca_df$Species), name = "Species") +  
  scale_color_manual(values = colors_fac, name = "Environment") +  
  scale_fill_manual(values = colors_fac, name = "Environment")


levels_sex <- levels(sex)

colors_sex <- c("salmon","steelblue")
names(colors_sex) <- levels_sex

compute_hull <- function(pca_df) pca_df[chull(pca_df$PC1, pca_df$PC2), ]
hulls <- pca_df %>%
  group_by(Sex) %>%
  do(compute_hull(.))

ggplot(pca_df, aes(x = PC1, y = PC2)) +
  geom_point(aes(shape = Species, color = Sex), size = 10) +  
  geom_text(aes(label = Individual), size = 2, vjust = 2) +  
  geom_polygon(data = hulls, aes(x = PC1, y = PC2, group = Sex, fill = Sex, color = Sex), alpha = 0.2) +  
  labs(
    title = "Shape's Principal Components", 
    x = paste("PC1 (", round(100 * PCA$sdev[1]^2 / sum(PCA$sdev^2), 1), "%)", sep = ""), 
    y = paste("PC2 (", round(100 * PCA$sdev[2]^2 / sum(PCA$sdev^2), 1), "%)", sep = ""), 
    shape = "Species", color = "Sex"
  ) +
 theme_minimal() + 
scale_shape_manual(values ​​= uPCH, labels = levels(pca_df$Species), name = "Species") + 
scale_color_manual(values ​​= colors_sex, name = "Sex") + # Define `colors_sex` with the colors corresponding to the sexes 
scale_fill_manual(values ​​= colors_sex, name = "Sex")

PC1 <- PCA$x[, 1] # Scores for PC1
PC2 <- PCA$x[, 2] # Scores for PC2
PC3 <- PCA$x[, 3]

preds_comb <- shape.predictor( 
gpa$coords, 
x = cbind(PC1, PC2), # Use PC1 and PC2 together 
Intercept = FALSE, 
pred1 = c(min(PC1), min(PC2)), # Min PC1 + Min PC2 
pred2 = c(max(PC1), min(PC2)), # Max PC1 + Min PC2
pred3 = c(min(PC1), max(PC2)), # Min PC1 + Max PC2
pred4 = c(max(PC1), max(PC2)), # Max PC1 + Max PC2
pred5 = c(min(PC2), min(PC3)),
pred6 = c(min(PC2), max(PC3)),
pred7 = c(max(PC2), min(PC3)),
pred8 = c(max(PC2), max(PC3))
)

M <- mshape(gpa$coords)

# Display the shapes associated with each combination of component extremes

par(mfrow = c(2, 2))

plotRefToTarget(M, preds_comb$pred3, main = "Min PC1 + Max PC2", mag = 2, outline = Sapajusoutline$outline,gridPars = GP, method = "points")
plotRefToTarget(M, preds_comb$pred4, main = "Max PC1 + Max PC2", mag = 2, outline = Sapajusoutline$outline, gridPars = GP, method = "points")
plotRefToTarget(M, preds_comb$pred1, main = "Min PC1 + Min PC2", mag = 2, outline = Sapajusoutline$outline, gridPars = GP, method = "points")
plotRefToTarget(M, preds_comb$pred2, main = "Max PC1 + Min PC2", mag = 2, outline = Sapajusoutline$outline, gridPars = GP, method = "points")

### For PC3
plotRefToTarget(M, preds_comb$pred6, main = "Min PC2 + Max PC3", mag = 2, outline = Sapajusoutline$outline, gridPars = GP, method = "points")
plotRefToTarget(M, preds_comb$pred8, main = "Max PC2 + Max PC3", mag = 2, outline = Sapajusoutline$outline, gridPars = GP, method = "points")
plotRefToTarget(M, preds_comb$pred5, main = "Min PC2 + Min PC3", mag = 2, outline = Sapajusoutline$outline, gridPars = GP, method = "points")
plotRefToTarget(M, preds_comb$pred7, main = "Max PC2 + Min PC3", mag = 2, outline = Sapajusoutline$outline, gridPars = GP, method = "points")

par(mfrow = c(1, 1))

####### Original PCA separated by sex ###
compute_hull <- function(pca_df) pca_df[chull(pca_df$PC1, pca_df$PC2), ]

hulls <- pca_df %>% 
group_by(Sex, Fac) %>% # Separate by `Sex` and `Fac` 
do(compute_hull(.)) %>% 
ungroup()

hulls <- hulls %>% 
mutate(Sex = factor(Sex), Fac = factor(Fac))

ggplot(pca_df, aes(x = PC1, y = PC2)) + 
geom_point(aes(shape = Species, color = Fac), size = 6) + 
geom_text(aes(label = Individual), size = 2, vjust = 2) + 
geom_polygon( 
data = hulls, 
aes(x = PC1, y = PC2, group = interaction(Sex, Fac), fill = Fac, color = Fac), 
alpha = 0.2, # Opacity for fill 
linewidth = 1 # Outline thickness for highlighting 
) + 
labs( 
title = "Shape's Principal Components (Separated by Sex)", 
x = paste("PC1 (", round(100 * PCA$sdev[1]^2 / sum(PCA$sdev^2), 1), "%)", sep = ""), 
y = paste("PC2 (", round(100 * PCA$sdev[2]^2 / sum(PCA$sdev^2), 1), "%)", sep = ""), 
shape = "Species", color = "Environment", fill = "Environment" 
) + 
theme_minimal() + 
scale_shape_manual(values ​​= uPCH, labels = levels(pca_df$Species), name = "Species") + 
scale_color_manual(values ​​= colors_fac, name = "Environment") + 
scale_fill_manual(values ​​= colors_fac, name = "Environment") + 
facet_wrap(~Sex) # Separate by sex

par(mfrow=c(1,2))

scorespc1<-PCA$x[,1] #scores for pc1
scorespc2<-PCA$x[,2] #scores for pc2
scorespc3<-PCA$x[,3] #scores for pc2

### View the shapes related to the extreme values ​​of each PC

predsPC1 <- shape.predictor(gpa$coords, scorespc1, Intercept = TRUE, predmin = min(scorespc1), predmax = max(scorespc1)) #estimating shape configurations based on pc scores
plotRefToTarget(mshape, predsPC1$predmin, mag=2, outline = Sapajusoutline$outline,gridPars = GP,  method = "points") #shape related to small pc scores
plotRefToTarget(mshape, predsPC1$predmax, mag=2, outline = Sapajusoutline$outline,gridPars = GP,  method = "points")#shape related to large pc scores

predsPC2 <- shape.predictor(gpa$coords, scorespc2, Intercept = TRUE, predmin = min(scorespc2), predmax = max(scorespc2)) #estimating shape configurations based on pc scores
plotRefToTarget(mshape, predsPC2$predmin, mag=2, outline = Sapajusoutline$outline,gridPars = GP,  method = "points") #shape related to small pc scores
plotRefToTarget(mshape, predsPC2$predmax, mag=2, outline = Sapajusoutline$outline,gridPars = GP,  method = "points")#shape related to large pc scores

predsPC3<- shape.predictor(gpa$coords, scorespc3, Intercept = TRUE, predmin = min(scorespc3), predmax = max(scorespc3)) #estimating shape configurations based on pc scores
plotRefToTarget(mshape, predsPC3$predmin, mag=2, outline = Sapajusoutline$outline,gridPars = GP,  method = "points") #shape related to small pc scores
plotRefToTarget(mshape, predsPC3$predmax, mag=2, outline = Sapajusoutline$outline,gridPars = GP,  method = "points")#shape related to large pc scores

par(mfrow=c(1,1))

dev.off()


####### PCA CORRELATION ##########

### Take the 3 PCs and correlate them to size ###

PC1 <- pca_df$PC1
PC2 <- pca_df$PC2
PC3 <- pca_df$PC3

correlation <- color(cbind(PC1, PC2, PC3, size))
print(correlation)

correlation <- color(cbind(PC1, PC2, PC3, biome))
print(correlation)

summary(PCA)
main_components <- PCA$x

bio_correlation <- color(cbind(main_components),as.numeric(biome))
print(bio_correlation)

correlation <- color(cbind(main_components),size)
print(correlation)


pca_data <- data.frame(pca_scores = main_components, Size = log(gpa$Csize))

size_correlations <- correlation
bio_cor <- correlation_bio

names(bio_cor) <- c("Comp1", "Comp2", "Comp3", "Comp4", "Comp5", "Comp6", 
"Comp7", "Comp8", "Comp9", "Comp10", "Comp11", "Comp12", 
"Comp13", "Comp14", "Comp15", "Comp16", "Comp17", "Comp18", 
"Comp19", "Comp20", "Comp21", "Comp22", "Comp23", "Comp24", 
"Comp25", "Comp26", "Comp27", "Comp28", "Comp29", "Comp30", 
"Comp31", "Comp32")
names(size_correlations) <- c("Comp1", "Comp2", "Comp3", "Comp4", "Comp5", "Comp6", 
"Comp7", "Comp8", "Comp9", "Comp10", "Comp11", "Comp12", 
"Comp13", "Comp14", "Comp15", "Comp16", "Comp17", "Comp18", 
"Comp19", "Comp20", "Comp21", "Comp22", "Comp23", "Comp24", 
"Comp25", "Comp26", "Comp27", "Comp28", "Comp29", "Comp30", 
"Comp31", "Comp32")


correlation_df <- data.frame(PC = names(bio_cor), Bio = bio_cor, Size = size_correlations)

correlation_df$PC <- factor(correlation_df$PC, 
levels = c("Comp1", "Comp2", "Comp3", "Comp4", "Comp5", 
"Comp6", "Comp7", "Comp8", "Comp9", "Comp10", 
"Comp11", "Comp12", "Comp13", "Comp14", 
"Comp15", "Comp16", "Comp17", "Comp18", 
"Comp19", "Comp20", "Comp21", "Comp22", 
"Comp23", "Comp24", "Comp25", "Comp26", 
"Comp27", "Comp28", "Comp29", "Comp30", 
"Comp31", "Comp32"))

ggplot(correlation_df, aes(x = PC, y = Bio, fill = Bio)) + 
geom_bar(stat = "identity") + 
scale_fill_gradientn(colors = c("red", "gray80", "steelblue"), 
values ​​= scales::rescale(c(-0.5, 0, 0.5)), 
name = "Color") + 
theme_classic() + 
labs(title = "Correlogram: Main Components vs. Biome", 
x = "Main Components", 
y = "Correlation") + 
theme(axis.text.x = element_text(angle = 45, hjust = 1))

correlation_df$PC <- factor(correlation_df$PC, 
levels = c("Comp1", "Comp2", "Comp3", "Comp4", "Comp5", 
"Comp6", "Comp7", "Comp8", "Comp9", "Comp10", 
"Comp11", "Comp12", "Comp13", "Comp14", 
"Comp15", "Comp16", "Comp17", "Comp18", 
"Comp19", "Comp20", "Comp21", "Comp22", 
"Comp23", "Comp24", "Comp25", "Comp26", 
"Comp27", "Comp28", "Comp29", "Comp30", 
"Comp31", "Comp32"))

ggplot(correlation_df, aes(x = PC, y = Size, fill = Size)) + 
geom_bar(stat = "identity") + 
scale_fill_gradientn(colors = c("red", "gray80", "steelblue"), 
values ​​= scales::rescale(c(-0.5, 0, 0.5)), 
name = "Color") + 
theme_classic() + 
labs(title = "Correlogram: Main Components vs. Size", 
x = "Main Components", 
y = "Correlation") + 
theme(axis.text.x = element_text(angle = 45, hjust = 1))


### this tells us if the correlation is likely

test_result <- color.test(PC1, size, method = "pearson")
print(test_result)

par(mfrow=c(3,1))

plot(pca_data$Size, PC1,
main = "Relationship between PC1 and logCS", xlab = "logCS", ylab = "PC1",
pch = 21, bg = colors_fac, col = colors_fac, cex = 2) +
abline(lm(PC1~Size, pca_data), col = "red")

test_result <- color.test(PC2, size, method = "pearson")
print(test_result)

plot(pca_data$Size, PC2,
main = "Relationship between PC2 and logCS", xlab = "logCS", ylab = "PC2",
pch = 21, bg = colors_fac, col = colors_fac, cex = 2) + 
abline(lm(PC2~Size,pca_data), col = "red")

test_result <- cor.test(PC3, size, method = "pearson")
print(test_result)

plot(pca_data$Size, PC3, 
main = "Relationship between PC3 and logCS", xlab = "logCS", ylab = "PC3", 
pch = 21, bg = colors_fac, col = colors_fac, cex = 2) + 
abline(lm(PC3~Size,pca_data), col = "blue")

#================================# CVA #================================####


# Load required packages
library(shapes)
library(geomorph)
library(car) # For dataEllipse

#================================================================
# 1. DATA PREPARATION
#================================================================
# This script assumes you have:
# - 'gpa': An object from geomorph::gpagen() containing Procrustes coordinates and Centroid Size.
# - 'Fac': A factor variable with group assignments (e.g., ecoregions) for each specimen.
# - 'colors_fac': A named vector of colors for your groups.
# - 'mshape': The mean shape, likely from gpa$consensus.
# - 'Sapajusoutline': An outline object for plotting shape changes.

# Calculate allometry-free shapes (residuals from a regression of shape on size)
allometry_model <- procD.lm(gpa$coords ~ log(gpa$Csize), iter = 999, RRPP = TRUE)
shape_residuals <- arrayspecs(allometry_model$residuals, p = dim(gpa$coords)[1], k = dim(gpa$coords)[2])
allometry_free_shape <- shape_residuals + array(gpa$consensus, dim(shape_residuals))

#================================================================
# 2. CANONICAL VARIATE ANALYSIS (CVA)
#================================================================

# Perform CVA on the allometry-free shape data
# Permutation test (10,000 rounds) assesses significance of group separation
cva_results <- CVA(allometry_free_shape, group = Fac, rounds = 10000, cv = TRUE)
print(cva_results)

dist_matrix <- cva_results$Dist$GroupdistMaha
pval_matrix <- cva_results$Dist$probsMaha

group_names <- levels(Fac)

results_table <- data.frame(
  Comparison = combn(group_names, 2, paste, collapse = " - "),
  Mahalanobis_Distance = as.vector(dist_matrix),
  p_value = as.vector(pval_matrix)
)

results_table$Mahalanobis_Distance <- round(results_table$Mahalanobis_Distance, 2)

cat("\n#================#================#===============#\n")
cat("--- PAIRWISE DISTANCES AND P-VALUES ---\n")
cat("#================#================#===============#\n")
print(results_table)

#================================================================
# CVA PLOT
#================================================================

# Set up plot parameters
par(mar = c(5, 5, 2, 2), pty = "s") # Set margins and square plot area

# Create the main plot of CVA scores
plot(cva_results$CVscores,
     asp = 1,
     pch = 21,
     bg = colors_fac,
     col = "black",
     cex = 2.5,
     xlab = paste0("Canonical Variate 1 (", round(cva_results$Var[1, 2], 1), "%)"),
     ylab = paste0("Canonical Variate 2 (", round(cva_results$Var[2, 2], 1), ")"),
     cex.lab = 1.2,
     cex.axis = 1.1
)

# Add 95% confidence ellipses for each group
for (i in 1:length(levels(Fac))) {
  dataEllipse(cva_results$CVscores[Fac == levels(Fac)[i], 1],
              cva_results$CVscores[Fac == levels(Fac)[i], 2],
              add = TRUE,
              levels = 0.95,
              col = colors_fac[levels(Fac)[i]],
              lwd = 2,
              plot.points = FALSE
  )
}

# Add a legend
legend("topleft",
       legend = levels(Fac),
       pch = 21,
       pt.bg = colors_fac[levels(Fac)],
       pt.cex = 2,
       cex = 1.1,
       bty = "n"
)

#================================================================
# 4. SUPPLEMENTARY VISUALIZATIONS
#================================================================

# Reset plot layout
par(mfrow = c(1, 1), mar = c(5, 4, 4, 2))

# --- Dendrogram of Mahalanobis Distances ---
maha_dist <- hclust(cva_results$Dist$GroupdistMaha)
maha_dist$labels <- levels(Fac)
plot(as.dendrogram(maha_dist),
     ylab = "Mahalanobis Distance",
     xlab = "Ecoregion Groups",
     main = "Morphological Distance Between Groups"
)

# --- Visualize Shape Changes Along CV Axes ---
# Get shape configurations at the extremes of the first two CV axes

cv1_min_shape <- min(cva_results$CVscores[, 1]) * matrix(cva_results$CVvis[, 1], nrow(mshape), ncol(mshape)) + mshape
cv1_max_shape <- max(cva_results$CVscores[, 1]) * matrix(cva_results$CVvis[, 1], nrow(mshape), ncol(mshape)) + mshape
cv2_min_shape <- min(cva_results$CVscores[, 2]) * matrix(cva_results$CVvis[, 2], nrow(mshape), ncol(mshape)) + mshape
cv2_max_shape <- max(cva_results$CVscores[, 2]) * matrix(cva_results$CVvis[, 2], nrow(mshape), ncol(mshape)) + mshape

# Plot the shape changes (requires 'geomorph' and a predefined outline)
# Example for CV1: shows the change from the negative to the positive end of the axis
par(mfrow = c(1, 2), mar = c(0, 0, 0, 0))

plotRefToTarget(cv1_min_shape, cv1_max_shape, 
                outline = Sapajusoutline$outline, method = "points", mag = 2, gridPars = GP)
plotRefToTarget(cv2_min_shape, cv2_max_shape, 
                outline = Sapajusoutline$outline, method = "points", mag = 2, gridPars = GP)

### CVA Corrigida ###
layout(matrix(c(4,5,0,1,1,2,1,1,3), 3), widths=c(1,1,1), heights=c(1,1,0.6))
par(mar=c(4, 4, 1, 1))

plot(cva_results$CVscores,
     asp = 1,
     pch = 21,
     bg = colors_fac,
     col = "black",
     cex = 2.5,
     xlab = paste0("Canonical Variate 1 (", round(cva_results$Var[1, 2], 1), "%)"),
     ylab = paste0("Canonical Variate 2 (", round(cva_results$Var[2, 2], 1), ")"),
     cex.lab = 1.2,
     cex.axis = 1.1
)

legend(+2.5,+3,legend=unique(Fac),pch=19,col=unique(colors_fac),title = "Ecoregions",
       cex = 1.2,
       pt.cex = 3,)

# add 95% ellipses
if (require(car)) {
  for(ii in 1:length(levels(Fac))){
    dataEllipse(cva.1$CVscores[Fac==levels(Fac)[ii],1],
                cva.1$CVscores[Fac==levels(Fac)[ii],2], 
                add=TRUE,levels=.95, col=unique(cores_CVA)[ii])}
}

#ordiellipse(cva$CVscores,group = fac,kind = "sd",conf = 0.95,col = legend_colors,lwd = 2)

text(
  cva$CVscores,
  labels = rownames(cva$CVscores),
  cex = 0.7,
  col = "black"
)

plotRefToTarget(mshape,cv1_min_shape,outline = Sapajusoutline$outline,gridPars = GP,  method = "points", mag = 2) 
plotRefToTarget(mshape,cv1_max_shape,outline = Sapajusoutline$outline,gridPars = GP,  method = "points", mag = 2)  
plotRefToTarget(mshape,cv2_min_shape,outline = Sapajusoutline$outline,gridPars = GP,  method = "points", mag = 2) 
plotRefToTarget(mshape,cv2_max_shape,outline = Sapajusoutline$outline,gridPars = GP,  method = "points", mag = 2) 
par(mfrow=c(1,1))


GP <- gridPar(
  pt.size = 0.8,
  pt.bg = "black",
  link.col = "darkgray",
  link.lwd = 1.5,
  out.col = "black",
  out.lwd = 2,
  tar.pt.bg = "black",
  tar.pt.size = 1,
  ref.pt.bg = "gray",
  ref.pt.size = 0.8
)
GP
# --- 2. Set up a 1x2 Plot Layout for side-by-side comparison ---
par(mfrow = c(1, 2), mar = c(1, 1, 2, 1))

# --- 3. Generate the Deformation Plots ---

# Plot 1: Savanna (Reference) vs. Amazon Forest (Target)
plotRefToTarget(
  cva_results$groupmeans[, , "SV"], # Reference shape
  cva_results$groupmeans[, , "AM"], # Target shape
  outline = Sapajusoutline$outline,
  gridPars = GP,
  method = "points",
  mag = 2
)
title("Savanna to Amazon Forest", cex.main = 1.5)

# Plot 2: Savanna (Reference) vs. Atlantic Forest (Target)
plotRefToTarget(
  cva_results$groupmeans[, , "SV"], # Reference shape
  cva_results$groupmeans[, , "AF"], # Target shape
  outline = Sapajusoutline$outline,
  gridPars = GP,
  method = "points",
  mag = 2
)
title("Savanna to Atlantic Forest", cex.main = 1.5)


# --- 4. Reset Plot Layout to Default ---
par(mfrow = c(1, 1))

### Shapiro - Wilk 
size <- gpa$Csize
head (size)

shapiro.test(size)
summary(size)

hist(size)
qqnorm(size) 
qqline(size)

dev.off()

# ANOVAs and ANCOVAs
#----------------------------------------------------------------

analysis_df <- data.frame(
  Species = species,
  Latitude = latitude,
  Biome = biome,
  Ecoregion = factor(factors$fac),
  Size = log(gpa$Csize),
  Sex = sex
)


# 2. UNIVARIATE ANALYSIS (CENTROID SIZE)
#----------------------------------------------------------------

# ANOVA: Size vs. Species
aov_species <- aov(Size ~ Species, data = analysis_df)
summary(aov_species)
TukeyHSD(aov_species)

# ANOVA: Size vs. Sex
aov_sex <- aov(Size ~ Sex, data = analysis_df)
summary(aov_sex)

# ANOVA: Size vs. Ecoregion
aov_ecoregion <- aov(Size ~ Ecoregion, data = analysis_df)
summary(aov_ecoregion)
TukeyHSD(aov_ecoregion)

# Visualization: Violin plot of Size by Ecoregion
# (Assuming 'colors_fac' is a predefined vector of colors)
ggplot(data = analysis_df, aes(x = Ecoregion, y = Size, fill = Ecoregion)) +
  geom_violin(trim = FALSE) +
  geom_jitter(width = 0.2, pch = 21, color = "black", bg = "cyan", size = 3) +
  geom_boxplot(width = 0.1, fill = "white", color = "black") +
  scale_fill_manual(
    values = colors_fac,
    labels = c("Amazon", "Atlantic Forest", "Savanna")
  ) +
  labs(x = "Ecoregion", y = "log(Centroid Size)", fill = NULL) +
  theme_minimal()


# 3. MULTIVARIATE ANALYSIS (SHAPE)
#----------------------------------------------------------------

# Create a geomorph data frame
gdf <- geomorph.data.frame(
  Shape = gpa$coords,
  Size = gpa$Csize,
  Species = species,
  Latitude = latitude,
  Biome = biome,
  Ecoregion = factor(factors$fac),
  Sex = sex
)

# Procrustes ANOVA (MANOVA) to test effects on shape
# Note: Using log(Size) to account for allometry
fit_allometry <- procD.lm(Shape ~ log(Size), data = gdf, iter = 999)
fit_biome <- procD.lm(Shape ~ log(Size) * Biome, data = gdf, iter = 999)
fit_species <- procD.lm(Shape ~ log(Size) * Species, data = gdf, iter = 999)
fit_ecoregion <- procD.lm(Shape ~ log(Size) * Ecoregion, data = gdf, iter = 999)
fit_sex <- procD.lm(Shape ~ log(Size) * Sex, data = gdf, iter = 999)
fit_interaction <- procD.lm(Shape ~ log(Size) * Ecoregion * Sex, data = gdf, iter = 999)

# Review model summaries
summary(fit_allometry)
summary(fit_biome)
summary(fit_species)
summary(fit_ecoregion)
summary(fit_sex)
summary(fit_interaction)

# Pairwise comparisons for significant factors (e.g., Species)
pairwise_results <- pairwise(fit_species, groups = gdf$Species)
summary(pairwise_results)

# Morphological distance tree (Phenogram) from pairwise distances
procrustes_distances <- summary(pairwise_results)$pairwise.tables$D
phenogram <- nj(procrustes_distances)

par(mfrow = c(1, 1))
plot(phenogram, main = "Morphological Phenogram (Neighbor-Joining)")

# 4. PCA-BASED MULTIVARIATE ANALYSIS (MANOVA WITH WILKS' LAMBDA)
#----------------------------------------------------------------
# This approach uses Principal Components as response variables.

# Assuming 'PCA' contains PCA results and 'PCA$x' are the scores.
# Using the first 16 PCs as response variables.
manova_sex <- manova(PCA$x[, 1:16] ~ Sex, data = analysis_df)
manova_ecoregion <- manova(PCA$x[, 1:16] ~ Ecoregion, data = analysis_df)
manova_species <- manova(PCA$x[, 1:16] ~ Species, data = analysis_df)

# Including size as a covariate
manova_sex_size <- manova(PCA$x[, 1:16] ~ Size + Sex, data = analysis_df)
manova_ecoregion_size <- manova(PCA$x[, 1:16] ~ Size + Ecoregion, data = analysis_df)

# Testing interactions
manova_species_interaction <- manova(PCA$x[, 1:16] ~ Size * Species, data = analysis_df)
manova_biome_interaction <- manova(PCA$x[, 1:16] ~ Size * Biome, data = analysis_df)

# Review MANOVA summaries with Wilks' test
summary(manova_sex, test = "Wilks")
summary(manova_ecoregion, test = "Wilks")
summary(manova_species, test = "Wilks")
summary(manova_sex_size, test = "Wilks")
summary(manova_ecoregion_size, test = "Wilks")
summary(manova_species_interaction, test = "Wilks")
summary(manova_biome_interaction, test = "Wilks")

# 5. ENVIRONMENTAL CORRELATION ANALYSIS (SIMPLE LINEAR MODELS)
#----------------------------------------------------------------

# Define environmental variables to test
env_variables <- c("lat", "long", "bio1", "bio2", "bio3", "bio4", "bio5", "bio6", "bio7", "bio8", "bio9", "bio10", "bio11", "bio12", "bio13", "bio14", "bio15", "bio16", "bio17", "bio18", "bio19", "npp", "humid", "soilmoist")

# Loop to run a linear model for each environmental variable
model_results <- list()
for (variable in env_variables) {
  formula <- as.formula(paste("CS ~", variable))
  model <- lm(formula, data = factors)
  model_results[[variable]] <- summary(model)
}

# Extract key results from models
p_values <- sapply(model_results, function(m) m$coefficients[2, "Pr(>|t|)"])
beta_coefficients <- sapply(model_results, function(m) m$coefficients[2, "Estimate"])
r_squared <- sapply(model_results, function(m) m$r.squared)

# Create a data frame for plotting
plot_data <- data.frame(
  Variable = names(beta_coefficients),
  Beta = beta_coefficients,
  P_Value = p_values,
  R_Squared = r_squared
) %>%
  mutate(Significant = P_Value < 0.05)

# Plot Beta coefficients
ggplot(plot_data, aes(x = Beta, y = reorder(Variable, Beta), color = Significant)) +
  geom_point(size = 4) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
  scale_color_manual(values = c("TRUE" = "deepskyblue3", "FALSE" = "rosybrown")) +
  labs(
    x = "Beta Coefficient",
    y = "Environmental Variable"
  ) +
  theme_minimal() +
  theme(legend.position = "none")

#=============== ALLOMETRY ####

# 6. ALLOMETRY ANALYSIS (SIZE-SHAPE RELATIONSHIP)
#----------------------------------------------------------------

# Model 1: Common allometry (Shape ~ Size)
fit_common_allometry <- procD.lm(Shape ~ log(Size), data = gdf, iter = 999)
summary(fit_common_allometry)

# Model 2: Allometry with interaction terms (testing for different slopes)
fit_interaction_allometry <- procD.lm(Shape ~ log(Size) * Sex * Ecoregion, data = gdf, iter = 999)
summary(fit_interaction_allometry)

# Visualize the common allometric pattern
plotAllometry(
  fit_common_allometry,
  size = gpa$Csize,
  logsz = TRUE,
  method = "PredLine",
  pch = as.numeric(pch),
  bg = colors_fac,
  cex = 3,
  col = colors_fac
)
legend(
  "bottomright",
  legend = levels(gdf$Ecoregion),
  pch = 22,
  pt.bg = unique(colors_fac),
  title = "Ecoregions",
  cex = 1.5
)


# 7. SIZE-SHAPE PCA & SHAPE PREDICTION
#----------------------------------------------------------------

# Perform a PCA on the size-shape space
size_shape_pca_plot <- plotAllometry(
  fit_common_allometry,
  size = gpa$Csize,
  logsz = TRUE,
  method = "size.shape",
  pch = as.numeric(pch),
  bg = colors_fac,
  cex = 5,
  col = colors_fac
)

# Extract PCA results and scores
size_shape_pca <- size_shape_pca_plot$size.shape.PCA
pc1_scores <- size_shape_pca$x[, 1]
pc2_scores <- size_shape_pca$x[, 2]

# Example test: Correlation between PC2 and sex
cor.test(pc2_scores, as.numeric(gdf$Sex))

# Predict shape changes along the primary axis of allometric variation (PC1)
min_score <- min(pc1_scores)
max_score <- max(pc1_scores)
shape_preds <- shape.predictor(
  gpa$coords,
  x = pc1_scores,
  Intercept = TRUE,
  predmin = min_score,
  predmax = max_score
)

# Visualize the predicted shapes at min and max PC1 scores
par(mfrow = c(1, 2))
plotRefToTarget(mshape(gpa$coords), shape_preds$predmin, mag = 2, outline = Sapajusoutline$outline, gridPars = GP, method = "points")
title("Shape at Minimum PC1 Score")
plotRefToTarget(mshape(gpa$coords), shape_preds$predmax, mag = 2, outline = Sapajusoutline$outline, gridPars = GP, method = "points")
title("Shape at Maximum PC1 Score")
par(mfrow = c(1, 1))


# 8. PAIRWISE COMPARISON OF ALLOMETRIC TRAJECTORIES
#----------------------------------------------------------------

# Fit models to compare between species
fit_species_interaction <- procD.lm(Shape ~ log(Size) * Species, data = gdf, print.progress = FALSE)
fit_no_interaction <- procD.lm(Shape ~ log(Size) + Species, data = gdf, print.progress = FALSE)

# Perform pairwise comparisons of slopes
pairwise_slopes <- pairwise(fit_species_interaction, fit.null = fit_no_interaction, groups = gdf$Species, print.progress = FALSE)

# Summarize pairwise results for different trajectory attributes
summary(pairwise_slopes, confidence = 0.95, test.type = "dist", stat.table = FALSE) # Distance between slopes
summary(pairwise_slopes, confidence = 0.95, test.type = "DL", stat.table = FALSE) # Difference in vector lengths
summary(pairwise_slopes, confidence = 0.95, test.type = "VC", angle.type = "deg", stat.table = FALSE) # Angle between vectors
summary(pairwise_slopes, confidence = 0.95, test.type = "var", stat.table = FALSE) # Difference in variance


# 9. PROCRUSTES DISTANCE HEATMAP
#----------------------------------------------------------------

# Assuming 'procrustes_distances' is a pre-computed distance matrix
species_names <- c(
  "Sapajus_apella", "Sapajus_cay", "Sapajus_libidinosus",
  "Sapajus_nigritus", "Sapajus_robustus", "Sapajus_xanthosternos"
)

procrustes_distances <- matrix(c(
  0.00000000, 0.02053338, 0.02600128, 0.03254954, 0.03394679, 0.04158019,
  0.02053338, 0.00000000, 0.02467312, 0.04457783, 0.03912553, 0.05200746,
  0.02600128, 0.02467312, 0.00000000, 0.03265863, 0.03809999, 0.04999526,
  0.03254954, 0.04457783, 0.03265863, 0.00000000, 0.04101024, 0.03830677,
  0.03394679, 0.03912553, 0.03809999, 0.04101024, 0.00000000, 0.03271396,
  0.04158019, 0.05200746, 0.04999526, 0.03830677, 0.03271396, 0.00000000
), nrow = 6, byrow = TRUE, dimnames = list(species_names, species_names))

# Melt the matrix for ggplot
heatmap_df <- melt(procrustes_distances)
colnames(heatmap_df) <- c("Species1", "Species2", "ProcrustesDist")

# Plot the heatmap
ggplot(heatmap_df, aes(x = Species1, y = Species2, fill = ProcrustesDist)) +
  geom_tile(color = "white") +
  scale_fill_viridis_c(option = "plasma", name = "Procrustes\nDistance") +
  coord_fixed() +
  labs(title = "Heatmap of Procrustes Distances Between Species") +
  theme_minimal(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.title = element_blank()
  )


# 10. MORPHOLOGICAL DISPARITY ANALYSIS
#----------------------------------------------------------------

# Calculate disparity among different groups
disparity_species <- morphol.disparity(Shape ~ 1, groups = ~ Species, data = gdf, iter = 999)
disparity_biome <- morphol.disparity(Shape ~ Size, groups = ~ Biome, data = gdf, iter = 999)
disparity_species_sex <- morphol.disparity(Shape ~ Size, groups = ~ Species + Sex, data = gdf, iter = 999)
disparity_ecoregion_sex <- morphol.disparity(Shape ~ Size, groups = ~ Ecoregion + Sex, data = gdf, iter = 999)

# Print results
print(disparity_species)
print(disparity_biome)
print(disparity_species_sex)
print(disparity_ecoregion_sex)

# ANALYSIS OF ALLOMETRY-FREE SHAPE (RESIDUALS)
#----------------------------------------------------------------

# Calculate allometry-free shapes by getting residuals
allometry_fit <- procD.lm(gpa$coords ~ log(gpa$Csize), iter = 999)
allometry_free_shape <- arrayspecs(allometry_fit$residuals, p = dim(gpa$coords)[1], k = dim(gpa$coords)[2])
allometry_free_shape <- allometry_free_shape + array(gpa$consensus, dim(allometry_free_shape)) # Add mean shape back

# Perform PCA on the allometry-free shapes
pca_residuals <- gm.prcomp(allometry_free_shape)
summary(pca_residuals)

# Create a data frame for plotting PCA results
df_pca_residuals <- data.frame(
  PC1 = pca_residuals$x[, "Comp1"],
  PC2 = pca_residuals$x[, "Comp2"],
  PC3 = pca_residuals$x[, "Comp3"],
  Species = gdf$Species,
  Sex = gdf$Sex,
  Ecoregion = gdf$Ecoregion
)

# Function to calculate convex hulls for plotting
compute_hull <- function(df) df[chull(df$PC1, df$PC2), ]
hulls <- df_pca_residuals %>%
  group_by(Ecoregion, Sex) %>%
  do(compute_hull(.)) %>%
  ungroup()

# Visualize PCA of residuals, faceted by sex
ggplot(df_pca_residuals, aes(x = PC1, y = PC2)) +
  geom_point(aes(shape = Species, color = Ecoregion), size = 5) +
  geom_polygon(data = hulls, aes(group = interaction(Sex, Ecoregion), fill = Ecoregion), alpha = 0.2) +
  labs(
    title = "PCA of Allometry-Free Shape",
    x = paste0("PC1 (", round(summary(pca_residuals)$importance[2, 1] * 100, 1), "%)"),
    y = paste0("PC2 (", round(summary(pca_residuals)$importance[2, 2] * 100, 1), "%)")
  ) +
  scale_shape_manual(values = uPCH, name = "Species") +
  scale_color_manual(values = colors_fac, name = "Ecoregion") +
  scale_fill_manual(values = colors_fac, name = "Ecoregion") +
  theme_minimal() +
  facet_wrap(~Sex)

# Test for shape differences in residuals
fit_residuals <- procD.lm(allometry_free_shape ~ Ecoregion * Sex, data = gdf)
summary(fit_residuals)

# Visualize shape changes along residual PCs
mean_residual_shape <- mshape(allometry_free_shape)
par(mfrow = c(2, 2))
# PC1
preds_PC1 <- shape.predictor(allometry_free_shape, x = df_pca_residuals$PC1, Intercept = TRUE)
plotRefToTarget(mean_residual_shape, preds_PC1$predmin, mag = 2, outline = Sapajusoutline$outline, gridPars = GP)
title("Min PC1")
plotRefToTarget(mean_residual_shape, preds_PC1$predmax, mag = 2, outline = Sapajusoutline$outline, gridPars = GP)
title("Max PC1")
# PC2
preds_PC2 <- shape.predictor(allometry_free_shape, x = df_pca_residuals$PC2, Intercept = TRUE)
plotRefToTarget(mean_residual_shape, preds_PC2$predmin, mag = 2, outline = Sapajusoutline$outline, gridPars = GP)
title("Min PC2")
plotRefToTarget(mean_residual_shape, preds_PC2$predmax, mag = 2, outline = Sapajusoutline$outline, gridPars = GP)
title("Max PC2")
par(mfrow = c(1, 1))


# SEX-SPECIFIC SHAPE ANALYSIS
#----------------------------------------------------------------

# Create a function to compare mean shapes of biomes to the overall mean
plot_biome_comparison <- function(gpa_data, factor_data, sex_label) {
  mean_shape_sex <- mshape(gpa_data$coords)
  coords_by_biome <- coords.subset(gpa_data$coords, factor_data$biome)
  
  par(mfrow = c(1, length(coords_by_biome)), mar = c(2, 2, 2, 2))
  for (biome_name in names(coords_by_biome)) {
    gpa_biome <- gpagen(coords_by_biome[[biome_name]], print.progress = FALSE)
    mean_shape_biome <- mshape(gpa_biome$coords)
    plotRefToTarget(mean_shape_sex, mean_shape_biome, mag = 3, outline = Sapajusoutline$outline, gridPars = GP)
    title(paste(sex_label, "-", biome_name))
  }
  par(mfrow = c(1, 1))
}

# Run comparisons for Males and Females
plot_biome_comparison(M_gpa, Males, "Males")
plot_biome_comparison(F_gpa, Females, "Females")

# Compare allometric models for each sex
fit_males <- procD.lm(M_gpa$coords ~ log(M_gpa$Csize) * as.factor(Males$fac), iter = 999)
fit_females <- procD.lm(F_gpa$coords ~ log(F_gpa$Csize) * as.factor(Females$fac), iter = 999)

summary(fit_males)
summary(fit_females)

# Pairwise comparison between sexes using the residual shape model
sex_comparison <- pairwise(fit_residuals, groups = gdf$Sex)
summary(sex_comparison)

# SHAPE-ENVIRONMENT COVARIATION (PARTIAL LEAST SQUARES)
#----------------------------------------------------------------

# Prepare environmental data
climatic_vars <- as.matrix(factors[, 11:31])
climatic_vars <- scale(climatic_vars, center = TRUE, scale = TRUE)

# Ensure row names match for analysis
rownames(climatic_vars) <- rownames(gdf$coords)

# --- PLS on Raw Shape Data ---

# Perform PLS between climatic variables and raw shape data
pls_raw_shape <- two.b.pls(A1 = climatic_vars, A2 = gpa$coords, iter = 9999)
summary(pls_raw_shape)

# --- PLS on Allometry-Free Shape Data ---

# Calculate allometry-free shapes by getting residuals from a PGLS
allometry_fit <- procD.lm(gpa$coords ~ log(gpa$Csize), iter = 999)
allometry_free_shape <- arrayspecs(allometry_fit$residuals, p = dim(gpa$coords)[1], k = dim(gpa$coords)[2])
allometry_free_shape <- allometry_free_shape + array(gpa$consensus, dim(allometry_free_shape)) # Add mean shape back

# Perform PLS between climatic variables and allometry-free shape
pls_allo_free <- two.b.pls(env = climatic_vars, Y = allometry_free_shape, iter = 9999)
summary(pls_allo_free)

# --- Visualization of PLS Results ---

# Helper function to create PLS score plots
create_pls_plot <- function(pls_result, title) {
  df_plot <- data.frame(
    PLS1 = pls_result$XScores[, 1],
    PLS2 = pls_result$YScores[, 1],
    Ecoregion = gdf$Ecoregion
  )
  
  # Define colors and shapes inside the function for encapsulation
  ecoregion_colors <- c("green3", "darkgreen", "goldenrod3")
  names(ecoregion_colors) <- levels(gdf$Ecoregion)
  ecoregion_shapes <- c(16, 15, 17) # Example shapes
  names(ecoregion_shapes) <- levels(gdf$Ecoregion)

  ggplot(df_plot, aes(x = PLS1, y = PLS2)) +
    geom_point(aes(color = Ecoregion, shape = Ecoregion), size = 6, alpha = 0.8) +
    scale_color_manual(values = ecoregion_colors) +
    scale_shape_manual(values = ecoregion_shapes) +
    labs(
      title = title,
      x = "PLS1 Scores (Climate)",
      y = "PLS1 Scores (Shape)"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
      legend.position = "bottom"
    )
}

# Generate plots for both analyses
print(create_pls_plot(pls_raw_shape, "PLS of Climate vs. Raw Shape"))
print(create_pls_plot(pls_allo_free, "PLS of Climate vs. Allometry-Free Shape"))

# Visualize shape deformations for the allometry-free PLS
mean_shape_resid <- mshape(allometry_free_shape)
pls_shape_scores <- pls_allo_free$YScores[, 1]
shape_preds <- shape.predictor(
  allometry_free_shape, 
  x = pls_shape_scores, 
  Intercept = TRUE
)

par(mfrow = c(1, 2))
plotRefToTarget(mean_shape_resid, shape_preds$predmin, mag = 2, outline = Sapajusoutline$outline, gridPars = GP)
title("Shape at Min PLS1 Score")
plotRefToTarget(mean_shape_resid, shape_preds$predmax, mag = 2, outline = Sapajusoutline$outline, gridPars = GP)
title("Shape at Max PLS1 Score")
par(mfrow = c(1, 1))

# Bar plot of climatic variable contributions to the first PLS axis
df_env_vectors <- data.frame(
  variable = rownames(pls_allo_free$left.pls.vectors),
  loading = pls_allo_free$left.pls.vectors[, 1]
)

ggplot(df_env_vectors, aes(x = reorder(variable, loading), y = loading, fill = loading > 0)) +
  geom_col() +
  coord_flip() +
  scale_fill_manual(values = c("steelblue", "firebrick"), guide = "none") +
  labs(
    title = "Climatic Loadings on PLS1 (Allometry-Free)",
    x = "Climatic Variable",
    y = "Loading"
  ) +
  theme_minimal()

# COMPARATIVE PHYLOGENETIC METHODS
#----------------------------------------------------------------

# --- Data and Tree Preparation ---

# Load phylogenetic tree
# tree <- read.tree("trees/AVGLOCSEX_Lima.tre")
# Or use the calibrated tree and prune tips
tree <- read.tree("trees/MCC_Lima_calibrated.tre")
tree <- drop.tip(tree, c("Sapajus_macrocephalus", "Sapajus_flavius", "Cebus_albifrons", "Cebus_capucinus"))
tree <- compute.brlen(tree, 1) # Standardize branch lengths

# Prepare shape data (using species means for species-level analysis)
shape_2d <- two.d.array(gpa$coords)
shape_means_2d <- rowsum(shape_2d, gdf$Species) / as.vector(table(gdf$Species))
shape_means <- arrayspecs(shape_means_2d, dim(gpa$coords)[1], dim(gpa$coords)[2])

# Prepare size data (using species means)
size_means <- rowsum(gpa$Csize, gdf$Species) / as.vector(table(gdf$Species))
size_means <- as.vector(size_means)
names(size_means) <- rownames(shape_means_2d)

# --- Phylogenetic Signal ---

# Test for phylogenetic signal in shape (using K-mult)
signal_shape <- physignal(shape_means, tree, iter = 999)
summary(signal_shape)

# Test for phylogenetic signal in size (using Pagel's Lambda and K)
signal_size_lambda <- phylosig(tree, size_means, method = "lambda", test = TRUE)
print(signal_size_lambda)

signal_size_k <- physignal(size_means, tree, iter = 999)
summary(signal_size_k)

# --- Phylogenetic ANOVA & PGLS ---

# Create a geomorph data frame with species means
# (Assuming 'temp.means', 'prec.means', 'bio.means' are pre-calculated)
gdf_means <- geomorph.data.frame(
  Shape = shape.means,
  Size = size_means,
  # Temp = temp.means,
  # Prec = prec.means,
  # Biome = bio.means,
  phy = tree
)

# Fit PGLS models to test hypotheses
# Example: Shape ~ Size * Biome
fit_pgls_biome <- procD.pgls(Shape ~ log(Size) * Biome, data = gdf_means, phy = tree, iter = 999)
summary(fit_pgls_biome)

# Example: Shape ~ Size * Precipitation
# fit_pgls_precip <- procD.pgls(Shape ~ log(Size) * Prec, data = gdf_means, phy = tree, iter = 999)
# summary(fit_pgls_precip)


# --- Visualization in Phylogenetic Context ---

# Perform phylogenetic PCA
phylo_pca <- gm.prcomp(shape_means, phy = tree)
summary(phylo_pca)

# Plot phylomorphospace
plot(phylo_pca, phylo = TRUE, main = "Phylomorphospace", pch = 21, bg = "lightblue", cex = 2)

# Visualize ancestral shape reconstructions
ancestral_shapes <- arrayspecs(phylo_pca$ancestors, dim(shape_means)[1], dim(shape_means)[2])
root_shape <- ancestral_shapes[, , 1] # Shape at the root

par(mfrow = c(1, 2))
plotRefToTarget(mshape(shape_means), root_shape, mag = 3, method = "points", outline = Sapajusoutline$outline, gridPars = GP)
title("Mean Shape vs. Root Ancestor")
# Compare root to a tip (e.g., Sapajus apella)
plotRefToTarget(root_shape, shape_means[, , "Sapajus_apella"], mag = 3, method = "points", outline = Sapajusoutline$outline, gridPars = GP)
title("Root vs. S. apella")
par(mfrow = c(1, 1))

# Map continuous trait (size) on the phylogeny
contMap(compute.brlen(tree, method = "Grafen"), size_means, fsize = 0.7)


# UTILITY SCRIPTS: DATA AVERAGING
#----------------------------------------------------------------

# --- Average by Species + Locality + Sex ---
average_by_group <- function(tps_file, factors_file, output_tps, output_env) {
  # This function can be expanded to perform the averaging logic
  # from the original script to keep the main analysis clean.
  # For now, it's a placeholder for that logic.
  print(paste("Averaging", tps_file, "based on", factors_file))
}

# Example usage:
# average_by_group("tps/jaw_18LM.tps", "global.csv", "tps/avglocsex.tps", "meanenv.txt")


# --- Average by Species ---
average_by_species <- function(tps_file, factors_file, output_tps, output_env) {
  # Placeholder for the averaging-by-species logic
  print(paste("Averaging", tps_file, "by species"))
}

# Example usage:
# average_by_species("tps/jaw_18LM.tps", "global.csv", "tps/avgspec.tps", "meanspec.txt")
