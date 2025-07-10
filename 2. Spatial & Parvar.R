# 1. SETUP
#----------------------------------------------------------------

# Load necessary libraries
library(geomorph)
library(ggplot2)
library(dplyr)
library(ape) 
library(phytools) 
library(spdep)
library(reshape2)
library(geosphere)
library(mpmcorrelogram)
library(vegan)
library(packfor)

# Load your data here (placeholders)
# gpa <- ...
# species <- ...
# sex <- ...
# biome <- ...
# factors <- ...
# Sapajusoutline <- ... # Outline data for plots
# GP <- ... # Grid parameters for plots
# locations <- read.csv("Planilhas/avglocsex.csv", sep = ";")
# Loc <- readland.tps("Datasets/tps/avglocsex.tps")
# tree <- read.tree("trees/AVGLOCSEX_Lima.tre")
# pls <- ... # Results from a previous PLS analysis

# Create the primary geomorph data frame
gdf <- geomorph.data.frame(
  Shape = gpa$coords,
  Size = gpa$Csize,
  Species = as.factor(species),
  Sex = as.factor(sex),
  Biome = as.factor(biome),
  Ecoregion = as.factor(factors$fac)
)


# SPATIAL & VARIATION PARTITIONING ANALYSIS
#----------------------------------------------------------------

# --- Data Preparation for Spatial Analysis ---

# Perform GPA on location-averaged data
L.gpa <- gpagen(Loc)

# Prepare geographic and trait distance matrices
coords <- as.matrix(locations[, c("long", "lat")])
geo_dist_m <- distm(coords, fun = distHaversine) / 1000 # Geographic distance in km
geo_dist <- as.dist(geo_dist_m)

logCS <- log(L.gpa$Csize)
size_dist <- dist(logCS) # Euclidean distance for size

shape_dist <- dist(two.d.array(L.gpa$coords)) # Procrustes distance for shape

# --- Spatial Autocorrelation (Mantel & Correlograms) ---

# Mantel-like test for size
mantel_size <- mantel.rtest(size_dist, geo_dist, nrepet = 999)
print(mantel_size)
plot(geo_dist, size_dist, pch = 16, col = "skyblue",
     xlab = "Geographic Distance (km)", ylab = "Size Distance (logCS)",
     main = "Geographic vs. Size Distance")
abline(lm(size_dist ~ geo_dist), col = "navy")

# Mantel-like test for shape
mantel_shape <- mantel.rtest(shape_dist, geo_dist, nrepet = 999)
print(mantel_shape)
plot(geo_dist, shape_dist, pch = 16, col = "salmon",
     xlab = "Geographic Distance (km)", ylab = "Shape Distance (Procrustes)",
     main = "Geographic vs. Shape Distance")
abline(lm(shape_dist ~ geo_dist), col = "firebrick")

# Spatial correlograms
mpmcorrelogram(dist(logCS), geo_dist, method = "pearson", alfa = 0.05, permutations = 999)
title(main = "Correlogram: Size vs. Geographic Distance")

mpmcorrelogram(shape_dist, geo_dist, method = "pearson", alfa = 0.05, permutations = 999)
title(main = "Correlogram: Shape vs. Geographic Distance")


# --- Variation Partitioning ---

# Prepare predictor matrices
# 1. Phylogeny
phylo_cov <- vcv.phylo(tree)
phylo_pca <- prcomp(phylo_cov)
phylo_vectors <- as.matrix(phylo_pca$x[, 1:3]) # Using first 3 phylo axes

# 2. Space (PCNM)
pcnm <- prcomp(geo_dist_m)
# Select significant spatial vectors
shape_pca_scores <- gm.prcomp(L.gpa$coords)$x
spatial_vectors_fwd <- forward.sel(shape_pca_scores, pcnm$x)
spatial_vectors <- as.matrix(pcnm$x[, spatial_vectors_fwd$order])

# 3. Size
size_vector <- as.matrix(log(L.gpa$Csize))

# 4. Environment (from PLS scores)
env_vectors <- as.matrix(pls$XScores)

# Run variation partitioning
var_partition <- varpart(
  shape_pca_scores,
  phylo_vectors,
  spatial_vectors,
  size_vector,
  env_vectors
)

# Plot the Venn diagram
plot(var_partition, cex = 1.2,
     Xnames = c("Phylogeny", "Space", "Size", "Environment"),
     bg = c("yellow", "deeppink", "navy", "green3"))

# Significance testing of unique contributions (conditional RDAs)
# Test unique effect of Phylogeny [a]
anova(rda(shape_pca_scores ~ phylo_vectors + Condition(spatial_vectors) + Condition(size_vector) + Condition(env_vectors)))

# Test unique effect of Space [b]
anova(rda(shape_pca_scores ~ spatial_vectors + Condition(phylo_vectors) + Condition(size_vector) + Condition(env_vectors)))

# Test unique effect of Size [c]
anova(rda(shape_pca_scores ~ size_vector + Condition(phylo_vectors) + Condition(spatial_vectors) + Condition(env_vectors)))

# Test unique effect of Environment [d]
anova(rda(shape_pca_scores ~ env_vectors + Condition(phylo_vectors) + Condition(spatial_vectors) + Condition(size_vector)))
