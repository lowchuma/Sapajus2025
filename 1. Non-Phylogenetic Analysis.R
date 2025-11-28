                        # ====  METADATA ==== 

#install_github("geomorphR/geomorph", ref = "Stable", build_vignettes = TRUE)

#install_github("mlcollyer/RRPP")
#install_github("fawda123/ggord")
#install.packages("geiger")
#install.packages("RColorBrewer")
#install.packages("vegan")
#install.packages("tidyverse")
#install.packages("rgl")

library(devtools)
library(tidyverse)
library(reshape2) 

# Statistics and Modeling

library(MASS) 
library(psych) 
library(klaR) 

# Data Visualization

library(ggplot2) # System for creating layer-based graphs
library(RColorBrewer) # Color palettes for visualization
library(ggord) # Creating ordination graphs with ggplot2
library(dotwhisker) # Coefficient graphs for statistical models
library(showtext)
library(sysfonts)
library(dplyr)
library(ggrepel)

# Ecology and Evolutionary Biology

library(vegan) # Ordination methods and diversity analysis
library(ape) # Phylogenetic and evolutionary analyses
library(phytools) # Tools for phylogenetic comparative biology
library(picante) # Integration of phylogenies and ecology
library(geiger) # Statistical methods for analyzing phylogenetic data
library(phangorn)

# Morphometrics

library(geomorph) # Geometric morphometric analysis
library(Morpho) # Tools for geometric morphometrics and mesh processing
library(shapes) # Routines for statistical analysis of landmark-based shapes
library(Rvcg) # 3D mesh processing and analysis
library(RRPP)
library(car)
library(reshape2)


# Spatial
library(cluster)
library(factoextra)
library(usdm) 
library(geosphere)
library(mpmcorrelogram)
library(spdep)
library(ade4)
library(adespatial)

              # ====  RAW DATA & COMMON PREDICTORS ====

tps <- readland.tps("tps/avglocsex.tps", specID =  "ID", readcurves = FALSE)
dim(tps)

factors <- read.csv("Plans/Factors.csv", sep=",")

names(factors)

species <- as.factor(factors$sp)
summary(species)

biome <- as.factor(factors$biome)
summary(biome)

sex <- as.factor(factors$sex)
summary (sex)

latitude <-as.numeric(factors$lat)

longitude <- as.numeric(factors$long)

environment <- as.factor(factors$fac)

# Run GPA - Aligns everything and takes the impact of raw data dimensionality

plot(tps)
gpa <- gpagen(tps)
link<-read.table("tps/link.txt") 
plot(gpa,link=link)

shp <- two.d.array(gpa$coords) 
cov_matrix <- cov(shp)  
print(cov_matrix)

## When necessary, find the mean specimen to draw the outline.

tps_0 <- readland.tps("tps/jaw_18LM.tps")
Y.gpa <- gpagen(tps_0)
findMeanSpec(Y.gpa$coords) # 90
global <- read.csv("Plans/global.csv", sep = ";")

plotOutliers(Y.gpa$coords, inspect.outliers = TRUE) ## verification


## Load the outline

drawinglandmark<-readland.tps("outline2/outline.tps")
outline<-read.table("outline2/outline.txt", header=FALSE)
summary(drawinglandmark)

mshape<-mshape(Y.gpa$coords)
Sapajusoutline<-warpRefOutline(file = "outline2/outline.txt", 
                               drawinglandmark[,,1],
                               mshape)    
#run dev.off() in case of error message 

## set outline configuration

grid.pars <- gridPar(grid.col = NULL) # if you want only the "invisible"


GP <- gridPar(n.col.cell = 100, pt.bg = "gray", pt.size = 0.8, tar.pt.bg = "cyan",
              tar.pt.size = 0.8, tar.out.col = "gray10", tar.out.cex = 0.5,
              grid.col = "white", grid.lwd = 0.5, txt.pos = 1, txt.col = "steelblue") # Custom grids

plotRefToTarget(mshape,Y.gpa$coords[,,90],outline = Sapajusoutline$outline,method="points", gridPars=GP)

dev.off()


            # ====  Parameters to standardize the plots ====

symbols <- c(
  Sapajus_apella = "\u25A0", # Square
  Sapajus_cay = "\u25BC", # Downward-pointing triangle
  Sapajus_libidinosus = "\u25B2", # Upward-pointing triangle
  Sapajus_nigritus = "\u25CF", # Circle
  Sapajus_robustus = "\u2666", # 
  Sapajus_xanthosternos = "\u2605" # 5-pointed star
)

symbols <- symbols[as.character(species)]

Fac <- as.factor(factors$fac) ## Environment factor

levels_fac <- levels(Fac)

colors_env <- c("AF" ="green3","AM" = "darkgreen", "SV" = "goldenrod1")
levels_env <- levels(as.factor(factors$env))
names(colors_env) <- levels_env


              # ====  Testing for confusing factors ====
              
               # if the interactions were not significant
       # we can move forward for the spatial autocorreletion tests

# Tests on Original data (all samples)

gdf_0 <- geomorph.data.frame(
  Shape = Y.gpa$coords,          # Procrustes coordinates
  Size = Y.gpa$Csize,            # Centroid Size
  Species = global$Species,        # Species factor
  Biome = global$Biome,       # Original Biome factor
  Sex = global$Sex            # Sex factor
)

print("Geomorph Data Frame created with factors")

# Run the Procrustes ANOVA for Shape
fit.shape <- procD.lm(Shape ~ Species * Sex,
                      data = gdf_0,
                      iter = 999,
                      RRPP = TRUE)
summary(fit.shape)

# Run the ANOVA for Size (SSD - Sexual Size Dimorphism)
fit.size <- procD.lm(Size ~ Species * Sex,
                     data = gdf_0,
                     iter = 999,
                     RRPP = TRUE)

summary(fit.size)


fit.full <- procD.lm(Shape ~ Size * Sex * Species,
                     data = gdf_0,
                     iter = 999,
                     RRPP = TRUE)

summary(fit.full)

print("Sex has no interaction, but explains a lot of the morphological variation")

                # ====  SPATIAL PARADIGM ==== 


library(raster)

# bio_files <- list.files("C:/Users/Lourenço/Downloads/wc2.1_2.5m_bio",
#                         pattern = "wc2.1_2.5m_bio", full.names = TRUE)
# elev_raster <- raster("C:/Users/Lourenço/Downloads/wc2.1_2.5m_elev/wc2.1_2.5m_elev.tif")
# npp_raster <- raster("C:/Users/Lourenço/Downloads/NPP/hdr.adf")
# humid_raster <- raster("C:/Users/Lourenço/Downloads/Humid/hdr.adf")
# soilph_raster <- raster("C:/Users/Lourenço/Downloads/SOILPH/hdr.adf")
# soilmoist_raster <- raster("C:/Users/Lourenço/Downloads/soilmoisture/hdr.adf")
# 
# # Stack all rasters
# bio_stack <- stack(bio_files)
# ref_raster <- bio_stack[[1]]  # use BIO1 as reference
# 
# elev_raster <- raster::resample(elev_raster, ref_raster, method = "bilinear")
# npp_raster <- raster::resample(npp_raster, ref_raster, method = "bilinear")
# humid_raster <- raster::resample(humid_raster, ref_raster, method = "bilinear")
# soilph_raster <- raster::resample(soilph_raster, ref_raster, method = "bilinear")
# soilmoist_raster <- raster::resample(soilmoist_raster, ref_raster, method = "bilinear")
# all_rasters <- stack(bio_stack, elev_raster, npp_raster, humid_raster, soilph_raster, soilmoist_raster)
# names(all_rasters) <- c(paste0("BIO", 1:19), "Elev", "NPP", "Humid", "Soil_P", "Soil_Moist")
# 
# library(sf)
# 
# coords_all <- data.frame(long = factors$long, lat = factors$lat)
# 
# # Assuming coords_all has columns "lon" and "lat"
# points_sf <- st_as_sf(coords_all, coords = c("long", "lat"), crs = 4326)  # WGS84
# 
# # Now transform to raster CRS
# points_sf <- st_transform(points_sf, crs = crs(all_rasters))
# 
# env_values <- raster::extract(all_rasters, points_sf)
# 
# # Combine with coordinates (optional)
# points_env <- bind_cols(coords_all, as_tibble(env_values))
# 
# env_values <- raster::extract(all_rasters, points_sf)
# 
# # Transform to tibble
# env_values_df <- as_tibble(env_values)
# 
# # Rename columns
# names(env_values_df) <- c(paste0("BIO", 1:length(bio_files)), "Elev", "NPP", "Humid", "Soil_P", "Soil_Moist")
# 
# # Combine with coordinates (optional)
# points_env <- bind_cols(coords_all, env_values_df)
# 
# # Remove old env_vars from factors
# old_vars <- c(paste0("bio",1:19), "npp", "soilmoist", "humid")
# factors <- factors %>% dplyr::select(-any_of(old_vars))
# 
# # Bind new environmental data
# factors <- bind_cols(factors, env_values_df)
# write.csv(factors, "Predictors.csv")

# === FULL ITERATIVE PIPELINE: Spatial MEMs -> Residual Clusters ===

# Purpose: produce clusters minimizing spatial autocorrelation

# ====  Robust MEM clustering script ====
utm_epsg        <- 32722
mem_thresh0     <- 0.1
mem_thresh_min  <- 0.01
max_iters       <- 6
K_candidates    <- 2:5
moran_I_cut     <- 0.15
top_MEM_fallback<- 5

library(dplyr)
library(sf)
library(spdep)
library(adespatial)
library(vegan)
library(cluster)
library(usdm)

# ====  1. Variable selection ====
env_vars <- c(paste0("BIO", 1:19), "Elev", "NPP", "Humid")
env_vars <- env_vars[env_vars %in% names(factors)]
env_mat <- dplyr::select(factors, dplyr::any_of(env_vars)) %>% 
  dplyr::select(where(is.numeric))

env_scaled <- scale(env_mat, center = TRUE)

# ====  2. Unique coordinates ====
coords_all <- data.frame(long = factors$long, lat = factors$lat)
coords_unique_df <- coords_all %>% distinct(long, lat) %>% mutate(loc_id = row_number())
map_loc <- coords_all %>% left_join(coords_unique_df, by = c("long","lat")) %>% pull(loc_id)

coords_unique_sf <- sf::st_as_sf(coords_unique_df, coords = c("long","lat"), crs = 4326, remove = FALSE)
coords_unique_utm <- st_transform(coords_unique_sf, crs = utm_epsg)
coords_unique_mat <- st_coordinates(coords_unique_utm)

# Preliminary spatial tests (unique localities)

library(vegan)
library(geosphere)
library(mpmcorrelogram)
library(dplyr)

# ==== 1. Coordinates ====
coords_unique <- coords_unique_mat  # n_loc x 2 matrix

# Detect if coordinates are lat/long or UTM
if(any(coords_unique[,1] < -180 | coords_unique[,1] > 180 |
       coords_unique[,2] < -90  | coords_unique[,2] > 90)){
  # Assume UTM in meters, Euclidean distances
  geo_dist <- dist(coords_unique) / 1000  # km
} else {
  # Latitude/Longitude in degrees, use geosphere
  geo_dist <- as.dist(distm(coords_unique, fun = distHaversine) / 1000)  # km
}

# ==== 2. Size distances (one representative per locality) ====
size_unique <- factors %>%
  mutate(loc_id = map_loc) %>%
  group_by(loc_id) %>%
  slice(1) %>% 
  ungroup() %>%
  pull(CS)
size_unique_log <- log(size_unique)

# ==== 3. Shape distances (Procrustes) ====
shape_array <- gpa$coords
n_lm <- dim(shape_array)[1]; n_dim <- dim(shape_array)[2]; n_ind <- dim(shape_array)[3]
shape_flat <- matrix(NA, nrow = n_ind, ncol = n_lm * n_dim)
for(i in 1:n_ind){ shape_flat[i, ] <- as.vector(t(shape_array[,,i])) }

shape_df <- as.data.frame(shape_flat) %>% mutate(loc_id = map_loc)
shape_unique_df <- shape_df %>%
  group_by(loc_id) %>%
  slice(1) %>%
  ungroup() %>%
  dplyr::select(-loc_id) %>%
  as.matrix()

# ====  4. Build Spatial Eigenvectors (MEMs) ====
cat("\n--- Building MEM spatial variables (para N=61) ====\n")
mem_unique <- dbmem(coords_unique_mat)
mem_full <- if(ncol(mem_unique) > 0) as.data.frame(mem_unique) else NULL
if(!is.null(mem_full)) colnames(mem_full) <- paste0("MEM", seq_len(ncol(mem_full)))

# ==== 5. Test Spatial Autocorrelation (RDA Method) ====
# Test if Shape is explained by spatial (MEM) vectors
rda_shape_space <- rda(shape_unique_df ~ ., data = mem_full)
anova_shape_space <- anova(rda_shape_space, permutations = 999)
cat("\n--- SAC Test (Shape vs MEMs) ====\n")
print(anova_shape_space) # Check p-value here
R2adj_shape <- RsquareAdj(rda_shape_space)$adj.r.squared
cat("Adjusted R² (Shape):", R2adj_shape, "\n")

# Test if Size is explained by spatial (MEM) vectors
rda_size_space <- rda(size_unique_log ~ ., data = mem_full)
anova_size_space <- anova(rda_size_space, permutations = 999)
cat("\n--- SAC Test (Size vs MEMs) ====\n")
print(anova_size_space) # Check p-value here
R2adj_size <- RsquareAdj(rda_size_space)$adj.r.squared
cat("Adjusted R² (Size):", R2adj_size, "\n")


# ==== 5. Correlograms (Visual Check) ====
cat("\n--- Generating Correlograms (Visual Check) ====\n")

# Calculate geographic distance matrix
if(any(coords_unique[,1] < -180 | coords_unique[,1] > 180 |
       coords_unique[,2] < -90  | coords_unique[,2] > 90)){
  geo_dist <- dist(coords_unique) / 1000  # Assumes UTM (km)
} else {
  geo_dist <- as.dist(distm(coords_unique, fun = distHaversine) / 1000) # Assumes Lat/Lon (km)
}

# Create morphological distance matrices
size_dist <- dist(size_unique_log)
shape_dist <- dist(shape_unique_df)

# Plot correlograms
mpmcorrelogram(size_dist, geo_dist, method = "pearson", alfa = 0.05, permutations = 999)
title(main = "Correlogram: Size vs Geographic Distance")

mpmcorrelogram(shape_dist, geo_dist, method = "pearson", alfa = 0.05, permutations = 999)
title(main = "Correlogram: Shape vs Geographic Distance")

# ==== Preliminary spatial tests for environmental variables ====

library(dplyr)
library(vegan)
library(geosphere)
library(mpmcorrelogram)

# Numeric environmental variables
env_mat <- factors %>% dplyr::select(any_of(env_vars)) %>% dplyr::select(where(is.numeric))

# Unique coordinates
coords_unique <- coords_unique_mat  # n_loc x 2

# Geographic distance (detects lat/long or UTM)
if(any(coords_unique[,1] < -180 | coords_unique[,1] > 180 |
       coords_unique[,2] < -90  | coords_unique[,2] > 90)){
  geo_dist <- dist(coords_unique) / 1000  # km, UTM
} else {
  geo_dist <- as.dist(distm(coords_unique, fun = distHaversine) / 1000)  # km
}

# Preliminary tests for each variable
for(varname in colnames(env_mat)){
  # Aggregate one observation per locality
  var_unique <- factors %>%
    mutate(loc_id = map_loc) %>%
    group_by(loc_id) %>%
    slice(1) %>%
    ungroup() %>%
    pull(varname)
  
  var_dist <- dist(var_unique)
  
  # Mantel test
  mantel_res <- mantel.rtest(var_dist, geo_dist, nrepet = 999)
  cat("Mantel test:", varname, "vs Geographic Distance -> r =",
      round(mantel_res$obs,3), "p =", mantel_res$pvalue, "\n")
  
  # Correlogram
  mpmcorrelogram(var_dist, geo_dist, method = "pearson", alfa = 0.05, permutations = 999)
  title(main = paste("Correlogram:", varname, "vs Geographic Distance"))
}

# ====  3. Connected KNN graph ====
knn_max <- 12
connected_k <- NA
for(k in 1:min(knn_max, nrow(coords_unique_mat)-1)){
  nb_try <- knn2nb(knearneigh(coords_unique_mat, k = k))
  if(n.comp.nb(nb_try)$nc == 1){ connected_k <- k; break }
}
if(is.na(connected_k)) connected_k <- min(4, nrow(coords_unique_mat)-1)
nb_unique <- knn2nb(knearneigh(coords_unique_mat, k = connected_k))
nb_unique <- make.sym.nb(nb_unique)
W_unique_bin <- nb2mat(nb_unique, style = "B", zero.policy = TRUE)
W_full_mat <- W_unique_bin[map_loc, map_loc]; diag(W_full_mat) <- 0
listw_full <- mat2listw(W_full_mat, style = "W", zero.policy = TRUE)

# ====  4. PCA on environment ====
pca_env <- rda(env_scaled[coords_unique_df$loc_id, , drop=FALSE])
pcs_env_sites <- vegan::scores(pca_env, display = "sites")

# ====  5. Build MEMs ====
mem_unique <- dbmem(coords_unique_mat)
mem_full <- if(ncol(mem_unique) > 0) as.data.frame(mem_unique) else NULL
if(!is.null(mem_full)) colnames(mem_full) <- paste0("MEM", seq_len(ncol(mem_full)))

# ====  6. Iterative MEM selection & clustering ====
iter <- 1
mem_thresh <- mem_thresh0
converged <- FALSE

while(iter <= max_iters && !converged){
  cat("\n--- Iteration", iter, "(mem_thresh =", mem_thresh, ") ====------------------\n")
  
  pcs_for_sel <- pcs_env_sites[, 1:min(3, ncol(pcs_env_sites)), drop = FALSE]
  selected_MEMs <- NULL
  
  # MEM selection
  if(!is.null(mem_full) && ncol(mem_full) > 0){
    fs_try <- try(forward.sel(Y = pcs_for_sel, X = mem_full, alpha = mem_thresh, nperm = 999), silent = TRUE)
    if(!inherits(fs_try, "try-error") && !is.null(fs_try) && nrow(fs_try) > 0){
      sel_idx <- fs_try$order
      selected_MEMs <- mem_full[, sel_idx, drop = FALSE]
      cat("Selected", ncol(selected_MEMs), "MEM(s) by forward.sel.\n")
    } else {
      mem_cor <- apply(mem_full, 2, function(x) abs(cor(x, pcs_for_sel[,1], use="pairwise.complete.obs")))
      top_idx <- order(mem_cor, decreasing = TRUE)[1:min(top_MEM_fallback, ncol(mem_full))]
      selected_MEMs <- mem_full[, top_idx, drop = FALSE]
      cat("Fallback: selected top", ncol(selected_MEMs), "MEMs by correlation.\n")
    }
  }
  
  # Residualize environment by MEMs only
  if(!is.null(selected_MEMs)){
    cov_df <- selected_MEMs
    env_resid <- sapply(seq_len(ncol(env_scaled)), function(j){
      resid(lm(env_scaled[coords_unique_df$loc_id,j] ~ ., data = cov_df))
    })
  } else {
    env_resid <- env_scaled[coords_unique_df$loc_id, , drop = FALSE]
  }
  env_resid <- scale(env_resid)
  
  # PCA on residuals
  pca_resid <- rda(env_resid)
  scores_resid_use <- vegan::scores(pca_resid, display = "sites")
  D_resid <- dist(scores_resid_use)
  
  # Test K candidates with detailed diagnostics
  sil_valid <- c(); valid_K <- c(); clusters_list <- list()
  
  for(K in K_candidates){
    cl <- cutree(hclust(D_resid, method = "ward.D2"), k = K)
    
    moran_obj <- moran.test(as.numeric(cl), nb2listw(nb_unique), zero.policy = TRUE)
    cat("K =", K, "-> Moran I =", round(moran_obj$estimate[1],3), "p =", round(moran_obj$p.value,3), "\n")
    if(abs(moran_obj$estimate[1]) > moran_I_cut | moran_obj$p.value < 0.05) next
    
    sil_val <- mean(silhouette(cl, D_resid)[, "sil_width"], na.rm = TRUE)
    cat("Silhouette for K =", K, ":", round(sil_val,3), "\n")
    
    sil_valid <- c(sil_valid, sil_val)
    valid_K <- c(valid_K, K)
    clusters_list[[as.character(K)]] <- cl
  }
  
  if(length(valid_K) == 0){
    cat("No valid K. Reducing MEM threshold.\n")
    mem_thresh <- max(mem_thresh / 2, mem_thresh_min)
    iter <- iter + 1
    next
  } else {
    K_best <- valid_K[which.max(sil_valid)]
    clusters_unique <- clusters_list[[as.character(K_best)]]
    cat("Selected K by silhouette:", K_best, "\n")
    
    # Moran final
    moran_obj <- moran.test(as.numeric(clusters_unique), nb2listw(nb_unique), zero.policy = TRUE)
    cat("Final Moran I (clusters) =", round(moran_obj$estimate[1],3), "p =", moran_obj$p.value, "\n")
    
    # Map back to all individuals
    factors$Cluster <- as.factor(clusters_unique[map_loc])
    
    converged <- TRUE
    cat("Converged at iteration", iter, ". Clusters assigned to factors$Cluster_final\n")
  }
  
  iter <- iter + 1
}

cat("\n=== Summary ===\n")
cat("Final number of MEMs used:", ifelse(!is.null(selected_MEMs), ncol(selected_MEMs), 0), "\n")
cat("Final number of clusters:", K_best, "\n")
cat("Cluster assignment preview:\n")

print(table(factors$Cluster, factors$sex))

# ====  7. Final spatial autocorrelation check ====
cat("\n=== Final Moran I check for all 83 individuals ===\n")

# Convert final cluster to numeric
clusters_numeric <- as.numeric(factors$Cluster)

# Moran I para os clusters mapeados aos 83 indivíduos
moran_final <- moran.test(clusters_numeric, listw_full, zero.policy = TRUE)

cat("Moran I statistic:", round(moran_final$estimate[1], 3), "\n")
cat("Expected I:", round(moran_final$estimate[2], 3), "\n")
cat("Variance:", round(moran_final$estimate[3], 3), "\n")
cat("p-value:", moran_final$p.value, "\n")

if(moran_final$p.value < 0.05) {
  cat("Warning: clusters still show significant spatial autocorrelation!\n")
} else {
  cat("Clusters are effectively spatially independent.\n")
}


# Map species to plotting symbols
species_factor_map <- factor(factors$sp)
pch_vector_map <- as.numeric(factors$pch)
species_pch_map <- tapply(pch_vector_map, species_factor_map, `[`, 1)

# Cluster summaries
cluster_env_means <- factors %>% group_by(Cluster) %>% summarise(across(all_of(env_vars), mean, na.rm = TRUE), .groups = "drop")
cluster_species_table <- table(factors$Cluster, factors$sp)
print(cluster_env_means, n = Inf, width = Inf)
print(cluster_species_table)
print(prop.table(cluster_species_table, margin = 1))

# ====  . Identify variables driving cluster separation ====
cat("\n=== Variables driving cluster separation ===\n")

env_selected <- env_scaled

# Map clusters for unique locations (61)
clusters_unique_factor <- as.factor(clusters_unique)

# Compute mean per cluster (centroid)
cluster_env_means <- t(sapply(levels(clusters_unique_factor), function(cl) {
  colMeans(env_selected[clusters_unique_factor == cl, , drop = FALSE])
}))
cluster_env_means <- as.data.frame(cluster_env_means)
cluster_env_means$Cluster <- levels(clusters_unique_factor)

# Compute ranges across clusters for each variable
names(cluster_env_means)
cluster_env_ranges <- apply(cluster_env_means[, 1:22], 2, function(x) max(x) - min(x))
cluster_env_ranking <- sort(cluster_env_ranges, decreasing = TRUE)
cat("Variables ranked by range across clusters:\n")
print(cluster_env_ranking)

# Optional: table of means per cluster for interpretation
cat("\nCluster centroids (means of VIF-selected variables):\n")
print(cluster_env_means[, c("Cluster", env_vars)])


# ====  Map cluster names back to all individuals (83) ====
# clusters_unique: vector of clusters for the 61 unique locations
# map_loc: indices mapping each individual to its unique location
clusters_all <- clusters_unique[map_loc]

cluster_names_short <- c(
  "1" = "Warm & Moderately Seasonal",
  "2" = "Thermally Stable & Dry-Seasonal",
  "3" = "Hot & Thermally Aseasonal",
  "4" = "Elevated & Thermally Seasonal"
)

factors$Cluster <- factor(clusters_all,
                          levels = names(cluster_names_short),
                          labels = cluster_names_short)

# Define intuitive colors for plotting (Climate-based)
cluster_colors <- c(
  "Hot & Thermally Aseasonal"        = "#FDE725FF", 
  "Warm & Moderately Seasonal"       = "#35B779FF",  
  "Thermally Stable & Dry-Seasonal"  = "#31688EFF",  
  "Elevated & Thermally Seasonal"    = "#440154FF"   
)
names(cluster_colors) <- levels(factors$Cluster)

# ====  PCA biplot for 61 unique locations ====
sites_scores_unique <- as.data.frame(vegan::scores(pca_env, display = "sites", choices = 1:2))
colnames(sites_scores_unique) <- c("PC1", "PC2")

# Map clusters for unique locations
sites_scores_unique$Cluster <- factor(clusters_unique,
                                      levels = names(cluster_names_short),
                                      labels = cluster_names_short)

# Variable scores for arrows
var_scores <- as.data.frame(vegan::scores(pca_env, display = "species", choices = 1:2))
colnames(var_scores) <- c("PC1", "PC2")
var_scores$Variable <- rownames(var_scores)

# Compute cluster centroids and hulls
centroids_unique <- sites_scores_unique %>%
  group_by(Cluster) %>%
  summarise(PC1 = mean(PC1), PC2 = mean(PC2), .groups = "drop")

find_hull <- function(df) df[chull(df$PC1, df$PC2), ]
hulls_unique <- sites_scores_unique %>%
  group_by(Cluster) %>%
  do(find_hull(.)) %>%
  ungroup()

# PCA summary for variance explained
pca_summary <- summary(pca_env)
prop_var <- pca_summary$cont$importance[2, 1:2] * 100

species_unique <- factors$sp[coords_unique_df$loc_id]

sites_scores_unique$Species <- as.factor(species_unique)
sites_scores_unique$species_symbol <- symbols[as.character(sites_scores_unique$Species)]

# Plot
pca_plot_unique <- ggplot() +
  geom_polygon(data = hulls_unique, aes(x = PC1, y = PC2, fill = Cluster, color = Cluster),
               alpha = 0.2, linewidth = 0.8) +
  geom_segment(data = var_scores, aes(x = 0, y = 0, xend = PC1, yend = PC2),
               arrow = arrow(length = unit(0.2, "cm")), color = "red", linewidth = 0.6) +
  geom_text_repel(data = var_scores, aes(x = PC1, y = PC2, label = Variable),
                  color = "red", size = 3, segment.alpha = 0.4, max.overlaps = Inf) +
  geom_text(data = sites_scores_unique,
            aes(x = PC1, y = PC2, label = species_symbol, color = Cluster),
            size = 8, vjust = 0.5) +
  geom_point(data = centroids_unique, aes(x = PC1, y = PC2, fill = Cluster),
             shape = 21, color = "black", size = 2, stroke = 1.2) +
  scale_fill_manual(name = "Climatic Cluster", values = cluster_colors) +
  scale_color_manual(values = cluster_colors, guide = "none") +
  labs(x = paste0("PC1 (", round(prop_var[1], 1), "%)"),
       y = paste0("PC2 (", round(prop_var[2], 1), "%)"),
       title = "PCA Biplot: Climatic Variables and Unique Locations") +
  coord_equal() +
  theme_bw(base_size = 14) +
  theme(panel.grid.minor = element_blank(),
        legend.position = "right",
        plot.title = element_text(hjust = 0.5))

print(pca_plot_unique)

# write.csv(factors,"Factors.csv")
          
        # ==== =========== MORPHOSPACE ================ 


                       # ====  PCA ====

# With Raw Shape Data

PCA <- gm.prcomp(gpa$coords)
summary(PCA)

PCs <- PCA$x
cor(PCs, log(gpa$Csize))
cor(PCs, as.numeric(sex))

cor <- matrix(NA,16,3, dimnames = list (paste("PC", 1:16), c("T-value", "r", "P-value")))

for(i in 1:16){
  
  r <- cor.test(PCs[,i], log(gpa$Csize))
  cor[i,1] = r$statistic
  cor[i,2] = r$estimate
  cor[i,3] = r$p.value
}

write.csv(cor, "corTest_Size.csv")

cor <- matrix(NA,16,3, dimnames = list (paste("PC", 1:16), c("T-value", "r", "P-value")))
for(i in 1:16){
  
  r <- cor.test(PCs[,i], as.numeric(sex))
  cor[i,1] = r$statistic
  cor[i,2] = r$estimate
  cor[i,3] = r$p.value
}
cor
write.csv(cor, "corTest_Sex.csv")

indv <- 1:nrow(PCA$x)

pca_df <- data.frame( 
  PC1 = PCA$x[, "Comp1"], 
  PC2 = PCA$x[, "Comp2"], 
  PC3 = PCA$x[, "Comp3"], 
  Species = species, 
  Individual = 1:nrow(PCA$x), Fac = factors$Cluster, Sex = sex)

species_colors <- c("red","navy","cyan3","saddlebrown","magenta4","black")
names(species_colors) <- levels(pca_df$Species)

pch <- as.numeric(factors$pch)
uPCH <- unique(pch)
names(uPCH) <- unique(pca_df$Species)

uPCH <- symbols

#============== PCA Plots ===========#

### Species PC1 vs. PC2

compute_hull <- function(pca_df) pca_df[chull(pca_df$PC1, pca_df$PC2), ]

hulls <- pca_df %>%
  group_by(Species) %>%
  do(compute_hull(.))

ggplot(pca_df, aes(x = PC1, y = PC2)) +
  geom_point(aes(color = Species, shape = Species), size = 8) +  
  geom_text(aes(label = Individual), size = 3, vjust = 1.5) +  
  geom_polygon(data = hulls, aes(group = Species, 
                                 fill = Species, 
                                 color = Species), alpha = 0.15) +  
  labs(title = "Shape's Principal Components", 
       x = paste("PC1 (", 
                 round(100 * PCA$sdev[1]^2 / sum(PCA$sdev^2), 1), "%)", 
                 sep = ""), 
       y = paste("PC2 (", 
                 round(100 * PCA$sdev[2]^2 / sum(PCA$sdev^2), 1), "%)",
                 sep = ""), color = "Species", shape = "Biomes") +
  theme_minimal() +
  scale_color_manual(values = species_colors, name = "Species") +  
  scale_fill_manual(values = species_colors, name = "Species") +  
  scale_shape_manual(values = uPCH, name = "Species")

### Bioclimatic Cluster PC1 vs. PC2


# Define colors for the clusters (adjust as needed)
cluster_colors <- cluster_colors

# Define shapes for species (optional, use if needed)
species_shapes <- symbols
species_shapes_mapped <- species_shapes[pca_df$Species]


# Function to calculate convex hulls
compute_hull <- function(df) df[chull(df$PC1, df$PC2), ]

# Calculate hulls for each climatic group
hulls <- pca_df %>%
  group_by(Fac) %>%
  do(compute_hull(.))

# Create the plot
pca_clim_plot <- ggplot(pca_df, aes(x = PC1, y = PC2)) +
  # Add convex hulls 
  geom_polygon(data = hulls, aes(fill = Fac, color = Fac), alpha = 0.2) +
  # Add points, colored by cluster, shaped by species (optional)
  geom_point(aes(color = Fac, shape = Species), size = 8, alpha = 0.8) + 
  # Or simpler: geom_point(aes(color = Fac), size = 3) +
  
  scale_color_manual(values = cluster_colors, name = "Climatic Group") +
  scale_fill_manual(values = cluster_colors, name = "Climatic Group") +
  scale_shape_manual(values = species_shapes, name = "Species") + # If using shapes
  
  labs(title = "Shape's Principal Components", 
       x = paste("PC1 (", 
                 round(100 * PCA$sdev[1]^2 / sum(PCA$sdev^2), 1), "%)", 
                 sep = ""), 
       y = paste("PC2 (", 
                 round(100 * PCA$sdev[2]^2 / sum(PCA$sdev^2), 1), "%)",
                 sep = ""), color = "Species", shape = "Biomes") +
  theme_bw(base_size = 14) +
  theme(plot.title = element_text(hjust = 0.5)) +
  coord_equal() # Keep aspect ratio 1:1 for morphospace

# Display the plot
print(pca_clim_plot)

### Clusters PC2 vs. PC3


# Define shapes for species (optional, use if needed)
species_shapes <- symbols
species_shapes_mapped <- species_shapes[pca_df$Species]


# Function to calculate convex hulls
compute_hull <- function(df) df[chull(df$PC2, df$PC3), ]

# Calculate hulls for each climatic group
hulls <- pca_df %>%
  group_by(Fac) %>%
  do(compute_hull(.))

# Create the plot
pca_clim_plot <- ggplot(pca_df, aes(x = PC2, y = PC3)) +
  # Add convex hulls 
  geom_polygon(data = hulls, aes(fill = Fac, color = Fac), alpha = 0.2) +
  # Add points, colored by cluster, shaped by species (optional)
  geom_point(aes(color = Fac, shape = Species), size = 8, alpha = 0.8) + 
  # Or simpler: geom_point(aes(color = Fac), size = 3) +
  
  scale_color_manual(values = cluster_colors, name = "Climatic Group") +
  scale_fill_manual(values = cluster_colors, name = "Climatic Group") +
  scale_shape_manual(values = species_shapes, name = "Species") + # If using shapes
  
  labs(
    title = "PCA Morphospace colored by Climatic Clusters",
    x = paste("PC2 (", round(summary(PCA)$importance[2,1]*100, 1), "%)", sep = ""),
    y = paste("PC3 (", round(summary(PCA)$importance[2,2]*100, 1), "%)", sep = "")
  ) +
  theme_bw(base_size = 10) +
  theme(plot.title = element_text(hjust = 0.5)) +
  coord_equal() # Keep aspect ratio 1:1 for morphospace

# Display the plot
print(pca_clim_plot)


### Sex PCA

levels_sex <- levels(sex)

colors_sex <- c("salmon","steelblue")
names(colors_sex) <- levels_sex

compute_hull <- function(pca_df) pca_df[chull(pca_df$PC1, pca_df$PC3), ]

hulls <- pca_df %>%
  group_by(Sex) %>%
  do(compute_hull(.))

ggplot(pca_df, aes(x = PC1, y = PC2)) +
  geom_point(aes(shape = Species, color = Sex), size = 10) +  
  geom_text(aes(label = Individual), size = 2, vjust = 2)  +  
  labs(
    title = "Shape's Principal Components", 
    x = paste("PC1 (", round(100 * PCA$sdev[1]^2 / sum(PCA$sdev^2), 1), "%)", sep = ""), 
    y = paste("PC2 (", round(100 * PCA$sdev[2]^2 / sum(PCA$sdev^2), 1), "%)", sep = ""), 
    shape = "Species", color = "Sex"
  ) +
  theme_minimal() + 
  scale_shape_manual(values = uPCH, labels = levels(pca_df$Species), name = "Species") + 
  scale_color_manual(values = colors_sex, name = "Sex") + # Define `colors_sex` with the colors corresponding to the sexes 
  scale_fill_manual(values = colors_sex, name = "Sex")


### Visualizing extreme PC scores shape deformations


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

plotRefToTarget(M, preds_comb$pred3, 
                main = "Min PC1 + Max PC2", mag = 2, 
              outline = Sapajusoutline$outline, gridPars = GP, method = "points")

plotRefToTarget(M, preds_comb$pred4, 
                main = "Max PC1 + Max PC2", mag = 2, 
              outline = Sapajusoutline$outline, gridPars = GP, method = "points")

plotRefToTarget(M, preds_comb$pred1, 
                main = "Min PC1 + Min PC2", mag = 2, 
              outline = Sapajusoutline$outline, gridPars = GP, method = "points")

plotRefToTarget(M, preds_comb$pred2, 
                main = "Max PC1 + Min PC2", mag = 2, 
              outline = Sapajusoutline$outline, gridPars = GP, method = "points")

### For PC3

#plotRefToTarget(M, preds_comb$pred6, main = "Min PC2 + Max PC3", mag = 2, outline = Sapajusoutline$outline, gridPars = GP, method = "points")

#plotRefToTarget(M, preds_comb$pred8, main = "Max PC2 + Max PC3", mag = 2, outline = Sapajusoutline$outline, gridPars = GP, method = "points")

#plotRefToTarget(M, preds_comb$pred5, main = "Min PC2 + Min PC3", mag = 2, outline = Sapajusoutline$outline, gridPars = GP, method = "points")

#plotRefToTarget(M, preds_comb$pred7, main = "Max PC2 + Min PC3", mag = 2, outline = Sapajusoutline$outline, gridPars = GP, method = "points")


par(mfrow=c(1,1))

dev.off()

# === === === === === === === === === === === === === === === === === === === #


                 # ====  Canonical Variate Analysis ====

cva_results <- CVA(gpa$coords, 
                   group = Fac, 
                   rounds = 10000, 
                   cv = TRUE)
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

cat("--- PAIRWISE DISTANCES AND P-VALUES ====\n")
print(results_table)

                    # ====  CVA PLOT ====

# Set up plot parameters
par(mar = c(5, 5, 2, 2), pty = "s") # Set margins and square plot area

# Create the main plot of CVA scores
plot(cva_results$CVscores,
     asp = 1,
     pch = 21,
     bg = colors_env[Fac],
     col = "black",
     cex = 2.5,
     xlab = paste0("Canonical Variate 1 (", round(cva_results$Var[1, 2], 1), "%)"),
     ylab = paste0("Canonical Variate 2 (", round(cva_results$Var[2, 2], 1), ")"),
     cex.lab = 1.2,
     cex.axis = 1.1
)

names(colors_env) <- levels(Fac)

# Add 95% confidence ellipses for each group
for (i in 1:length(levels(Fac))) {
  dataEllipse(cva_results$CVscores[Fac == levels(Fac)[i], 1],
              cva_results$CVscores[Fac == levels(Fac)[i], 2],
              add = TRUE,
              levels = 0.95,
              col = colors_env[levels(Fac)[i]],
              lwd = 2,
              plot.points = FALSE
  )
}

# Add a legend
legend("topleft",
       legend = levels(Fac),
       pch = 21,
       pt.bg = colors_env[levels(Fac)],
       pt.cex = 2,
       cex = 1.1,
       bty = "n"
)


           # ====  Calculate allometry-sex-free shapes ====
           # (residuals from a regression of shape on size) 

# mems <- dbmem(coords_all)
# space <- as.matrix(mems)
# shape <- two.d.array(gpa$coords)
# nrow(shape)
# shape_pcnm <- forward.sel(shape,space) 
# shape_pcnm$variables
# 
# space <- as.matrix(mems[,shape_pcnm$variables])

allometry_model <- procD.lm(gpa$coords ~ log(gpa$Csize) + sex,
                            iter = 999, RRPP = TRUE)

shape_residuals <- arrayspecs(allometry_model$residuals, 
                              p = dim(gpa$coords)[1], 
                              k = dim(gpa$coords)[2])

allometry_free_shape <- shape_residuals + array(gpa$consensus, 
                                            dim(shape_residuals))
AlloPCA <- gm.prcomp(allometry_free_shape)
summary(AlloPCA)
PCs <- AlloPCA$x 
cor(PCs, log(gpa$Csize))
cor(PCs, as.numeric(sex))

cor <- matrix(NA,16,3, dimnames = list (paste("PC", 1:16), c("T-value", "r", "P-value")))

for(i in 1:16){
  
  r <- cor.test(PCs[,i], log(gpa$Csize))
  cor[i,1] = r$statistic
  cor[i,2] = r$estimate
  cor[i,3] = r$p.value
}
cor
write.csv(cor, "corTest_AlloSize.csv")

cor <- matrix(NA,16,3, dimnames = list (paste("PC", 1:16), c("T-value", "r", "P-value")))
for(i in 1:16){
  
  r <- cor.test(PCs[,i], as.numeric(sex))
  cor[i,1] = r$statistic
  cor[i,2] = r$estimate
  cor[i,3] = r$p.value
}
cor
write.csv(cor, "corTest_AlloSex.csv")

# Perform CVA on the allometry-free shape data
# Permutation test (10,000 rounds) assesses significance of group separation

cva_results <- CVA(allometry_free_shape, 
                   group = Fac, 
                   rounds = 10000, 
                   cv = TRUE)
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

cat("--- PAIRWISE DISTANCES AND P-VALUES ====\n")
print(results_table)

# Set up plot parameters
par(mar = c(5, 5, 2, 2), pty = "s") # Set margins and square plot area

# Create the main plot of CVA scores
plot(cva_results$CVscores,
     asp = 1,
     pch = 21,
     bg = colors_env[Fac],
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
              col = colors_env[levels(Fac)[i]],
              lwd = 2,
              plot.points = FALSE
  )
}

# Add a legend
legend("topleft",
       legend = levels(Fac),
       pch = 21,
       pt.bg = colors_env[levels(Fac)],
       pt.cex = 2,
       cex = 1.1,
       bty = "n"
)


             # ====  CVA of Shape into Climatic Clusters ====

Env <- as.factor(factors$Cluster)

cva_results <- CVA(gpa$coords, 
                   group = Env, 
                   rounds = 10000, 
                   cv = TRUE)
print(cva_results)

dist_matrix <- cva_results$Dist$GroupdistMaha
pval_matrix <- cva_results$Dist$probsMaha

group_names <- levels(Env)

results_table <- data.frame(
  Comparison = combn(group_names, 2, paste, collapse = " - "),
  Mahalanobis_Distance = as.vector(dist_matrix),
  p_value = as.vector(pval_matrix)
)

results_table$Mahalanobis_Distance <- round(results_table$Mahalanobis_Distance, 2)

cat("--- PAIRWISE DISTANCES AND P-VALUES ====\n")
print(results_table)

                          # ====  CVA PLOT ====

cva_scores <- as.data.frame(cva_results$CVscores)
colnames(cva_scores) <- c("CV1", "CV2", "CV3")

cva_scores$Species <- factors$sp
cva_scores$Cluster <- factors$Cluster

species_shapes <- uPCH

names(cluster_colors) <- levels(cva_scores$Cluster)

ggplot(cva_scores, aes(x = CV1, y = CV2)) +
  
  # 95% Confidence Ellipses per Cluster
  stat_ellipse(aes(color = Cluster, group = Cluster),
               type = "norm", level = 0.95,
               linewidth = 1) +
  
  # Points → species shape + cluster color
  geom_point(aes(shape = Species, fill = Cluster, color = Cluster),
             size = 10, stroke = 0.8) +
  
  # Scales
  scale_shape_manual(values = species_shapes) +
  scale_color_manual(values = cluster_colors) +
  scale_fill_manual(values = cluster_colors) +
  
  # Labels with explained variance
  labs(
    title = "CVA Morphospace — Species Shape + Climatic Clusters",
    x = paste0("CV1 (", round(cva_results$Var[1,2], 1), "%)"),
    y = paste0("CV2 (", round(cva_results$Var[2,2], 1), "%)")
  ) +
  
  coord_equal() +
  theme_bw(base_size = 14) +
  theme(
    panel.grid.minor = element_blank(),
    legend.position = "right",
    plot.title = element_text(hjust = 0.5)
  )


#### With Non-Allometric Shape


cva_results <- CVA(allometry_free_shape, 
                   group = Env, 
                   rounds = 10000, 
                   cv = TRUE)
print(cva_results) # diminui sem a alometria e o sexo

dist_matrix <- cva_results$Dist$GroupdistMaha
pval_matrix <- cva_results$Dist$probsMaha

group_names <- levels(Env)

results_table <- data.frame(
  Comparison = combn(group_names, 2, paste, collapse = " - "),
  Mahalanobis_Distance = as.vector(dist_matrix),
  p_value = as.vector(pval_matrix)
)

results_table$Mahalanobis_Distance <- round(results_table$Mahalanobis_Distance, 2)

cat("--- PAIRWISE DISTANCES AND P-VALUES ====\n")
print(results_table)

library(viridis)

# ==== Step 1: Create an empty square matrix ====
group_names <- levels(Env)  # or your cluster names
n <- length(group_names)
dist_mat <- matrix(NA, nrow = n, ncol = n,
                   dimnames = list(group_names, group_names))

# ==== Step 2: Fill in the upper triangle with your distances ====
dist_mat[upper.tri(dist_mat)] <- as.vector(dist_matrix)

# ==== Step 3: Mirror to the lower triangle ====
dist_mat[lower.tri(dist_mat)] <- t(dist_mat)[lower.tri(dist_mat)]

# ==== Step 4: Diagonal = 0 ====
diag(dist_mat) <- 0

# ==== Step 5: Melt into long format for ggplot2 ====
dist_long <- melt(dist_mat, varnames = c("Cluster1", "Cluster2"), value.name = "Mahalanobis_Distance")

# ==== Step 6: Plot heatmap ====
ggplot(dist_long, aes(x = Cluster1, y = Cluster2, fill = Mahalanobis_Distance)) +
  geom_tile(color = "white") +
  geom_text(aes(label = round(Mahalanobis_Distance, 2)), color = "black", size = 4) +
  scale_fill_viridis(option = "magma", name = "Mahalanobis D") +
  theme_minimal(base_size = 14) +
  labs(title = "Pairwise Mahalanobis Distances Between Clusters",
       x = "Cluster", y = "Cluster") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

scores <- cva_results$CVscores[, 2] 
names(scores) <- rownames(cva_biome$CVscores)

idx_min <- which.min(scores)
idx_max <- which.max(scores)

shape_min <- allometry_free_shape[,, idx_min]
shape_max <- allometry_free_shape[,, idx_max]


par(mfrow = c(1,2))
# par(mar = c(0, 0, 0, 0))

plotRefToTarget(
  M1 = mshape(gpa$coords),
  M2 = shape_min,
  mag = 1.2,
  outline = Sapajusoutline$outline,
  method = "points",
  gridPars = GP
)

plotRefToTarget(
  M1 = mshape(gpa$coords),
  M2 = shape_max,
  mag = 1.2,
  outline = Sapajusoutline$outline,
  method = "points",
  gridPars = GP
)

par(mar = c(5,4,4,2)+0.1)


# ====  CVA PLOT ====

cva_scores <- as.data.frame(cva_results$CVscores)
colnames(cva_scores) <- c("CV1", "CV2", "CV3")

cva_scores$Species <- factors$sp
cva_scores$Cluster <- factors$Cluster

species_shapes <- uPCH

names(cluster_colors) <- levels(cva_scores$Cluster)

ggplot(cva_scores, aes(x = CV1, y = CV2)) +
  
  # 95% Confidence Ellipses per Cluster
  stat_ellipse(aes(color = Cluster, group = Cluster),
               type = "norm", level = 0.95,
               linewidth = 1) +
  
  # Points → species shape + cluster color
  geom_point(aes(shape = Species, fill = Cluster, color = Cluster),
             size = 10, stroke = 0.8) +
  
  # Scales
  scale_shape_manual(values = species_shapes) +
  scale_color_manual(values = cluster_colors) +
  scale_fill_manual(values = cluster_colors) +
  
  # Labels with explained variance
  labs(
    title = "CVA Morphospace — Species Shape + Climatic Clusters",
    x = paste0("CV1 (", round(cva_results$Var[1,2], 1), "%)"),
    y = paste0("CV2 (", round(cva_results$Var[2,2], 1), "%)")
  ) +
  
  coord_equal() +
  theme_bw(base_size = 14) +
  theme(
    panel.grid.minor = element_blank(),
    legend.position = "right",
    plot.title = element_text(hjust = 0.5)
  )


              # ====  LDA for two-way factor ====

# Define the two-level Biome factor
# Ensure this factor is correctly ordered/synchronized
biome_factor <- as.factor(factors$biome)

# Check levels (should show 'Forest' and 'Savanna')
cat("Levels of Biome factor:\n")
print(levels(biome_factor))
cat("Table of Biome factor:\n")
print(table(biome_factor))


# ==== 2. PERFORM LDA/CVA (2 Groups) ====
cat("\n--- Running CVA/LDA for Biome (Forest vs. Savanna) ====\n")

# Use the CVA function (it handles 2 groups as LDA)
# Using allometry-free shape
cva_biome <- CVA(allometry_free_shape, 
                 group = biome_factor, 
                 rounds = 10000, # Number of permutations for significance
                 cv = TRUE)     # Perform cross-validation

# ==== 3. PRINT RESULTS ====

# Print the main results summary
cat("\n--- CVA/LDA Summary ====\n")
print(cva_biome) 

# Print Mahalanobis distance and p-value between the two groups
cat("\n--- Pairwise Distance & Significance ====\n")
# Accessing might differ slightly depending on package version, check print(cva_biome) output
maha_dist_biome <- cva_biome$Dist$GroupdistMaha # Distance between group 1 and 2
p_val_biome <- cva_biome$Dist$probsMaha      # P-value for the distance
cat(paste("Mahalanobis Distance (Forest vs Savanna):", round(maha_dist_biome, 3), "\n"))
cat(paste("Permutation P-value:", p_val_biome, "\n"))


# Print Cross-Validation results
cat("\n--- Cross-Validation Results ====\n")
print(cva_biome$CV) 
cv_table_biome <- cva_biome$CVcv
overall_accuracy_biome <- sum(diag(cv_table_biome)) / sum(cv_table_biome) * 100
cat(paste("\nOverall Cross-Validated Accuracy:", round(overall_accuracy_biome, 2), "%\n"))
print("Cross-Validation Confusion Matrix (%):")
print(prop.table(cv_table_biome, margin = 1) * 100)


# ==== 4. (OPTIONAL) PLOT SCORES ====
# Since there are only 2 groups, there's only one Canonical Variate (CV1)
# A density plot or histogram is more informative than a scatter plot

cat("\n--- Plotting Scores (Density Plot) ====\n")
cva_scores_df <- data.frame(CV1 = cva_biome$CVscores[,1], Biome = biome_factor)

ggplot(cva_scores_df, aes(x = CV1, fill = Biome)) +
  geom_density(alpha = 0.6) +
  scale_fill_manual(values = c("Forest" = "darkgreen", "Savanna" = "goldenrod1")) + # Adjust colors if needed
  labs(title = "LDA Scores Distribution by Biome",
       x = "Canonical Variate 1 (LDA Axis)",
       y = "Density") +
  theme_bw()


#### Plot Mean Shapes of CVA

# Get the mean shapes for each cluster
mean_shapes <- cva_results$groupmeans 

cluster_names <- dimnames(mean_shapes)[[3]] 

# Generate all unique pairwise combinations
comparison_pairs <- combn(cluster_names, 2, simplify = FALSE)

# Create the list in the required format (Ref vs Tar)
comparisons <- lapply(comparison_pairs, function(pair) {
  c(Ref = pair[1], Tar = pair[2]) 
})

cat("Defined pairwise comparisons for the 4 clusters:\n")
print(comparisons)

outline_ref_points <- drawinglandmark[,,1] # Get coordinates from the TPS file

# Graphic parameters
GP_pairwise <- gridPar( 
  pt.bg = "grey60", 
  pt.size = 0.8, 
  tar.pt.bg = "black", 
  tar.pt.size = 0.8, 
  link.col = "red", 
  link.lwd = 1, 
  link.lty = 1, 
  out.col = "grey60"
)

# Layout
num_comparisons <- length(comparisons)
plot_rows <- ifelse(num_comparisons <= 3, 1, 2)
plot_cols <- ceiling(num_comparisons / plot_rows)
par(mfrow = c(plot_rows, plot_cols), mar = c(1, 1, 3, 1))

mag_factor <- 3

for (comp in comparisons) { 
  
  ref_group <- comp["Ref"] 
  tar_group <- comp["Tar"] 
  
  M1_coords <- mean_shapes[, , ref_group] # reference 
  M2_coords <- mean_shapes[, , tar_group] # target 
  
  plotRefToTarget( 
    M1 = M1_coords, 
    M2 = M2_coords, 
    outline = Sapajusoutline$outline, 
    gridPars = GP, 
    method = "points", 
    mag = mag_factor 
  ) 
  
  title(paste(ref_group, "\nvs\n", tar_group))
}

par(mfrow = c(1,1))
dev.off()

forest_idx <- which(biome_factor == "Forest")
savanna_idx <- which(biome_factor == "Savanna")

Forest_shape <- mshape(allometry_free_shape[,,forest_idx])
Savanna_shape <- mshape(allometry_free_shape[,,savanna_idx])

par(mfrow=c(1,1))

par(mar = c(0, 0, 0, 0))  

plotRefToTarget(
  M1 = Forest_shape,
  M2 = Savanna_shape,
  mag = 3,  # magnificação
  outline = Sapajusoutline$outline,
  method = "points",
  gridPars = GP
)

title("Forest → Savanna\nMean Shape Difference (magnified x3)")
par(mar = c(5, 4, 4, 2) + 0.1)

        # ====  END of Morphological Ordination Analysis ====


                  # ====  Linear Models ====

## Test for normal distribution

### Shapiro - Wilk 

size <- gpa$Csize
head (size)

shapiro.test(size)
summary(size)

size <- log(gpa$Csize)

shapiro.test(size)
summary(size)

# hist(size)
# qqnorm(size) 
# qqline(size)

# dev.off()

                      # ====  ANOVAs and ANCOVAs ====

size_mat <- as.matrix(size)
space <- as.matrix(mems)
size_pcnm <- forward.sel(size_mat,space) 
size_pcnm$variables

CSspace <- as.matrix(mems[,size_pcnm$variables])

analysis_df <- data.frame(
  Species = species,
  Latitude = latitude,
  Biome = biome, Ecorregion = factor(factors$fac),
  Climatic = factor(factors$Cluster),
  Size = log(gpa$Csize),
  Sex = sex,
  MEMs = CSspace
)


           # ====  UNIVARIATE ANALYSIS (CENTROID SIZE) ====

# ==== 1. Check for Normality (Shapiro-Wilk Test by Group) ====
# ANOVA assumes that the *residuals* are normal, but checking the
# normality of the data within each group is good practice.

cat("--- Normality Test (Shapiro-Wilk) by Groups ====\n")

names(factors)

groups <- factors[,c(4,7,10,11,37)]
dim (groups)

for (i in 1:5) { normality_groups <- analysis_df %>%
  group_by(groups[,i]) %>%
  summarise(
    W_statistic = shapiro.test(Size)$statistic,
    p_value = shapiro.test(Size)$p.value
  )

print(normality_groups)
  Sys.sleep(0.1) 
}


# Rule: If any p_value is < 0.05, the group is NOT normal.
# If most groups are not normal, prefer the Kruskal-Wallis test.

# ==== 2. Check for Homogeneity of Variances (Levene's Test) ====
# This is the most important premise for ANOVA.

cat("\n--- Homogeneity of Variances Test (Levene) ====\n")

# The 'center = median' makes the test more robust to non-normal data.

for (i in 1:5) {levene_test <- leveneTest(Size ~ groups[,i], 
                           data = analysis_df, 
                           center = median)
                  
print(levene_test)
Sys.sleep(1) 

}

# Rule:
# If p > 0.05: The variances ARE homogeneous. You CAN use ANOVA.
# If p < 0.05: The variances are NOT homogeneous. Use Kruskal-Wallis.

# ==== Decision: Which test to run? ====
#
# PATH A (Parametric):
# - If Levene's Test gave p > 0.05 AND the data are (reasonably) normal.
#
# PATH B (Non-Parametric):
# - If Levene's Test gave p < 0.05 OR the data are not normal.

# ==== PATH A: Parametric Test (ANOVA + Tukey) ====

cat("\n--- Parametric ANOVA and Tukey Post-hoc ====\n")

library(purrr)
library(dplyr)

# Define models in a named list
models <- list(
  "Sex" = lm(Size ~ Sex, data = analysis_df),
  "Sex + Ecorregion" = lm(Size ~ Sex + Ecorregion, data = analysis_df),
  "Sex + Climatic" = lm(Size ~ Sex + Climatic, data = analysis_df),
  "Sex + Biome" = lm(Size ~ Sex + Biome, data = analysis_df)
)
summary(models$`Sex + Climatic`)
summary(models$`Sex + Biome`)
summary(models$`Sex + Ecorregion`)

# Extract ANOVA results safely
results_df <- map_df(models, function(model) {
  an <- as.data.frame(anova(model))
  an$term <- rownames(an)
  rownames(an) <- NULL
  an
}, .id = "Model")
summary(lm(Size ~ Sex, data= analysis_df))

# Rename the columns consistently
colnames(results_df) <- c(
  "Model", "Df", "SumSq", "MeanSq", "Fvalue", "Pvalue", "Term"
)

# Optional: reorder for readability
results_df <- results_df %>%
  relocate(Model, Term, .before = Df)

# Save as CSV
write.csv(results_df, "ANOVA_results.csv", row.names = FALSE)
model <- lm(Size ~ Sex + Climatic, data = analysis_df)
summary(model)
TukeyHSD(model)

library(broom)

posthoc_df <- map_df(models[-1], function(model) {
  tk <- TukeyHSD(model)
  
  bind_rows(
    lapply(names(tk), function(f) {
      tidy(tk[[f]]) %>%
        mutate(Effect = f)
    })
  )
}, .id = "Model") %>%
  relocate(Model, Effect)

# Save posthoc results
write.csv(posthoc_df, "ANOVA_posthoc.csv", row.names = FALSE)

# Visualization: Violin plot of Size by Ecorregion

ggplot(data = analysis_df, aes(x = Ecorregion, y = Size, fill = Ecorregion)) +
  geom_violin(trim = FALSE) +
  geom_jitter(width = 0.25, pch = 21, color = "black", bg = "gray", size = 3) +
  geom_boxplot(width = 0.1, fill = "white", color = "black") +
  scale_fill_manual(
    values = colors_env,
    labels = c("Amazon", "Atlantic Forest", "Savanna")
  ) +
  labs(x = "Ecorregion", y = "log(Centroid Size)", fill = NULL) +
  theme_minimal()


library(ggnewscale)

# Create the plot
ggplot(data = analysis_df, aes(x = Climatic, y = Size)) +
  
  # Layer 1: Violin plot (filled by Climatic cluster)
  geom_violin(aes(fill = Climatic), trim = FALSE, alpha = 0.6) +
  scale_fill_manual(
    name = "Climatic Cluster",
    values = cluster_colors # Use cluster colors
  ) +
  
  # Introduce a new scale for fill before the points layer
  new_scale_fill() +
  
  # Layer 2: Jittered points (shape & fill by Species)
  geom_jitter(aes(shape = Species, fill = Species), # <-- fill mapped to Species
              width = 0.25,
              color = "black",     # Outline for shapes 21-25
              size = 4) +          # <-- Increased size
  
  # Layer 3: Box plot (overlaid)
  geom_boxplot(width = 0.1, fill = "white", color = "black", alpha = 0.5) +
  
  # Layer 4: Scales for points
  scale_shape_manual(
    name = "Species",
    values = uPCH # Use the species PCH codes
  ) +
  scale_fill_manual(
    name = "Species",         # Use "Species" as legend title
    values = species_colors # <-- Use species colors
  ) +
  
  # Layer 5: Labels and Theme
  labs(x = "Climatic Cluster", y = "log(Centroid Size)") +
  theme_minimal() +
  theme(legend.position = "right")


# ==== 6. Linear Models: log(Centroid Size) vs Environmental Variables ====
# (This section analyzes SIZE, not shape)
cat("\n--- 6. Linear Models: log(Size) vs Environmental Variables ====\n")
if (!"CS" %in% names(factors)) factors$CS <- log(gdf$Size)

all_climatic_variables <- c("BIO1", "BIO2", "BIO3", "BIO4", "BIO5", "BIO6",
                            "BIO7", "BIO8", "BIO9", "BIO10", "BIO11", "BIO12",
                            "BIO13", "BIO14", "BIO15", "BIO16", "BIO17",
                            "BIO18", "BIO19" , "NPP", "Humid", "Elev")
env_variables_to_test <- all_climatic_variables[all_climatic_variables %in% names(factors)]

factors_std_env <- cbind(analysis_df,scale(factors[, env_variables_to_test]))
names(factors_std_env)
factors_std_env$Size

model_results_size <- list()
conf_intervals_size <- list() # To store confidence intervals separately

for (variable in env_variables_to_test) {
  # Check if variable has zero variance after scaling
  if (sd(factors_std_env[[variable]], na.rm = TRUE) == 0) {
    cat(paste("Skipping", variable, "- zero variance after scaling.\n"))
    next # Skip to next variable
  }
  formula <- as.formula(paste("Size ~ Sex + MEMs.MEM1 + MEMs.MEM5 + MEMs.MEM13 +", variable))
  model <- lm(formula, data = factors_std_env)
  model_results_size[[variable]] <- summary(model)
  # Try to calculate confidence intervals, handle potential errors
  tryCatch({
    conf_intervals_size[[variable]] <- confint(model)[6, ] # Get CI for the predictor
  }, error = function(e) {
    cat(paste("Could not calculate CI for", variable, ":", e$message, "\n"))
    conf_intervals_size[[variable]] <- c(NA, NA) # Assign NA on error
  })
}

# Prepare data only for models that ran successfully and have CI
valid_vars <- names(model_results_size)
p_values_size <- sapply(model_results_size[valid_vars], function(m) coef(m)[6, "Pr(>|t|)"])
beta_coeffs_size <- sapply(model_results_size[valid_vars], function(m) coef(m)[6, "Estimate"])

r_coeffs_size <- sapply(model_results_size[valid_vars], function(m) (m)[9])

model_results_size$BIO1[9]

# Extract CIs safely
ci_matrix <- do.call(rbind, conf_intervals_size[valid_vars])
ci_low_size <- ci_matrix[, 1]
ci_high_size <- ci_matrix[, 2]


plot_data_size <- data.frame(
  term = valid_vars,
  estimate = beta_coeffs_size,
  conf.low = ci_low_size,
  conf.high = ci_high_size,
  p.value = p_values_size
) %>%
  filter(!is.na(conf.low)) %>% # Remove rows where CI calculation failed
  mutate(Significant = p.value < 0.05)

# Plot using dotwhisker (ensure package is loaded)
# library(dotwhisker)
dwplot_size <- dwplot(plot_data_size, dot_args = list(aes(color = Significant)), whisker_args = list(aes(color=Significant))) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
  scale_color_manual(values = c("TRUE" = "deepskyblue", "FALSE" = "grey60"), guide = "none") +
  labs(
    title = "Effect of Environmental Variables on log(Centroid Size)",
    x = "Standardized Beta Coefficient",
    y = "Environmental Variable"
  ) +
  theme_minimal()
print(dwplot_size)

# Convert the model results into a clean dataframe for export
lm_results_df <- data.frame(
  Predictor = valid_vars,
  Estimate = as.numeric(beta_coeffs_size),
  R2 = as.numeric(r_coeffs_size),
  CI_low = as.numeric(ci_low_size),
  CI_high = as.numeric(ci_high_size),
  P_value = as.numeric(p_values_size)
) %>%
  mutate(Significant = P_value < 0.05) %>%
  arrange(P_value)

write.csv(lm_results_df, "LinearModels_Size_Environmental.csv", row.names = FALSE)


             # ====  MULTIVARIATE ANALYSIS FOR SHAPE ==== 


gdf <- geomorph.data.frame(
  Shape = gpa$coords,
  NonAllo = allometry_free_shape,
  Size = gpa$Csize,
  Sex = factor(factors$sex),
  Species = factor(factors$sp),
  Biome = factor(factors$biome),
  Ecoregion = factor(factors$fac),
  Cluster = factor(factors$Cluster)
)

for(var in colnames(space)) {
  gdf[[var]] <- space[,var]
}

run_procD <- function(factor_name) {
  formula_text <- paste0("Shape ~ log(Size) + Sex + ", factor_name)
  fit <- procD.lm(as.formula(formula_text), data = gdf, iter = 999, RRPP = TRUE)
  
  cat("\n--- ANOVA for", factor_name, "---\n")
  print(summary(fit))
  
  # Pairwise post-hoc
  pval <- summary(fit)$table[factor_name, "Pr(>F)"]
  if (!is.na(pval) && pval < 0.05) {
    cat("\n--- Pairwise Comparisons for", factor_name, "---\n")
    posthoc <- pairwise(fit, groups = gdf[[factor_name]], covariate = NULL)
    print(summary(posthoc, test.type = "dist"))
  } else {
    cat("Factor", factor_name, "not significant (p =", round(pval, 4), ")\n")
  }
  return(fit)
}
fit_null <- run_procD("1")
fit_species    <- run_procD("Species")

fit_cluster    <- run_procD("Cluster")
fit_ecoregion  <- run_procD("Ecoregion")
fit_biome      <- run_procD("Biome")

morphol.disparity(fit_cluster, groups = gdf$Cluster)
morphol.disparity(fit_ecoregion, groups = gdf$Ecoregion)
morphol.disparity(fit_biome, groups = gdf$Biome)
morphol.disparity(fit_species, groups = gdf$Species)

species_levels <- unique(species)

species_means <- list()

for (sp in species_levels) {
  inds <- which(species == sp)           
  species_means[[sp]] <- mshape(allometry_free_shape[,,inds])  
}


species_matrix <- t(sapply(species_means, function(x) as.vector(x)))

distances_species <- dist(species_matrix)

# NJ tree
phenogram_species <- nj(as.matrix(distances_species))
plot(phenogram_species, main="Morphological Phenogram (Species - NJ)")

# Heatmap
dist_matrix_species <- as.matrix(distances_species)
heatmap_df_species <- melt(dist_matrix_species)
colnames(heatmap_df_species) <- c("Species1", "Species2", "ProcrustesDist")

ggplot(heatmap_df_species, aes(x=Species1, y=Species2, fill=ProcrustesDist)) +
  geom_tile(color="white") +
  scale_fill_viridis_c(option="plasma", name="Procrustes\nDistance") +
  coord_fixed() +
  theme_minimal(base_size = 12) +
  theme(axis.text.x = element_text(angle=45, hjust=1)) +
  labs(title="Heatmap of Procrustes Distances Between Species")

interaction.plot(
  x.factor = gdf$Cluster,   
  trace.factor = gdf$Sex,   
  response = gdf$Size,      
  fun = mean,
  col = c("blue","red"),    
  lty = 1,                  
  pch = 16,                 
  type = "b",               
  xlab = "Cluster",
  ylab = "Mean Centroid Size",
  main = "Interaction Plot: Cluster x Sex"
)

# ====  If wanted, save the results ====

library(geomorph)
library(dplyr)
library(openxlsx)

factors_to_run <- c("Cluster", "Ecoregion", "Biome", "Species")

# ==== Lists to store results ====
factor_anova <- list()
factor_posthoc <- list()
factor_disparity <- list()

for(fac in factors_to_run){
  
  # Run global model
  fit <- procD.lm(as.formula(paste0("Shape ~ log(Size) + Sex + MEM5 + MEM7 +", fac)),
                  
                  data = gdf, iter = 999, RRPP = TRUE)
  
  # ==== ANOVA: only the factor ====
  anova_tbl <- as.data.frame(summary(fit)$table)
  
  if(fac %in% rownames(anova_tbl)){
    fac_row <- anova_tbl[fac,, drop=FALSE] 
    fac_row$Term <- fac 
    fac_row$Model <- fac 
    rownames(fac_row) <- NULL 
    fac_row <- fac_row[, c("Model", "Term", setdiff(names(fac_row), c("Model","Term")))] 
    factor_anova[[fac]] <- fac_row 
  }
  
}

# ==== Save Factors ANOVA ====
all_factor_anova <- bind_rows(factor_anova)
write.csv(all_factor_anova, "procD_lm_models_Nsummary.csv", row.names=FALSE)

        # ====  SHAPE-ENVIRONMENT COVARIATION (PARTIAL LEAST SQUARES) ====

# Prepare environmental data
names(factors)
climatic_vars <- as.matrix(factors[, 13:34])
climatic_vars <- scale(climatic_vars, center = TRUE, scale = TRUE)

# Ensure row names match for analysis
nrow(climatic_vars)
rownames(climatic_vars) <- dimnames(gpa$coords)[[3]]

# ==== PLS on Raw Shape Data ====

# Perform PLS between climatic variables and raw shape data
pls_raw_shape <- two.b.pls(A1 = climatic_vars, A2 = gpa$coords, iter = 999)
summary(pls_raw_shape)

# ==== PLS on Allometry-Free Shape Data ====

# Calculate allometry-free shapes by getting residuals from a PGLS
allometry_fit <- procD.lm(gpa$coords ~ log(gpa$Csize) + sex, iter = 999)
allometry_free_shape <- arrayspecs(allometry_fit$residuals, p = dim(gpa$coords)[1], k = dim(gpa$coords)[2])
allometry_free_shape <- allometry_free_shape + array(gpa$consensus, dim(allometry_free_shape)) # Add mean shape back

# Perform PLS between climatic variables and allometry-free shape
pls_allo_free <- two.b.pls(A1 = climatic_vars, 
                           A2 = allometry_free_shape, iter = 999)
summary(pls_allo_free)

# ====  PLS plots ====

scores_df <- as.data.frame(pls_allo_free$YScores[, 1:2])
colnames(scores_df) <- c("PLS1", "PLS2")

eig <- pls_allo_free$svd$d^2
var1 <- round((eig[1] / sum(eig)) * 100, 1)
var2 <- round((eig[2] / sum(eig)) * 100, 1)

scores_df$Cluster <- factor(gdf$Cluster)
scores_df$Species <- factor(gdf$Species)

pch_species_vector <- as.numeric(factors$pch)
scores_df$PCH <- pch_species_vector

species_names_factor <- factor(gdf$Species)
unique_species_names <- levels(species_names_factor)
species_pch_map <- tapply(pch_species_vector, species_names_factor, `[`, 1)
names(species_pch_map) <- unique_species_names

cluster_plot_colors <- cluster_colors

hull_data <- scores_df %>%
  group_by(Cluster) %>%
  filter(n() > 2) %>%
  slice(chull(PLS1, PLS2)) %>%
  ungroup()

plot_scores_gg <- ggplot(scores_df, aes(x = PLS1, y = PLS2)) +
  geom_polygon(
    data = hull_data,
    aes(fill = Cluster, color = Cluster),
    alpha = 0.20,
    linewidth = 1
  ) +
  geom_point(
    aes(shape = Species, fill = Cluster, color = Cluster),
    size = 4,
    stroke = 0.8
  ) +
  scale_shape_manual(name = "Species", values = species_pch_map) +
  scale_fill_manual(name = "Climate Cluster", values = cluster_plot_colors) +
  scale_color_manual(name = "Climate Cluster", values = cluster_plot_colors) +
  labs(
    x = paste0("PLS1 (", var1, "%)"),
    y = paste0("PLS2 (", var2, "%)")
  ) +
  theme_bw(base_size = 14) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(linetype = "dashed", color = "grey90"),
    legend.position = "right",
    plot.title = element_text(hjust = 0.5)
  )

print(plot_scores_gg)

# ====  PLS Loadings ====
loadings_df <- as.data.frame(pls_allo_free$left.pls.vectors[, 1:2])
colnames(loadings_df) <- c("PLS1", "PLS2")
rownames(loadings_df) <- colnames(climatic_vars)

loadings_plot_df <- data.frame(
  Variable = rownames(loadings_df),
  Loading = loadings_df$PLS1
)
loadings_plot_df <- loadings_plot_df[order(loadings_plot_df$Loading, decreasing = TRUE), ]
loadings_plot_df$Color <- ifelse(loadings_plot_df$Loading > 0, "steelblue", "firebrick")

# Plot
par(mar = c(12, 4, 2, 1))
barplot(
  loadings_plot_df$Loading,
  names.arg = loadings_plot_df$Variable,
  col = loadings_plot_df$Color,
  border = NA,
  las = 2,
  ylab = "Loading PLS1"
)
abline(h = 0, col = "black", lwd = 1)

# ====  PLS Loadings dataframe ====
loadings_df <- as.data.frame(pls_allo_free$left.pls.vectors[, 1:2])
colnames(loadings_df) <- c("PLS1", "PLS2")
loadings_df$Variable <- colnames(climatic_vars)

score_range <- apply(scores_df[, c("PLS1", "PLS2")], 2, function(x) max(x) - min(x))
loading_range <- apply(loadings_df[, c("PLS1", "PLS2")], 2, function(x) max(x) - min(x))
scale_factor <- min(score_range / loading_range) * 1  # 0.5 deixa um pouco de folga

loadings_df$PLS1_scaled <- loadings_df$PLS1 * scale_factor
loadings_df$PLS2_scaled <- loadings_df$PLS2 * scale_factor

# ====  Plot Scores + Loadings  ====
plot_scores_loadings <- ggplot(scores_df, aes(x = PLS1, y = PLS2)) +
  geom_polygon(
    data = hull_data,
    aes(fill = Cluster, color = Cluster),
    alpha = 0.20,
    linewidth = 1
  ) +
  geom_point(
    aes(shape = Species, fill = Cluster, color = Cluster),
    size = 4,
    stroke = 0.8
  ) +
  geom_segment(
    data = loadings_df,
    aes(x = 0, y = 0, xend = PLS1_scaled, yend = PLS2_scaled),
    arrow = arrow(length = unit(0.3, "cm")),
    color = "red",
    linewidth = 0.8
  ) +
  geom_text(
    data = loadings_df,
    aes(x = PLS1_scaled, y = PLS2_scaled, label = Variable),
    hjust = 0.5, vjust = -0.5,
    size = 3
  ) +
  scale_shape_manual(name = "Species", values = species_pch_map) +
  scale_fill_manual(name = "Climate Cluster", values = cluster_plot_colors) +
  scale_color_manual(name = "Climate Cluster", values = cluster_plot_colors) +
  labs(
    x = paste0("PLS1 (", var1, "%)"),
    y = paste0("PLS2 (", var2, "%)")
  ) +
  theme_bw(base_size = 14) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(linetype = "dashed", color = "grey90"),
    legend.position = "right",
    plot.title = element_text(hjust = 0.5)
  )

print(plot_scores_loadings)

#============= Average specimens by locality & sex ==============#
  
# tps_0 <- readland.tps("tps/jaw_18LM.TPS", specID="imageID")
# Y.gpa <- gpagen(tps_0)
# factors <- read.csv("Plans/global.csv", sep = ";")
# 
# data <- global
# sapajus <- tps_0
# 
# avgterm <- as.factor(paste(data$Species, data$Lat, data$LongZ))
# x <- two.d.array(sapajus)
# means <- rowsum(x, avgterm) / as.vector(table(avgterm))
# 
# Y <- arrayspecs(means, dim(sapajus)[1], dim(sapajus)[2], sep = NULL)
# 
# writeland.tps(Y, file = "tps/avglocsex_2.tps")
# 
# names(data)
# meanenv <- rowsum(data[, c("Lat", "Long")], avgterm) / as.vector(table(avgterm))
#   
# write.table(meanenv, file = "meanenv.txt")
