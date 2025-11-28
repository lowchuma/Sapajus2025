                # ==== Phylogenetic Paradigm ====

    # ==== Creating Polytomous Tree to Increase Sample Size ====

library(ggplot2)
library(geomorph)
library(phytools)
library(vegan)
library(geiger)
library(phangorn)
library(tidyverse)
library(RColorBrewer)
library(usdm)
library(spdep)
library(adespatial)

# === INPUTS ===
tree <- read.tree("trees/Cebinae_Lima_Ultrametric.tre")
data <- read.csv("Plans/Factors.csv", sep = ",")
tps <- readland.tps("tps/avglocsex.tps")
Y.gpa <- gpagen(tps)

# Check tree labels
tree$tip.label
tree$edge.length

# Define outgroups exactly as in the tree
outgroups <- c("Cebus_albifrons", "Cebus_capucinus", "Sapajus_flavius", "Sapajus_macrocephalus")
ingroups <- c(paste(unique(as.character(data$sp))))

# Remove outgroups (only if they actually exist in the tree)
tree_work <- drop.tip(tree, intersect(tree$tip.label, outgroups))
tree_work <- keep.tip(tree, intersect(tree$tip.label, ingroups))

# Check remaining tips in the tree
tree_work$tip.label
length(tree_work$tip)  # number of remaining species
tree_work$edge.length

# Ensure the species column is named "sp"
data$sp <- as.character(data$sp)

# Create unique IDs per individual (current data order == GPA order)
names(data)
data$ID <- NA
data <- data %>%
  group_by(sp) %>%
  mutate(ID = paste0(sp, "_", row_number())) %>%
  ungroup()
head(data[,c(1:3, 37)])

# Synchronize GPA labels with data IDs
dimnames(Y.gpa$coords)[[3]] <- data$ID

# Count number of individuals per species
count_ind <- data %>%
  count(sp, name = "n_ind")

count_ind

# Expand tree to represent individuals (polytomies)
branch_length <- 0.1
tree_expanded <- tree_work

for(i in seq_len(nrow(count_ind))) {
  
  sp <- count_ind$sp[i]
  n <- count_ind$n_ind[i]
  
  if(sp %in% tree_expanded$tip.label && n > 1) {
    
    new_tips <- paste0(sp, "_", seq_len(n))
    
    # Locate original species tip
    node_target <- which(tree_expanded$tip.label == sp)
    
    # Add new individuals as a polytomy
    for(j in seq_len(n)) {
      tree_expanded <- bind.tip(
        tree_expanded,
        tip.label = new_tips[j],
        where = node_target,
        position = branch_length
      )
    }
    
    # Remove the original single-species tip
    tree_expanded <- drop.tip(tree_expanded, sp)
  }
}

# Inspect updated tip labels
tree_expanded$tip.label
length(tree_expanded$tip.label)

# Reorder data according to tree
data_sync <- data[match(tree_expanded$tip.label, data$ID), ]
# write.csv(data_sync, "Factors_Sync.csv")

# Reorder GPA coordinates according to the reordered data
coords_sync <- Y.gpa$coords[,, match(tree_expanded$tip.label, data$ID)]
dimnames(coords_sync)[[3]] <- tree_expanded$tip.label

size_sync <- Y.gpa$Csize[match(tree_expanded$tip.label, data$ID)]
names(size_sync) <- tree_expanded$tip.label

# Plotting options
# plot(tree_expanded, cex = 0.6)
# title("Expanded tree with all individuals")
# 
# plot(tree_expanded, cex = 0.7, no.margin = TRUE)
# axisPhylo()
# 

plot(tree_expanded, cex = 0.5, no.margin = TRUE, label.offset = 0.02)

plot(tree_expanded, type = "fan", cex = 0.5)

# write.tree(tree_expanded, "trees/Polytomic_Lima.tre")

# ==== Select best Phylo ====

library(phytools)
library(ape)
library(phangorn)
library(ape)
library(RRPP)

tree1 <- read.tree("trees/Polytomic_Lima.tre")
tree2 <- read.tree("trees/Polytomic_Upham.tre")
tree3 <- read.tree("trees/Polytomic_Wright.tre")

trees <- list(tree1, tree2, tree3)
tree_names <- c("Lima", "Upham", "Wright")

results <- data.frame(Tree = tree_names,
                      R2 = NA,
                      AICc = NA,
                      stringsAsFactors = FALSE)

fits <- list()

gdf <- geomorph.data.frame(
  Shape   = coords_sync,
  Size    = log(size_sync),
  Sex     = factor(data_sync$sex),
  Cluster = factor(data_sync$Cluster), Ecoregion = factor(data_sync$fac),
  Biome   = factor(data_sync$biome)
)

dimnames(gdf$Shape)[[3]] <- tree_expanded$tip.label
names(gdf$Size) <- tree_expanded$tip.label
names(gdf$Sex) <- tree_expanded$tip.label
names(gdf$Cluster) <- tree_expanded$tip.label
names(gdf$Biome) <- tree_expanded$tip.label

geiger::name.check(data_sync,tree1)
geiger::name.check(data_sync,tree2)
geiger::name.check(data_sync,tree3)

PCA <- gm.prcomp(gdf$Shape)

# ==== Best model and Tree for Shape ====

library(mvMORPH)

trees_list <- list(
  "Lima" = tree1,
  "Upham" = tree2,
  "Wright" = tree3
)

formulas_list <- list(
  "Null" = PCA$x[,1:16] ~ 1,
  "Full" = PCA$x[,1:16] ~ log(size_sync) + data_sync$sex
)

models <- c("BM", "OU", "EB", "lambda")

raw_results <- list()
counter <- 1

print("Starting iteractions...")

# for(t_name in names(trees_list)){
#   current_tree <- trees_list[[t_name]]
#   
#   for(f_name in names(formulas_list)){
#     current_formula <- formulas_list[[f_name]]
#     
#     for(m in models){
#       
#       cat(paste(t_name, "|", f_name, "|", m, "\n"))
#       
#       fit <- try(mvgls(
#         current_formula,
#         tree = current_tree,
#         model = m,
#         method = "PL-LOOCV",
#         REML = TRUE
#       ), silent = TRUE)
#       
#       if(inherits(fit, "try-error")){
#         raw_results[[counter]] <- data.frame(
#           Tree = t_name,
#           Type = f_name,
#           Model = m,
#           logLik = NA,
#           GIC = NA,
#           EIC = NA,
#           Param = NA,
#           AIC = NA,
#           Convergence = "Error"
#         )
#       } else {
#         gic_val <- try(as.numeric(GIC(fit)$GIC), silent = TRUE)
#         if(inherits(gic_val, "try-error")) gic_val <- NA
#         
#         eic_val <- try(as.numeric(EIC(fit)$EIC), silent = TRUE)
#         if(inherits(eic_val, "try-error")) eic_val <- NA
#         
#         param_val <- if(is.null(fit$param) || all(is.na(fit$param))) NA else as.numeric(fit$param[1])
#         aic_calc <- -2 * as.numeric(logLik(fit)) + 2 * (ncol(PCA$x[,1:16]) + length(fit$coefficients))
#         
#         raw_results[[counter]] <- data.frame(
#           Tree = t_name,
#           Type = f_name,
#           Model = m,
#           logLik = as.numeric(logLik(fit)),
#           GIC = gic_val,
#           EIC = eic_val,
#           Param = param_val,
#           AIC = aic_calc,
#           Convergence = "Success"
#         )
#       }
#       counter <- counter + 1
#     }
#   }
# }
# 
# final_table <- do.call(rbind, raw_results)
# 
# best_models <- final_table[order(final_table$GIC), ]
# 
# best_full <- subset(final_table, Type == "Full")
# best_full <- best_full[order(best_full$GIC), ]
# 
# best_null <- subset(final_table, Type == "Null")
# best_null <- best_null[order(best_full$GIC), ]
# 
# print("--- TOP 5 MODELS (Only Null) ---")
# print(head(best_null, 5))
# 
# print("--- TOP 5 MODELS (Only Full) ---")
# print(head(best_full, 5))
# 
# write.csv(final_table, "Shape_BestModels.csv")

# ==== Best Model and Tree for Size ====

library(nlme)

# 1. PREPARATION
# Tree list
trees_list <- list(
  "Lima" = tree1,
  "Upham" = tree2,
  "Wright" = tree3
)

# Dataframe for gls (ensuring alignment)
df_size <- data.frame(
  logSize = log(size_sync),
  Sex = factor(data_sync$sex)
)
rownames(df_size) <- tree1$tip.label 

# 2. DEFINE FORMULAS AND MODELS
formulas_list <- list(
  "Null" = logSize ~ 1,
  "Full" = logSize ~ Sex
)

models <- c("BM", "OU", "EB", "lambda")

results_size <- list()
counter <- 1

print("Starting Model Selection for SIZE (gls)...")

# # 3. NESTED LOOP
# for(t_name in names(trees_list)){
#   current_tree <- trees_list[[t_name]]
#   
#   for(f_name in names(formulas_list)){
#     current_form <- formulas_list[[f_name]]
#     
#     for(m in models){
#       
#       cat(paste("Running:", t_name, "|", f_name, "|", m, "\n"))
#       
#       # Define correlation structure based on model
#       cor_struct <- switch(m,
#                            "BM" = corBrownian(1, phy = current_tree),
#                            "OU" = corMartins(1, phy = current_tree), 
#                            "EB" = corBlomberg(1, phy = current_tree), # ACDC / Early Burst
#                            "lambda" = corPagel(1, phy = current_tree)
#       )
#       
#       # Run gls model using ML for comparison
#       fit <- try(gls(
#         current_form,
#         data = df_size,
#         correlation = cor_struct,
#         method = "REML" 
#       ), silent = TRUE)
#       
#       if(inherits(fit, "try-error")){
#         results_size[[counter]] <- data.frame(
#           Tree = t_name,
#           Type = f_name,
#           Model = m,
#           logLik = NA,
#           AIC = NA,
#           AICc = NA,
#           Param = NA,
#           Convergence = "Error"
#         )
#       } else {
#         # Parameter extraction (alpha, lambda, etc.)
#         coef_struct <- coef(fit$modelStruct$corStruct, unconstrained = FALSE)
#         param_val <- if(length(coef_struct) > 0) as.numeric(coef_struct) else NA
#         
#         # AICc Calculation (Small sample correction)
#         k <- length(coef(fit)) + length(coef_struct) + 1 # +1 for residual variance
#         n <- nrow(df_size)
#         aic_val <- AIC(fit)
#         aicc_val <- aic_val + (2 * k * (k + 1)) / (n - k - 1)
#         
#         results_size[[counter]] <- data.frame(
#           Tree = t_name,
#           Type = f_name,
#           Model = m,
#           logLik = as.numeric(logLik(fit)),
#           AIC = aic_val,
#           AICc = aicc_val,
#           Param = param_val,
#           Convergence = "Success"
#         )
#       }
#       counter <- counter + 1
#     }
#   }
# }
# 
# # 4. CONSOLIDATE RESULTS
# size_table <- do.call(rbind, results_size)
# 
# # Separate and Sort by AICc
# best_size_null <- subset(size_table, Type == "Null" & Param > 0)
# best_size_null <- best_size_null[order(best_size_null$AICc), ]
# 
# best_size_full <- subset(size_table, Type == "Full" & Param > 0)
# best_size_full <- best_size_full[order(best_size_full$AICc), ]
# 
# print("--- TOP 5 SIZE MODELS (Null) ---")
# print(head(best_size_null, 5))
# 
# print("--- TOP 5 SIZE MODELS (Full - Sex) ---")
# print(head(best_size_full, 5))
# 
# write.csv(size_table, "Size_BestModels.csv")

# ==== The best tree ====

tree <- tree1 # Lima
# model = "lambda"
dev.off()

# 1. Phylogenetic Signal ====

# Shape

physignal_form <- physignal(A = coords_sync, phy = tree, iter = 999)
summary(physignal_form)
plot(physignal_form)
plot(physignal_form$PACA, phylo = FALSE)
physignal_form$K.by.p

# Decomposes shape's signal

PSe.shape <- physignal.eigen(coords_sync, phy = tree)
summary(PSe.shape)
PSe.shape$eig.obs$values # Signal vanishes in multidimensional data
plot(PSe.shape$Kmult)

plot(PSe.shape)
plot(PSe.shape, type = "vectors")
KC.plot <- plot(PSe.shape$KC)
add.tree(KC.plot, tree, edge.col = 6)

# Final PhySignal
PS.shape <- physignal.z(A = coords_sync,
            lambda = "front", PAC.no = 5,
            phy = tree1, iter = 999)
summary(PS.shape)
names(PS.shape)
plot(PS.shape)

# Shape free of allometry and sex

allometry_model <- procD.lm(Shape ~ Size + Sex, data = gdf, 
                              iter = 999, RRPP = TRUE)

summary(allometry_model)

shape_residuals <- arrayspecs(allometry_model$residuals, 
                              p = dim(coords_sync)[1], 
                              k = dim(coords_sync)[2])

allometry_free_shape <- shape_residuals + array(Y.gpa$consensus, 
                                                dim(shape_residuals))

nonAllo_sync <- allometry_free_shape[,,match(tree_expanded$tip.label, 
                                     data_sync$ID)]

dimnames(nonAllo_sync)[[3]] <- tree_expanded$tip.label  

physignal_form <- physignal(A = nonAllo_sync, phy = tree, iter = 999)
summary(physignal_form)
plot(physignal_form)
# plot(physignal_form$PACA, phylo = FALSE)
physignal_form$K.by.p

# Decomposes shape's signal

# PSe.shape <- physignal.eigen(nonAllo_sync, phy = tree) # does not run with this array

PS.shape <- physignal.z(A = nonAllo_sync,
            lambda = "front", PAC.no = 5,
            phy = tree, iter = 999)
summary(PS.shape)
plot(PS.shape)

# Centroid Size

names(size_sync) <- data_sync$ID

physignal_size <- physignal(A = log(size_sync), phy = tree, iter = 999)

summary(physignal_size)
plot(physignal_size)

PS_size <- physignal.z(A = log(size_sync), phy = tree, 
                       iter = 999, lambda = "front")

summary(PS_size)
plot(PS_size)

# Select best lambdas

best_lambda_shp <- PS.shape$lambda
best_lambda_sz <- PS_size$lambda

# Final modeling with Best tree (Lima) and Best Model (lambda)
library(mvMORPH)


# 2. procD.pgls for clusters and biomes ====

# This analysis uses lambda = 1 by default
# but the estimated lambda for the different independent variables:
# Y.gpa$coords = 0.2962
# NonAlloSex = 0.3224
# logCS = 0.1673

fit_procD_cov <- procD.pgls(coords_sync ~ log(Size) + Sex, 
                                 data = gdf, phy = tree, iter = 999,
                                 SS.type = "II", lambda = best_lambda_shp)

summary(fit_procD_cov)

fit_procD_clusters <- procD.pgls(coords_sync ~ log(Size) + Sex + Cluster, 
                                         data = gdf, phy = tree, iter = 999,
                                 SS.type = "II", lambda = best_lambda_shp)
summary(fit_procD_clusters)

fit_procD_biomes <- procD.pgls(coords_sync ~ log(Size) + Sex + Biome, 
                               data = gdf, phy = tree, iter = 999,
                               SS.type = "II", lambda = best_lambda_shp)
summary(fit_procD_biomes) # biomes are full attached to phylogenetic relations

fit_procD_ecoregions <- procD.pgls(coords_sync ~ log(Size) + Sex + Ecoregion, 
                               data = gdf, phy = tree, iter = 999,
                               SS.type = "II", lambda = best_lambda_shp)
summary(fit_procD_ecoregions) # ecoregions are too


# ===========================================#
# 3. Run procD.pgls for each climatic variable
# ===========================================#

# Assign each column from env_matrix_scaled to gdf

bioclim <- data_sync[,c(13:34)]
bioclim <- scale(bioclim,center=T,scale=T)
rownames(bioclim) <-data_sync$ID

for (var in colnames(bioclim)) {
  gdf[[var]] <- bioclim[, var]
}

# Space

coords <- data_sync[,c("long", "lat")]

pcnm <- pcnm(dist(coords))
pcnm <- as.matrix(pcnm$vectors)
rownames(pcnm) <-data_sync$ID

shape <- two.d.array(coords_sync)
shape_pcnm <- forward.sel(shape,pcnm)
shape_pcnm$variables

space <- as.matrix(pcnm[,shape_pcnm$variables])
colnames(space) <- shape_pcnm$variables
rownames(space) <-data_sync$ID
name.check (tree,space)

for (var in colnames(space)) {
  gdf[[var]] <- space[, var]
}

names(gdf)

Space <- cbind(pcnm[,shape_pcnm$variables])

# example 

fit <- procD.pgls(Shape ~ log(Size) + Sex + Space + BIO1, 
                  data = gdf, phy = tree, iter = 999,
                  SS.type = "II", lambda = best_lambda_shp)
summary(fit)

# # Initialize list to store results
# results_list <- list()
# 
# # Loop through each environmental variable
# for (var_name in colnames(bioclim)) {
#   
#   # Create the formula dynamically
#   fmla <- as.formula(paste("Shape ~ Size + Sex + Space +", var_name))
#   
#   # Fit the procD.pgls model
#   fit <- try(procD.pgls(fmla, data = gdf, phy = tree,
#                         iter = 999,SS.type = "II", 
#                         lambda = best_lambda_shp), silent = TRUE)
#   
#   if (inherits(fit, "try-error")) {
#     message(paste("Model failed for", var_name, "- skipping"))
#     next
#   }
#   
#   # Extract ANOVA table
#   aov_table <- fit$aov.table
#   
#   # Get row index corresponding to the environmental variable
#   if (!(var_name %in% rownames(aov_table))) {
#     message(paste("Variable", var_name, "not found in ANOVA table - skipping"))
#     next
#   }
#   
#   row_idx <- var_name
#   
#   # Save key values
#   F_value <- aov_table[row_idx, "F"]
#   Z_value <- aov_table[row_idx, "Z"]
#   P_value <- aov_table[row_idx, "Pr(>F)"]
#   Rsq_value <- aov_table[row_idx, "Rsq"]
#   
#   n <- fit$ANOVA$n
#   K_model <- fit$ANOVA$p
#   RSS_obs <- fit$ANOVA$RSS.model[1]
#   
#   AICc_PGLS <- n * log(RSS_obs / n) + 2 * K_model + (2 * K_model * (K_model + 1)) / (n - K_model - 1)
# 
#   results_list[[var_name]] <- data.frame(
#     Variable = var_name,
#     AICc = AICc_PGLS, 
#     F_Value = F_value,
#     Z_Value = Z_value,
#     P_Value = P_value,
#     Rsq = Rsq_value
#   )
# }
# 
# # Combine all results into one dataframe
# final_table <- do.call(rbind, results_list)
# 
# # Print and save
# 
# print(final_table)
# 
# # write.csv(final_table, "procD_pgls_results.csv", row.names = FALSE)
# 
# best_row <- final_table[which.min(final_table$AICc), ]
# best_var <- best_row$Variable
# best_var

# ==================================================================#
#  PARTIAL EFFECTS PLOT (BEST PGLS MODEL VS. OLS)
# ==================================================================#

cat("\n--- Generating Partial Effects Plot for Best PGLS Model ---\n")

# Load required libraries for new blocks
library(nlme)      # For gls
library(car)       # For Anova
library(emmeans)   # For post-hoc tests

# --- 1.1 Identify Best Variable & Shape Scores (As before) ---

shape_resid_pca <- gm.prcomp(allometry_free_shape)
shape_scores_y <- shape_resid_pca$x[, "Comp1"]
predictor_x <- gdf[[best_var]]

# --- 1.2 Create Dataframe for Plotting (NOW INCLUDES FAC) ---
Fac <- as.factor(data_sync$fac) ## Environment factor
plot_df <- data.frame(
  Predictor = predictor_x,
  Shape_Score_PC1 = shape_scores_y,
  Species = data_sync$sp,  # Use Species factor from gdf
  Fac = Fac,                # Environment factor
  Space = Space,
  Size = gdf$Size, Sex = gdf$Sex
)
rownames(plot_df) <- tree$tip.label

library(nlme)
library(ggplot2)
library(tidyr)

# --- 1. Fit Models ---
# OLS (Standard Linear Regression)
fit_ols <- lm(Shape_Score_PC1 ~ Size + Sex + Space + Predictor, data = plot_df)

# PGLS (Phylogenetic Generalized Least Squares)
# Using the optimized Lambda (~0.32)
phy_cor <- corPagel(0.32, phy = tree)
fit_pgls <- gls(
  Shape_Score_PC1 ~ Size + Sex + Space + Predictor,
  data = plot_df,
  correlation = phy_cor
)

# --- 2. Extract Coefficients ---
coef_ols <- coef(fit_ols)
coef_pgls <- coef(fit_pgls)

cat(paste("Best PGLS predictor identified:", best_var, "\n"))
cat(paste("OLS (non-phylo) Slope:", round(coef_ols[2], 4), "\n"))
cat(paste("PGLS (phylo) Slope:   ", round(coef_pgls[2], 4), "\n"))

# --- 3. Calculate Adjusted Data for Plotting ---
# We remove the effects of confounding variables (Size, Sex, Space)
# to visualize the partial effect of the Predictor.
# Formula: Adjusted_Y = Residuals + Intercept + (Slope * Predictor)

# For OLS
plot_df$Adjusted_OLS <- residuals(fit_ols) + 
  coef_ols["(Intercept)"] + 
  coef_ols["Predictor"] * plot_df$Predictor

# For PGLS
plot_df$Adjusted_PGLS <- residuals(fit_pgls) + 
  coef_pgls["(Intercept)"] + 
  coef_pgls["Predictor"] * plot_df$Predictor

# --- 4. Generate Regression Lines ---
x_range <- seq(min(plot_df$Predictor), max(plot_df$Predictor), length.out = 100)

pred_lines <- data.frame(Predictor = x_range)
pred_lines$OLS <- coef_ols["(Intercept)"] + coef_ols["Predictor"] * x_range
pred_lines$PGLS <- coef_pgls["(Intercept)"] + coef_pgls["Predictor"] * x_range

# Pivot for ggplot
lines_long <- pivot_longer(pred_lines, 
                           cols = c("OLS", "PGLS"), 
                           names_to = "Model", 
                           values_to = "Shape_Predicted")

# --- 5. Plot ---
species_pch_map <- c(
  "Sapajus_apella" = "\u25A0",        # Square
  "Sapajus_cay" = "\u25BC",           # Downward-pointing triangle
  "Sapajus_libidinosus" = "\u25B2",  # Upward-pointing triangle
  "Sapajus_nigritus" = "\u25CF",      # Circle
  "Sapajus_robustus" = "\u2666",      # Diamond
  "Sapajus_xanthosternos" = "\u2605" # 5-pointed star
)

colors_env <- c("AF" ="green3","AM" = "darkgreen", "SV" = "goldenrod1")


# Using Adjusted_PGLS for points to match the main statistical test
ggplot() +
  geom_line(data = lines_long, aes(x = Predictor, y = Shape_Predicted, linetype = Model, color = Model), size = 1.2) +
  
  geom_point(data = plot_df, 
             aes(x = Predictor, y = Adjusted_PGLS, shape = Species, color = Fac), 
             size = 10, alpha = 0.8) +
  
  scale_color_manual(values = c(colors_env, "OLS"="black", "PGLS"="blue")) +
  scale_linetype_manual(values = c("PGLS" = "solid", "OLS" = "dashed")) +
  scale_shape_manual(values = species_pch_map) +
  
  labs(title = paste("Partial Effect of", best_var, "on Shape"),
       x = best_var,
       y = "Shape Score (Adjusted for Size, Sex, Space)", 
       caption = "Points represents shape residuals + predictor effect (PGLS adjusted)") +
  
  guides(
    shape = guide_legend(order = 1),
    fill = guide_legend(order = 2),
    color = guide_legend(order = 3),
    linetype = guide_legend(order = 3))
    +
  theme_bw()


pca <- gm.prcomp(allometry_free_shape)
pc1 <- pca$x[,1]

min1 <- min(pc1)
max1 <- max(pc1)
M <- mshape(coords_sync)

pred_min <- shape.predictor( 
  allometry_free_shape, 
  x = pc1, 
  pred1 = min(pc1), min(pc1))

pred_max <- shape.predictor( 
  allometry_free_shape, 
  x = pc1, 
  pred1 = max(pc1), max(pc1))

drawinglandmark <- readland.tps("outline2/outline.tps")
Sapajusoutline <- warpRefOutline(file="outline2/outline.txt", drawinglandmark[,,1], M)
dev.off()

par(mfrow=c(1,2), xpd=FALSE)

GP <- gridPar(n.col.cell = 100, pt.bg = "gray", pt.size = 0.8, tar.pt.bg = "cyan",
              tar.pt.size = 0.8, tar.out.col = "gray10", tar.out.cex = 0.5,
              grid.col = "white", grid.lwd = 0.5, txt.pos = 1, txt.col = "steelblue") # Custom grids


plotRefToTarget(M, pred_min$pred1,
                main="Min PC1", mag=2,
                outline=Sapajusoutline$outline, gridPars=GP,
                method="points")

plotRefToTarget(M, pred_max$pred1,
                main="Min PC1", mag=2,
                outline=Sapajusoutline$outline, gridPars=GP,
                method="points")

par(mfrow=c(1,1), xpd=FALSE)


# ==== Phylogenetic Partial Least Squares ====

# --- 10.1 Prepare Allometry- and Sex-Free Shape (N=83) ---
# Residualize shape against confounders (Size and Sex)
# We use procD.lm because PGLS (K=0.052) is negligible and this is faster.

names(data_sync)
climatic_vars <- as.matrix(data_sync[, 13:34])
climatic_vars <- scale(climatic_vars, center = TRUE, scale = TRUE)
rownames(climatic_vars) <- data_sync$ID

fit_confounders <- procD.pgls(Shape ~ Size + Sex + Space, 
                              data = gdf, 
                              phy = tree,
                              iter = 999, lambda = best_lambda_shp)
summary(fit_confounders)

shape_residuals <- arrayspecs(fit_confounders$pgls.residuals, 
                              p = dim(Y.gpa$coords)[1], 
                              k = dim(Y.gpa$coords)[2])

allometry_free_shape <- shape_residuals + array(Y.gpa$consensus, 
                                                dim(shape_residuals))

# ==================================================================#
# 1. PHYLOGENETIC MANUAL PLS (with Optimized Lambda) - FINAL VERSION
# ==================================================================#
# OBJECTIVE: Run the PLS with Lambda=0.32
# (Transforming BOTH blocks)

library(geomorph)
library(ape)

# --- 1. DEFINE THE CORRECT LAMBDA (from your diagnoses) ---
lambda <- best_lambda_shp # (from Free-Form, Report 26)
n_species <- length(tree$tip.label)

# --- 2. PREPARE THE DATA BLOCKS (As you did) ---
# (Assuming 'allometry_free_shape' and 'climatic_vars' exist)

n_spec <- dim(allometry_free_shape)[3]
p <- dim(allometry_free_shape)[1]
k <- dim(allometry_free_shape)[2]

# Block Y (Shape)
shape_mat <- matrix(aperm(allometry_free_shape, c(3,1,2)),
                    
                    nrow = n_spec, ncol = p*k)
rownames(shape_mat) <- data_sync$ID

# Block X (Climate)
clim_mat <- climatic_vars
rownames(clim_mat) <- data_sync$ID

# --- 3. CREATE THE TRANSFORMATION MATRIX (Its Logic) ---
# Create the VCV matrix scaled by Lambda
tree_mv <- tree # Use the tree from Lima
tree_mv$edge.length[tree_mv$edge.length == 0] <- 1e-8 # Avoid errors
tree_mv$root.edge <- 0

# Create the C_lambda matrix (Pagel)
C_lambda <- (1 - lambda) * diag(n_species) + lambda * vcv(tree_mv)
rownames(C_lambda) <- tree_mv$tip.label
colnames(C_lambda) <- tree_mv$tip.label

# Ensure that the data is in the same order as the VCV matrix
shape_mat <- shape_mat[tree_mv$tip.label, ]
clim_mat <- clim_mat[tree_mv$tip.label, ]

# Create the Cholesky matrix
L <- t(chol(C_lambda))

# --- 4. Transform BOTH blocks with the CORRECTED VCV ---
# (This is the correction)
corrected_phy_shape <- solve(L) %*% shape_mat
corrected_phy_climb <- solve(L) %*% clim_mat

# --- 5. Run the PLS (OLS) on the transformed data ---
cat("\n--- PHYLOGENETIC PLS (WITH LAMBDA = 0.37) ---\n")
pls_results <- two.b.pls(A1 = corrected_phy_climb,
                               A2 = corrected_phy_shape,
                               iter = 999,
                               print.progress = TRUE)

summary(pls_results)

# ==== 10.3 PLS Plots ====

print("Generating PLS plots...")

scores_df <- as.data.frame(pls_results$YScores[, 1:2]) # YScores = Shape Scores
colnames(scores_df) <- c("PLS1", "PLS2")

eig <- pls_results$svd$d^2
var1 <- round((eig[1] / sum(eig)) * 100, 1)
var2 <- round((eig[2] / sum(eig)) * 100, 1)

# Add factors from gdf and data_sync
scores_df$Cluster <- gdf$Cluster
scores_df$Species <- factor(data_sync$sp)

# Prepare PCH (plotting symbols)
pch_species_vector <- as.numeric(data_sync$pch)
scores_df$PCH <- pch_species_vector

species_names_factor <- factor(data_sync$sp)
unique_species_names <- levels(species_names_factor)
species_pch_map <- c(
  "Sapajus_apella" = "\u25A0",        # Square
  "Sapajus_cay" = "\u25BC",           # Downward-pointing triangle
  "Sapajus_libidinosus" = "\u25B2",  # Upward-pointing triangle
  "Sapajus_nigritus" = "\u25CF",      # Circle
  "Sapajus_robustus" = "\u2666",      # Diamond
  "Sapajus_xanthosternos" = "\u2605" # 5-pointed star
)
names(species_pch_map) <- unique_species_names

# Colors
cluster_colors <- c(
  "Hot & Thermally Aseasonal"        = "#FDE725FF", 
  "Warm & Moderately Seasonal"       = "#35B779FF",  
  "Thermally Stable & Dry-Seasonal"  = "#31688EFF",  
  "Elevated & Thermally Seasonal"    = "#440154FF"   
)

names(cluster_colors) <- levels(scores_df$Cluster)

# if (!exists("cluster_colors")) {
#   print("Warning: 'cluster_colors' not defined. Using default palette.")
#   cluster_colors <- scales::hue_pal()(length(unique(scores_df$Cluster)))
# }

scores_df$Cluster <- factor(scores_df$Cluster, levels = names(cluster_colors))
scores_df$Eco <- factor(data_sync$fac, levels = names(colors_env))

# Convex Hulls
hull_data <- scores_df %>%
  group_by(Eco) %>%
  filter(n() > 2) %>%
  slice(chull(PLS1, PLS2)) %>%
  ungroup()

# ggplot Scores Plot
plot_scores_gg <- ggplot(scores_df, aes(x = PLS1, y = PLS2)) +
  geom_polygon(
    data = hull_data,
    aes(fill = Eco, color = Eco),
    alpha = 0.20,
    linewidth = 1
  ) +
  geom_point(
    aes(shape = Species, fill = Eco, color = Eco),
    size = 10,
    stroke = 0.8
  ) +
  scale_shape_manual(name = "Species", values = species_pch_map) +
  scale_fill_manual(name = "Ecoregions / Biomes", values = colors_env) +
  scale_color_manual(name = "Ecoregions / Biomes", values = colors_env) +
  labs(
    x = paste0("PLS1 (", var1, "%)"),
    y = paste0("PLS2 (", var2, "%)"),
    title = "PLS: Shape (Pure) vs. Climate"
  ) +
  theme_bw(base_size = 14) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(linetype = "dashed", color = "grey90"),
    legend.position = "right",
    plot.title = element_text(hjust = 0.5)
  )

print(plot_scores_gg)


# ==== 10.4 PLS Loadings (Environmental Vectors) ====
loadings_df <- as.data.frame(pls_results$left.pls.vectors[, 1:2]) # Left = A1 = Climate
colnames(loadings_df) <- c("PLS1", "PLS2")
rownames(loadings_df) <- colnames(bioclim)

loadings_plot_df <- data.frame(
  Variable = rownames(loadings_df),
  Loading = loadings_df$PLS1
)
loadings_plot_df <- loadings_plot_df[order(loadings_plot_df$Loading, decreasing = TRUE), ]
loadings_plot_df$Color <- ifelse(loadings_plot_df$Loading > 0, "steelblue", "firebrick")

# Loadings Bar Plot
print("Generating PLS1 Loadings bar plot...")
par(mar = c(12, 4, 4, 1)) # Increase bottom margin
barplot(
  loadings_plot_df$Loading,
  names.arg = loadings_plot_df$Variable,
  col = loadings_plot_df$Color,
  border = NA,
  las = 2, # Perpendicular labels
  ylab = "Loading PLS1",
  main = "Contribution of Climatic Variables to PLS1"
)
abline(h = 0, col = "black", lwd = 1)
par(mar = c(5, 4, 4, 2)) # Reset default margin

cat("\n--- PLS Loadings Plot finished ---\n")


# ==== 10.5 PLS Scores + Loadings Combined Plot ====
# (This is the new block you provided)

print("Generating combined PLS Scores & Loadings plot...")

# Re-create loadings_df for scaling
loadings_df <- as.data.frame(pls_results$left.pls.vectors[, 1:2])
colnames(loadings_df) <- c("PLS1", "PLS2")
loadings_df$Variable <- colnames(climatic_vars)

# Calculate scaling factor
score_range <- apply(scores_df[, c("PLS1", "PLS2")], 2, function(x) max(x) - min(x))
loading_range <- apply(loadings_df[, c("PLS1", "PLS2")], 2, function(x) max(x) - min(x))
scale_factor <- min(score_range / loading_range) * 1

loadings_df$PLS1_scaled <- loadings_df$PLS1 * scale_factor
loadings_df$PLS2_scaled <- loadings_df$PLS2 * scale_factor

# Combined plot
ggplot(scores_df, aes(x = PLS1, y = PLS2)) +
  geom_polygon(
    data = hull_data,
    aes(fill = Eco, color = Eco),
    alpha = 0.20,
    linewidth = 1
  ) +
  geom_point(
    aes(shape = Species, fill = Eco, color = Eco),
    size = 10,
    stroke = 0.8
  ) +
  geom_segment(
    data = loadings_df,
    aes(x = 0, y = 0, xend = PLS1_scaled, yend = PLS2_scaled),
    arrow = arrow(length = unit(0.3, "cm")),
    color = "blue",
    linewidth = 0.8
  ) +
  geom_text(
    data = loadings_df,
    aes(x = PLS1_scaled, y = PLS2_scaled, label = Variable),
    hjust = 0.5, vjust = -0.5,
    size = 6
  ) +
  scale_shape_manual(name = "Species", values = species_pch_map) +
  scale_fill_manual(name = "Ecoregions / Biomes", values = colors_env) +
  scale_color_manual(name = "Ecoregions / Biomes", values = colors_env) +
  labs(
    x = paste0("PLS1 (", var1, "%)"),
    y = paste0("PLS2 (", var2, "%)"),
    title = "Phylogenetic PLS: Shape (Pure) vs. Climate"
  ) +
  theme_bw(base_size = 14) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(linetype = "dashed", color = "grey90"),
    legend.position = "right",
    plot.title = element_text(hjust = 0.5)
  )


cat("\n--- Combined PLS Scores & Loadings Plot finished ---\n")

# Clean Plot

scale_factor <- (max(scores_df$PLS1) - min(scores_df$PLS1)) /
  (max(loadings_df$PLS1) - min(loadings_df$PLS1))

loadings_df$PLS1_scaled <- loadings_df$PLS1 * scale_factor
loadings_df$PLS2_scaled <- loadings_df$PLS2 * scale_factor * 0.05

x_buffer <- (max(scores_df$PLS1) - min(scores_df$PLS1)) * 0.20
y_buffer <- (max(scores_df$PLS2) - min(scores_df$PLS2)) * 0.20

text_offset_x <- (max(scores_df$PLS1) - min(scores_df$PLS1)) * 0.03
text_offset_y <- (max(scores_df$PLS2) - min(scores_df$PLS2)) * 0.03

loadings_df$LabelX <- loadings_df$PLS1_scaled + text_offset_x
loadings_df$LabelY <- loadings_df$PLS2_scaled + text_offset_y

max_right <- max(scores_df$PLS1, loadings_df$PLS1_scaled, loadings_df$LabelX)
min_left  <- min(scores_df$PLS1, loadings_df$PLS1_scaled, loadings_df$LabelX)

range_x <- max_right - min_left
hard_left  <- min_left  - range_x * 0.30
hard_right <- max_right + range_x * 0.30

ggplot(scores_df, aes(x = PLS1, y = PLS2)) +
  geom_polygon(
    data = hull_data,
    aes(fill = Eco, color = Eco),
    alpha = 0.20,
    linewidth = 1
  ) +
  geom_point(
    aes(shape = Species, fill = Eco, color = Eco),
    size = 10, stroke = 0.8
  ) +
  geom_segment(
    data = loadings_df,
    aes(x = 0, y = 0, xend = PLS1_scaled, yend = PLS2_scaled),
    arrow = arrow(length = unit(0.3, "cm")),
    color = "blue",
    linewidth = 0.8
  ) +
  geom_text(
    data = loadings_df,
    aes(x = LabelX, y = LabelY, label = Variable),
    size = 6
  ) +
  scale_shape_manual(values = species_pch_map) +
  scale_fill_manual(values = colors_env) +
  scale_color_manual(values = colors_env) +
  labs(
    x = paste0("PLS1 (", var1, "%)"),
    y = paste0("PLS2 (", var2, "%)"),
    title = "Phylogenetic PLS: Shape (Pure) vs. Climate"
  ) +
  theme_bw(base_size = 14) +
  theme(
    plot.margin = margin(0, 0, 0, 0),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(linetype = "dashed", color = "grey90"),
    legend.position = "none",
    plot.title = element_text(hjust = 0.5)
  ) +
  coord_cartesian(
    xlim = c(hard_left, hard_right),
    ylim = c(min(scores_df$PLS2) - y_buffer,
             max(scores_df$PLS2) + y_buffer),
    clip = "off"
  )


# ==================================================================#
# EXTREMES BY EACH PLS AXIS SEPARATELY (NO COMBOS)
# ==================================================================#
library(geomorph)

M <- mshape(allometry_free_shape)

preds_pls1 <- shape.predictor(
  allometry_free_shape, 
  x = scores_df$PLS1,
  min_pls1 = min(scores_df$PLS1),
  max_pls1 = max(scores_df$PLS1)
)

preds_pls2 <- shape.predictor(
  allometry_free_shape, 
  x = scores_df$PLS2,
  min_pls2 = min(scores_df$PLS2),
  max_pls2 = max(scores_df$PLS2)
)

par(mfrow=c(2,2), mar=c(1,1,3,1)) 

plotRefToTarget(M, preds_pls1$min_pls1, 
                main="Min PLS1 (Negative)", mag=2, 
                outline=Sapajusoutline$outline, method="points", 
                gridPars = GP)

plotRefToTarget(M, preds_pls1$max_pls1, 
                main="Max PLS1 (Positive)", mag=2, 
                outline=Sapajusoutline$outline, method="points", 
                gridPars = GP)

plotRefToTarget(M, preds_pls2$min_pls2, 
                main="Min PLS2 (Negative)", mag=2, 
                outline=Sapajusoutline$outline, method="points", 
                gridPars = GP)

plotRefToTarget(M, preds_pls2$max_pls2, 
                main="Max PLS2 (Positive)", mag=2, 
                outline=Sapajusoutline$outline, method="points", 
                gridPars = GP)

par(mfrow=c(1,1))

# ==================================================================#
# 1. Selecting the Best model for log(Size)
# ==================================================================#

library(nlme)
library(ape)
library(dplyr)

plot_df_base <- data.frame(
  Size = log(size_sync),
  Sex = factor(data_sync$sex)
)
rownames(plot_df_base) <- tree$tip.label

env_vars <- colnames(bioclim)
models <- c("lambda")
results_list <- list()
Size_Space <- forward.sel(size_sync,pcnm)
SSpace <- cbind(pcnm[,Size_Space$variables])
colnames(SSpace) <- Size_Space$variables

# for(var_name in env_vars){
#   
#   # Prepare data for current variable
#   plot_df <- plot_df_base
#   plot_df$EnvVar <- bioclim[, var_name]
#   
#   # Bind spatial vectors to dataframe to use in formula
#   # We use cbind to ensure all PCNMs are available
#   plot_df_full <- cbind(plot_df, SSpace)
#   
#   # Construct formula dynamically to include all PCNMs
#   spatial_terms <- paste(colnames(SSpace), collapse = " + ")
#   fmla <- as.formula(paste("Size ~ Sex +", spatial_terms, "+ EnvVar"))
#   
#   model_results <- list()
#   
#   for(model in models){
#     
#     # Define correlation structure
#     cor_obj <- switch(model,
#                       BM = corBrownian(1, phy = tree),
#                       OU = corMartins(1, phy = tree, fixed = FALSE),
#                       lambda = corPagel(0.16, phy = tree, fixed = TRUE) # Fixed lambda as requested
#     )
#     
#     # Fit model
#     fit <- try(gls(
#       fmla,
#       data = plot_df_full,
#       correlation = cor_obj,
#       method = "ML"
#     ), silent = TRUE)
#     
#     if(inherits(fit, "try-error")){
#       model_results[[model]] <- data.frame(
#         Model = model, logLik = NA, AICc = NA, Coef = NA, P_value = NA
#       )
#     } else {
#       summ <- summary(fit)
#       coef_val <- summ$tTable["EnvVar","Value"]
#       p_val <- summ$tTable["EnvVar","p-value"]
#       
#       n <- nrow(plot_df_full)
#       k <- length(coef(fit))
#       RSS <- sum(resid(fit)^2) # Approximate for GLS context
#       AICc <- AIC(fit) + (2*k*(k+1))/(n-k-1) # Use AIC(fit) directly
#       
#       model_results[[model]] <- data.frame(
#         Model = model,
#         logLik = as.numeric(logLik(fit)),
#         AICc = AICc,
#         Coef = coef_val,
#         P_value = p_val
#       )
#     }
#   }
#   results_list[[var_name]] <- do.call(rbind, model_results)
# }  
#   
# final_table <- do.call(rbind, lapply(names(results_list), function(v){
#   df <- results_list[[v]]
#   df$Variable <- v
#   df
# }))
# 
# best_models <- final_table %>% group_by(Variable) %>% slice_min(AICc, n=1)
# print(best_models)
# 
# best_models_ordered <- best_models %>% arrange(AICc)
# print(best_models_ordered)
# write.csv(best_models_ordered, "Size_PGLS.csv")
# 
# best_row <- best_models_ordered[1, ]
# best_var <- best_row$Variable
# best_model <- best_row$Model
# cat("Best variable:", best_var, "with:", best_model, "\n")

Fac <- factor(data_sync$fac)          
colors_env <- c("AF"="green3","AM"="darkgreen","SV"="goldenrod1")
names(colors_env) <- levels(Fac)

plot_df <- data.frame(
  Predictor = predictor_x,
  Shape_Score_PC1 = shape_scores_y,
  Species = data_sync$sp,  # Use Species factor from gdf
  Fac = Fac,                # Environment factor
  Space = Space,
  Size = gdf$Size, Sex = gdf$Sex,
  EnvVar = bioclim[,best_var]
)

rownames(plot_df) <- tree$tip.label

library(nlme)
library(ggplot2)
library(tidyr)

# --- 1. Fit Models (Including all covariates) ---
# We assume 'SSpace' represents the spatial eigenvectors (PCNMs) selected for Size
# Ensure 'SSpace' is bound to plot_df or available in the environment

fit_ols <- lm(Size ~ Sex + SSpace + EnvVar, data = plot_df)

# PGLS with fixed Lambda for Size (~0.16)
fit_pgls <- gls(Size ~ Sex + SSpace + EnvVar, 
                data = plot_df,
                correlation = corPagel(0.16, phy = tree, fixed = TRUE),
                method = "ML")

# --- 2. Extract Coefficients ---
coef_ols <- coef(fit_ols)
coef_pgls <- coef(fit_pgls)

# --- 3. Calculate Adjusted Data ---
# We remove the effects of Sex and Spatial Structure to visualize 
# the partial effect of the environmental variable (EnvVar).
# Formula: Adjusted_Y = Residuals + Intercept + (Slope * EnvVar)

# For OLS
plot_df$Adjusted_Size_OLS <- residuals(fit_ols) + 
  coef_ols["(Intercept)"] + 
  coef_ols["EnvVar"] * plot_df$EnvVar

# For PGLS (Use this for plotting points)
plot_df$Adjusted_Size_PGLS <- residuals(fit_pgls) + 
  coef_pgls["(Intercept)"] + 
  coef_pgls["EnvVar"] * plot_df$EnvVar

# --- 4. Generate Regression Lines ---
x_range <- seq(min(plot_df$EnvVar), max(plot_df$EnvVar), length.out = 100)

pred_lines <- data.frame(EnvVar = x_range)
pred_lines$OLS <- coef_ols["(Intercept)"] + coef_ols["EnvVar"] * x_range
pred_lines$PGLS <- coef_pgls["(Intercept)"] + coef_pgls["EnvVar"] * x_range

# Pivot for ggplot
lines_long <- pivot_longer(pred_lines, 
                           cols = c("OLS", "PGLS"), 
                           names_to = "Model", 
                           values_to = "Predicted_Size")

library(ggplot2)

# --- 1. Define Mappings ---
# We convert the requested Unicode shapes to their R 'pch' equivalents
# (21-25) to allow separate control of fill (Environment) and border (Black).
species_pch_map_filled <- c(
  "Sapajus_apella" = 22,        # Square (Filled)
  "Sapajus_cay" = 25,           # Downward Triangle (Filled)
  "Sapajus_libidinosus" = 24,   # Upward Triangle (Filled)
  "Sapajus_nigritus" = 21,      # Circle (Filled)
  "Sapajus_robustus" = 23,      # Diamond (Filled)
  "Sapajus_xanthosternos" = 21  # Circle (Filled) - Reusing 21 as 21-25 are the only filled options
)

colors_env <- c("AF" ="green3", "AM" = "darkgreen", "SV" = "goldenrod1")

# --- 2. Plot (corrigido) ---
ggplot() +
  # A. Regression Lines
  geom_line(
    data = lines_long,
    aes(x = EnvVar, y = Predicted_Size, linetype = Model, color = Model),
    size = 1.2
  ) +
  
  # B. Points (PGLS)
  geom_point(
    data = plot_df,
    aes(x = EnvVar, y = Adjusted_Size_PGLS, shape = Species, fill = Fac),
    size = 8, alpha = 0.8, color = "black", stroke = 0.8
  ) +
  
  # --- 3. Scales ---
  scale_color_manual(
    name = "Regression Model",
    values = c("OLS" = "black", "PGLS" = "blue")
  ) +
  scale_linetype_manual(
    name = "Regression Model",
    values = c("PGLS" = "solid", "OLS" = "dashed")
  ) +
  scale_fill_manual(
    name = "Environment",
    values = colors_env,
    guide = guide_legend(override.aes = list(shape = 21, size = 4))
  ) +
  scale_shape_manual(
    name = "Species",
    values = species_pch_map_filled
  ) +
  
  # --- 4. Tema & Labels ---
  labs(
    title = paste("Partial Effect of", best_var, "on Mandibular Size"),
    x = "Mean Temperature of Driest Quarter (BIO9)",
    y = "log(Centroid Size) (Adjusted for Sex & Space)",
    caption = "PGLS (solid blue) vs. OLS (dashed black)"
  ) +
  theme_bw(base_size = 12) +
  theme(
    legend.position = "right",
    plot.title = element_text(hjust = 0.5, size = 14),
    axis.title = element_text(size = 12)
  ) +
  guides(
    shape = guide_legend(order = 1),
    fill = guide_legend(order = 2),
    color = guide_legend(order = 3),
    linetype = guide_legend(order = 3)
  )

# ==================================================================#
# 2. PHYLOGENETIC vs. NON-PHYLOGENETIC ANOVA FOR log(Size)
# ==================================================================#

cat("\n--- Running OLS vs. PGLS comparisons for log(Size) ---\n")

# Load required libraries
library(nlme)
library(car)

# --- 2.1 Prepare Dataframe (as before) ---
size_df <- data.frame(
  logSize = gdf$Size,
  Sex = gdf$Sex,
  Cluster = gdf$Cluster,
  Biome = gdf$Biome,
  Ecoregion = factor(data_sync$fac)
)

rownames(size_df) <- tree$tip.label 

# Define the PGLS correlation structure (BM, lambda = 1)
phy_cor <- corPagel(best_lambda_sz, phy = tree)

# --- 2.2 Test 1: Cluster ---
cat("\n--- [Cluster] OLS (non-phylo) Model ---\n")
fit_ols_cluster <- gls(logSize ~ Sex + Cluster, data = size_df, method = "ML")
print(car::Anova(fit_ols_cluster, type="II"))

cat("\n--- [Cluster] PGLS (phylo) Model ---\n")
fit_pgls_cluster <- gls(logSize ~ Sex + Cluster, data = size_df, correlation = phy_cor, method = "ML")
print(car::Anova(fit_pgls_cluster, type="II"))


# --- 2.3 Test 2: Ecoregion ---
cat("\n--- [Ecoregion] OLS (non-phylo) Model ---\n")
fit_ols_ecoregion <- gls(logSize ~ Sex + Ecoregion, data = size_df, method = "ML")
print(car::Anova(fit_ols_ecoregion, type="II"))

cat("\n--- [Ecoregion] PGLS (phylo) Model ---\n")
fit_pgls_ecoregion <- gls(logSize ~ Sex + Ecoregion, data = size_df, correlation = phy_cor, method = "ML")
print(car::Anova(fit_pgls_ecoregion, type="II"))

# --- 2.4 Test 3: Biome ---
cat("\n--- [Biome] OLS (non-phylo) Model ---\n")
fit_ols_biome <- gls(logSize ~ Sex + Biome, data = size_df, method = "ML")
print(car::Anova(fit_ols_biome, type="II"))

cat("\n--- [Biome] PGLS (phylo) Model ---\n")
fit_pgls_biome <- gls(logSize ~ Sex + Biome, data = size_df, correlation = phy_cor, method = "ML")
print(car::Anova(fit_pgls_biome, type="II")) # PERFECT COLINEARITY

cat("\n--- Size ANOVAs (OLS vs PGLS) Complete ---\n")
?car::Anova
# --- Visualization in Phylogenetic Context ---

# Perform phylogenetic PCA

avgterm <- paste(as.factor(data_sync$sp))
x <- two.d.array(coords_sync)
means <- rowsum(x, avgterm) / as.vector(table(avgterm))

shape_means <- arrayspecs(means, dim(coords_sync)[1], dim(coords_sync)[2], sep = NULL)

phylo_pca <- gm.prcomp(shape_means, phy = tree_work,
                       align.to.phy = FALSE, GLS = FALSE, transform = FALSE)
summary(phylo_pca)

# Plot phylomorphospace
plot(phylo_pca, phylo = TRUE, main = "Phylomorphospace", pch = 21, bg = "lightblue", cex = 2)

# Visualize ancestral shape reconstructions
ancestral_shapes <- arrayspecs(phylo_pca$ancestors, dim(shape_means)[1], dim(shape_means)[2])
root_shape <- ancestral_shapes[, , 1] # Shape at the root

par(mfrow = c(3, 2))
plotRefToTarget(mshape(shape_means), root_shape, mag = 3, method = "points", outline = Sapajusoutline$outline, gridPars = GP)
title("Mean Shape vs. Root Ancestor")

# Compare root to a tip (e.g., Sapajus apella)
plotRefToTarget(root_shape, shape_means[, , "Sapajus_xanthosternos"], mag = 3, method = "points", outline = Sapajusoutline$outline, gridPars = GP)
title("Root vs. S. xanthosternos")

par(mfrow = c(1, 1))

# Map continuous trait (size) on the phylogeny
names(data)

size_means <- rowsum(size_sync, avgterm) / as.vector(table(avgterm))
names(size_means) <- rownames(size_means)
contMap(compute.brlen(tree_work, method = "Grafen"), size_means, fsize = 1)

library(phytools)

obj <- contMap(
  compute.brlen(tree_work, method = "Grafen"),
  size_means,
  fsize = 1,
  plot = FALSE
)

obj$lims  
obj$tree$edge.length  

obj$lwd <- 25

plot(obj, lwd = obj$lwd, fsize = 1)
