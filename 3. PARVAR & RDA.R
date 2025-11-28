            # ==== PHYLOGENETIC COMPARATIVE ANALYSIS ====

library(ggplot2)
library(geomorph)
library(phytools)
library(vegan)
library(geiger)
library(phangorn)
library(tidyverse)
library(RColorBrewer)
library(usdm)
library(spaMM)
library(foreach)
library(doParallel)
library(spdep)
library(adespatial)

# gm analyses ====

sapajus<-readland.tps("tps/avglocsex.tps", specID="ID")
gpa<-gpagen(sapajus)
dim(sapajus)

data<-read.csv("Plans/Factors.csv", sep=",")
dim(data)
names(data)

#link<-define.links(gpa$consensus, ptsize = 2, links = NULL)
#write.table(link,file="link.txt")

link <- read.table("tps/link.txt")


drawinglandmark<-readland.tps("outline2/outline.tps")
outline<-read.table("outline2/outline.txt", header=FALSE)

summary(drawinglandmark)
mshape<-mshape(gpa$coords)
summary(mshape)

summary(outline)
Sapajusoutline<-warpRefOutline(file = "outline2/outline.txt", drawinglandmark[,,1],mshape)#run dev.off() in case of error message ##set outline configuration

dev.off()

names(data)

### Phylogeny

library(ape)
library(dplyr)

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
dimnames(gpa$coords)[[3]] <- data$ID

tree <- read.tree("trees/Polytomic_Lima.tre")

# Inspect updated tip labels

tree$tip.label
length(tree$tip.label)
tree$edge.length

plot(tree, type = "phylogram", cex = 0.5)
axisPhylo()

# Reorder data according to tree
data_sync <- data[match(tree$tip.label, data$ID), ]
rownames(data_sync) <- data_sync$ID

# Reorder GPA coordinates according to the reordered data

coords_sync <- gpa$coords[,, match(tree$tip.label, data$ID)]
dimnames(coords_sync)[[3]] <- tree$tip.label

size_sync <- gpa$Csize[match(tree$tip.label, data$ID)]
names(size_sync) <- tree$tip.label

sex_free <- procD.lm(coords_sync ~ data$sex,
                     iter = 999, RRPP = TRUE)

shape_residuals <- arrayspecs(sex_free$residuals, 
                              p = dim(coords_sync)[1], 
                              k = dim(coords_sync)[2])

sex_free_shape <- shape_residuals + array(coords_sync, 
                                                dim(shape_residuals))
gdf <- geomorph.data.frame(
  Shape = coords_sync,
  Free_Shape = sex_free_shape,
  Size = size_sync,
  Sex = factor(data_sync$sex),
  Species = factor(data_sync$sp),
  Biome = factor(data_sync$biome),
  Ecoregion = factor(data_sync$fac),
  Cluster = factor(data_sync$Cluster)
)

dimnames(gdf$Shape)[[3]] <- data_sync$ID
dimnames(gdf$Free_Shape)[[3]] <- data_sync$ID

name.check(tree,data)

cov <- vcv.phylo(tree)
phylo <- prcomp(cov)
summary(phylo)

phy <- as.matrix(phylo$x)
name.check(tree, phy)

### Shape (Y)

pcamorph <- gm.prcomp(coords_sync)
summary(pcamorph)

plot(pcamorph, pch = as.numeric(data_sync$pch), col = as.factor(data_sync$fac))

pcsgm <- as.matrix(pcamorph$x)
rownames(pcsgm) <- data_sync$ID

name.check(tree,pcsgm)

# Space

coords <- data_sync[,c("long", "lat")]

# geographicaldistance <- as.matrix(dist(coords))
# pcnmcoord<-prcomp(geographicaldistance)
# summary(pcnmcoord)

# pcnm<-pcnmcoord$x
# rownames(pcnm) <-data$sp_ives

#or

pcnmcoord <- pcnm(dist(coords))
pcnm <- as.matrix(pcnmcoord$vectors)
summary(pcnm)

rownames(pcnm) <-data_sync$ID
name.check (tree,pcnm)

#write.table(pcnm, "pcnm.xls", row.names= TRUE, col.names = TRUE, sep = " ")

### Environment

names(data)
bioclim <- data_sync[,c(13:34)]
bioclim <- scale(bioclim,center=T,scale=T)
rownames(bioclim) <-data_sync$ID
head(bioclim)
name.check (tree,bioclim)

data_std <- data_sync
data_std[,c(13:34)] <- scale(data_sync[,c(13:34)])

envpca <- prcomp(bioclim)
summary(envpca)

env <- as.matrix(envpca$x)
rownames(env) <- data_sync$ID
name.check(tree,env)


# Covariables (Size + Sex)

shapiro.test(data_sync$CS)
shapiro.test(log(data_sync$CS))

size <- as.matrix(log(size_sync))
names(size) <- data_sync$ID
name.check (tree,size)

sex <- as.matrix(as.numeric(as.factor(data_sync$sex)))
rownames(sex) <- data_sync$ID
name.check (tree,sex)

### Selection


#install.packages("packfor", repos="http://R-Forge.R-project.org")

shape_pcnm <- forward.sel(pcsgm,pcnm) 
shape_pcnm$variables

size_pcnm <- forward.sel(size,pcnm) ## 31, 5, 2, 7 /// 51,80,4,54,79
size_pcnm$variables

phy_clean <- forward.sel(pcsgm,phy) 
phy_clean$variables

env_clean <- forward.sel(pcsgm,env)
env_clean$variables

sex_matrix <- model.matrix(~ sex - 1, data = data_sync) 
sex_sel <- forward.sel(pcsgm, sex_matrix, nperm = 9999, alpha = 0.05)

data_sync$LogCS <- log(size_sync)

size_matrix <- model.matrix(~ LogCS - 1, data = data_sync) 
size_sel <- forward.sel(pcsgm, size_matrix, nperm = 9999, alpha = 0.05)

size_sel_data <- size_matrix[, size_sel$variables, drop = FALSE]
sex_sel_data <- sex_matrix[, sex_sel$variables, drop = FALSE]

covariables_sel <- cbind(size_sel_data, sex_sel_data) 


## Elements

spatial <- pcnm[,shape_pcnm$variables]


# environmental vars

env_ <- env[,env_clean$variables]

# Phylogenetic vectors

phy_ <- phy[,phy_clean$variables]


# =================== Variation Partition (PARVAR) ===================

parvar <- varpart(pcsgm,
                  ~ phy_,
                  ~ covariables_sel,
                  ~ env_,
                  ~ spatial)

plot(parvar,cex=1.2,Xnames=c("Phylogeny","Size + Sex", "Environment", "Spatial Structure"),bg=c("yellow", "deeppink","navy","green3"))


### for all parvar's variables #=====================#

y <- pcsgm
ST <- spatial
Env <- env_
Phy <- phy_
Cov <- covariables_sel

p_values <- c()

rda.result <- rda(y ~ Phy)
anova_result <- anova(rda.result, step = 1000, perm.max = 2000)
p_values <- c(p_values, anova_result$`Pr(>F)`[1])

rda.result_alo <- rda(y ~ Cov)
anova_result <- anova(rda.result_alo, step = 1000, perm.max = 2000)
p_values <- c(p_values, anova_result$`Pr(>F)`[1])

rda.result <- rda(y ~ Env)
anova_result <- anova(rda.result, step = 1000, perm.max = 2000)
p_values <- c(p_values, anova_result$`Pr(>F)`[1])

rda.result <- rda(y ~ ST)
anova_result <- anova(rda.result, step = 1000, perm.max = 2000)
p_values <- c(p_values, anova_result$`Pr(>F)`[1])

rda.result <- rda(y ~ Phy + Cov)
anova_result <- anova(rda.result, step = 1000, perm.max = 2000)
p_values <- c(p_values, anova_result$`Pr(>F)`[1])

rda.result <- rda(y ~ Phy + Env)
anova_result <- anova(rda.result, step = 1000, perm.max = 2000)
p_values <- c(p_values, anova_result$`Pr(>F)`[1])

rda.result <- rda(y ~ Phy + ST)
anova_result <- anova(rda.result, step = 1000, perm.max = 2000)
p_values <- c(p_values, anova_result$`Pr(>F)`[1])

rda.result <- rda(y ~ Cov + Env)
anova_result <- anova(rda.result, step = 1000, perm.max = 2000)
p_values <- c(p_values, anova_result$`Pr(>F)`[1])

rda.result <- rda(y ~ Cov + ST)
anova_result <- anova(rda.result, step = 1000, perm.max = 2000)
p_values <- c(p_values, anova_result$`Pr(>F)`[1])

rda.result <- rda(y ~ Env + ST)
anova_result <- anova(rda.result, step = 1000, perm.max = 2000)
p_values <- c(p_values, anova_result$`Pr(>F)`[1])

rda.result <- rda(y ~ Phy + Cov + Env)
anova_result <- anova(rda.result, step = 1000, perm.max = 2000)
p_values <- c(p_values, anova_result$`Pr(>F)`[1])

rda.result <- rda(y ~ Phy + Cov + ST)
anova_result <- anova(rda.result, step = 1000, perm.max = 2000)
p_values <- c(p_values, anova_result$`Pr(>F)`[1])

rda.result <- rda(y ~ Phy + Env + ST)
anova_result <- anova(rda.result, step = 1000, perm.max = 2000)
p_values <- c(p_values, anova_result$`Pr(>F)`[1])

rda.result <- rda(y ~ Cov + Env + ST)
anova_result <- anova(rda.result, step = 1000, perm.max = 2000)
p_values <- c(p_values, anova_result$`Pr(>F)`[1])

rda.result <- rda(y ~ Env + Cov + Phy + ST)
anova_result <- anova(rda.result, step = 1000, perm.max = 2000)
p_values <- c(p_values, anova_result$`Pr(>F)`[1])


### RDA condition

rda.result <- rda(y ~ Phy + Condition(Cov) + Condition(Env) + Condition(ST))
anova_result <- anova(rda.result, step = 1000, perm.max = 2000)
p_values <- c(p_values, anova_result$`Pr(>F)`[1])

rda.result <- rda(y ~ Cov + Condition(Phy) + Condition(Env) + Condition(ST))
anova_result <- anova(rda.result, step = 1000, perm.max = 2000)
p_values <- c(p_values, anova_result$`Pr(>F)`[1])

rda.result <- rda(y ~ Env + Condition(Phy) + Condition(Cov) + Condition(ST))
anova_result <- anova(rda.result, step = 1000, perm.max = 2000)
p_values <- c(p_values, anova_result$`Pr(>F)`[1])
plot(rda.result)

rda.result <- rda(y ~ ST + Condition(Phy) + Condition(Cov) + Condition(Env))
anova_result <- anova(rda.result, step = 1000, perm.max = 2000)
p_values <- c(p_values, anova_result$`Pr(>F)`[1])

results_df <- data.frame(
  p_value = p_values
)

print(results_df)

parvar$p_values <- results_df

fract <- parvar$part$fract
indfract <- parvar$part$indfract

results_fract <- data.frame(
  Partition = rownames(fract),
  Df = fract$Df,
  R_square = fract$R.square,
  Adj_R_square = fract$Adj.R.square,
  Testable = fract$Testable
)

results_indfract <- data.frame(
  Partition = rownames(indfract),
  Df = indfract$Df,
  R_square = indfract$R.square,
  Adj_R_square = indfract$Adj.R.square,
  Testable = indfract$Testable
)

final_results <- rbind(results_fract, results_indfract)

mapping <- c("X1" = "Phylogeny", 
             "X2" = "Size & Sex", 
             "X3" = "Environment", 
             "X4" = "Spatial Structure")

final_results$Partition <- gsub("X1", mapping["X1"], final_results$Partition)
final_results$Partition <- gsub("X2", mapping["X2"], final_results$Partition)
final_results$Partition <- gsub("X3", mapping["X3"], final_results$Partition)
final_results$Partition <- gsub("X4", mapping["X4"], final_results$Partition)
print(final_results)

p_values_complete <- c(p_values, rep(NA, nrow(final_results) - length(p_values)))
final_results$p_value <- p_values_complete
print(final_results)

write.csv(final_results, "parvar_correct.csv", row.names = FALSE)


# =============== PGLS best models for CS ================== #

phy_models <- c( "BM", "Pagel", "Martins", "Grafen")
size_pcnm

results_df <- data.frame()

for (bio in colnames(bioclim)) {
  for (phyM in phy_models) {
    
    formula_base <- as.formula(paste("size ~ sex + pcnm[,1] + pcnm[,5] +", bio))
    
    cor_struct <- switch(
      phyM,
      "BM"      = corBrownian(1, phy = tree),
      "Pagel"   = corPagel(1, phy = tree, fixed = FALSE),
      "Martins" = corMartins(1, phy = tree, fixed = FALSE),
      "Grafen"  = corGrafen(1, phy = tree, fixed = FALSE)
    )
    
    fit <- try(
      gls(
        formula_base,
        data = as.data.frame(cbind(size_sync, bioclim, pcnm)),
        correlation = cor_struct,
        method = "REML",
        control = glsControl(opt = "optim", msMaxIter = 200)
      ),
      silent = TRUE
    )
    
    # skip if model failed
    if (inherits(fit, "try-error")) next
    if (is.na(logLik(fit)[1])) next
    
    # --- EXTRAI PARÂMETRO DA ESTRUTURA AJUSTADA ---
    param_val <- NA
    cor_fitted <- try(fit$modelStruct$corStruct, silent = TRUE)
    if (!inherits(cor_fitted, "try-error") && !is.null(cor_fitted)) {
      # coef() sobre o corStruct ajustado devolve o(s) parâmetro(s)
      cor_coef_try <- try(coef(cor_fitted, unconstrained = FALSE), silent = TRUE)
      if (!inherits(cor_coef_try, "try-error") && length(cor_coef_try) > 0) {
        # pega o primeiro valor (p.ex. lambda, kappa, alpha, rho...)
        param_val <- as.numeric(cor_coef_try[1])
      } else {
        # corBrownian normalmente não tem parâmetro estimável (interpretação: 1)
        if (inherits(cor_fitted, "corBrownian")) param_val <- 1
      }
    }
    
    # Para Pagel, garante que o lambda esteja em [0,1]
    if (phyM == "Pagel" && (!is.finite(param_val) || param_val < 0 || param_val > 1)) next
    
    coefs <- summary(fit)$tTable
    if (!(bio %in% rownames(coefs))) next
    
    coef_row <- coefs[bio, ]
    
    results_df <- rbind(
      results_df,
      data.frame(
        Variable = bio,
        PhyModel = phyM,
        Estimate = coef_row[1],
        StdError = coef_row[2],
        tValue = coef_row[3],
        pValue = coef_row[4],
        AIC = AIC(fit),
        logLik = as.numeric(logLik(fit)),
        Parameter = param_val,
        converged = TRUE,
        stringsAsFactors = FALSE
      )
    )
  }
}

rownames(results_df) <- NULL
write.csv(results_df, "gls_bio_results.csv", row.names = FALSE)
print(head(results_df, 20))

best_models <- results_df %>%
  group_by(PhyModel) %>%
  filter(AIC == min(AIC, na.rm = TRUE)) %>%
  arrange(PhyModel)

write.csv(best_models, "best_models_by_phy.csv", row.names = FALSE)
print(best_models)

which.min(best_models$AIC)

# ---- Scale BIO9 ----
names(data_sync)
dat <- cbind(data_sync[,c(1:34,37,38)], pcnm[,size_pcnm$variables])
dat$BIO9_scaled <- as.numeric((scale(data_sync$BIO9)))
dat$sex <- factor(data_sync$sex)
dat$Cluster <- as.factor(data_sync$Cluster)
dat$sp <- as.factor(data_sync$sp)

# ---- Prediction sequence ----
bio_seq <- seq(
  min(dat$BIO9_scaled, na.rm = TRUE),
  max(dat$BIO9_scaled, na.rm = TRUE),
  length.out = 200
)

# ---- Mean covariates ----
pcnm1_m <- mean(dat$PCNM1, na.rm = TRUE)
pcnm5_m  <- mean(dat$PCNM5,  na.rm = TRUE)

# ---- Reference sex (controle) ----
sex_ref <- levels(dat$sex)[1]

newdat <- data.frame(
  BIO9_scaled = bio_seq,
  sex = factor(rep(sex_ref, length(bio_seq)), levels = levels(dat$sex)),
  PCNM1 = pcnm1_m,
  PCNM5  = pcnm5_m
)

# ---- Models ----
lm_fit <- lm(size ~ sex + PCNM1 + PCNM5 + BIO9_scaled, data = dat)

cor_struct <- corGrafen(1, phy = tree)
pgls_fit <- gls(
  size ~ sex + PCNM1 + PCNM5 + BIO9_scaled,
  data = dat,
  correlation = cor_struct,
  method = "REML"
)

# ---- Predictions ----
newdat$lm_pred   <- predict(lm_fit, newdata = newdat)
newdat$pgls_pred <- predict(pgls_fit, newdata = newdat)

# ---- Cluster colors ----
cluster_colors <- c(
  "Hot & Thermally Aseasonal"       = "#FDE725FF",
  "Warm & Moderately Seasonal"      = "#35B779FF",
  "Thermally Stable & Dry-Seasonal" = "#31688EFF",
  "Elevated & Thermally Seasonal"   = "#440154FF"
)
cluster_colors <- cluster_colors[levels(as.factor(data_sync$Cluster))]

# ---- Unicode species symbols ----
symbols <- c(
  Sapajus_apella        = "\u25A0", # square
  Sapajus_cay           = "\u25BC", # triangle down
  Sapajus_libidinosus   = "\u25B2", # triangle up
  Sapajus_nigritus      = "\u25CF", # circle
  Sapajus_robustus      = "\u2666", # diamond
  Sapajus_xanthosternos = "\u2605"  # star
)

symbols <- symbols[levels(dat$sp)]

# ---- Plot ----
ggplot() +
  # PGLS curve
  geom_line(
    data = newdat,
    aes(x = BIO9_scaled, y = pgls_pred),
    color = "red", size = 1.4
  ) +
  # LM curve
  geom_line(
    data = newdat,
    aes(x = BIO9_scaled, y = lm_pred),
    color = "blue", linetype = "dashed", size = 1.1
  ) +
  # Points using unicode shapes
  geom_text(
    data = dat,
    aes(
      x = BIO9_scaled,
      y = size,
      color = Cluster,
      label = symbols[sp]
    ),
    size = 8, alpha = 0.85
  ) +
  scale_color_manual(values = cluster_colors, name = "Climate Cluster") +
  labs(
    title = "Size vs BIO9 (scaled)",
    subtitle = "LM (dashed) & PGLS (solid)",
    x = "BIO9 (Mean Temperature of Driest Quarter)",
    y = "Size"
  ) +
  theme_minimal(base_size = 14)
