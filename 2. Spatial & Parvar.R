
### Spatial dependence - Test for locality bias
#install.packages("spdep")
library(spdep)
library(dplyr)
#citation("dplyr")
#citation("spdep")
#citation("mpmcorrelogram")

locations <- read.csv("Planilhas/avglocsex.csv", sep = ";") ## samples averaged by locations, no pseudoreplication
Loc <- readland.tps("Datasets/tps/avglocsex.tps")
L.gpa <- gpagen(Loc)

coords <- cbind(locations$lat,locations$long)
# Required packages
#install.packages(c("reshape2", "geosphere", "mpmcorrelogram"))
library(reshape2)
library(geosphere)
library(mpmcorrelogram)

# --- 1. Prepare data ---

# Extract and log-transform centroid size
CS <- as.numeric(locations$CS)
logCS <- log(CS)

# Euclidean distance matrix for log(CS)
sizeresdist <- as.matrix(dist(logCS))

# Geographic coordinates matrix
coords <- as.matrix(locations[, c("long", "lat")])

# Compute geographic distance matrix in kilometers
geo_km <- distm(coords, fun = distHaversine) / 1000  # Convert meters to km
geographicaldistance <- as.dist(geo_km)

# --- 2. Correlation: Geographic distance vs. Size distance ---

# Melt to vector form
geodistx <- melt(as.matrix(geographicaldistance))
sizeresdisty <- melt(sizeresdist)

# Filter valid data
valid_data <- geodistx[, 3] != 0 & sizeresdisty[, 3] != 0
x2 <- geodistx[, 3][valid_data]
y2 <- sizeresdisty[, 3][valid_data]

# Plot
plot(x2, y2, pch = 19, col = "blue",
     xlab = "Geographic distance (km)",
     ylab = "Size distance (logCS)",
     main = "Correlation between geographic distance and size")
abline(lm(y2 ~ x2), col = "black", lty = 2)

# Pearson correlation test
cor.test(x2, y2, method = "pearson")

# --- 3. Correlogram: Size vs Geographic distance ---
mpm <- mpmcorrelogram(
  dist(logCS), geographicaldistance,
  method = "pearson",
  alfa = 0.05, permutations = 999,
  simil = FALSE, plot = TRUE, print = TRUE
)
title(main = "Correlogram: Size (logCS) vs Geographic distance (km)")

# --- 4. Correlation: Geographic distance vs. Shape distance ---

# Procrustes shape distance
shapeDist <- dist(two.d.array(L.gpa$coords))
shapeDist_mat <- as.matrix(shapeDist)
shapedisty <- melt(shapeDist_mat)

# Filter valid data
valid_data <- geodistx[, 3] != 0 & shapedisty[, 3] != 0
x2 <- geodistx[, 3][valid_data]
y2 <- shapedisty[, 3][valid_data]

# Plot
plot(x2, y2, pch = 19, col = "red",
     xlab = "Geographic distance (km)",
     ylab = "Shape distance (Procrustes)",
     main = "Correlation between geographic distance and shape")
abline(lm(y2 ~ x2), col = "black", lty = 2)

# Pearson correlation test
cor.test(x2, y2, method = "pearson")

# --- 5. Correlogram: Shape vs Geographic distance ---
mpm <- mpmcorrelogram(
  shapeDist, geographicaldistance,
  method = "pearson",
  alfa = 0.05, permutations = 999,
  simil = FALSE, plot = TRUE, print = TRUE
)
title(main = "Correlogram: Shape vs Geographic distance (km)")

# Instale os pacotes se necessário
# install.packages("ggplot2")
# install.packages("dplyr")

library(ggplot2)
library(dplyr)

# ---- Dados do Correlograma: Size vs. Geographic Distance ----
correlog_size <- data.frame(
  Class = 1:13,
  Distance_Range = c("0–272", "273–546", "546–819", "819–1091", "1091–1364",
                     "1364–1637", "1637–1910", "1910–2183", "2183–2456",
                     "2456–2729", "2729–3002", "3002–3274", "3274–3548"),
  rM = c(0.0146, 0.0362, 0.0339, 0.0366, -0.0479, 0.0143, -0.0099, 0.0356,
         -0.0368, -0.0516, -0.0046, -0.0543, -0.0337),
  p = c(0.256, 0.061, 0.119, 0.103, 0.034, 0.311, 0.343, 0.078,
        0.089, 0.027, 0.394, 0.041, 0.173)
)

correlog_size$Significant <- correlog_size$p < 0.05

# ---- Dados do Correlograma: Shape vs. Geographic Distance ----
correlog_shape <- data.frame(
  Class = 1:13,
  Distance_Range = c("0–272", "273–546", "546–819", "819–1091", "1091–1364",
                     "1364–1637", "1637–1910", "1910–2183", "2183–2456",
                     "2456–2729", "2729–3002", "3002–3274", "3274–3548"),
  rM = c(0.0281, 0.0360, 0.0464, 0.0049, -0.0197, -0.0265, -0.0125, 0.0037,
         -0.0276, 0.0068, 0.0195, -0.0363, -0.0531),
  p = c(0.112, 0.063, 0.043, 0.431, 0.213, 0.155, 0.317, 0.452,
        0.157, 0.419, 0.263, 0.120, 0.097)
)

correlog_shape$Significant <- correlog_shape$p < 0.05

# ---- Função para plotar correlograma ----
plot_correlogram <- function(data, title_text) {
  ggplot(data, aes(x = factor(Class), y = rM, fill = Significant)) +
    geom_bar(stat = "identity", color = "black") +
    scale_fill_manual(values = c("TRUE" = "red", "FALSE" = "gray")) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    scale_x_discrete(labels = data$Distance_Range) +
    labs(
      x = "Distance Class (km)",
      y = "Moran's r",
      fill = "Significant (p < 0.05)",
      title = title_text
    ) +
    theme_minimal(base_size = 14) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
}

# ---- Gerar os plots ----
plot_correlogram(correlog_size, "Spatial Correlogram: Size vs. Geographic Distance")
plot_correlogram(correlog_shape, "Spatial Correlogram: Shape vs. Geographic Distance")


pcnm<-prcomp(geographicaldistance) # Já cria os elementos espacialmente independentes
summary(pcnm)


### Testes t ###

t_test_result <- t.test(LogCS ~ biome, data = locations)
print(t_test_result) ## p < 0,001
t_test_result <- t.test(LogCS ~ sex, data = locations)
print(t_test_result) ## p < 0,001


#####variation partitioning#########
coordsavg <- data.frame(lat = locations$lat, lng = locations$long)
geographicaldistance <- as.matrix(dist(coordsavg))

pcnmS<-prcomp(geographicaldistance)
pcaS<-gm.prcomp(L.gpa$coords)
summary(pcnmS)
require(packfor)
forward.sel(pcaS$x, pcnmS$x)
forward.sel(pcaS$x, pls$XScores)

tree<-read.tree("trees/AVGLOCSEX_Lima.tre")
rownames(data) <- data$sp_ives
plot(tree, cex = 0.5)
axisPhylo()

cov<-vcv.phylo(tree)
phylo<-prcomp(cov)
summary(phylo)
phy<-as.matrix(phylo$x[,1:3])

require(vegan)
parvar<-varpart(pcaS$x,~phy, ~pcnmS$x[,c(51,80,4,54,79,42,8,21)],~log(L.gpa$Csize),~pls$XScores)
plot(parvar,cex=1.2,Xnames=c("Phylogeny","Spatial Structure","Size","Environment"),bg=c("yellow", "deeppink","navy","green3"))

parvar
y<-pcaS$x
size<-as.matrix(log(L.gpa$Csize))
pcnm<-model.matrix(~pcnmS$x[,51]+pcnmS$x[,80]+pcnmS$x[,4]+pcnmS$x[,54]+pcnmS$x[,79]+pcnmS$x[,42]+pcnmS$x[,8]+pcnmS$x[,21])
env<-model.matrix(~pls$XScores)
rda.result <- rda(y ~ size + pcnm + env)
anova(rda.result, step = 10000, perm.max = 10000)

y<-pcaS$x
sf<-as.matrix(pcnmS$x[,c(51,80,4,54,79,42,8,21)])
env <-pls$XScores

### for all parvar's variables #####

p_values <- c()

rda.result <- rda(y ~ phy)
anova_result <- anova(rda.result, step = 1000, perm.max = 2000)
p_values <- c(p_values, anova_result$`Pr(>F)`[1])

rda.result_alo <- rda(y ~ size)
anova_result <- anova(rda.result_alo, step = 1000, perm.max = 2000)
p_values <- c(p_values, anova_result$`Pr(>F)`[1])

rda.result <- rda(y ~ env)
anova_result <- anova(rda.result, step = 1000, perm.max = 2000)
p_values <- c(p_values, anova_result$`Pr(>F)`[1])

rda.result <- rda(y ~ sf)
anova_result <- anova(rda.result, step = 1000, perm.max = 2000)
p_values <- c(p_values, anova_result$`Pr(>F)`[1])

rda.result <- rda(y ~ phy + size)
anova_result <- anova(rda.result, step = 1000, perm.max = 2000)
p_values <- c(p_values, anova_result$`Pr(>F)`[1])

rda.result <- rda(y ~ phy + env)
anova_result <- anova(rda.result, step = 1000, perm.max = 2000)
p_values <- c(p_values, anova_result$`Pr(>F)`[1])

rda.result <- rda(y ~ phy + sf)
anova_result <- anova(rda.result, step = 1000, perm.max = 2000)
p_values <- c(p_values, anova_result$`Pr(>F)`[1])

rda.result <- rda(y ~ size + env)
anova_result <- anova(rda.result, step = 1000, perm.max = 2000)
p_values <- c(p_values, anova_result$`Pr(>F)`[1])

rda.result <- rda(y ~ size + sf)
anova_result <- anova(rda.result, step = 1000, perm.max = 2000)
p_values <- c(p_values, anova_result$`Pr(>F)`[1])

rda.result <- rda(y ~ env + sf)
anova_result <- anova(rda.result, step = 1000, perm.max = 2000)
p_values <- c(p_values, anova_result$`Pr(>F)`[1])

rda.result <- rda(y ~ phy + size + env)
anova_result <- anova(rda.result, step = 1000, perm.max = 2000)
p_values <- c(p_values, anova_result$`Pr(>F)`[1])

rda.result <- rda(y ~ phy + size + sf)
anova_result <- anova(rda.result, step = 1000, perm.max = 2000)
p_values <- c(p_values, anova_result$`Pr(>F)`[1])

rda.result <- rda(y ~ phy + env + sf)
anova_result <- anova(rda.result, step = 1000, perm.max = 2000)
p_values <- c(p_values, anova_result$`Pr(>F)`[1])

rda.result <- rda(y ~ size + env + sf)
anova_result <- anova(rda.result, step = 1000, perm.max = 2000)
p_values <- c(p_values, anova_result$`Pr(>F)`[1])

rda.result <- rda(y ~ env + size + phy + sf)
anova_result <- anova(rda.result, step = 1000, perm.max = 2000)
p_values <- c(p_values, anova_result$`Pr(>F)`[1])


### RDA condicionada

rda.result <- rda(y ~ phy + Condition(size) + Condition(env) + Condition(sf))
anova_result <- anova(rda.result, step = 1000, perm.max = 2000)
p_values <- c(p_values, anova_result$`Pr(>F)`[1])

rda.result <- rda(y ~ size + Condition(phy) + Condition(env) + Condition(sf))
anova_result <- anova(rda.result, step = 1000, perm.max = 2000)
p_values <- c(p_values, anova_result$`Pr(>F)`[1])

rda.result <- rda(y ~ env + Condition(phy) + Condition(size) + Condition(sf))
anova_result <- anova(rda.result, step = 1000, perm.max = 2000)
p_values <- c(p_values, anova_result$`Pr(>F)`[1])
plot(rda.result)

rda.result <- rda(y ~ sf + Condition(phy) + Condition(size) + Condition(env))
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
             "X2" = "Size", 
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

#write.csv(final_results, "parvar_wright.csv", row.names = FALSE)

# Ajuste de modelos de RDA (Redundancy Analysis) para partição da variação

y <- pcaS$x
size <- log(L.gpa$Csize)
pcnm2 <- model.matrix(y ~ pcnm)
env2<-model.matrix(~pls$XScores[,1]+pls$XScores[,2]+pls$XScores[,3])
#?model.matrix
rda.result <- rda(y ~ pcnm2 + env2)
anova(rda.result, step = 10000, perm.max = 10000)
