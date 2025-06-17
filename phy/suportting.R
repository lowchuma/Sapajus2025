sapajus <- dados
data <- fatores
link <- read.table("TPS/links.txt")
bioclim <- data
# simple shape x bio10, bio13, bio15, bio16
# alometry x bio6, bio8, bio10, bio13, bio14, bio15, bio17, bio19
# pgls.SEy x bio4
# pncm3, 4 e 7
# bio4:PCNM7 = 0.97293    0.4790  2.031230  0.0475
# DREDGE: PCNM 20 (?)

plot(data$bio4, data$LogCS)
abline(fit, col = "red")
points(data$bio4, data$LogCS, pch = data$pch, col = "blue")
text(data$bio4, data$LogCS, data$label, col = "black")
title(main = "Relação entre bio4 e LogCS")
xlab("bio4")
ylab("LogCS")
legend("topright", c("Linha de regressão"), lty = 1, col = "red")
xlim(c(min(data$bio4), max(data$bio4)))
ylim(c(min(data$LogCS), max(data$LogCS)))
grid(TRUE)
dev.off()
coefs <- fit$coefficients
coefs.se <- fit$coefficients.standard.error
pvals <- fit$p.values
plot(coefs)

# Criando o data frame
df_pch<- data.frame(biome = biome)

# Vetores de referência para os sexos e os valores pch correspondentes
biomas <- c("Savanna", "Forest")
pch_values <- c(18, 15)

# Adicionando a nova coluna 'pch' usando match
df_pch$pch <- pch_values[match(df_pch$biome, biomas)]

# Visualizando o data frame atualizado
print(df_pch)

# Plotar o modelo com ggplot2 e diferenciar biomas por cor
ggplot(data, aes(x = bio4, y = LogCS, color = sex)) +
  geom_point(pch = as.numeric(df_pch$pch)) +
  geom_smooth(method = "lm", se = TRUE) +
  labs(title = "Phylogeny-Generalized Least Squares", x = "Climate Variability", y = "LogCS") +  scale_color_manual(breaks = unique(data$sex), values = c("darkblue", "red")) +  # Exemplo de paleta de cores
  theme_minimal()


# interaction models
fit <- procD.lm(gpa$coords ~ data$LogCS * data$biome * data$bio1, iter=999)
summary(fit)

fit <- procD.lm(gpa$coords ~ data$LogCS * data$biome * data$bio4, iter=999)
summary(fit)

fit <- procD.lm(gpa$coords ~ data$LogCS * data$biome * data$bio2, iter=999)
summary(fit)

fit <- procD.lm(gpa$coords ~ data$LogCS * data$biome * data$bio6, iter=999)
summary(fit)

fit <- procD.lm(gpa$coords ~ data$LogCS * data$biome * data$bio7, iter=999)
summary(fit)

fit <- procD.lm(gpa$coords ~ data$LogCS * data$biome * data$bio8, iter=999) # *
summary(fit)

fit <- procD.lm(gpa$coords ~ data$LogCS * data$biome * data$bio9, iter=999)
summary(fit)

fit <- procD.lm(gpa$coords ~ data$LogCS * data$biome * data$bio10, iter=999) # *
summary(fit)

fit <- procD.lm(gpa$coords ~ data$LogCS * data$biome * data$bio11, iter=999) 
summary(fit)

fit <- procD.lm(gpa$coords ~ data$LogCS * data$biome * data$bio12, iter=999)
summary(fit)

fit <- procD.lm(gpa$coords ~ data$LogCS * data$biome * data$bio13, iter=999)
summary(fit)

fit <- procD.lm(gpa$coords ~ data$LogCS * data$biome * data$bio14, iter=999)
summary(fit)

fit <- procD.lm(gpa$coords ~ data$LogCS * data$biome * data$bio15, iter=999)
summary(fit)

fit <- procD.lm(gpa$coords ~ data$LogCS * data$biome * data$bio16, iter=999)
summary(fit)

fit <- procD.lm(gpa$coords ~ data$LogCS * data$biome * data$bio17, iter=999)
summary(fit)

fit <- procD.lm(gpa$coords ~ data$LogCS * data$biome * data$bio18, iter=999)
summary(fit)

fit <- procD.lm(gpa$coords ~ data$LogCS * data$biome * data$bio19, iter=999)
summary(fit)

pcnm3<-pcnm[,3]
pcnm4<-pcnm[,4]
pcnm7<-pcnm[,7]

full.model <- lm(data$LogCS~pcnm3+pcnm4+pcnm7)
summary(full.model)

fit<-pgls.SEy(LogCS~bio4+pcnm3+pcnm4+pcnm7,data=data,tree=tree)
summary(fit)
