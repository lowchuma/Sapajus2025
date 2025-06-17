Lima <- readxl::read_xlsx("Lima et al. (2018)/AIC - Lima.xlsx")
head(Lima)

Upham <- readxl::read_xlsx("Upham et al. (2019)/AIC - Upham.xlsx")

Wright <- readxl::read_xlsx("Wright et al. (2015)/AIC - Wright.xlsx")
head(Wright)

L <- Lima$AIC
W <- Wright$AIC
U <- Upham$AIC

# --- 1. Instalar e carregar pacotes necessários ---
# Se você não tiver os pacotes, instale-os primeiro (remova o '#' da linha abaixo)
# install.packages(c("readxl", "dplyr", "car"))

library(readxl)
library(dplyr)
library(car)

# --- 2. Carregar seus dados (como no seu script) ---
# O código abaixo assume que os arquivos estão no seu diretório de trabalho.
# Lima <- readxl::read_xlsx("Lima et al. (2018)/AIC - Lima.xlsx")
# Upham <- readxl::read_xlsx("Upham et al. (2019)/AIC - Upham.xlsx")
# Wright <- readxl::read_xlsx("Wright et al. (2015)/AIC - Wright.xlsx")
#
# L <- Lima$AIC
# W <- Wright$AIC
# U <- Upham$AIC


# --- 3. Preparar o data frame para a análise ---
# A ANOVA precisa de um data frame com os valores em uma coluna e os grupos em outra.
dados_aov <- data.frame(
  AIC = c(L, W, U),
  Estudo = factor(rep(c("Lima", "Wright", "Upham"), times = c(length(L), length(W), length(U))))
)

# Visualizar o data frame criado
print(head(dados_aov))
print(summary(dados_aov))


# --- 4. Verificar a premissa de homogeneidade das variâncias ---
# A ANOVA assume que a variância é igual entre os grupos. Usamos o Teste de Levene.
# Hipótese nula (H0): As variâncias são iguais.
# Se p > 0.05, assumimos que as variâncias são homogêneas.
teste_levene <- leveneTest(AIC ~ Estudo, data = dados_aov)
cat("--- Teste de Levene para Homogeneidade das Variâncias ---\n")
print(teste_levene)


# --- 5. Realizar a Análise de Variância (ANOVA) ---
# A ANOVA testa se há uma diferença estatisticamente significativa entre as médias dos grupos.
# Hipótese nula (H0): As médias dos grupos são iguais.
# Se p < 0.05, rejeitamos H0 e concluímos que pelo menos um grupo é diferente.
modelo_aov <- aov(AIC ~ Estudo, data = dados_aov)
cat("\n\n--- Resultados da ANOVA ---\n")
print(summary(modelo_aov))


# --- 6. Teste Post-Hoc de Tukey (Tukey's HSD) ---
# Se a ANOVA for significativa (p < 0.05), este teste compara cada par de grupos (Lima vs. Wright, etc.)
# para ver quais são significativamente diferentes um do outro.
# Se o "p adj" (p-valor ajustado) for < 0.05, a diferença entre aquele par é significativa.
teste_tukey <- TukeyHSD(modelo_aov)
cat("\n\n--- Teste Post-Hoc de Tukey HSD ---\n")
print(teste_tukey)


# --- 7. Visualização com Boxplot ---
boxplot(AIC ~ Estudo,
        data = dados_aov,
        main = "Comparação dos Valores de AIC por Estudo",
        xlab = "Estudos",
        ylab = "Valor de AIC",
        col = c("#D7E8D4", "#BDE0FE", "#FFDAB9"),
        border = "black"
)

# Carregue o pacote se ainda não o fez
library(dplyr)

# Carregue suas planilhas (ajuste o caminho se necessário)
lima_shape <- readxl::read_xlsx("Lima et al. (2018)/AIC - Lima.xlsx", sheet = 1)
wright_shape <- readxl::read_xlsx("Wright et al. (2015)/AIC - Wright.xlsx", sheet = 1)
upham_shape <- readxl::read_xlsx("Upham et al. (2019)/AIC - Upham.xlsx", sheet = 1)

# Carregue suas planilhas (ajuste o caminho se necessário)
lima_size <- readxl::read_xlsx("Lima et al. (2018)/AIC - Lima.xlsx", sheet = 2)
wright_size <- readxl::read_xlsx("Wright et al. (2015)/AIC - Wright.xlsx", sheet = 2)
upham_size <- readxl::read_xlsx("Upham et al. (2019)/AIC - Upham.xlsx", sheet = 2)

# Melhor modelo SHAPE de Lima
lima_shape[which.min(lima_shape$AIC), ]

# Melhor modelo SHAPE de Wright
wright_shape[which.min(wright_shape$AIC), ]

# Melhor modelo SIZE de Lima
lima_size[which.min(lima_size$AIC), ]

# Melhor modelo SIZE de Wright
wright_size[which.min(wright_size$AIC), ]

# Lista completa de AICs para Wright "Shape"
print(wright_shape$AIC)

# Lista completa de AICs para Wright "Size"
print(wright_size$AIC)
