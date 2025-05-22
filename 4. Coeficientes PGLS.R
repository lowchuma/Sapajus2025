setwd("C:/Users/Lourenço/OneDrive/UFSM/LabMastozoo/Sapajus/Data/PGLS.SEy")


#BROWNIAN

# Pacotes necessários
library(tidyverse)  # Para manipulação de dados

# Lista de variáveis (bio1 até bio19)
bio_vars <- paste0("bio", 1:19)

#TAMANHO

run_pgls <- function(var, data, tree) {
  # Cria a fórmula dinamicamente
  formula <- as.formula(paste("CS ~", var))
  
  # Tenta rodar o modelo
  fit <- try(pgls.SEy(formula, data = data, tree = tree), silent = TRUE)
  
  # Verifica se o modelo foi rodado corretamente
  if (inherits(fit, "try-error")) {
    return(data.frame(Variable = var, Beta = NA, T_value = NA, P_value = NA, AIC = NA, BIC = NA, Sigma = NA))
  }
  
  # Extrai o sumário do modelo
  fited <- summary(fit)
  
  # Extrai os coeficientes da tabela tTable
  coef_summary <- fited$tTable
  
  # Inicializa as variáveis
  beta <- NA
  t_value <- NA
  p_value <- NA
  
  # Verifica se a tabela está estruturada corretamente
  if (!is.null(coef_summary) && nrow(coef_summary) >= 2) {
    beta <- coef_summary[2, "Value"]
    t_value <- coef_summary[2, "t-value"]
    p_value <- coef_summary[2, "p-value"]
  }
  
  # Extrai AIC, BIC e sigma
  aic_value <- fited$AIC
  bic_value <- fited$BIC
  sigma_value <- fited$sigma
  
  # Retorna um data frame com os resultados
  return(data.frame(Variable = var, Beta = beta, T_value = t_value, P_value = p_value, AIC = aic_value, BIC = bic_value, Sigma = sigma_value))
}

# Aplicar a função para cada variável bio
results <- map_df(bio_vars, ~run_pgls(.x, data, tree))

# Também pode salvar em CSV se preferir
#write.csv(results, "wright_simple_pgls_B.csv", row.names = FALSE)


### FORMA 

# Função para rodar o PGLS e extrair os resultados
run_pgls <- function(var, data, tree) {
  # Cria a fórmula dinamicamente
  formula <- as.formula(paste("y ~", var))
  
  # Tenta rodar o modelo
  fit <- try(pgls.SEy(formula, data = data, tree = tree), silent = TRUE)
  
  # Verifica se o modelo foi rodado corretamente
  if (inherits(fit, "try-error")) {
    return(data.frame(Variable = var, Beta = NA, T_value = NA, P_value = NA, AIC = NA, BIC = NA, Sigma = NA))
  }
  
  # Extrai o sumário do modelo
  fited <- summary(fit)
  
  # Extrai os coeficientes da tabela tTable
  coef_summary <- fited$tTable
  
  # Inicializa as variáveis
  beta <- NA
  t_value <- NA
  p_value <- NA
  
  # Verifica se a tabela está estruturada corretamente
  if (!is.null(coef_summary) && nrow(coef_summary) >= 2) {
    beta <- coef_summary[2, "Value"]
    t_value <- coef_summary[2, "t-value"]
    p_value <- coef_summary[2, "p-value"]
  }
  
  # Extrai AIC, BIC e sigma
  aic_value <- fited$AIC
  bic_value <- fited$BIC
  sigma_value <- fited$sigma
  
  # Retorna um data frame com os resultados
  return(data.frame(Variable = var, Beta = beta, T_value = t_value, P_value = p_value, AIC = aic_value, BIC = bic_value, Sigma = sigma_value))
}

# Aplicar a função para cada variável bio
resultsY <- map_df(bio_vars, ~run_pgls(.x, data, tree))

# Também pode salvar em CSV se preferir
#write.csv(resultsY, "wrightY_simple_pgls_B.csv", row.names = FALSE)

library(ggplot2)
library(ggrepel)  # Para rótulos que não se sobrepõem

par(mfrow=c(1,2))

# Gráfico de coeficientes com personalizações
results %>%
  mutate(Variable = factor(Variable, levels = Variable[order(Beta)])) %>%
  ggplot(aes(x = Variable, y = Beta)) +
  geom_col(aes(fill = P_value < 0.05), color = "black") +  # Preenchimento e contorno
  scale_fill_manual(values = c("TRUE" = "dodgerblue", "FALSE" = "salmon"), 
                    labels = c("TRUE" = "Significativa (p < 0.05)", "FALSE" = "Não Significativa")) +
  geom_errorbar(aes(ymin = Beta - 1.96 * Sigma, ymax = Beta + 1.96 * Sigma), width = 0.2, color = "darkgrey") +
  coord_flip() +
  labs(title = "Coeficientes Beta dos Modelos PGLS (Brownian)",
       subtitle = "Análise das variáveis preditoras em relação ao tamanho",
       x = "Variáveis",
       y = "Coeficiente Beta",
       fill = "Significância") +
  theme_minimal(base_size = 15) +  # Aumenta o tamanho base da fonte
  theme(legend.position = "top") +  # Coloca a legenda no topo
  geom_text_repel(aes(label = ifelse(P_value < 0.05, round(Beta, 2), "")),  # Rótulos para coeficientes significativos
                  box.padding = 0.5, point.padding = 0.5, show.legend = FALSE) +
  theme(plot.title = element_text(hjust = 0.5, size = 20, face = "bold"),  # Alinhamento e estilo do título
        plot.subtitle = element_text(hjust = 0.5, size = 16, face = "italic"))  # Estilo do subtítulo


# Gráfico de coeficientes com personalizações
resultsY %>%
  mutate(Variable = factor(Variable, levels = Variable[order(Beta)])) %>%
  ggplot(aes(x = Variable, y = Beta)) +
  geom_col(aes(fill = P_value < 0.05), color = "black") +  # Preenchimento e contorno
  scale_fill_manual(values = c("TRUE" = "dodgerblue", "FALSE" = "salmon"), 
                    labels = c("TRUE" = "Significativa (p < 0.05)", "FALSE" = "Não Significativa")) +
  geom_errorbar(aes(ymin = Beta - 1.96 * Sigma, ymax = Beta + 1.96 * Sigma), width = 0.2, color = "darkgrey") +
  coord_flip() +
  labs(title = "Coeficientes Beta dos Modelos PGLS (Brownian)",
       subtitle = "Análise das variáveis preditoras em relação à forma",
       x = "Variáveis",
       y = "Coeficiente Beta",
       fill = "Significância") +
  theme_minimal(base_size = 15) +  # Aumenta o tamanho base da fonte
  theme(legend.position = "top") +  # Coloca a legenda no topo
  geom_text_repel(aes(label = ifelse(P_value < 0.05, round(Beta, 2), "")),  # Rótulos para coeficientes significativos
                  box.padding = 0.5, point.padding = 0.5, show.legend = FALSE) +
  theme(plot.title = element_text(hjust = 0.5, size = 20, face = "bold"),  # Alinhamento e estilo do título
        plot.subtitle = element_text(hjust = 0.5, size = 16, face = "italic"))  # Estilo do subtítulo

par(mfrow=c(1,1))


library(dplyr)
library(purrr)

# Lista de variáveis (bio1 até bio19)
bio_vars <- paste0("bio", 1:19)

## TAMANHO + PCNM

run_pgls <- function(var, data, tree) {
  # Cria a fórmula dinamicamente
  formula <- as.formula(paste("CS ~", var, "+ pcnm1 + pcnm5"))
  
  # Tenta rodar o modelo
  fit <- try(pgls.SEy(formula, data = data, tree = tree), silent = TRUE)
  
  # Verifica se o modelo foi rodado corretamente
  if (inherits(fit, "try-error")) {
    return(data.frame(Variable = var, Beta = NA, T_value = NA, P_value = NA, AIC = NA, BIC = NA, Sigma = NA))
  }
  
  # Extrai o sumário do modelo
  fited <- summary(fit)
  
  # Extrai os coeficientes da tabela tTable
  coef_summary <- fited$tTable
  
  # Inicializa as variáveis
  beta <- NA
  t_value <- NA
  p_value <- NA
  
  # Verifica se a tabela está estruturada corretamente
  if (!is.null(coef_summary) && nrow(coef_summary) >= 2) {
    beta <- coef_summary[2, "Value"]
    t_value <- coef_summary[2, "t-value"]
    p_value <- coef_summary[2, "p-value"]
  }
  
  # Extrai AIC, BIC e sigma
  aic_value <- fited$AIC
  bic_value <- fited$BIC
  sigma_value <- fited$sigma
  
  # Retorna um data frame com os resultados
  return(data.frame(Variable = var, Beta = beta, T_value = t_value, P_value = p_value, AIC = aic_value, BIC = bic_value, Sigma = sigma_value))
}

# Aplicar a função para cada variável bio
results <- map_df(bio_vars, ~run_pgls(.x, data, tree))

# Salvar em CSV
#write.csv(results, "wright_pcnm_pgls.csv", row.names = FALSE)

### FORMA + PCNM
run_pgls <- function(var, data, tree) {
  # Cria a fórmula dinamicamente
  formula <- as.formula(paste("CS ~", var, "+ pcnm1 + pcnm5"))
  
  # Tenta rodar o modelo
  fit <- try(pgls.SEy(formula, data = data, tree = tree), silent = TRUE)
  
  # Verifica se o modelo foi rodado corretamente
  if (inherits(fit, "try-error")) {
    return(data.frame(Variable = var, Beta = NA, T_value = NA, P_value = NA, AIC = NA, BIC = NA, Sigma = NA))
  }
  
  # Extrai o sumário do modelo
  fited <- summary(fit)
  
  # Extrai os coeficientes da tabela tTable
  coef_summary <- fited$tTable
  
  # Inicializa as variáveis
  beta <- NA
  t_value <- NA
  p_value <- NA
  
  # Verifica se a tabela está estruturada corretamente
  if (!is.null(coef_summary) && nrow(coef_summary) >= 2) {
    beta <- coef_summary[2, "Value"]
    t_value <- coef_summary[2, "t-value"]
    p_value <- coef_summary[2, "p-value"]
  }
  
  # Extrai AIC, BIC e sigma
  aic_value <- fited$AIC
  bic_value <- fited$BIC
  sigma_value <- fited$sigma
  
  # Retorna um data frame com os resultados
  return(data.frame(Variable = var, Beta = beta, T_value = t_value, P_value = p_value, AIC = aic_value, BIC = bic_value, Sigma = sigma_value))
}

# Aplicar a função para cada variável bio
results <- map_df(bio_vars, ~run_pgls(.x, data, tree))
#write.csv(resultsY, "wrightY_pcnm_pgls_B.csv", row.names = FALSE)

## FORMA + PCNM

run_pgls <- function(var, data, tree) {
  # Cria a fórmula dinamicamente
  formula <- as.formula(paste("CS ~", var, "+ pcnm1 + pcnm5"))
  
  # Tenta rodar o modelo
  fit <- try(pgls.SEy(formula, data = data, tree = tree), silent = TRUE)
  
  # Verifica se o modelo foi rodado corretamente
  if (inherits(fit, "try-error")) {
    return(data.frame(Variable = var, Beta = NA, T_value = NA, P_value = NA, AIC = NA, BIC = NA, Sigma = NA))
  }
  
  # Extrai o sumário do modelo
  fited <- summary(fit)
  
  # Extrai os coeficientes da tabela tTable
  coef_summary <- fited$tTable
  
  # Inicializa as variáveis
  beta <- NA
  t_value <- NA
  p_value <- NA
  
  # Verifica se a tabela está estruturada corretamente
  if (!is.null(coef_summary) && nrow(coef_summary) >= 2) {
    beta <- coef_summary[2, "Value"]
    t_value <- coef_summary[2, "t-value"]
    p_value <- coef_summary[2, "p-value"]
  }
  
  # Extrai AIC, BIC e sigma
  aic_value <- fited$AIC
  bic_value <- fited$BIC
  sigma_value <- fited$sigma
  
  # Retorna um data frame com os resultados
  return(data.frame(Variable = var, Beta = beta, T_value = t_value, P_value = p_value, AIC = aic_value, BIC = bic_value, Sigma = sigma_value))
}

# Aplicar a função para cada variável bio
results <- map_df(bio_vars, ~run_pgls(.x, data, tree))

# Salvar em CSV
#write.csv(results, "wright_pcnm_pgls.csv", row.names = FALSE)

### FORMA + PCNM
run_pgls <- function(var, data, tree) {
  # Cria a fórmula dinamicamente
  formula <- as.formula(paste("CS ~", var, "+ pcnm1 + pcnm5"))
  
  # Tenta rodar o modelo
  fit <- try(pgls.SEy(formula, data = data, tree = tree), silent = TRUE)
  
  # Verifica se o modelo foi rodado corretamente
  if (inherits(fit, "try-error")) {
    return(data.frame(Variable = var, Beta = NA, T_value = NA, P_value = NA, AIC = NA, BIC = NA, Sigma = NA))
  }
  
  # Extrai o sumário do modelo
  fited <- summary(fit)
  
  # Extrai os coeficientes da tabela tTable
  coef_summary <- fited$tTable
  
  # Inicializa as variáveis
  beta <- NA
  t_value <- NA
  p_value <- NA
  
  # Verifica se a tabela está estruturada corretamente
  if (!is.null(coef_summary) && nrow(coef_summary) >= 2) {
    beta <- coef_summary[2, "Value"]
    t_value <- coef_summary[2, "t-value"]
    p_value <- coef_summary[2, "p-value"]
  }
  
  # Extrai AIC, BIC e sigma
  aic_value <- fited$AIC
  bic_value <- fited$BIC
  sigma_value <- fited$sigma
  
  # Retorna um data frame com os resultados
  return(data.frame(Variable = var, Beta = beta, T_value = t_value, P_value = p_value, AIC = aic_value, BIC = bic_value, Sigma = sigma_value))
}

# Aplicar a função para cada variável bio
resultsY <- map_df(bio_vars, ~run_pgls(.x, data, tree))
#write.csv(resultsY, "wrightY_pcnm_pgls_B.csv", row.names = FALSE)

# Gráfico de coeficientes com personalizações
results %>%
  mutate(Variable = factor(Variable, levels = Variable[order(Beta)])) %>%
  ggplot(aes(x = Variable, y = Beta)) +
  geom_col(aes(fill = P_value < 0.05), color = "black") +  # Preenchimento e contorno
  scale_fill_manual(values = c("TRUE" = "dodgerblue", "FALSE" = "salmon"), 
                    labels = c("TRUE" = "Significativa (p < 0.05)", "FALSE" = "Não Significativa")) +
  geom_errorbar(aes(ymin = Beta - 1.96 * Sigma, ymax = Beta + 1.96 * Sigma), width = 0.2, color = "darkgrey") +
  coord_flip() +
  labs(title = "Coeficientes Beta dos Modelos PGLS + PCNM (Brownian)",
       subtitle = "Análise das variáveis preditoras em relação ao tamanho",
       x = "Variáveis",
       y = "Coeficiente Beta",
       fill = "Significância") +
  theme_minimal(base_size = 15) +  # Aumenta o tamanho base da fonte
  theme(legend.position = "top") +  # Coloca a legenda no topo
  geom_text_repel(aes(label = ifelse(P_value < 0.05, round(Beta, 2), "")),  # Rótulos para coeficientes significativos
                  box.padding = 0.5, point.padding = 0.5, show.legend = FALSE) +
  theme(plot.title = element_text(hjust = 0.5, size = 20, face = "bold"),  # Alinhamento e estilo do título
        plot.subtitle = element_text(hjust = 0.5, size = 16, face = "italic"))  # Estilo do subtítulo


# Gráfico de coeficientes com personalizações
resultsY %>%
  mutate(Variable = factor(Variable, levels = Variable[order(Beta)])) %>%
  ggplot(aes(x = Variable, y = Beta)) +
  geom_col(aes(fill = P_value < 0.05), color = "black") +  # Preenchimento e contorno
  scale_fill_manual(values = c("TRUE" = "dodgerblue", "FALSE" = "salmon"), 
                    labels = c("TRUE" = "Significativa (p < 0.05)", "FALSE" = "Não Significativa")) +
  geom_errorbar(aes(ymin = Beta - 1.96 * Sigma, ymax = Beta + 1.96 * Sigma), width = 0.2, color = "darkgrey") +
  coord_flip() +
  labs(title = "Coeficientes Beta dos Modelos PGLS + PCNM (Brownian)",
       subtitle = "Análise das variáveis preditoras em relação à forma",
       x = "Variáveis",
       y = "Coeficiente Beta",
       fill = "Significância") +
  theme_minimal(base_size = 15) +  # Aumenta o tamanho base da fonte
  theme(legend.position = "top") +  # Coloca a legenda no topo
  geom_text_repel(aes(label = ifelse(P_value < 0.05, round(Beta, 2), "")),  # Rótulos para coeficientes significativos
                  box.padding = 0.5, point.padding = 0.5, show.legend = FALSE) +
  theme(plot.title = element_text(hjust = 0.5, size = 20, face = "bold"),  # Alinhamento e estilo do título
        plot.subtitle = element_text(hjust = 0.5, size = 16, face = "italic"))  # Estilo do subtítulo

par(mfrow=c(1,1))

################################################ PAGEL's ####################################################

library(dplyr)
library(purrr)


bio_vars <- paste0("bio", 1:19)


### TAMANHO
run_pgls <- function(var, data, tree) {
  # Cria a fórmula dinamicamente
  formula <- as.formula(paste("CS ~", var))
  
  # Tenta rodar o modelo com correlação de Pagel
  fit <- try(pgls.SEy(formula, data = data, tree = tree, corClass = corPagel), silent = TRUE)
  
  # Verifica se o modelo foi rodado corretamente
  if (inherits(fit, "try-error")) {
    return(data.frame(Variable = var, Beta = NA, T_value = NA, P_value = NA, AIC = NA, BIC = NA, Sigma = NA))
  }
  
  # Extrai o sumário do modelo
  fited <- summary(fit)
  
  # Extrai os coeficientes da tabela tTable
  coef_summary <- fited$tTable
  
  # Inicializa as variáveis
  beta <- NA
  t_value <- NA
  p_value <- NA
  
  # Verifica se a tabela está estruturada corretamente
  if (!is.null(coef_summary) && nrow(coef_summary) >= 2) {
    beta <- coef_summary[2, "Value"]
    t_value <- coef_summary[2, "t-value"]
    p_value <- coef_summary[2, "p-value"]
  }
  
  # Extrai AIC, BIC e sigma
  aic_value <- fited$AIC
  bic_value <- fited$BIC
  sigma_value <- fited$sigma
  
  # Retorna um data frame com os resultados
  return(data.frame(Variable = var, Beta = beta, T_value = t_value, P_value = p_value, AIC = aic_value, BIC = bic_value, Sigma = sigma_value))
}

# Aplicar a função para cada variável bio
results <- map_df(bio_vars, ~run_pgls(.x, data, tree))

# Salvar em CSV
write.csv(results, "Pagel/wright_pgls_pagel.csv", row.names = FALSE)


### FORMA 
run_pgls <- function(var, data, tree) {
  # Cria a fórmula dinamicamente
  formula <- as.formula(paste("y ~", var))
  
  # Tenta rodar o modelo com correlação de Pagel
  fit <- try(pgls.SEy(formula, data = data, tree = tree, corClass = corPagel), silent = TRUE)
  
  # Verifica se o modelo foi rodado corretamente
  if (inherits(fit, "try-error")) {
    return(data.frame(Variable = var, Beta = NA, T_value = NA, P_value = NA, AIC = NA, BIC = NA, Sigma = NA))
  }
  
  # Extrai o sumário do modelo
  fited <- summary(fit)
  
  # Extrai os coeficientes da tabela tTable
  coef_summary <- fited$tTable
  
  # Inicializa as variáveis
  beta <- NA
  t_value <- NA
  p_value <- NA
  
  # Verifica se a tabela está estruturada corretamente
  if (!is.null(coef_summary) && nrow(coef_summary) >= 2) {
    beta <- coef_summary[2, "Value"]
    t_value <- coef_summary[2, "t-value"]
    p_value <- coef_summary[2, "p-value"]
  }
  
  # Extrai AIC, BIC e sigma
  aic_value <- fited$AIC
  bic_value <- fited$BIC
  sigma_value <- fited$sigma
  
  # Retorna um data frame com os resultados
  return(data.frame(Variable = var, Beta = beta, T_value = t_value, P_value = p_value, AIC = aic_value, BIC = bic_value, Sigma = sigma_value))
}

# Aplicar a função para cada variável bio
resultsY <- map_df(bio_vars, ~run_pgls(.x, data, tree))

# Salvar em CSV
write.csv(resultsY, "Pagel/wrightY_pgls_pagel.csv", row.names = FALSE)


# Gráfico de coeficientes com personalizações
results %>%
  mutate(Variable = factor(Variable, levels = Variable[order(Beta)])) %>%
  ggplot(aes(x = Variable, y = Beta)) +
  geom_col(aes(fill = P_value < 0.05), color = "black") +  # Preenchimento e contorno
  scale_fill_manual(values = c("TRUE" = "dodgerblue", "FALSE" = "salmon"), 
                    labels = c("TRUE" = "Significativa (p < 0.05)", "FALSE" = "Não Significativa")) +
  geom_errorbar(aes(ymin = Beta - 1.96 * Sigma, ymax = Beta + 1.96 * Sigma), width = 0.2, color = "darkgrey") +
  coord_flip() +
  labs(title = "Coeficientes Beta dos Modelos PGLS (Pagel)",
       subtitle = "Análise das variáveis preditoras em relação ao tamanho",
       x = "Variáveis",
       y = "Coeficiente Beta",
       fill = "Significância") +
  theme_minimal(base_size = 15) +  # Aumenta o tamanho base da fonte
  theme(legend.position = "top") +  # Coloca a legenda no topo
  geom_text_repel(aes(label = ifelse(P_value < 0.05, round(Beta, 2), "")),  # Rótulos para coeficientes significativos
                  box.padding = 0.5, point.padding = 0.5, show.legend = FALSE) +
  theme(plot.title = element_text(hjust = 0.5, size = 20, face = "bold"),  # Alinhamento e estilo do título
        plot.subtitle = element_text(hjust = 0.5, size = 16, face = "italic"))  # Estilo do subtítulo


# Gráfico de coeficientes com personalizações
resultsY %>%
  mutate(Variable = factor(Variable, levels = Variable[order(Beta)])) %>%
  ggplot(aes(x = Variable, y = Beta)) +
  geom_col(aes(fill = P_value < 0.05), color = "black") +  # Preenchimento e contorno
  scale_fill_manual(values = c("TRUE" = "dodgerblue", "FALSE" = "salmon"), 
                    labels = c("TRUE" = "Significativa (p < 0.05)", "FALSE" = "Não Significativa")) +
  geom_errorbar(aes(ymin = Beta - 1.96 * Sigma, ymax = Beta + 1.96 * Sigma), width = 0.2, color = "darkgrey") +
  coord_flip() +
  labs(title = "Coeficientes Beta dos Modelos PGLS (Pagel)",
       subtitle = "Análise das variáveis preditoras em relação à forma",
       x = "Variáveis",
       y = "Coeficiente Beta",
       fill = "Significância") +
  theme_minimal(base_size = 15) +  # Aumenta o tamanho base da fonte
  theme(legend.position = "top") +  # Coloca a legenda no topo
  geom_text_repel(aes(label = ifelse(P_value < 0.05, round(Beta, 2), "")),  # Rótulos para coeficientes significativos
                  box.padding = 0.5, point.padding = 0.5, show.legend = FALSE) +
  theme(plot.title = element_text(hjust = 0.5, size = 20, face = "bold"),  # Alinhamento e estilo do título
        plot.subtitle = element_text(hjust = 0.5, size = 16, face = "italic"))  # Estilo do subtítulo

par(mfrow=c(1,1))




### TAMANHO + PCNM
run_pgls <- function(var, data, tree) {
  # Cria a fórmula dinamicamente
  formula <- as.formula(paste("CS ~", var, "+ pcnm1 + pcnm5"))
  
  # Tenta rodar o modelo com correlação de Pagel
  fit <- try(pgls.SEy(formula, data = data, tree = tree, corClass = corPagel), silent = TRUE)
  
  # Verifica se o modelo foi rodado corretamente
  if (inherits(fit, "try-error")) {
    return(data.frame(Variable = var, Beta = NA, T_value = NA, P_value = NA, AIC = NA, BIC = NA, Sigma = NA))
  }
  
  # Extrai o sumário do modelo
  fited <- summary(fit)
  
  # Extrai os coeficientes da tabela tTable
  coef_summary <- fited$tTable
  
  # Inicializa as variáveis
  beta <- NA
  t_value <- NA
  p_value <- NA
  
  # Verifica se a tabela está estruturada corretamente
  if (!is.null(coef_summary) && nrow(coef_summary) >= 2) {
    beta <- coef_summary[2, "Value"]
    t_value <- coef_summary[2, "t-value"]
    p_value <- coef_summary[2, "p-value"]
  }
  
  # Extrai AIC, BIC e sigma
  aic_value <- fited$AIC
  bic_value <- fited$BIC
  sigma_value <- fited$sigma
  
  # Retorna um data frame com os resultados
  return(data.frame(Variable = var, Beta = beta, T_value = t_value, P_value = p_value, AIC = aic_value, BIC = bic_value, Sigma = sigma_value))
}

# Aplicar a função para cada variável bio
results <- map_df(bio_vars, ~run_pgls(.x, data, tree))

# Salvar em CSV
write.csv(results, "wright_pcnm_pgls_pagel.csv", row.names = FALSE)



### FORMA + PCNM
run_pgls <- function(var, data, tree) {
  # Cria a fórmula dinamicamente
  formula <- as.formula(paste("y ~", var, "+ pcnm1 + pcnm5"))
  
  # Tenta rodar o modelo com correlação de Pagel
  fit <- try(pgls.SEy(formula, data = data, tree = tree, corClass = corPagel), silent = TRUE)
  
  # Verifica se o modelo foi rodado corretamente
  if (inherits(fit, "try-error")) {
    return(data.frame(Variable = var, Beta = NA, T_value = NA, P_value = NA, AIC = NA, BIC = NA, Sigma = NA))
  }
  
  # Extrai o sumário do modelo
  fited <- summary(fit)
  
  # Extrai os coeficientes da tabela tTable
  coef_summary <- fited$tTable
  
  # Inicializa as variáveis
  beta <- NA
  t_value <- NA
  p_value <- NA
  
  # Verifica se a tabela está estruturada corretamente
  if (!is.null(coef_summary) && nrow(coef_summary) >= 2) {
    beta <- coef_summary[2, "Value"]
    t_value <- coef_summary[2, "t-value"]
    p_value <- coef_summary[2, "p-value"]
  }
  
  # Extrai AIC, BIC e sigma
  aic_value <- fited$AIC
  bic_value <- fited$BIC
  sigma_value <- fited$sigma
  
  # Retorna um data frame com os resultados
  return(data.frame(Variable = var, Beta = beta, T_value = t_value, P_value = p_value, AIC = aic_value, BIC = bic_value, Sigma = sigma_value))
}

# Aplicar a função para cada variável bio
resultsY <- map_df(bio_vars, ~run_pgls(.x, data, tree))

# Salvar em CSV
write.csv(resultsY, "wrightY_pcnm_pgls_pagel.csv", row.names = FALSE)


# Gráfico de coeficientes com personalizações
results %>%
  mutate(Variable = factor(Variable, levels = Variable[order(Beta)])) %>%
  ggplot(aes(x = Variable, y = Beta)) +
  geom_col(aes(fill = P_value < 0.05), color = "black") +  # Preenchimento e contorno
  scale_fill_manual(values = c("TRUE" = "dodgerblue", "FALSE" = "salmon"), 
                    labels = c("TRUE" = "Significativa (p < 0.05)", "FALSE" = "Não Significativa")) +
  geom_errorbar(aes(ymin = Beta - 1.96 * Sigma, ymax = Beta + 1.96 * Sigma), width = 0.2, color = "darkgrey") +
  coord_flip() +
  labs(title = "Coeficientes Beta dos Modelos PGLS + PCNM (Pagel)",
       subtitle = "Análise das variáveis preditoras em relação ao tamanho",
       x = "Variáveis",
       y = "Coeficiente Beta",
       fill = "Significância") +
  theme_minimal(base_size = 15) +  # Aumenta o tamanho base da fonte
  theme(legend.position = "top") +  # Coloca a legenda no topo
  geom_text_repel(aes(label = ifelse(P_value < 0.05, round(Beta, 2), "")),  # Rótulos para coeficientes significativos
                  box.padding = 0.5, point.padding = 0.5, show.legend = FALSE) +
  theme(plot.title = element_text(hjust = 0.5, size = 20, face = "bold"),  # Alinhamento e estilo do título
        plot.subtitle = element_text(hjust = 0.5, size = 16, face = "italic"))  # Estilo do subtítulo


# Gráfico de coeficientes com personalizações
resultsY %>%
  mutate(Variable = factor(Variable, levels = Variable[order(Beta)])) %>%
  ggplot(aes(x = Variable, y = Beta)) +
  geom_col(aes(fill = P_value < 0.05), color = "black") +  # Preenchimento e contorno
  scale_fill_manual(values = c("TRUE" = "dodgerblue", "FALSE" = "salmon"), 
                    labels = c("TRUE" = "Significativa (p < 0.05)", "FALSE" = "Não Significativa")) +
  geom_errorbar(aes(ymin = Beta - 1.96 * Sigma, ymax = Beta + 1.96 * Sigma), width = 0.2, color = "darkgrey") +
  coord_flip() +
  labs(title = "Coeficientes Beta dos Modelos PGLS + PCNM (Pagel)",
       subtitle = "Análise das variáveis preditoras em relação à forma",
       x = "Variáveis",
       y = "Coeficiente Beta",
       fill = "Significância") +
  theme_minimal(base_size = 15) +  # Aumenta o tamanho base da fonte
  theme(legend.position = "top") +  # Coloca a legenda no topo
  geom_text_repel(aes(label = ifelse(P_value < 0.05, round(Beta, 2), "")),  # Rótulos para coeficientes significativos
                  box.padding = 0.5, point.padding = 0.5, show.legend = FALSE) +
  theme(plot.title = element_text(hjust = 0.5, size = 20, face = "bold"),  # Alinhamento e estilo do título
        plot.subtitle = element_text(hjust = 0.5, size = 16, face = "italic"))  # Estilo do subtítulo

par(mfrow=c(1,1))

