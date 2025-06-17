#--------------------------------------------------------------------------
# FINAL SCRIPT FOR SIZE (logCS) PGLS ANALYSIS
#--------------------------------------------------------------------------
# This script runs multiple evolutionary models for body size (logCS)
# against environmental data, controlled for spatial structure (PCNM). It
# consolidates all results, performs model selection, and visualizes the best model.
#--------------------------------------------------------------------------

# --- 1. REQUIRED PACKAGES ---
# Make sure all are installed
library(tidyverse) # For data manipulation and plots (dplyr, purrr)
library(nlme)      # Essential for the gls function
library(ape)       # For correlation structures (corPagel, etc.)
library(progress)  # For a visual progress bar

#--------------------------------------------------------------------------
# --- 2. DATA PREPARATION AND STANDARDIZATION ---
#--------------------------------------------------------------------------
# CRITICAL: Assuming 'data', 'tree', and 'pcnm' are loaded.
# 'data' must contain a column 'CS' (Centroid Size).

# --- 2a. Variable Selection and Data Merging ---
# List of environmental predictor variables
bio_vars <- c(paste0("bio", 1:19), "humid", "npp", "soilmoist")

# Definition of spatial covariates (PCNMs)
pcnm_covariates <- c("pcnm1", "pcnm5") 
colnames(pcnm) <- paste0("pcnm", 1:ncol(pcnm))

# To prevent duplication on re-runs, remove any existing pcnm columns first.
data <- data %>% dplyr::select(-starts_with("pcnm"))
data <- bind_cols(data, as.data.frame(pcnm))


# --- 2b. CRITICAL: Ensure Species Column for PGLS ---
# The PGLS functions need a column with species names to match data to the tree.
if("sp_ives" %in% names(data) && !"species" %in% names(data)) {
  data$species <- data$sp_ives
} else if (!"species" %in% names(data)) {
  stop("A column named 'species' or 'sp_ives' with names matching the tree tips is required in your 'data' frame.")
}
rownames(data) <- data$species

# --- 2c. Data Standardization ---
# Standardizes the new response variable (logCS) and all predictors.
if(!"CS" %in% names(data)) stop("Data frame must contain a 'CS' (Centroid Size) column.")

data_scaled <- data %>%
  mutate(
    logCS = as.numeric(scale(log(CS))), # NEW RESPONSE VARIABLE
    across(all_of(bio_vars), ~ as.numeric(scale(.)))
  )

#--------------------------------------------------------------------------
# --- 3. MASTER PGLS FUNCTION ---
#--------------------------------------------------------------------------
run_pgls_master <- function(predictor_var, evo_model_name, cor_structure, data, covariates = NULL) {
  
  # NEW RESPONSE VARIABLE: logCS
  if (is.null(covariates)) {
    formula_str <- paste("logCS ~", predictor_var)
  } else {
    covariates_str <- paste(covariates, collapse = " + ")
    formula_str <- paste("logCS ~", predictor_var, "+", covariates_str)
  }
  formula <- as.formula(formula_str)
  
  fit_control <- glsControl(maxIter = 200, msMaxIter = 200)
  
  fit <- try(
    gls(formula, data = data, correlation = cor_structure, method = "REML", control = fit_control),
    silent = TRUE
  )
  
  if (inherits(fit, "try-error")) {
    return(data.frame(
      Variable = predictor_var, Beta = NA, T_value = NA, P_value = NA, AIC = NA, BIC = NA, Std_Error = NA,
      Evo_Param_Value = NA, Evolutionary_Model = evo_model_name, Model_Type = ifelse(is.null(covariates), "Pure", "PCNM")
    ))
  }
  
  fited <- summary(fit)
  coef_summary <- fited$tTable
  beta <- coef_summary[2, "Value"]; t_value <- coef_summary[2, "t-value"]; p_value <- coef_summary[2, "p-value"]; std_error <- coef_summary[2, "Std.Error"]
  evo_param <- try(coef(fit$modelStruct$corStruct, unconstrained = FALSE), silent = TRUE)
  if (inherits(evo_param, "try-error") || length(evo_param) == 0) evo_param <- NA
  
  data.frame(
    Variable = predictor_var, Beta = beta, T_value = t_value, P_value = p_value, AIC = fited$AIC,
    BIC = fited$BIC, Std_Error = std_error, Evo_Param_Value = evo_param,
    Evolutionary_Model = evo_model_name, Model_Type = ifelse(is.null(covariates), "Pure", "PCNM")
  )
}

#--------------------------------------------------------------------------
# --- 4. EXECUTION OF ALL ANALYSES ---
#--------------------------------------------------------------------------
evo_models_list <- list(
  "Brownian" = corBrownian(1, tree, form = ~species),
  "Pagel"    = corPagel(1, tree, fixed = FALSE, form = ~species),
  "Martins"  = corMartins(1, tree, fixed = FALSE, form = ~species),
  "Grafen"   = corGrafen(0.5, tree, fixed = FALSE, form = ~species)
)

analysis_grid <- expand_grid(predictor = bio_vars, evo_model_name = names(evo_models_list))

pb <- progress_bar$new(format = " Running models [:bar] :percent in :elapsed | ETA: :eta", total = nrow(analysis_grid) * 2, clear = FALSE, width = 80)

message("Starting execution of all PGLS models for Size (logCS)...")

all_results <- pmap_df(analysis_grid, function(predictor, evo_model_name) {
  pb$tick(); result_pure <- run_pgls_master(predictor, evo_model_name, evo_models_list[[evo_model_name]], data_scaled, NULL)
  pb$tick(); result_pcnm <- run_pgls_master(predictor, evo_model_name, evo_models_list[[evo_model_name]], data_scaled, pcnm_covariates)
  bind_rows(result_pure, result_pcnm)
})

#--------------------------------------------------------------------------
# --- 5. FINAL CONSOLIDATION AND EXPORT ---
#--------------------------------------------------------------------------
final_results_table <- all_results %>% arrange(Variable, Evolutionary_Model, Model_Type)
#write.csv(final_results_table, "PGLS_Size_Lima.csv", row.names = FALSE)
message("\nMain analysis completed successfully!")
message("All results have been saved to 'PGLS_Size_Lima.csv'")

#--------------------------------------------------------------------------
# --- 6. MODEL COMPARISON AND SELECTION ---
#--------------------------------------------------------------------------
message("\nCalculating AIC-based model ranks...")
aic_ranking <- final_results_table %>%
  filter(!is.na(AIC)) %>% group_by(Variable) %>%
  mutate(delta_AIC = AIC - min(AIC, na.rm = TRUE), 
         likelihood = exp(-0.5 * delta_AIC),
         akaike_weight = likelihood / sum(likelihood)) %>%
  arrange(Variable, delta_AIC) %>% ungroup()
#write.csv(aic_ranking, "AIC_Size_Lima.csv", row.names = FALSE)

best_models_summary <- aic_ranking %>%
  group_by(Variable) %>% slice_min(order_by = delta_AIC, n = 1, with_ties = FALSE) %>%
  dplyr::select(Variable, Best_Evo_Model = Evolutionary_Model, Best_Model_Type = Model_Type, AIC, akaike_weight)
message("Summary of the best models based on AIC:"); print(as.data.frame(best_models_summary))


#--------------------------------------------------------------------------
# --- 7. VISUALIZING THE BEST MODEL ---
#--------------------------------------------------------------------------
message("\nIdentifying best model and generating a regression plot...")

best_model_info <- aic_ranking %>%
  filter(Model_Type == "PCNM") %>%
  slice_min(order_by = AIC, n = 1, with_ties = FALSE)

if (nrow(best_model_info) > 0) {
  
  best_predictor_name <- best_model_info$Variable
  best_evo_model_name <- best_model_info$Evolutionary_Model
  
  message(paste("Best model found: logCS ~", best_predictor_name, "with PCNMs and", best_evo_model_name, "evolution."))
  
  # Fit the final, best model again to get predictions for plotting
  best_cor_structure <- evo_models_list[[best_evo_model_name]]
  best_formula_str <- paste("logCS ~", best_predictor_name, "+", paste(pcnm_covariates, collapse = " + "))
  best_formula <- as.formula(best_formula_str)
  
  final_model_fit <- gls(best_formula, data = data_scaled, correlation = best_cor_structure, method = "REML")
  
  # Add predicted values to the data frame
  data_scaled$predicted_logCS <- predict(final_model_fit)
  
  # --- Generate and save the regression plot ---
  # Create a data frame for ggplot, including the 'fac' column for coloring
  plot_data <- data.frame(
    predictor = data_scaled[[best_predictor_name]],
    response = data_scaled$logCS,
    predicted = data_scaled$predicted_logCS,
    ecoregion = as.factor(data_scaled$fac)
  )
  
  cores_fac <- c("AF" ="green3","AM" = "darkgreen", "SV" = "#FFDB58")
  names(cores_fac) <- levels(plot_data$ecoregion)
  
  regression_plot <- ggplot(plot_data, aes(x = predictor, y = response)) +
    geom_point(aes(color = ecoregion), size = 3, alpha = 0.8) + 
    geom_line(aes(y = predicted), color = "black", linewidth = 1.2) + 
    scale_color_manual(name = "Ecoregion", values = cores_fac) +
    labs(
      title = "PGLS Regression of Size vs. Environment",
      subtitle = paste("Best Model: log(CS) ~", best_predictor_name, "+ PCNMs (", best_evo_model_name, ")"),
      x = paste("Standardized", best_predictor_name),
      y = "Standardized log(Centroid Size)"
    ) +
    theme_minimal(base_size = 15) +
    theme(legend.position = "bottom")
  
  ggsave("PGLS_Plot_Size_Lima.pdf", plot = regression_plot, width = 8, height = 6)
  
  message("Regression plot saved to 'PGLS_Plot_Size_Lima.pdf'")
  
} else {
  message("Visualization skipped: Could not identify a best PCNM model from the results.")
}

message("\nAll analyses are complete.")

# --- Generate and save the regression plot ---
# Create a data frame for ggplot, including the 'fac' column for coloring
plot_data <- data.frame(
  predictor = data_scaled[[best_predictor_name]],
  response = data_scaled$logCS,
  ecoregion = as.factor(data_scaled$fac)
)

cores_fac <- c("AF" = "green3", "AM" = "darkgreen", "SV" = "#FFDB58")
names(cores_fac) <- levels(plot_data$ecoregion)

# Extrai os coeficientes do modelo para a linha de regressão
model_coeffs <- coef(final_model_fit)
intercept <- model_coeffs["(Intercept)"]
slope <- model_coeffs[best_predictor_name]

# Bloco de código revisado e limpo
# Shapes com preenchimento (21–25)
bioma_shapes <- c("AF" = 21, "AM" = 22, "SV" = 24)

regression_plot <- ggplot(plot_data, aes(x = predictor, y = response)) +
  geom_point(
    aes(fill = ecoregion, shape = ecoregion),
    size = 8,
    alpha = 0.85,
    stroke = 0.7,
    color = "black"  # contorno
  ) +
  geom_abline(intercept = intercept, slope = slope, color = "black", linewidth = 1.2) +
  scale_fill_manual(name = "Ecoregion", values = cores_fac) +
  scale_shape_manual(name = "Ecoregion", values = bioma_shapes) +
  labs(
    title = "PGLS Regression of Size vs. Environment",
    subtitle = paste("Best Model: log(CS) ~", best_predictor_name, "+ PCNMs (", best_evo_model_name, ")"),
    x = paste("Standardized", best_predictor_name),
    y = "Standardized log(Centroid Size)"
  ) +
  theme_classic(base_size = 15) +
  theme(
    legend.position = "bottom",
    legend.title = element_text(size = 13),
    legend.text = element_text(size = 11),
    plot.title = element_text(face = "bold", size = 18),
    plot.subtitle = element_text(size = 13),
    axis.title = element_text(size = 14)
  ) +
  guides(
    fill = guide_legend(override.aes = list(shape = 21)),
    shape = guide_legend(override.aes = list(fill = "white"))
  )

print(regression_plot)
ggsave("PGLS_Plot_Size_Lima.pdf", plot = regression_plot, width = 8, height = 6)

message("Regression plot saved to 'PGLS_Plot_Size_Lima.pdf'")
