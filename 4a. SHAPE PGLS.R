#--------------------------------------------------------------------------
# FINAL SCRIPT FOR COMPLETE COMPARATIVE PGLS ANALYSIS
#--------------------------------------------------------------------------
# This script runs multiple evolutionary models for allometry-free shape data
# and data controlled for spatial structure (PCNM), and consolidates all
# results into a single table. It also performs model comparison
# using ANOVA (LRT), AIC-based ranking, and visualizes the best model.
#--------------------------------------------------------------------------

# --- 1. REQUIRED PACKAGES ---
# Make sure all are installed
library(tidyverse) # For data manipulation and plots (dplyr, purrr)
library(nlme)      # Essential for the gls function
library(ape)       # For correlation structures (corPagel, etc.)
library(geomorph)  # For morphometric functions and plotting
library(packfor)   # Only if you need forward.sel
library(progress)  # For a visual progress bar

#--------------------------------------------------------------------------
# --- 2. DATA PREPARATION AND STANDARDIZATION ---
#--------------------------------------------------------------------------
# CRITICAL: Assuming 'data', 'tree', and 'pcnm' are loaded.
# 'data' must contain a column with species names matching the tree tips.

# --- 2a. Variable Selection and Data Merging ---
# List of environmental predictor variables
bio_vars <- c(paste0("bio", 1:19), "humid", "npp", "soilmoist")

# Definition of spatial covariates (PCNMs)
pcnm_covariates <- c("pcnm2", "pcnm31", "pcnm1") 
colnames(pcnm) <- paste0("pcnm", 1:ncol(pcnm))
data <- bind_cols(data, as.data.frame(pcnm))

# --- 2b. CRITICAL: Ensure Species Column for PGLS ---
# The PGLS functions need a column with species names to match data to the tree.
# We assume your species names are in 'sp_ives' and create a 'species' column.
if("sp_ives" %in% names(data) && !"species" %in% names(data)) {
  data$species <- data$sp_ives
} else if (!"species" %in% names(data)) {
  stop("A column named 'species' or 'sp_ives' with names matching the tree tips is required in your 'data' frame.")
}
rownames(data) <- data$species

# --- 2c. Data Standardization ---
# Standardizes the response variable and all predictors.
data_scaled <- data %>%
  mutate(
    shape.2d = scale(shape.2d),
    across(all_of(bio_vars), scale)
  )

#--------------------------------------------------------------------------
# --- 3. MASTER PGLS FUNCTION ---
#--------------------------------------------------------------------------
run_pgls_master <- function(predictor_var, evo_model_name, cor_structure, data, covariates = NULL) {
  
  if (is.null(covariates)) {
    formula_str <- paste("shape.2d ~", predictor_var)
  } else {
    covariates_str <- paste(covariates, collapse = " + ")
    formula_str <- paste("shape.2d ~", predictor_var, "+", covariates_str)
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

message("Starting execution of all PGLS models...")

all_results <- pmap_df(analysis_grid, function(predictor, evo_model_name) {
  pb$tick(); result_pure <- run_pgls_master(predictor, evo_model_name, evo_models_list[[evo_model_name]], data_scaled, NULL)
  pb$tick(); result_pcnm <- run_pgls_master(predictor, evo_model_name, evo_models_list[[evo_model_name]], data_scaled, pcnm_covariates)
  bind_rows(result_pure, result_pcnm)
})

#--------------------------------------------------------------------------
# --- 5. FINAL CONSOLIDATION AND EXPORT ---
#--------------------------------------------------------------------------
final_results_table <- all_results %>% arrange(Variable, Evolutionary_Model, Model_Type)
#write.csv(final_results_table, "PGLS_Lima.csv", row.names = FALSE)
message("\nMain analysis completed successfully!")
message("All results have been saved to 'PGLS_Final_Compiled_Results.csv'")

#--------------------------------------------------------------------------
# --- 6. MODEL COMPARISON AND SELECTION ---
#--------------------------------------------------------------------------
message("\nCalculating AIC-based model ranks...")
aic_ranking <- final_results_table %>%
  filter(!is.na(AIC)) %>% group_by(Variable) %>%
  mutate(delta_AIC = AIC - min(AIC), likelihood = exp(-0.5 * delta_AIC), akaike_weight = likelihood / sum(likelihood)) %>%
  arrange(Variable, delta_AIC) %>% ungroup()
#write.csv(aic_ranking, "AIC - Lima.csv", row.names = FALSE)

best_models_summary <- aic_ranking %>%
  group_by(Variable) %>% slice_min(order_by = delta_AIC, n = 1, with_ties = FALSE) %>%
  dplyr::select(Variable, Best_Evo_Model = Evolutionary_Model, Best_Model_Type = Model_Type, AIC, akaike_weight)
message("Summary of the best models based on AIC:"); print(as.data.frame(best_models_summary))


#--------------------------------------------------------------------------
# --- 7. VISUALIZING THE BEST MODEL ---
#--------------------------------------------------------------------------
message("\nIdentifying best model and generating prediction plot...")

best_model_info <- aic_ranking %>%
  filter(Model_Type == "PCNM") %>%
  slice_min(order_by = AIC, n = 1, with_ties = FALSE)

  best_predictor_name <- best_model_info$Variable
  best_evo_model_name <- best_model_info$Evolutionary_Model
  
  message(paste("Best model found: shape ~", best_predictor_name, "with PCNMs and", best_evo_model_name, "evolution."))
  
  # Fit the final, best multivariate PGLS model again for plotting
  best_cor_structure <- evo_models_list[[best_evo_model_name]]
  best_formula <- as.formula(paste("shape.2d ~", best_predictor_name, "+", paste(pcnm_covariates, collapse = " + ")))
  # Verifique se sua linha está assim, especialmente a parte "phy = tree"
  gdf <- geomorph.data.frame(shape = adj.shape, phy = tree, data)
  final_model_fit <- procD.pgls(best_formula, phy = phy, data = gdf, corStruct = best_cor_structure, iter = 999)
  
  # --- 7a. Generate and save the shape change plot (TPS grid) ---
  
  # Define a color palette for consistency
  cores_fac <- c("AF" ="green3","AM" = "darkgreen", "SV" = "#FFDB58")
  # Create a vector of colors corresponding to each specimen's ecoregion
  plot_colors <- cores_fac[as.factor(gdf$fac)]
  
  #pdf("PGLS_Plot_Lima_ShapeChange.pdf", width = 8, height = 8)
  
  # CORRECTED: Use the 'plot_colors' vector for the 'bg' argument to color points.
  plot(final_model_fit, type = "regression", reg.type = "RegScore", 
       predictor = gdf[[best_predictor_name]],
       main = paste("Shape Change Along", best_predictor_name),
       pch = 21, cex = 3, bg = plot_colors, col = "gray20") # 'col' for point border
  
  # Add a legend manually
  legend("topleft", legend = names(cores_fac), 
         pch = 21, pt.bg = cores_fac, 
         title = "Ecoregion", bty = "n")
  
  dev.off()
  
  message("Shape change visualization plot saved to 'PGLS_Plot_Lima_ShapeChange.pdf'")
  
  # --- 7b. Generate and save the regression score plot ---
  
  # Extract regression scores
  regression_scores <- final_model_fit$pgls.residuals
  
  # Create a data frame for ggplot, including the 'fac' column for coloring
  plot_data <- data.frame(
    predictor = gdf[[best_predictor_name]],
    reg_score = regression_scores[,1], # Use the first column of regression scores
    ecoregion = as.factor(gdf$fac)     # Add the ecoregion factor
  )
  
  # Define a color palette for consistency
  cores_fac <- c("AF" ="green3","AM" = "darkgreen", "SV" = "goldenrod1")
  names(cores_fac) <- levels(plot_data$ecoregion)
  
  # Create the plot
  regression_plot <- ggplot(plot_data, aes(x = predictor, y = reg_score)) +
    geom_point(aes(color = ecoregion), size = 8, alpha = 0.8) + 
    geom_smooth(method = "lm", se = FALSE, color = "black", linewidth = 1.2, aes(group = 1)) + 
    scale_color_manual(name = "Ecoregion", values = cores_fac) +
    labs(
      title = "PGLS Regression Plot",
      subtitle = paste("Relationship between Shape and", best_predictor_name, "(controlled for phylogeny & space)"),
      x = paste("Predictor values for", best_predictor_name),
      y = "Regression Score (Shape variation along the predictor axis)"
    ) +
    theme_minimal(base_size = 15) +
    theme(legend.position = "bottom")
  
 #ggsave("PGLS_Plot_Lima_Regression.pdf", plot = regression_plot, width = 8, height = 6)
  message("Regression score plot saved to 'PGLS_Plot_Lima_Regression.pdf'")
  
  names(cores_fac)
  names(bioma_shapes)
  
  bioma_shapes <- c("AF" = 16, "AM" = 15, "SV" = 17)  # Shapes por bioma
  
  regression_plot <- ggplot(plot_data, aes(x = predictor, y = reg_score)) +
    geom_point(aes(color = ecoregion, shape = ecoregion), size = 8, alpha = 0.9) +
    geom_smooth(method = "lm", se = FALSE, color = "black", linewidth = 1.2, aes(group = 1)) +
    scale_color_manual(name = "Ecoregion", values = cores_fac) +
    scale_shape_manual(name = "Ecoregion", values = bioma_shapes) +
    labs(
      title = "PGLS Regression Plot",
      subtitle = paste("Shape variation along", best_predictor_name, "\n(controlled for phylogeny and spatial autocorrelation)"),
      x = paste("Predictor:", best_predictor_name),
      y = "Regression Score (shape axis)"
    ) +
    theme_classic(base_size = 15) +
    theme(
      legend.position = "bottom",
      legend.title = element_text(size = 13),
      legend.text = element_text(size = 11),
      plot.title = element_text(face = "bold", size = 18),
      plot.subtitle = element_text(size = 13),
      axis.title = element_text(size = 14)
    )
  
  print(regression_plot)
  
  # --- 7b. Predict shapes at the extremes of the best predictor ---
  # This section now uses the logic from your example script.
  
  # Get the vector of the best predictor variable from the geomorph data frame
  predictor_vec <- gdf[[best_predictor_name]]
  
  # Use shape.predictor to calculate the shapes at min and max values
  # Note: PGLS is not used here, this is a visualization of the raw relationship
  # as shown in your example script.
  preds <- shape.predictor(
    adj.shape, 
    x = predictor_vec,
    Intercept = FALSE, # Use FALSE as it's a single predictor
    predmin = min(predictor_vec), 
    predmax = max(predictor_vec)
  )
  
  # Calculate the mean shape (mshape)
  mshape <- mshape(adj.shape)
  
  # --- 7c. Generate and save the visualization plot ---
  # This requires your visualization objects ('Sapajusoutline', 'GP') to be loaded.
  
  if (exists("Sapajusoutline") && exists("GP")) {
    
    pdf("PGLS_Plot_Lima_ShapeExtremes.pdf", width = 12, height = 6)
    
    # Set up a 1x2 plot layout
    par(mfrow = c(1, 2))
    
    # Plot 1: Shape at the minimum predictor value
    plotRefToTarget(mshape, preds$predmin, 
                    mag = 2, outline = Sapajusoutline$outline, 
                    method = "points", gridPars = GP)
    title(paste("Shape at Min", best_predictor_name))
    
    # Plot 2: Shape at the maximum predictor value
    plotRefToTarget(mshape, preds$predmax, 
                    mag = 2, outline = Sapajusoutline$outline, 
                    method = "points", gridPars = GP)
    title(paste("Shape at Max", best_predictor_name))
    
    # Close the PDF device
    dev.off()
    
    # Reset plot layout
    par(mfrow = c(1, 1))
    
    message("Plot of shape extremes saved to 'PGLS_Plot_Lima_ShapeExtremes.pdf'")
    
  } else {
    message("Visualization of shape extremes skipped: 'Sapajusoutline' and/or 'GP' objects not found.")
  }
  
message("\nAll analyses are complete.")
  
