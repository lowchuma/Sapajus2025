#################################################################
# STEP 1 (CORRECTED VERSION): LOAD AND SYNCHRONIZE DATA
#################################################################

# Load required packages
library(ape)
library(phylolm)

# Load your 3 phylogenies
tree_lima <- read.tree("Datasets/trees/AVGLOCSEX_Lima.tre")
tree_upham <- read.tree("Datasets/trees/AVGLOCSEX_Upham.tre")
tree_wright <- read.tree("Datasets/trees/AVGLOCSEX_Wright.tre")

# Load your data WITHOUT setting row names yet
original_data <- read.csv("Planilhas/avglocsex.csv", sep = ";")


# --- CRITICAL CORRECTION STEP ---

# Inspect your data frame to find the column name with the IDs
# Run this command to see all column names:
# names(original_data) 

# Assuming the column with names like 'Sapajus_robustus_8', etc., is called "sp_ives".
# **REPLACE 'sp_ives' WITH THE CORRECT COLUMN NAME IF DIFFERENT.**
# This line will set the row names of your data frame from that column.
rownames(original_data) <- original_data$sp_ives

# As a safety measure, ensure there are no spaces (replace with "_")
rownames(original_data) <- gsub(" ", "_", rownames(original_data))


# --- SYNCHRONIZATION (should work now) ---

# Put the original trees in a named list
original_tree_list <- list(
  Lima2018 = tree_lima,
  Upham2019 = tree_upham,
  Wright2015 = tree_wright
)

# Find the common species between ALL trees and the DATA
common_species <- rownames(original_data)
for (tree_name in names(original_tree_list)) {
  common_species <- intersect(common_species, original_tree_list[[tree_name]]$tip.label)
}

# Verification check
if (length(common_species) == 0) {
  stop("ERROR: No common species found. Check if the column name used for row names is correct.")
} else {
  cat(paste("\nSuccess! Found", length(common_species), "common samples.\n"))
}

# Filter the final data frame
final_data <- original_data[common_species, ]

# Prune the trees and create the final list
tree_list <- list()
for (tree_name in names(original_tree_list)) {
  original_tree <- original_tree_list[[tree_name]]
  tree_list[[tree_name]] <- drop.tip(original_tree, setdiff(original_tree$tip.label, common_species))
}

cat("Synchronization complete. You can now proceed to the AIC analysis.\n")

#################################################################
# STEP 4: COMPARING PHYLOGENIES WITH AIC AND AKAIKE WEIGHTS (CORRECTED VERSION)
#################################################################

# Create a data frame to store the comparison results
comparison_results <- data.frame(
  Phylogeny = names(tree_list),
  AIC = numeric(length(tree_list))
)

cat("Starting phylogeny comparison via AIC...\n")

# Loop to fit the model on each tree from your clean list of trees
for (i in 1:nrow(comparison_results)) {
  
  tree_name <- comparison_results$Phylogeny[i]
  current_tree <- tree_list[[tree_name]]
  
  # --- START OF CORRECTION ---
  # Remove the "root edge" from the tree, if it exists.
  current_tree$root.edge <- NULL
  # --- END OF CORRECTION ---
  
  # Ensure the data is in the same order as the current tree (good practice)
  ordered_data <- final_data[current_tree$tip.label, ]
  
  cat(paste("...Fitting model for tree:", tree_name, "\n"))
  
  # Fit the phylogenetic linear model with 'phylolm'
  # Replace 'log(CS) ~ bio12' with your main relationship of interest.
  fit <- try(phylolm(log(CS) ~ bio12, data = ordered_data, phy = current_tree, model = "lambda"))
  
  # Store the AIC value if the model ran successfully
  if (!inherits(fit, "try-error")) {
    comparison_results$AIC[i] <- AIC(fit)
  } else {
    comparison_results$AIC[i] <- NA # Mark as NA if the model failed
    cat(paste("!!! WARNING: The model for tree", tree_name, "failed to converge. !!!\n"))
  }
}

# --- Calculation of Akaike Weights ---

# Remove any rows where the model may have failed
comparison_results <- na.omit(comparison_results)

# Find the lowest AIC value (the best model in the set)
min_aic <- min(comparison_results$AIC)

# Calculate Delta AIC (ΔAIC): the difference of each model from the best one
comparison_results$delta_AIC <- comparison_results$AIC - min_aic

# Calculate the Akaike Weight
comparison_results$akaike_weight <- exp(-0.5 * comparison_results$delta_AIC) / sum(exp(-0.5 * comparison_results$delta_AIC))

# Round the results for better visualization
comparison_results <- transform(comparison_results,
                                AIC = round(AIC, 2),
                                delta_AIC = round(delta_AIC, 2),
                                akaike_weight = round(akaike_weight, 4)
)


cat("\n#################################################\n")
cat("--- FINAL COMPARISON RESULT ---\n")
cat("#################################################\n")
print(comparison_results)

# Identify the winning phylogeny
if(nrow(comparison_results) > 0) {
  winning_phylogeny <- comparison_results[which.min(comparison_results$AIC), ]
  cat("\nConclusion: The phylogeny with the most statistical support is '", winning_phylogeny$Phylogeny, "'\n", sep = "")
} else {
  cat("\nConclusion: No model converged successfully. Check your data and trees.\n")
}