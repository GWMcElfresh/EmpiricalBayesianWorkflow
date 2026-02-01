#!/usr/bin/env Rscript
# Example script demonstrating the zero-inflated beta mixed model workflow

# Source the main workflow functions
source("zero_inflated_beta_workflow.R")

# Example 1: Run the full workflow with split conformal prediction
cat("\n========================================\n")
cat("Example 1: Split Conformal Prediction\n")
cat("========================================\n\n")

results_split <- run_zibeta_workflow(use_jackknife = FALSE)

# Access results
cat("\nModel Summary:\n")
cat("- Frequentist model (glmmTMB): ", class(results_split$model_freq)[1], "\n")
cat("- Bayesian model (brms): ", class(results_split$model_brms)[1], "\n")
cat("- Conformal coverage: ", sprintf("%.2f%%", results_split$conformal_results$coverage * 100), "\n")

# Example 2: Run with Jackknife+ (more computationally intensive)
# Uncomment to run:
# cat("\n\n========================================\n")
# cat("Example 2: Jackknife+ Conformal Prediction\n")
# cat("========================================\n\n")
# 
# results_jackknife <- run_zibeta_workflow(use_jackknife = TRUE)
# 
# cat("\nComparison:\n")
# cat("- Split conformal coverage: ", sprintf("%.2f%%", results_split$conformal_results$coverage * 100), "\n")
# cat("- Jackknife+ coverage: ", sprintf("%.2f%%", results_jackknife$conformal_results$coverage * 100), "\n")

cat("\n\nWorkflow completed successfully!\n")
