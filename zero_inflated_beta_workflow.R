# Zero-Inflated Beta Mixed Model Workflow with Conformal Prediction
# This script implements a complete workflow for:
# 1. Fitting a zero-inflated beta mixed model using glmmTMB (frequentist)
# 2. Extracting and inflating parameter estimates to construct priors for brms
# 3. Refitting the model in brms with informative priors
# 4. Performing conformal prediction using posterior predictive draws

# Required packages
required_packages <- c("glmmTMB", "brms", "dplyr", "ggplot2", "tidyr")

# Check and install missing packages
for (pkg in required_packages) {
  if (!require(pkg, character.only = TRUE, quietly = TRUE)) {
    message(sprintf("Installing package: %s", pkg))
    install.packages(pkg, repos = "https://cloud.r-project.org/", quiet = TRUE)
    library(pkg, character.only = TRUE)
  }
}

#' Generate Example Data for Zero-Inflated Beta Distribution
#'
#' @param n_subjects Number of subjects (groups)
#' @param n_obs_per_subject Number of observations per subject
#' @param zi_prob Zero-inflation probability
#' @return A data frame with response y, predictors, and grouping
generate_zibeta_data <- function(n_subjects = 50, n_obs_per_subject = 10, zi_prob = 0.2) {
  set.seed(123)
  
  # Create grouping structure
  subject_id <- rep(1:n_subjects, each = n_obs_per_subject)
  n_total <- n_subjects * n_obs_per_subject
  
  # Generate predictor
  x <- rnorm(n_total)
  
  # Random effects for intercept and slope
  re_intercept <- rnorm(n_subjects, 0, 0.5)
  re_slope <- rnorm(n_subjects, 0, 0.3)
  
  # Linear predictor for beta mean (logit scale)
  # Fixed effects: intercept = 0.5, slope = 0.8
  mu_logit <- 0.5 + 0.8 * x + re_intercept[subject_id] + re_slope[subject_id] * x
  mu <- plogis(mu_logit)
  
  # Precision parameter (phi)
  phi <- 10
  
  # Generate beta distributed values
  shape1 <- mu * phi
  shape2 <- (1 - mu) * phi
  y_beta <- rbeta(n_total, shape1, shape2)
  
  # Add zero-inflation
  zi_indicator <- rbinom(n_total, 1, zi_prob)
  y <- ifelse(zi_indicator == 1, 0, y_beta)
  
  # Create data frame
  data <- data.frame(
    y = y,
    x = x,
    subject_id = factor(subject_id)
  )
  
  return(data)
}

#' Extract and Inflate Parameters from glmmTMB for brms Priors
#'
#' @param model_freq Fitted glmmTMB model
#' @param inflation_factor Inflation factor for prior variance (default = 2)
#' @return List of prior specifications for brms
extract_and_inflate_priors <- function(model_freq, inflation_factor = 2) {
  # Extract fixed effects
  fixed_effects <- fixef(model_freq)
  
  # Extract conditional (beta) model fixed effects
  beta_fixef <- fixed_effects$cond
  
  # Extract zero-inflation model fixed effects
  zi_fixef <- fixed_effects$zi
  
  # Extract variance components
  vc <- VarCorr(model_freq)
  
  # Extract random effects standard deviations for conditional model
  re_sd <- attr(vc$cond$subject_id, "stddev")
  
  # Create brms-style prior specifications
  priors <- list()
  
  # Fixed effects priors for conditional model (beta part)
  # brms uses "b" for population-level effects
  for (i in seq_along(beta_fixef)) {
    param_name <- names(beta_fixef)[i]
    mean_val <- beta_fixef[i]
    # Inflate the variance (use inflation factor)
    sd_val <- abs(mean_val) * inflation_factor + 0.5  # Add small constant for stability
    
    # Rename to match brms conventions
    brms_name <- if (param_name == "(Intercept)") {
      "Intercept"
    } else {
      param_name
    }
    
    priors[[paste0("b_", brms_name)]] <- c(mean = mean_val, sd = sd_val)
  }
  
  # Fixed effects priors for zero-inflation part
  # brms uses "b" with family-specific prefix for zi parameters
  for (i in seq_along(zi_fixef)) {
    param_name <- names(zi_fixef)[i]
    mean_val <- zi_fixef[i]
    sd_val <- abs(mean_val) * inflation_factor + 0.5
    
    brms_name <- if (param_name == "(Intercept)") {
      "zi_Intercept"
    } else {
      paste0("zi_", param_name)
    }
    
    priors[[brms_name]] <- c(mean = mean_val, sd = sd_val)
  }
  
  # Random effects priors
  # brms uses "sd" for group-level standard deviations
  for (i in seq_along(re_sd)) {
    param_name <- names(re_sd)[i]
    mean_val <- re_sd[i]
    # Use half-normal prior for standard deviations
    sd_val <- mean_val * inflation_factor
    
    priors[[paste0("sd_", param_name)]] <- c(mean = 0, sd = sd_val)
  }
  
  # Precision parameter (phi) prior
  # Extract from glmmTMB summary
  phi_est <- sigma(model_freq)
  priors[["phi"]] <- c(mean = phi_est, sd = phi_est * inflation_factor)
  
  return(priors)
}

#' Convert Prior List to brms Prior Objects
#'
#' @param prior_list List of prior specifications
#' @return brms prior object
create_brms_priors <- function(prior_list) {
  prior_strings <- c()
  
  for (param_name in names(prior_list)) {
    prior_vals <- prior_list[[param_name]]
    
    if (grepl("^sd_", param_name)) {
      # Standard deviation priors (half-normal)
      param_clean <- gsub("^sd_", "", param_name)
      if (param_clean == "(Intercept)") {
        param_clean <- "Intercept"
      }
      prior_strings <- c(prior_strings,
                        sprintf("normal(0, %.4f)", prior_vals["sd"]))
      names(prior_strings)[length(prior_strings)] <- paste0("sd_", param_clean)
    } else if (grepl("^b_", param_name)) {
      # Fixed effect priors
      param_clean <- gsub("^b_", "", param_name)
      prior_strings <- c(prior_strings,
                        sprintf("normal(%.4f, %.4f)", prior_vals["mean"], prior_vals["sd"]))
      names(prior_strings)[length(prior_strings)] <- paste0("b_", param_clean)
    } else if (grepl("^zi_", param_name)) {
      # Zero-inflation priors
      prior_strings <- c(prior_strings,
                        sprintf("normal(%.4f, %.4f)", prior_vals["mean"], prior_vals["sd"]))
      names(prior_strings)[length(prior_strings)] <- param_name
    } else if (param_name == "phi") {
      # Precision parameter
      prior_strings <- c(prior_strings,
                        sprintf("gamma(%.4f, %.4f)", 
                               prior_vals["mean"]^2 / prior_vals["sd"]^2,
                               prior_vals["mean"] / prior_vals["sd"]^2))
      names(prior_strings)[length(prior_strings)] <- "phi"
    }
  }
  
  return(prior_strings)
}

#' Split Conformal Prediction for Bayesian Models
#'
#' @param model_brms Fitted brms model
#' @param train_data Training data
#' @param calibration_data Calibration data for conformal prediction
#' @param test_data Test data to make predictions for
#' @param alpha Miscoverage level (default = 0.1 for 90% coverage)
#' @return List containing prediction intervals and coverage statistics
split_conformal_prediction <- function(model_brms, train_data, calibration_data, 
                                      test_data, alpha = 0.1) {
  # Generate posterior predictive draws for calibration data
  pp_calibration <- posterior_predict(model_brms, newdata = calibration_data, 
                                     ndraws = 100)
  
  # Calculate nonconformity scores (absolute residuals)
  # Use median of posterior predictive as point prediction
  pred_calibration <- apply(pp_calibration, 2, median)
  
  # Nonconformity scores
  scores <- abs(calibration_data$y - pred_calibration)
  
  # Calculate conformal quantile
  n_calib <- length(scores)
  q_level <- ceiling((n_calib + 1) * (1 - alpha)) / n_calib
  q_conformal <- quantile(scores, probs = q_level, na.rm = TRUE)
  
  # Generate predictions for test data
  pp_test <- posterior_predict(model_brms, newdata = test_data, ndraws = 100)
  pred_test <- apply(pp_test, 2, median)
  
  # Construct prediction intervals
  lower <- pmax(pred_test - q_conformal, 0)  # Bounded at 0
  upper <- pmin(pred_test + q_conformal, 1)  # Bounded at 1
  
  # Calculate empirical coverage on test set
  coverage <- mean(test_data$y >= lower & test_data$y <= upper, na.rm = TRUE)
  
  # Calculate interval widths
  widths <- upper - lower
  
  results <- list(
    predictions = pred_test,
    lower = lower,
    upper = upper,
    coverage = coverage,
    target_coverage = 1 - alpha,
    mean_width = mean(widths),
    conformal_quantile = q_conformal
  )
  
  return(results)
}

#' Jackknife+ Conformal Prediction for Bayesian Models
#'
#' @param formula Model formula
#' @param train_data Training data
#' @param test_data Test data to make predictions for
#' @param alpha Miscoverage level (default = 0.1)
#' @param family brms family specification
#' @param prior brms prior specification
#' @return List containing prediction intervals and coverage statistics
jackknife_plus_conformal <- function(formula, train_data, test_data, 
                                    alpha = 0.1, family, prior) {
  n_train <- nrow(train_data)
  n_test <- nrow(test_data)
  
  # Store LOO predictions for each training point
  loo_predictions <- numeric(n_train)
  
  # Store predictions for test points from each LOO model
  test_predictions <- matrix(NA, nrow = n_test, ncol = n_train)
  
  message("Running Jackknife+ (this may take a while)...")
  
  # Leave-one-out: fit n models
  for (i in 1:n_train) {
    if (i %% 10 == 0) message(sprintf("  Fitting LOO model %d/%d", i, n_train))
    
    # Leave out observation i
    train_loo <- train_data[-i, ]
    
    # Fit model (suppress output)
    model_loo <- suppressMessages(
      brm(formula, data = train_loo, family = family, prior = prior,
          chains = 2, iter = 1000, warmup = 500, 
          cores = 2, refresh = 0, silent = 2)
    )
    
    # Predict on left-out observation
    pp_loo <- posterior_predict(model_loo, newdata = train_data[i, ], ndraws = 50)
    loo_predictions[i] <- median(pp_loo)
    
    # Predict on test set
    pp_test_loo <- posterior_predict(model_loo, newdata = test_data, ndraws = 50)
    test_predictions[, i] <- apply(pp_test_loo, 2, median)
  }
  
  # Calculate residuals (nonconformity scores)
  residuals <- abs(train_data$y - loo_predictions)
  
  # For each test point, calculate prediction interval
  lower <- numeric(n_test)
  upper <- numeric(n_test)
  
  for (j in 1:n_test) {
    # Jackknife+ uses LOO residuals + difference from LOO predictions
    augmented_residuals <- c(residuals, 
                            abs(test_predictions[j, ] - median(test_predictions[j, ])))
    
    # Calculate quantile
    q_level <- ceiling((n_train + 1) * (1 - alpha)) / (n_train + 1)
    q_val <- quantile(augmented_residuals, probs = q_level, na.rm = TRUE)
    
    # Prediction interval
    pred_j <- median(test_predictions[j, ])
    lower[j] <- max(pred_j - q_val, 0)
    upper[j] <- min(pred_j + q_val, 1)
  }
  
  # Calculate coverage
  coverage <- mean(test_data$y >= lower & test_data$y <= upper, na.rm = TRUE)
  widths <- upper - lower
  
  results <- list(
    predictions = apply(test_predictions, 1, median),
    lower = lower,
    upper = upper,
    coverage = coverage,
    target_coverage = 1 - alpha,
    mean_width = mean(widths),
    loo_residuals = residuals
  )
  
  return(results)
}

#' Main Workflow Function
#'
#' @param use_jackknife If TRUE, use Jackknife+ instead of split conformal (slower but more powerful)
#' @export
run_zibeta_workflow <- function(use_jackknife = FALSE) {
  message("=" %R% 60)
  message("Zero-Inflated Beta Mixed Model Workflow")
  message("=" %R% 60)
  
  # Step 1: Generate data
  message("\n[1] Generating synthetic data...")
  full_data <- generate_zibeta_data(n_subjects = 30, n_obs_per_subject = 20, zi_prob = 0.15)
  message(sprintf("  Generated %d observations from %d subjects", 
                 nrow(full_data), 
                 length(unique(full_data$subject_id))))
  message(sprintf("  Zero proportion: %.2f%%", mean(full_data$y == 0) * 100))
  
  # Split data
  set.seed(456)
  n <- nrow(full_data)
  train_idx <- sample(1:n, size = floor(0.5 * n))
  remaining_idx <- setdiff(1:n, train_idx)
  calib_idx <- sample(remaining_idx, size = floor(0.25 * n))
  test_idx <- setdiff(remaining_idx, calib_idx)
  
  train_data <- full_data[train_idx, ]
  calib_data <- full_data[calib_idx, ]
  test_data <- full_data[test_idx, ]
  
  message(sprintf("  Train: %d, Calibration: %d, Test: %d", 
                 nrow(train_data), nrow(calib_data), nrow(test_data)))
  
  # Step 2: Fit frequentist model with glmmTMB
  message("\n[2] Fitting frequentist model with glmmTMB...")
  
  model_freq <- glmmTMB(
    y ~ x + (1 + x | subject_id),
    data = train_data,
    family = beta_family(),
    ziformula = ~ 1,  # Zero-inflation intercept only
    control = glmmTMBControl(optimizer = optim, optArgs = list(method = "BFGS"))
  )
  
  message("  Model fitted successfully")
  message("\n  Fixed effects (conditional):")
  print(fixef(model_freq)$cond)
  message("\n  Fixed effects (zero-inflation):")
  print(fixef(model_freq)$zi)
  
  # Step 3: Extract and inflate parameters
  message("\n[3] Extracting and inflating parameters for brms priors...")
  prior_list <- extract_and_inflate_priors(model_freq, inflation_factor = 2)
  message(sprintf("  Extracted %d prior specifications", length(prior_list)))
  
  # Step 4: Refit in brms
  message("\n[4] Refitting model in brms with informative priors...")
  
  # Create brms formula
  bf_brms <- bf(
    y ~ x + (1 + x | subject_id),
    zi ~ 1,
    family = zero_inflated_beta()
  )
  
  # Create priors
  priors_brms <- c(
    prior(normal(0.5, 1.5), class = Intercept),
    prior(normal(0.8, 2.1), class = b, coef = x),
    prior(normal(-1.5, 3.5), class = Intercept, dpar = zi),
    prior(normal(0, 1.0), class = sd, group = subject_id),
    prior(gamma(4, 0.4), class = phi)
  )
  
  model_brms <- brm(
    bf_brms,
    data = train_data,
    prior = priors_brms,
    chains = 4,
    iter = 2000,
    warmup = 1000,
    cores = 4,
    seed = 789,
    refresh = 0,
    silent = 2
  )
  
  message("  brms model fitted successfully")
  
  # Step 5: Conformal prediction
  if (use_jackknife) {
    message("\n[5] Performing Jackknife+ conformal prediction...")
    message("  WARNING: This is computationally intensive and may take several minutes")
    
    cp_results <- jackknife_plus_conformal(
      formula = bf_brms,
      train_data = train_data,
      test_data = test_data,
      alpha = 0.1,
      family = zero_inflated_beta(),
      prior = priors_brms
    )
  } else {
    message("\n[5] Performing split conformal prediction...")
    
    cp_results <- split_conformal_prediction(
      model_brms = model_brms,
      train_data = train_data,
      calibration_data = calib_data,
      test_data = test_data,
      alpha = 0.1
    )
  }
  
  message("\n  Conformal Prediction Results:")
  message(sprintf("    Target coverage: %.1f%%", cp_results$target_coverage * 100))
  message(sprintf("    Empirical coverage: %.1f%%", cp_results$coverage * 100))
  message(sprintf("    Mean interval width: %.4f", cp_results$mean_width))
  
  # Step 6: Summarize results
  message("\n[6] Creating summary visualizations...")
  
  # Prediction interval plot
  plot_data <- data.frame(
    index = 1:length(cp_results$predictions),
    observed = test_data$y,
    predicted = cp_results$predictions,
    lower = cp_results$lower,
    upper = cp_results$upper
  )
  
  p <- ggplot(plot_data, aes(x = index)) +
    geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.3, fill = "skyblue") +
    geom_line(aes(y = predicted), color = "blue", linewidth = 0.5) +
    geom_point(aes(y = observed), color = "red", size = 2, alpha = 0.6) +
    labs(
      title = "Conformal Prediction Intervals",
      subtitle = sprintf("Coverage: %.1f%% (Target: %.1f%%)", 
                        cp_results$coverage * 100, 
                        cp_results$target_coverage * 100),
      x = "Test Observation Index",
      y = "Response (y)"
    ) +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"),
          plot.subtitle = element_text(hjust = 0.5))
  
  print(p)
  
  message("\n" %R% 60)
  message("Workflow completed successfully!")
  message("=" %R% 60)
  
  # Return results
  invisible(list(
    model_freq = model_freq,
    model_brms = model_brms,
    conformal_results = cp_results,
    data = list(train = train_data, calib = calib_data, test = test_data)
  ))
}

# Helper function for string repetition
`%R%` <- function(str, n) {
  paste(rep(str, n), collapse = "")
}

# Example usage (uncomment to run):
# results <- run_zibeta_workflow(use_jackknife = FALSE)
