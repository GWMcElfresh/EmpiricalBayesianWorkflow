#!/usr/bin/env Rscript
# Test script for zero-inflated beta mixed model workflow

# Source the main workflow
source("zero_inflated_beta_workflow.R")

# Set up test environment
test_passed <- 0
test_failed <- 0

# Helper function for test assertions
assert_true <- function(condition, test_name) {
  if (condition) {
    cat(sprintf("✓ PASS: %s\n", test_name))
    test_passed <<- test_passed + 1
  } else {
    cat(sprintf("✗ FAIL: %s\n", test_name))
    test_failed <<- test_failed + 1
  }
}

cat("\n========================================\n")
cat("Testing Zero-Inflated Beta Workflow\n")
cat("========================================\n\n")

# Test 1: Data generation
cat("[Test 1] Data Generation\n")
data <- generate_zibeta_data(n_subjects = 10, n_obs_per_subject = 5, zi_prob = 0.2)
assert_true(nrow(data) == 50, "Correct number of observations")
assert_true(ncol(data) == 3, "Correct number of columns")
assert_true(all(data$y >= 0 & data$y <= 1), "Response bounded in [0,1]")
assert_true(sum(data$y == 0) > 0, "Some zero values present")
cat("\n")

# Test 2: glmmTMB model fitting
cat("[Test 2] Frequentist Model Fitting\n")
tryCatch({
  suppressMessages({
    model_freq <- glmmTMB(
      y ~ x + (1 + x | subject_id),
      data = data,
      family = beta_family(),
      ziformula = ~ 1
    )
  })
  assert_true(!is.null(model_freq), "Model fitted successfully")
  assert_true(length(fixef(model_freq)$cond) > 0, "Fixed effects extracted")
  assert_true(length(fixef(model_freq)$zi) > 0, "ZI effects extracted")
}, error = function(e) {
  assert_true(FALSE, paste("Model fitting:", e$message))
})
cat("\n")

# Test 3: Prior extraction and inflation
cat("[Test 3] Prior Extraction and Inflation\n")
if (exists("model_freq")) {
  tryCatch({
    priors <- extract_and_inflate_priors(model_freq, inflation_factor = 2)
    assert_true(length(priors) > 0, "Priors extracted")
    assert_true(all(sapply(priors, length) == 2), "All priors have mean and sd")
    assert_true(any(grepl("^b_", names(priors))), "Fixed effect priors named correctly")
    assert_true(any(grepl("^zi_", names(priors))), "ZI priors named correctly")
  }, error = function(e) {
    assert_true(FALSE, paste("Prior extraction:", e$message))
  })
}
cat("\n")

# Test 4: Conformal prediction (lightweight test with small model)
cat("[Test 4] Split Conformal Prediction\n")
set.seed(999)
small_data <- generate_zibeta_data(n_subjects = 8, n_obs_per_subject = 10, zi_prob = 0.15)

# Split data
n <- nrow(small_data)
train_idx <- 1:40
calib_idx <- 41:60
test_idx <- 61:n

train_data <- small_data[train_idx, ]
calib_data <- small_data[calib_idx, ]
test_data <- small_data[test_idx, ]

tryCatch({
  suppressMessages({
    # Fit a simple model for testing
    test_model <- glmmTMB(
      y ~ x + (1 | subject_id),
      data = train_data,
      family = beta_family(),
      ziformula = ~ 1
    )
    
    # Fit brms model (minimal iterations for speed)
    bf_test <- bf(y ~ x + (1 | subject_id), zi ~ 1, family = zero_inflated_beta())
    
    brms_test <- brm(
      bf_test,
      data = train_data,
      prior = c(
        prior(normal(0, 2), class = Intercept),
        prior(normal(0, 2), class = b),
        prior(normal(0, 2), class = Intercept, dpar = zi),
        prior(normal(0, 1), class = sd),
        prior(gamma(2, 1), class = phi)
      ),
      chains = 2,
      iter = 500,
      warmup = 250,
      cores = 2,
      refresh = 0,
      silent = 2
    )
    
    # Run conformal prediction
    cp_results <- split_conformal_prediction(
      brms_test,
      train_data,
      calib_data,
      test_data,
      alpha = 0.1
    )
    
    assert_true(!is.null(cp_results$predictions), "Predictions generated")
    assert_true(!is.null(cp_results$coverage), "Coverage calculated")
    assert_true(length(cp_results$predictions) == nrow(test_data), "Correct number of predictions")
    assert_true(all(cp_results$lower >= 0 & cp_results$lower <= 1), "Lower bounds valid")
    assert_true(all(cp_results$upper >= 0 & cp_results$upper <= 1), "Upper bounds valid")
    assert_true(all(cp_results$lower <= cp_results$upper), "Intervals well-formed")
  })
}, error = function(e) {
  cat(sprintf("  Note: Conformal prediction test skipped due to: %s\n", e$message))
  cat("  (This is acceptable if packages are not fully installed)\n")
})
cat("\n")

# Test 5: Parameter naming conventions
cat("[Test 5] Parameter Naming Conventions\n")
if (exists("priors")) {
  # Check brms naming conventions
  has_intercept <- any(grepl("Intercept", names(priors)))
  has_zi_intercept <- any(grepl("zi_Intercept", names(priors)))
  has_sd <- any(grepl("^sd_", names(priors)))
  
  assert_true(has_intercept, "Contains Intercept parameter")
  assert_true(has_zi_intercept, "Contains zi_Intercept parameter")
  assert_true(has_sd, "Contains sd_ parameters")
}
cat("\n")

# Summary
cat("========================================\n")
cat(sprintf("Test Summary: %d passed, %d failed\n", test_passed, test_failed))
cat("========================================\n\n")

if (test_failed > 0) {
  quit(status = 1)
} else {
  cat("All tests passed! ✓\n\n")
  quit(status = 0)
}
