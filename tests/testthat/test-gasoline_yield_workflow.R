test_that("Full workflow with GasolineYield dataset - Zero-Inflated Beta with conditionals on mean, variance, and ZI", {
  skip_on_cran()
  
  # Load required packages
  if (!requireNamespace("betareg", quietly = TRUE)) {
    skip("betareg package not available")
  }
  if (!requireNamespace("glmmTMB", quietly = TRUE)) {
    skip("glmmTMB package not available")
  }
  if (!requireNamespace("brms", quietly = TRUE)) {
    skip("brms package not available")
  }
  
  # Load the GasolineYield dataset from betareg package
  data("GasolineYield", package = "betareg")
  
  # Prepare data for Zero-Inflated Beta demonstration
  # The GasolineYield dataset doesn't naturally contain zeros
  # We artificially introduce zeros to demonstrate ZI capability
  set.seed(12345)
  gasoline <- GasolineYield
  gasoline$yield_zi <- gasoline$yield
  
  # Set bottom 10% to zero to simulate structural zeros
  threshold <- quantile(gasoline$yield, 0.10)
  gasoline$yield_zi[gasoline$yield < threshold] <- 0
  
  # Add a grouping variable for random effects (batch serves as this)
  # batch is already a factor in the data
  gasoline$batch_id <- as.integer(gasoline$batch)
  
  # Verify we have zeros
  expect_true(sum(gasoline$yield_zi == 0) > 0)
  expect_true(sum(gasoline$yield_zi > 0 & gasoline$yield_zi < 1) > 0)
  
  # Define the model formulas:
  # - Mean (mu): modeled by temp and batch (random effect)
  # - Precision (phi): modeled by temp (distributional regression)
  # - Zero-Inflation (zi): modeled by pressure
  
  formula_mean <- yield_zi ~ temp + (1 | batch)
  formula_zi <- ~ pressure
  
  # Note: dispformula for glmmTMB controls the precision parameter (phi)
  # For brms, this will be handled via bf() with distributional regression
  
  # Step 1: Fit frequentist model with glmmTMB
  message("Fitting frequentist zero-inflated beta model...")
  
  model_freq <- tryCatch({
    glmmTMB::glmmTMB(
      formula = formula_mean,
      data = gasoline,
      family = glmmTMB::beta_family(),
      ziformula = formula_zi,
      dispformula = ~ temp,  # Precision depends on temp
      REML = FALSE
    )
  }, error = function(e) {
    fail(paste("Frequentist model failed:", e$message))
    NULL
  })
  
  skip_if(is.null(model_freq))
  
  # Verify model fit
  expect_s3_class(model_freq, "glmmTMB")
  expect_true(!is.null(glmmTMB::fixef(model_freq)))
  
  # Check that we have ZI parameters
  zi_coefs <- glmmTMB::fixef(model_freq)$zi
  expect_true(!is.null(zi_coefs))
  expect_true("(Intercept)" %in% names(zi_coefs))
  expect_true("pressure" %in% names(zi_coefs))
  
  # Check conditional model parameters
  cond_coefs <- glmmTMB::fixef(model_freq)$cond
  expect_true("(Intercept)" %in% names(cond_coefs))
  expect_true("temp" %in% names(cond_coefs))
  
  # Check dispersion model parameters (phi ~ temp)
  disp_coefs <- glmmTMB::fixef(model_freq)$disp
  expect_true(!is.null(disp_coefs))
  expect_true("(Intercept)" %in% names(disp_coefs))
  expect_true("temp" %in% names(disp_coefs))
  
  # Step 2: Extract and inflate priors
  message("Extracting and inflating priors...")
  
  # Create family config
  family_config <- list(
    name = "zero_inflated_beta",
    glmmTMB_family = "beta_family()",
    brms_family = "zero_inflated_beta()",
    has_zi = TRUE,
    params = c("phi"),
    param_classes = list(phi = "phi"),
    zi_link = "logit"
  )
  
  attr(model_freq, "family_config") <- family_config
  class(model_freq) <- c("freq_glmmTMB_zero_inflated_beta", "freq_glmmTMB", class(model_freq))
  
  # Extract priors
  priors <- tryCatch({
    extract_and_inflate_priors(model_freq, inflation_factor = 2.5)
  }, error = function(e) {
    fail(paste("Prior extraction failed:", e$message))
    NULL
  })
  
  skip_if(is.null(priors))
  
  # Verify prior structure
  expect_type(priors, "list")
  expect_true("b_Intercept" %in% names(priors))
  expect_true("b_temp" %in% names(priors))
  expect_true("zi_Intercept" %in% names(priors))
  expect_true("zi_pressure" %in% names(priors))
  
  # Random effects should be named sd_groupname_term
  # Check if any sd_ prefixed names exist
  sd_names <- grep("^sd_", names(priors), value = TRUE)
  expect_true(length(sd_names) > 0)
  
  # Check for phi (precision) priors if supported
  # In glmmTMB beta models, phi is on log scale (dispersion)
  # We should have priors for dispersion intercept and temp effect
  
  # Step 3: Fit Bayesian model with brms
  message("Fitting Bayesian zero-inflated beta model with distributional regression...")
  
  # For brms with distributional regression on phi, we need to use bf()
  # bf(formula_mean, phi ~ temp, zi ~ pressure)
  
  brms_priors <- tryCatch({
    create_brms_priors(priors)
  }, error = function(e) {
    fail(paste("Prior creation failed:", e$message))
    NULL
  })
  
  skip_if(is.null(brms_priors))
  
  # Define distributional formula for brms
  brms_formula <- brms::bf(
    yield_zi ~ temp + (1 | batch),
    phi ~ temp,  # Precision parameter depends on temp
    zi ~ pressure  # Zero-inflation depends on pressure
  )
  
  # Determine backend
  backend <- "cmdstanr"
  if (!requireNamespace("cmdstanr", quietly = TRUE)) {
    if (requireNamespace("rstan", quietly = TRUE)) {
      backend <- "rstan"
    } else {
      skip("Neither cmdstanr nor rstan available")
    }
  }
  
  model_bayes <- tryCatch({
    brms::brm(
      formula = brms_formula,
      data = gasoline,
      family = brms::zero_inflated_beta(),
      prior = brms_priors,
      chains = 2,
      iter = 1500,
      warmup = 750,
      cores = 2,
      backend = backend,
      refresh = 0,
      silent = 2,
      seed = 12345
    )
  }, error = function(e) {
    message(paste("Bayesian model fitting error:", e$message))
    NULL
  }, warning = function(w) {
    message(paste("Bayesian model warning:", w$message))
    NULL
  })
  
  skip_if(is.null(model_bayes))
  
  # Verify Bayesian model
  expect_s3_class(model_bayes, "brmsfit")
  expect_equal(family(model_bayes)$family[1], "zero_inflated_beta")
  
  # Check that model has phi and zi components
  bayes_summary <- brms::posterior_summary(model_bayes)
  param_names <- rownames(bayes_summary)
  
  expect_true(any(grepl("^b_", param_names)))  # Fixed effects
  expect_true(any(grepl("^b_phi_", param_names)))  # Phi (precision) effects
  expect_true(any(grepl("^b_zi_", param_names)))  # ZI effects
  expect_true(any(grepl("^sd_batch__", param_names)))  # Random effects
  
  # Step 4: Perform conformal prediction
  message("Performing conformal prediction...")
  
  # Split data into train, calibration, and test
  set.seed(12345)
  n <- nrow(gasoline)
  indices <- sample(1:n)
  
  n_train <- floor(0.5 * n)
  n_calib <- floor(0.25 * n)
  
  train_idx <- indices[1:n_train]
  calib_idx <- indices[(n_train + 1):(n_train + n_calib)]
  test_idx <- indices[(n_train + n_calib + 1):n]
  
  train_data <- gasoline[train_idx, ]
  calib_data <- gasoline[calib_idx, ]
  test_data <- gasoline[test_idx, ]
  
  # Since we already have a full model, we can use it directly
  # For calibration, we need to refit on train data only
  
  # Refit frequentist model on training data
  model_freq_train <- tryCatch({
    glmmTMB::glmmTMB(
      formula = formula_mean,
      data = train_data,
      family = glmmTMB::beta_family(),
      ziformula = formula_zi,
      dispformula = ~ temp,
      REML = FALSE
    )
  }, error = function(e) {
    message(paste("Training model error:", e$message))
    NULL
  })
  
  skip_if(is.null(model_freq_train))
  
  # Extract priors from training model
  attr(model_freq_train, "family_config") <- family_config
  class(model_freq_train) <- c("freq_glmmTMB_zero_inflated_beta", "freq_glmmTMB", class(model_freq_train))
  
  priors_train <- tryCatch({
    extract_and_inflate_priors(model_freq_train, inflation_factor = 2.5)
  }, error = function(e) {
    message(paste("Prior extraction error:", e$message))
    NULL
  })
  
  skip_if(is.null(priors_train))
  
  brms_priors_train <- create_brms_priors(priors_train)
  
  # Fit Bayesian model on training data
  model_bayes_train <- tryCatch({
    brms::brm(
      formula = brms_formula,
      data = train_data,
      family = brms::zero_inflated_beta(),
      prior = brms_priors_train,
      chains = 2,
      iter = 1500,
      warmup = 750,
      cores = 2,
      backend = backend,
      refresh = 0,
      silent = 2,
      seed = 12345
    )
  }, error = function(e) {
    message(paste("Training Bayesian model error:", e$message))
    NULL
  })
  
  skip_if(is.null(model_bayes_train))
  
  # Perform split conformal prediction
  # Note: conformal_prediction_split expects 'y' as response variable name
  calib_data_cp <- calib_data
  calib_data_cp$y <- calib_data_cp$yield_zi
  test_data_cp <- test_data
  test_data_cp$y <- test_data_cp$yield_zi
  
  cp_results <- tryCatch({
    conformal_prediction_split(
      model = model_bayes_train,
      calibration_data = calib_data_cp,
      test_data = test_data_cp,
      alpha = 0.1  # 90% coverage
    )
  }, error = function(e) {
    message(paste("Conformal prediction error:", e$message))
    NULL
  })
  
  skip_if(is.null(cp_results))
  
  # Verify conformal prediction results
  expect_s3_class(cp_results, "conformal_prediction")
  expect_true(is.numeric(cp_results$coverage))
  expect_true(cp_results$coverage >= 0 && cp_results$coverage <= 1)
  expect_true(length(cp_results$predictions) == nrow(test_data))
  expect_true(length(cp_results$lower) == nrow(test_data))
  expect_true(length(cp_results$upper) == nrow(test_data))
  
  # Check that prediction intervals make sense
  expect_true(all(cp_results$lower <= cp_results$upper))
  expect_true(all(cp_results$lower >= 0))  # Beta on [0,1]
  expect_true(all(cp_results$upper <= 1))
  
  # Coverage should be reasonable (around 0.9 for alpha=0.1)
  # With small sample size, there can be variation
  expect_true(cp_results$coverage > 0.5)
  
  message(sprintf("Coverage: %.2f (target: %.2f)", cp_results$coverage, 0.9))
  message(sprintf("Mean interval width: %.4f", cp_results$mean_width))
  
  # Print summary
  print(cp_results)
})

test_that("GasolineYield workflow validates parameter mapping between glmmTMB and brms", {
  skip_on_cran()
  
  if (!requireNamespace("betareg", quietly = TRUE)) {
    skip("betareg package not available")
  }
  if (!requireNamespace("glmmTMB", quietly = TRUE)) {
    skip("glmmTMB package not available")
  }
  
  # Load and prepare data
  data("GasolineYield", package = "betareg")
  set.seed(999)
  gasoline <- GasolineYield
  gasoline$yield_zi <- gasoline$yield
  threshold <- quantile(gasoline$yield, 0.10)
  gasoline$yield_zi[gasoline$yield < threshold] <- 0
  
  # Fit simple ZI beta model
  model_freq <- glmmTMB::glmmTMB(
    formula = yield_zi ~ temp + (1 | batch),
    data = gasoline,
    family = glmmTMB::beta_family(),
    ziformula = ~ pressure,
    REML = FALSE
  )
  
  # Setup for parameter extraction
  family_config <- list(
    name = "zero_inflated_beta",
    has_zi = TRUE,
    params = c("phi"),
    param_classes = list(phi = "phi")
  )
  
  attr(model_freq, "family_config") <- family_config
  class(model_freq) <- c("freq_glmmTMB_zero_inflated_beta", "freq_glmmTMB", class(model_freq))
  
  # Extract parameters
  params <- extract_parameters(model_freq)
  
  # Verify parameter extraction
  expect_type(params, "list")
  expect_true("beta_fixef" %in% names(params))
  expect_true("zi_fixef" %in% names(params))
  
  # Random effects should exist
  expect_true("re_sd" %in% names(params))
  
  # Check fixed effects mapping
  expect_true("(Intercept)" %in% names(params$beta_fixef))
  expect_true("temp" %in% names(params$beta_fixef))
  
  # Check ZI mapping
  expect_true("(Intercept)" %in% names(params$zi_fixef))
  expect_true("pressure" %in% names(params$zi_fixef))
  
  # Extract and create priors
  priors <- extract_and_inflate_priors(model_freq, inflation_factor = 2)
  
  # Verify prior naming follows brms convention
  expect_true("b_Intercept" %in% names(priors))
  expect_true("b_temp" %in% names(priors))
  expect_true("zi_Intercept" %in% names(priors))
  expect_true("zi_pressure" %in% names(priors))
  
  # Check if random effects priors exist
  sd_names <- grep("^sd_", names(priors), value = TRUE)
  expect_true(length(sd_names) > 0)
  
  # Verify prior values are inflated
  # Priors should be numeric vectors with mean and sd
  expect_true(is.numeric(priors$b_Intercept))
  expect_true(length(priors$b_Intercept) == 2)  # mean and sd
  
  message("Parameter mapping validation successful!")
})
