# Tests for family mapping functionality

test_that("get_supported_families returns correct structure", {
  families <- get_supported_families()
  
  expect_type(families, "list")
  expect_true(length(families) > 0)
  expect_true("zero_inflated_beta" %in% names(families))
  expect_true("beta" %in% names(families))
  expect_true("gamma" %in% names(families))
})

test_that("get_family_config retrieves correct configuration", {
  config <- get_family_config("zero_inflated_beta")
  
  expect_type(config, "list")
  expect_equal(config$name, "zero_inflated_beta")
  expect_true(config$has_zi)
  expect_false(config$has_hurdle)
  expect_true("phi" %in% config$params)
})

test_that("get_family_config fails for unsupported family", {
  expect_error(get_family_config("unsupported_family"))
})

test_that("map_parameter_name works correctly", {
  config <- get_family_config("zero_inflated_beta")
  
  # Test intercept mapping
  expect_equal(map_parameter_name("(Intercept)", config, "cond"), "Intercept")
  expect_equal(map_parameter_name("(Intercept)", config, "zi"), "zi_Intercept")
  
  # Test regular parameter mapping
  expect_equal(map_parameter_name("x", config, "cond"), "x")
  expect_equal(map_parameter_name("x", config, "zi"), "zi_x")
})

test_that("create_brms_family works for different families", {
  skip_if_not_installed("brms")
  
  # Zero-inflated beta
  config_zib <- get_family_config("zero_inflated_beta")
  family_zib <- create_brms_family(config_zib)
  expect_s3_class(family_zib, "brmsfamily")
  
  # Beta
  config_beta <- get_family_config("beta")
  family_beta <- create_brms_family(config_beta)
  expect_s3_class(family_beta, "brmsfamily")
  
  # Gamma
  config_gamma <- get_family_config("gamma")
  family_gamma <- create_brms_family(config_gamma)
  expect_s3_class(family_gamma, "brmsfamily")
})

test_that("create_glmmTMB_family works for different families", {
  skip_if_not_installed("glmmTMB")
  
  # Zero-inflated beta (uses beta_family in glmmTMB)
  config_zib <- get_family_config("zero_inflated_beta")
  family_zib <- create_glmmTMB_family(config_zib)
  expect_s3_class(family_zib, "family")
  
  # Gamma
  config_gamma <- get_family_config("gamma")
  family_gamma <- create_glmmTMB_family(config_gamma)
  expect_s3_class(family_gamma, "family")
})
