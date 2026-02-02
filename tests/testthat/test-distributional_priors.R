

library(testthat)
devtools::load_all(".")
# library(EmpiricalBayesianWorkflow)
library(glmmTMB)
library(brms)

test_that("Distributional priors are correctly constructed", {
  skip_if_not_installed("glmmTMB")
  skip_if_not_installed("brms")
  
  set.seed(123)
  N <- 50
  x <- rnorm(N)
  data <- data.frame(y=rnorm(N), x=x)
  
  # 1. Gaussian with dispformula (sigma ~ x)
  fit_gauss <- glmmTMB(y ~ x, dispformula = ~ x, data = data, family = gaussian())
  # Ensure attributes for method dispatch if not automatically handled by package in test env
  # (When running checks, S3 dispatch usually works if package is loaded)
  if (!inherits(fit_gauss, "freq_glmmTMB")) {
      class(fit_gauss) <- c("freq_glmmTMB", class(fit_gauss))
      attr(fit_gauss, "family_config") <- get_family_config("gaussian")
  }

  priors <- extract_and_inflate_priors(fit_gauss)
  brms_priors <- create_brms_priors(priors)
  
  # Check structure
  expect_true(any(brms_priors$dpar == "sigma" & brms_priors$class == "Intercept"))
  expect_true(any(brms_priors$dpar == "sigma" & brms_priors$class == "b" & brms_priors$coef == "x"))
  
  # 2. Beta with dispformula (phi ~ x)
  y_beta <- rbeta(N, 0.5, 0.5)
  data_beta <- data.frame(y=y_beta, x=x)
  fit_beta <- glmmTMB(y ~ x, dispformula = ~ x, data = data_beta, family = beta_family())
  
  if (!inherits(fit_beta, "freq_glmmTMB")) {
      class(fit_beta) <- c("freq_glmmTMB", class(fit_beta))
      attr(fit_beta, "family_config") <- get_family_config("beta")
  }
  
  priors_beta <- extract_and_inflate_priors(fit_beta)
  brms_priors_beta <- create_brms_priors(priors_beta)
  
  expect_true(any(brms_priors_beta$dpar == "phi" & brms_priors_beta$class == "Intercept"))
  expect_true(any(brms_priors_beta$dpar == "phi" & brms_priors_beta$class == "b" & brms_priors_beta$coef == "x"))
  
  # 3. ZI Beta with dispformula (phi ~ x)
  # The mapping should be identical to regular Beta
  # We test to ensure the family config allows it
  fit_zibeta <- glmmTMB(y ~ x, dispformula = ~ x, ziformula = ~ 1, data = data_beta, family = beta_family())
  if (!inherits(fit_zibeta, "freq_glmmTMB")) {
      class(fit_zibeta) <- c("freq_glmmTMB", class(fit_zibeta))
      attr(fit_zibeta, "family_config") <- get_family_config("zero_inflated_beta")
  }

  priors_zibeta <- extract_and_inflate_priors(fit_zibeta)
  brms_priors_zibeta <- create_brms_priors(priors_zibeta)

  expect_true(any(brms_priors_zibeta$dpar == "phi" & brms_priors_zibeta$class == "Intercept"))
  expect_true(any(brms_priors_zibeta$dpar == "phi" & brms_priors_zibeta$class == "b" & brms_priors_zibeta$coef == "x"))

  # 4. ZI NBinom2 with dispformula (shape ~ x)
  y_nb <- rnbinom(N, size=2, mu=exp(x))
  data_nb <- data.frame(y=y_nb, x=x)
  fit_zinb <- glmmTMB(y ~ x, dispformula = ~ x, ziformula = ~ 1, data = data_nb, family = nbinom2())
  if (!inherits(fit_zinb, "freq_glmmTMB")) {
      class(fit_zinb) <- c("freq_glmmTMB", class(fit_zinb))
      attr(fit_zinb, "family_config") <- get_family_config("zero_inflated_nbinom2")
  }
  
  priors_zinb <- extract_and_inflate_priors(fit_zinb)
  brms_priors_zinb <- create_brms_priors(priors_zinb)
  
  expect_true(any(brms_priors_zinb$dpar == "shape" & brms_priors_zinb$class == "Intercept"))
  expect_true(any(brms_priors_zinb$dpar == "shape" & brms_priors_zinb$class == "b" & brms_priors_zinb$coef == "x"))
})

