

# library(EmpiricalBayesianWorkflow)
devtools::load_all(".")
library(glmmTMB)
library(brms)

# Mock model object to avoid fitting time
mock_glmmTMB <- function(family_name, disp_fixef) {
    model <- list()
    class(model) <- "glmmTMB"
    
    # Mock extract_parameters result structure
    # We can test extract_and_inflate_priors by mocking the input it expects 'model' to provide
    # BUT extract_and_inflate_priors calls extract_parameters(model).
    # extract_parameters.freq_glmmTMB calls glmmTMB::fixef(model) etc.
    # So we need a real model or we need to mock extract_parameters explicitly?
    # It's cleaner to just fit a simple model.
}

cat("Fitting reference distributional models...\n")
set.seed(123)
N <- 50
x <- rnorm(N)
data <- data.frame(y=rnorm(N), x=x)

# Gaussian with dispformula
# sigma ~ x
fit_gauss <- glmmTMB(y ~ x, dispformula = ~ x, data = data, family = gaussian())
# Manually add attributes expected by extract_parameters.freq_glmmTMB
class(fit_gauss) <- c("freq_glmmTMB", class(fit_gauss))
attr(fit_gauss, "family_config") <- get_family_config("gaussian")

# Generate priors
cat("Constructing priors for Gaussian distributional...\n")
priors <- extract_and_inflate_priors(fit_gauss)
cat("Printing raw priors list:\n")
print(priors)
brms_priors <- create_brms_priors(priors)
print(brms_priors)

# Check for sigma parameters
cat("\nChecking for sigma_Intercept and sigma_x matches...\n")
has_sigma_int <- any(brms_priors$dpar == "sigma" & brms_priors$class == "Intercept")
has_sigma_x <- any(brms_priors$dpar == "sigma" & brms_priors$class == "b" & brms_priors$coef == "x")

if (has_sigma_int && has_sigma_x) {
    cat("Gaussian distributional (sigma) PASS\n")
} else {
    cat("Gaussian distributional (sigma) FAIL\n")
    print(brms_priors)
}

# Beta with dispformula (phi)
# phi ~ x
y_beta <- rbeta(N, 0.5, 0.5)
data_beta <- data.frame(y=y_beta, x=x)
fit_beta <- glmmTMB(y ~ x, dispformula = ~ x, data = data_beta, family = beta_family())
class(fit_beta) <- c("freq_glmmTMB", class(fit_beta))
attr(fit_beta, "family_config") <- get_family_config("beta")

cat("\nConstructing priors for Beta distributional...\n")
priors_beta <- extract_and_inflate_priors(fit_beta)
brms_priors_beta <- create_brms_priors(priors_beta)
# print(brms_priors_beta)

has_phi_int <- any(brms_priors_beta$dpar == "phi" & brms_priors_beta$class == "Intercept")
has_phi_x <- any(brms_priors_beta$dpar == "phi" & brms_priors_beta$class == "b" & brms_priors_beta$coef == "x")

if (has_phi_int && has_phi_x) {
    cat("Beta distributional (phi) PASS\n")
} else {
    cat("Beta distributional (phi) FAIL\n")
    print(brms_priors_beta)
}
