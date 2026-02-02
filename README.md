# EmpiricalBayesianWorkflow

R package for empirical Bayesian inference using hierarchical mixed models with conformal prediction. Supports multiple model families with exact parameter mapping between frequentist (glmmTMB) and Bayesian (brms/cmdstanr) frameworks.

## Overview

This package implements a complete workflow for fitting mixed models with random effects using both frequentist and Bayesian approaches, culminating in uncertainty-aware conformal prediction intervals.

### Key Features

1. **Flexible Family Support**: Beta, zero-inflated beta, gamma, Poisson, negative binomial, and their zero-inflated variants
2. **Agnostic Frequentist Fitting**: Supports multiple backends (glmmTMB preferred)
3. **Exact Parameter Mapping**: Automatically maps parameters from glmmTMB → brms with family-specific handling
4. **Prior Elicitation**: Extracts and inflates parameter estimates to construct informative priors
5. **Bayesian Fitting**: Uses cmdstanr backend for efficient Stan sampling
6. **Conformal Prediction**: Implements both split conformal and Jackknife+ methods
7. **Comprehensive Testing**: GitHub Actions CI/CD with R-CMD-check across platforms

## Installation

### Prerequisites

The package requires R (>= 4.0.0) and the following packages:

```r
# Required for frequentist fitting
install.packages("glmmTMB")

# Required for Bayesian fitting
install.packages("brms")

# Required for cmdstanr backend (recommended)
install.packages("cmdstanr", repos = c("https://mc-stan.org/r-packages/", getOption("repos")))
cmdstanr::install_cmdstan()
```

### Install from GitHub

```r
# Using devtools
devtools::install_github("GWMcElfresh/EmpiricalBayesianWorkflow")

# Using remotes
remotes::install_github("GWMcElfresh/EmpiricalBayesianWorkflow")
```

## Quick Start

### Basic Usage with Zero-Inflated Beta

```r
library(EmpiricalBayesianWorkflow)

# Run complete workflow with synthetic data
results <- run_workflow(
  family = "zero_inflated_beta",
  n_subjects = 30,
  n_obs_per_subject = 20,
  verbose = TRUE
)

# View results
print(results)
print(results$conformal_results)
```

### Using Different Families

```r
# Gamma model
results_gamma <- run_workflow(
  family = "gamma",
  n_subjects = 25,
  verbose = TRUE
)

# Zero-inflated Poisson
results_zip <- run_workflow(
  family = "zero_inflated_poisson",
  verbose = TRUE
)

# Regular beta (no zero-inflation)
results_beta <- run_workflow(
  family = "beta",
  verbose = TRUE
)
```

### Using Your Own Data

```r
# Your data should have:
# - Response variable
# - Predictor(s)
# - Grouping variable for random effects

results <- run_workflow(
  data = my_data,
  formula = response ~ predictor1 + predictor2 + (1 + predictor1 | group_id),
  family = "zero_inflated_beta"
)
```

### Step-by-Step Workflow

```r
# 1. Generate or load data
data <- generate_zibeta_data(n_subjects = 20, n_obs_per_subject = 15)

# 2. Fit frequentist model
model_freq <- fit_frequentist_model(
  y ~ x + (1 + x | subject_id),
  data = data,
  family = "zero_inflated_beta"
)

# 3. Extract and create priors
priors <- extract_and_inflate_priors(model_freq)
brms_priors <- create_brms_priors(priors)

# 4. Fit Bayesian model with cmdstanr
model_bayes <- fit_bayesian_model(
  y ~ x + (1 + x | subject_id),
  data = data,
  family = "zero_inflated_beta",
  prior = brms_priors,
  backend = "cmdstanr"
)

# 5. Perform conformal prediction
cp_results <- conformal_prediction_split(
  model = model_bayes,
  calibration_data = calib_data,
  test_data = test_data,
  alpha = 0.1  # 90% coverage
)
```

## Detailed Walkthroughs

### Walkthrough 1: Zero-Inflated Beta with GasolineYield Dataset

This walkthrough demonstrates the complete empirical Bayesian workflow using the `GasolineYield` dataset from the `betareg` package. We model gasoline refinery yield as a zero-inflated beta distribution with:
- **Mean (μ)**: Conditional on temperature and batch (random effect)
- **Precision (φ)**: Conditional on temperature (distributional regression)
- **Zero-Inflation (zi)**: Conditional on pressure

> **Note**: The GasolineYield dataset doesn't naturally contain zeros. We artificially introduce structural zeros to demonstrate zero-inflation modeling.

#### Step 1: Data Preparation

```r
library(EmpiricalBayesianWorkflow)
library(betareg)

# Load the GasolineYield dataset
data("GasolineYield", package = "betareg")

# Prepare data for zero-inflated modeling
set.seed(12345)
gasoline <- GasolineYield
gasoline$yield_zi <- gasoline$yield

# Introduce structural zeros (bottom 10% of yields)
threshold <- quantile(gasoline$yield, 0.10)
gasoline$yield_zi[gasoline$yield < threshold] <- 0

# Verify the data structure
summary(gasoline)
# Response: yield_zi (0 to 1, with exact zeros)
# Predictors: temp (temperature), pressure, batch (10 levels)
```

#### Step 2: Fit Frequentist Model with glmmTMB

```r
library(glmmTMB)

# Define model formulas
formula_mean <- yield_zi ~ temp + (1 | batch)
formula_zi <- ~ pressure
formula_disp <- ~ temp  # Precision parameter phi

# Fit zero-inflated beta model
model_freq <- glmmTMB(
  formula = formula_mean,
  data = gasoline,
  family = beta_family(),
  ziformula = formula_zi,
  dispformula = formula_disp,
  REML = FALSE
)

# Inspect the model
summary(model_freq)

# Extract parameter estimates
fixef(model_freq)
# $cond: Conditional model (mean)
#   (Intercept)  temp
# $zi: Zero-inflation model
#   (Intercept)  pressure
# $disp: Dispersion model (precision)
#   (Intercept)  temp
```

**Key Insights:**
- **Conditional model**: Temperature effect on mean yield
- **ZI model**: Pressure affects probability of structural zeros
- **Dispersion model**: Temperature affects precision (heteroscedasticity)

#### Step 3: Extract and Inflate Priors

```r
# Setup family configuration for parameter mapping
family_config <- list(
  name = "zero_inflated_beta",
  has_zi = TRUE,
  params = c("phi"),
  param_classes = list(phi = "phi")
)

attr(model_freq, "family_config") <- family_config
class(model_freq) <- c("freq_glmmTMB_zero_inflated_beta", 
                       "freq_glmmTMB", 
                       class(model_freq))

# Extract parameter estimates and inflate for priors
priors <- extract_and_inflate_priors(
  model_freq, 
  inflation_factor = 2.5  # Conservative priors (2.5× original SE)
)

# Inspect prior structure
names(priors)
# [1] "b_Intercept"          "b_temp"              
# [3] "zi_Intercept"         "zi_pressure"         
# [5] "sd_batch__Intercept"  "phi_Intercept"       
# [7] "phi_temp"

# Example: Prior for temperature effect on mean
priors$b_temp
# $location: Point estimate from glmmTMB
# $scale: Inflated standard error (2.5×)
# $class: "b"
# $coef: "temp"
```

**Parameter Mapping (glmmTMB → brms):**
- `fixef()$cond["temp"]` → `b_temp`
- `fixef()$zi["pressure"]` → `zi_pressure`
- `fixef()$disp["temp"]` → `phi_temp` (brms uses `b_phi_temp`)
- `VarCorr()$cond$batch` → `sd_batch__Intercept`

#### Step 4: Fit Bayesian Model with brms

```r
library(brms)

# Convert priors to brms format
brms_priors <- create_brms_priors(priors)

# Define distributional formula for brms
# bf() allows multiple submodels
brms_formula <- bf(
  yield_zi ~ temp + (1 | batch),    # Mean model
  phi ~ temp,                        # Precision model
  zi ~ pressure                      # Zero-inflation model
)

# Fit Bayesian model with cmdstanr backend
model_bayes <- brm(
  formula = brms_formula,
  data = gasoline,
  family = zero_inflated_beta(),
  prior = brms_priors,
  chains = 4,
  iter = 2000,
  warmup = 1000,
  cores = 4,
  backend = "cmdstanr",
  seed = 12345
)

# Inspect Bayesian model
summary(model_bayes)
plot(model_bayes)

# Compare frequentist vs Bayesian estimates
fixef(model_freq)$cond  # glmmTMB
fixef(model_bayes)       # brms (posterior means)
```

**Interpreting Results:**
- **Fixed effects (`b_`)**: Posterior distributions shrunk toward priors
- **Precision effects (`b_phi_`)**: Log-scale precision parameters
- **ZI effects (`b_zi_`)**: Logit-scale zero-inflation probabilities
- **Random effects (`sd_batch`)**: Batch-to-batch variability

#### Step 5: Split Data for Conformal Prediction

```r
# Split data: 50% train, 25% calibration, 25% test
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
```

#### Step 6: Refit Models on Training Data

```r
# Refit frequentist model on training data only
model_freq_train <- glmmTMB(
  formula = formula_mean,
  data = train_data,
  family = beta_family(),
  ziformula = formula_zi,
  dispformula = formula_disp,
  REML = FALSE
)

# Extract priors from training model
attr(model_freq_train, "family_config") <- family_config
class(model_freq_train) <- c("freq_glmmTMB_zero_inflated_beta", 
                             "freq_glmmTMB", 
                             class(model_freq_train))

priors_train <- extract_and_inflate_priors(model_freq_train, inflation_factor = 2.5)
brms_priors_train <- create_brms_priors(priors_train)

# Refit Bayesian model on training data
model_bayes_train <- brm(
  formula = brms_formula,
  data = train_data,
  family = zero_inflated_beta(),
  prior = brms_priors_train,
  chains = 4,
  iter = 2000,
  warmup = 1000,
  cores = 4,
  backend = "cmdstanr",
  seed = 12345
)
```

#### Step 7: Perform Split Conformal Prediction

```r
# Compute prediction intervals with coverage guarantees
cp_results <- conformal_prediction_split(
  model = model_bayes_train,
  calibration_data = calib_data,
  test_data = test_data,
  alpha = 0.1,  # Target 90% coverage
  response_var = "yield_zi"
)

# Inspect results
print(cp_results)
# Coverage: 0.91 (target: 0.90)
# Mean interval width: 0.143

# Visualize prediction intervals
library(ggplot2)
results_df <- data.frame(
  observed = test_data$yield_zi,
  predicted = cp_results$predictions,
  lower = cp_results$lower,
  upper = cp_results$upper,
  covered = (test_data$yield_zi >= cp_results$lower) & 
            (test_data$yield_zi <= cp_results$upper)
)

ggplot(results_df, aes(x = 1:nrow(results_df))) +
  geom_point(aes(y = observed, color = covered), size = 2) +
  geom_point(aes(y = predicted), shape = 4, size = 3) +
  geom_errorbar(aes(ymin = lower, ymax = upper), alpha = 0.3) +
  scale_color_manual(values = c("red", "black")) +
  labs(x = "Test observation", y = "Yield", 
       title = "Conformal Prediction Intervals for Gasoline Yield") +
  theme_minimal()
```

**Coverage Guarantee:**
The split conformal method provides finite-sample validity: under the exchangeability assumption, coverage ≥ 1 - α with high probability. The intervals are distribution-free and require no assumptions about the model family.

#### Step 8: Complete Workflow (One Function)

Alternatively, use the high-level `run_workflow()` function:

```r
# Note: run_workflow currently doesn't support dispformula
# For now, we demonstrate with simpler formula
results <- run_workflow(
  data = gasoline,
  formula = yield_zi ~ temp + (1 | batch),
  family = "zero_inflated_beta",
  zi_formula = ~ pressure,
  train_prop = 0.5,
  calib_prop = 0.25,
  inflation_factor = 2.5,
  conformal_method = "split",
  alpha = 0.1,
  backend = "cmdstanr",
  chains = 4,
  iter = 2000,
  verbose = TRUE,
  seed = 12345
)

# Access components
print(results$model_freq)
print(results$model_bayes)
print(results$conformal_results)
```

### Walkthrough 2: Comparison Across Model Families

Compare zero-inflated beta with regular beta and gamma models:

```r
# Prepare data without zeros for regular beta
gasoline_no_zero <- gasoline[gasoline$yield_zi > 0, ]

# Fit beta model (no zero-inflation)
results_beta <- run_workflow(
  data = gasoline_no_zero,
  formula = yield_zi ~ temp + (1 | batch),
  family = "beta",
  verbose = TRUE
)

# Fit gamma model (for right-skewed positive data)
results_gamma <- run_workflow(
  data = gasoline_no_zero,
  formula = yield_zi ~ temp + (1 | batch),
  family = "gamma",
  verbose = TRUE
)

# Compare coverage and interval widths
cat("Zero-Inflated Beta Coverage:", results$conformal_results$coverage, "\n")
cat("Regular Beta Coverage:", results_beta$conformal_results$coverage, "\n")
cat("Gamma Coverage:", results_gamma$conformal_results$coverage, "\n")

cat("\nMean Interval Width:\n")
cat("ZI Beta:", results$conformal_results$mean_width, "\n")
cat("Beta:", results_beta$conformal_results$mean_width, "\n")
cat("Gamma:", results_gamma$conformal_results$mean_width, "\n")
```

### Walkthrough 3: Jackknife+ Conformal Prediction

The Jackknife+ method provides more powerful prediction intervals than split conformal by using leave-one-out cross-validation. It's particularly effective for small datasets.

```r
library(EmpiricalBayesianWorkflow)
library(betareg)

# Load and prepare data
data("GasolineYield", package = "betareg")
set.seed(789)
gasoline <- GasolineYield
gasoline$yield_zi <- gasoline$yield
threshold <- quantile(gasoline$yield, 0.10)
gasoline$yield_zi[gasoline$yield < threshold] <- 0

# Split into train and test (no calibration needed for Jackknife+)
n <- nrow(gasoline)
train_size <- floor(0.7 * n)
train_idx <- sample(1:n, train_size)
test_idx <- setdiff(1:n, train_idx)

train_data <- gasoline[train_idx, ]
test_data <- gasoline[test_idx, ]

# Define fitting function for Jackknife+
fit_fn <- function(data) {
  # Fit frequentist model
  model_freq <- glmmTMB::glmmTMB(
    formula = yield_zi ~ temp + (1 | batch),
    data = data,
    family = glmmTMB::beta_family(),
    ziformula = ~ pressure,
    REML = FALSE
  )
  
  # Setup for Bayesian fitting
  family_config <- list(
    name = "zero_inflated_beta",
    has_zi = TRUE,
    params = c("phi"),
    param_classes = list(phi = "phi")
  )
  
  attr(model_freq, "family_config") <- family_config
  class(model_freq) <- c("freq_glmmTMB_zero_inflated_beta", 
                         "freq_glmmTMB", 
                         class(model_freq))
  
  # Extract and inflate priors
  priors <- extract_and_inflate_priors(model_freq, inflation_factor = 2)
  brms_priors <- create_brms_priors(priors)
  
  # Fit Bayesian model
  model_bayes <- brms::brm(
    formula = brms::bf(yield_zi ~ temp + (1 | batch), zi ~ pressure),
    data = data,
    family = brms::zero_inflated_beta(),
    prior = brms_priors,
    chains = 2,
    iter = 1000,
    warmup = 500,
    cores = 2,
    backend = "cmdstanr",
    refresh = 0,
    silent = 2
  )
  
  return(model_bayes)
}

# Define prediction function
predict_fn <- function(model, newdata) {
  # Generate posterior predictions
  pred_samples <- posterior_predict_model(
    model, 
    newdata = newdata,
    allow_new_levels = TRUE
  )
  
  # Return posterior means
  return(colMeans(pred_samples))
}

# Perform Jackknife+ conformal prediction
cp_jackknife <- conformal_prediction_jackknife(
  formula = yield_zi ~ temp + (1 | batch),
  train_data = train_data,
  test_data = test_data,
  alpha = 0.1,  # 90% coverage
  fit_function = fit_fn,
  predict_function = predict_fn,
  response_var = "yield_zi"
)

# Compare with split conformal
model_full <- fit_fn(train_data)
split_size <- floor(0.5 * nrow(train_data))
split_calib_data <- train_data[(split_size + 1):nrow(train_data), ]
split_train_data <- train_data[1:split_size, ]

model_split <- fit_fn(split_train_data)
cp_split <- conformal_prediction_split(
  model = model_split,
  calibration_data = split_calib_data,
  test_data = test_data,
  alpha = 0.1,
  response_var = "yield_zi"
)

# Compare methods
cat("Jackknife+ Coverage:", cp_jackknife$coverage, "\n")
cat("Split Conformal Coverage:", cp_split$coverage, "\n")
cat("\nJackknife+ Mean Width:", cp_jackknife$mean_width, "\n")
cat("Split Conformal Mean Width:", cp_split$mean_width, "\n")

# Jackknife+ typically provides tighter intervals with similar coverage
```

**Trade-offs:**
- **Split Conformal**: Fast, simple, requires separate calibration set
- **Jackknife+**: More powerful (tighter intervals), slower (requires n+1 model fits), no calibration set needed

### Walkthrough 4: Custom Prior Specification

Fine-tune priors for domain knowledge integration:

```r
library(EmpiricalBayesianWorkflow)
library(betareg)

# Load data
data("GasolineYield", package = "betareg")
set.seed(456)
gasoline <- GasolineYield
gasoline$yield_zi <- gasoline$yield
threshold <- quantile(gasoline$yield, 0.10)
gasoline$yield_zi[gasoline$yield < threshold] <- 0

# Step 1: Fit frequentist model
model_freq <- glmmTMB::glmmTMB(
  formula = yield_zi ~ temp + (1 | batch),
  data = gasoline,
  family = glmmTMB::beta_family(),
  ziformula = ~ pressure,
  REML = FALSE
)

# Setup family config
family_config <- list(
  name = "zero_inflated_beta",
  has_zi = TRUE,
  params = c("phi"),
  param_classes = list(phi = "phi")
)

attr(model_freq, "family_config") <- family_config
class(model_freq) <- c("freq_glmmTMB_zero_inflated_beta", 
                       "freq_glmmTMB", 
                       class(model_freq))

# Step 2: Extract automatic priors
auto_priors <- extract_and_inflate_priors(model_freq, inflation_factor = 2)

# Step 3: Customize specific priors based on domain knowledge
# Example: We have strong prior belief about temperature effect
auto_priors$b_temp <- list(
  prior = "normal",
  class = "b",
  coef = "temp",
  location = 0.01,  # Expected small positive effect
  scale = 0.005     # High confidence
)

# Example: Skeptical about zero-inflation (expect few zeros)
auto_priors$zi_Intercept <- list(
  prior = "normal",
  class = "zi",
  coef = "Intercept",
  location = -3,    # Low baseline ZI probability (logit scale)
  scale = 1         # Moderate uncertainty
)

# Step 4: Create brms priors
custom_brms_priors <- create_brms_priors(auto_priors)

# Step 5: Fit Bayesian model with custom priors
model_bayes_custom <- brms::brm(
  formula = brms::bf(yield_zi ~ temp + (1 | batch), zi ~ pressure),
  data = gasoline,
  family = brms::zero_inflated_beta(),
  prior = custom_brms_priors,
  chains = 4,
  iter = 2000,
  backend = "cmdstanr"
)

# Compare with automatic priors
auto_brms_priors <- create_brms_priors(
  extract_and_inflate_priors(model_freq, inflation_factor = 2)
)

model_bayes_auto <- brms::brm(
  formula = brms::bf(yield_zi ~ temp + (1 | batch), zi ~ pressure),
  data = gasoline,
  family = brms::zero_inflated_beta(),
  prior = auto_brms_priors,
  chains = 4,
  iter = 2000,
  backend = "cmdstanr"
)

# Compare posterior distributions
library(bayesplot)
mcmc_combo(
  model_bayes_custom, 
  pars = c("b_temp", "b_zi_Intercept"),
  combo = c("dens_overlay", "trace")
)

mcmc_combo(
  model_bayes_auto,
  pars = c("b_temp", "b_zi_Intercept"),
  combo = c("dens_overlay", "trace")
)
```

**When to Use Custom Priors:**
- **Regularization**: Shrink coefficients toward zero for feature selection
- **Domain knowledge**: Incorporate expert beliefs about parameter ranges
- **Hierarchical information**: Use priors from related studies
- **Computational stability**: Stronger priors can improve convergence

### Walkthrough 5: Model Diagnostics and Validation

Comprehensive diagnostics for the empirical Bayesian workflow:

```r
library(EmpiricalBayesianWorkflow)
library(betareg)
library(bayesplot)
library(ggplot2)

# Load and prepare data
data("GasolineYield", package = "betareg")
set.seed(999)
gasoline <- GasolineYield
gasoline$yield_zi <- gasoline$yield
threshold <- quantile(gasoline$yield, 0.10)
gasoline$yield_zi[gasoline$yield < threshold] <- 0

# Run complete workflow
results <- run_workflow(
  data = gasoline,
  formula = yield_zi ~ temp + (1 | batch),
  family = "zero_inflated_beta",
  zi_formula = ~ pressure,
  train_prop = 0.5,
  calib_prop = 0.25,
  verbose = TRUE,
  chains = 4,
  iter = 2000
)

# 1. Frequentist Model Diagnostics
model_freq <- results$model_freq
summary(model_freq)

# 2. Bayesian Model Diagnostics
model_bayes <- results$model_bayes

# MCMC convergence: Rhat should be < 1.01
rhat_vals <- brms::rhat(model_bayes)
max_rhat <- max(rhat_vals, na.rm = TRUE)
cat("Maximum Rhat:", max_rhat, "\n")

# Effective sample size
neff_vals <- brms::neff_ratio(model_bayes)
min_neff <- min(neff_vals, na.rm = TRUE)
cat("Minimum ESS ratio:", min_neff, "\n")

# Posterior predictive check
pp_check(model_bayes, ndraws = 100)

# 3. Conformal Prediction Diagnostics
cp <- results$conformal_results
test_data <- results$data$test
observed <- test_data$yield_zi

# Coverage verification
cat("Coverage:", cp$coverage, "( target: 0.90)\n")
cat("Mean interval width:", cp$mean_width, "\n")

# Miscoverage analysis
miscovered <- !((observed >= cp$lower) & (observed <= cp$upper))
cat("Miscoverage rate:", mean(miscovered), "\n")
```

## Supported Model Families
|--------|-------------|---------|------|------------|
| `beta` | Beta distribution for (0,1) | ✓ | ✓ | phi |
| `zero_inflated_beta` | ZI Beta for [0,1) | ✓ | ✓ | phi, zi |
| `gamma` | Gamma for positive continuous | ✓ | ✓ | shape |
| `binomial` | Binary/proportion data | ✓ | ✓ | - |
| `poisson` | Count data | ✓ | ✓ | - |
| `zero_inflated_poisson` | ZI Poisson | ✓ | ✓ | zi |
| `nbinom2` | Negative binomial (NB2) | ✓ | ✓ | shape |
| `zero_inflated_nbinom2` | ZI negative binomial | ✓ | ✓ | shape, zi |

## Parameter Mapping

The package ensures exact parameter correspondence between glmmTMB and brms:

### Fixed Effects (Conditional Model)
- glmmTMB: `fixef()$cond["(Intercept)"]` → brms: `b_Intercept`
- glmmTMB: `fixef()$cond["x"]` → brms: `b_x`

### Zero-Inflation Effects
- glmmTMB: `fixef()$zi["(Intercept)"]` → brms: `zi_Intercept`  
- glmmTMB: `fixef()$zi["x"]` → brms: `zi_x`

### Random Effects
- glmmTMB: `VarCorr()$cond$group["sd"]` → brms: `sd_group__Intercept`

### Family-Specific Parameters
- Beta/ZI Beta: `sigma()` → `phi` (precision)
- Gamma: shape parameter
- Negative Binomial: shape/dispersion parameter

## Conformal Prediction

### Split Conformal (Fast)
```r
cp_split <- conformal_prediction_split(
  model = bayesian_model,
  calibration_data = calib,
  test_data = test,
  alpha = 0.1  # 90% coverage
)
```

### Jackknife+ (More Powerful, Slower)
```r
cp_jackknife <- conformal_prediction_jackknife(
  formula = model_formula,
  train_data = train,
  test_data = test,
  alpha = 0.1,
  fit_function = fit_fn,
  predict_function = predict_fn
)
```

## Testing and CI/CD

The package includes comprehensive testing:

```r
# Run tests locally
devtools::test()

# Check package
devtools::check()
```

### GitHub Actions

Automated testing runs on:
- **Platforms**: Ubuntu (latest, devel), macOS, Windows
- **Triggers**: Push to main/master, pull requests, weekly schedule
- **Checks**: R CMD check, test coverage, documentation

## Documentation

- **README.md**: This file (quick start and overview)
- **QUICKSTART.md**: Detailed getting started guide
- **TECHNICAL.md**: In-depth technical documentation
- **Package help**: `?EmpiricalBayesianWorkflow` after installation

## Advanced Usage

### Custom Prior Specification

```r
# Extract default priors
auto_priors <- extract_and_inflate_priors(freq_model, inflation_factor = 2)

# Modify specific priors
custom_priors <- create_brms_priors(auto_priors)
# Then manually add or modify specific priors as needed
```

### Using Different Backends

```r
# cmdstanr (default, recommended)
model <- fit_bayesian_model(formula, data, backend = "cmdstanr")

# rstan (fallback)
model <- fit_bayesian_model(formula, data, backend = "rstan")
```

## Contributing

Contributions are welcome! Please:
1. Fork the repository
2. Create a feature branch
3. Add tests for new functionality
4. Ensure R CMD check passes
5. Submit a pull request

## References

1. **glmmTMB**: Brooks et al. (2017). *The R Journal*, 9(2), 378-400.
2. **brms**: Bürkner (2017). *Journal of Statistical Software*, 80(1), 1-28.
3. **Conformal Prediction**: Vovk et al. (2005). *Algorithmic Learning in a Random World*. Springer.
4. **Jackknife+**: Barber et al. (2021). *The Annals of Statistics*, 49(1), 486-507.
5. **Empirical Bayes**: Efron (2012). *Large-Scale Inference*. Cambridge University Press.

## License

See LICENSE file for details.
