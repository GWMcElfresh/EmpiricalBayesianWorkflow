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

## Supported Families

| Family | Description | glmmTMB | brms | Parameters |
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
