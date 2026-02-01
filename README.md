# EmpiricalBayesianWorkflow

General pipeline for empirical Bayesian inference using hierarchical models, with a focus on zero-inflated beta mixed models and conformal prediction.

## Overview

This repository implements a complete workflow for fitting zero-inflated beta mixed models with random effects using both frequentist and Bayesian approaches, culminating in uncertainty-aware conformal prediction intervals.

### Key Features

1. **Frequentist Model Fitting**: Uses `glmmTMB` to fit zero-inflated beta mixed models with random intercepts and slopes
2. **Prior Elicitation**: Automatically extracts and inflates parameter estimates to construct informative priors
3. **Parameter Renaming**: Carefully maps parameters to match `brms` conventions
4. **Bayesian Refitting**: Refits the model in `brms` using Stan with informative priors
5. **Conformal Prediction**: Implements both split conformal and Jackknife+ methods for calibrated prediction intervals
6. **Uncertainty Quantification**: Assesses whether new observations fall within prediction intervals

## Installation

### Prerequisites

The workflow requires R (>= 4.0.0) and the following packages:

```r
install.packages(c("glmmTMB", "brms", "dplyr", "ggplot2", "tidyr"))
```

Note: `brms` requires a working Stan installation. See [brms installation guide](https://paul-buerkner.github.io/brms/) for details.

## Usage

### Basic Example

```r
# Source the workflow functions
source("zero_inflated_beta_workflow.R")

# Run the complete workflow with split conformal prediction
results <- run_zibeta_workflow(use_jackknife = FALSE)

# Access components
frequentist_model <- results$model_freq
bayesian_model <- results$model_brms
conformal_results <- results$conformal_results

# View coverage statistics
print(conformal_results$coverage)
```

### Running the Example Script

```bash
Rscript example.R
```

### Advanced Usage

For Jackknife+ conformal prediction (more computationally intensive but potentially more powerful):

```r
results <- run_zibeta_workflow(use_jackknife = TRUE)
```

## Workflow Details

### 1. Data Generation

The workflow includes synthetic data generation for demonstration:
- Zero-inflated beta distributed responses
- Random effects for subjects (intercepts and slopes)
- Continuous predictor variable

### 2. Frequentist Model (glmmTMB)

```r
model <- glmmTMB(
  y ~ x + (1 + x | subject_id),
  data = train_data,
  family = beta_family(),
  ziformula = ~ 1
)
```

### 3. Prior Construction

Parameters are extracted and inflated:
- Fixed effects: Normal priors with mean = estimate, SD = inflated
- Random effects: Half-normal priors for standard deviations
- Precision parameter: Gamma prior
- All parameters renamed to match brms conventions

### 4. Bayesian Model (brms)

```r
model_brms <- brm(
  bf(y ~ x + (1 + x | subject_id), zi ~ 1),
  data = train_data,
  family = zero_inflated_beta(),
  prior = informative_priors
)
```

### 5. Conformal Prediction

Two methods available:

#### Split Conformal
- Fast and efficient
- Uses separate calibration set
- Provides valid prediction intervals

#### Jackknife+
- More computationally intensive
- Leave-one-out cross-validation
- Can provide tighter intervals

## Output

The workflow produces:

1. **Fitted Models**: Both frequentist and Bayesian model objects
2. **Prediction Intervals**: Lower and upper bounds for test observations
3. **Coverage Statistics**: Empirical coverage vs. target coverage
4. **Visualizations**: Plots showing predictions with uncertainty intervals

## Mathematical Details

### Zero-Inflated Beta Distribution

The response $y$ follows:

$$
y \sim 
\begin{cases}
0 & \text{with probability } \pi \\
\text{Beta}(\mu\phi, (1-\mu)\phi) & \text{with probability } 1-\pi
\end{cases}
$$

### Mixed Model Structure

- **Conditional mean**: $\text{logit}(\mu_{ij}) = \beta_0 + \beta_1 x_{ij} + b_{0i} + b_{1i}x_{ij}$
- **Zero-inflation**: $\text{logit}(\pi) = \gamma_0$
- **Random effects**: $(b_{0i}, b_{1i})^T \sim N(0, \Sigma)$

### Conformal Prediction

For split conformal with miscoverage level $\alpha$:

1. Compute nonconformity scores on calibration set: $s_i = |y_i - \hat{y}_i|$
2. Find quantile: $q = \text{quantile}(s, (n+1)(1-\alpha)/n)$
3. Prediction interval: $[\hat{y} - q, \hat{y} + q]$ (bounded to $[0,1]$)

## License

See LICENSE file for details.

## References

- Shafer, G., & Vovk, V. (2008). A tutorial on conformal prediction. *Journal of Machine Learning Research*, 9, 371-421.
- Brooks-Bartlett, J. (2018). *Probabilistic programming & Bayesian methods for hackers*. Addison-Wesley.
- Bürkner, P. C. (2017). brms: An R package for Bayesian multilevel models using Stan. *Journal of Statistical Software*, 80(1), 1-28.
