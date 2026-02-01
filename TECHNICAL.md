# Technical Documentation

## Zero-Inflated Beta Mixed Model Workflow

This document provides technical details about the implementation of the zero-inflated beta mixed model workflow with conformal prediction.

## Table of Contents

1. [Overview](#overview)
2. [Statistical Background](#statistical-background)
3. [Implementation Details](#implementation-details)
4. [Parameter Mapping](#parameter-mapping)
5. [Conformal Prediction](#conformal-prediction)
6. [API Reference](#api-reference)

## Overview

The workflow implements a complete pipeline for:
- **Frequentist estimation** using `glmmTMB`
- **Prior elicitation** from frequentist estimates
- **Bayesian refitting** using `brms/Stan`
- **Conformal prediction** for uncertainty quantification

## Statistical Background

### Zero-Inflated Beta Distribution

The zero-inflated beta distribution is appropriate for response variables that:
- Are bounded in the interval [0, 1]
- Can take the exact value of 0 (or 1) with non-negligible probability
- Examples: proportions, rates, probabilities with structural zeros

The model combines:
1. **Point mass at zero**: P(Y = 0) = π
2. **Beta distribution**: Y | Y > 0 ~ Beta(μφ, (1-μ)φ)

Where:
- π is the zero-inflation probability
- μ is the mean of the beta component
- φ is the precision parameter

### Mixed Effects Structure

The workflow supports random intercepts and slopes:

```
logit(μᵢⱼ) = (β₀ + b₀ᵢ) + (β₁ + b₁ᵢ)xᵢⱼ
logit(πᵢⱼ) = γ₀

(b₀ᵢ, b₁ᵢ)ᵀ ~ N(0, Σ)
```

Where:
- i indexes subjects/groups
- j indexes observations within subjects
- β are fixed effects
- b are random effects
- Σ is the random effects covariance matrix

## Implementation Details

### 1. Data Generation (`generate_zibeta_data`)

Creates synthetic data with:
- Specified number of subjects and observations per subject
- Random effects for intercepts and slopes
- Zero-inflation mechanism
- Beta-distributed non-zero responses

**Key parameters:**
- `n_subjects`: Number of groups/subjects
- `n_obs_per_subject`: Observations per group
- `zi_prob`: Probability of structural zeros

### 2. Frequentist Model Fitting

Uses `glmmTMB` with:
- `family = beta_family()`: Beta distribution for (0,1) responses
- `ziformula = ~1`: Zero-inflation model specification
- Random effects: `(1 + x | subject_id)`

**Model structure:**
```r
glmmTMB(
  y ~ x + (1 + x | subject_id),
  data = train_data,
  family = beta_family(),
  ziformula = ~ 1
)
```

### 3. Prior Extraction (`extract_and_inflate_priors`)

The function:
1. Extracts fixed effects from conditional model
2. Extracts fixed effects from zero-inflation model
3. Extracts variance components for random effects
4. Inflates variances by a specified factor (default: 2×)
5. Renames parameters to match brms conventions

**Inflation rationale:**
- Accounts for estimation uncertainty in frequentist estimates
- Provides weakly informative priors
- Prevents overconfident priors from dominating the data

### 4. Parameter Mapping

#### glmmTMB → brms Naming Convention

| glmmTMB | brms | Description |
|---------|------|-------------|
| `fixef()$cond["(Intercept)"]` | `Intercept` | Conditional model intercept |
| `fixef()$cond["x"]` | `b_x` | Conditional model coefficient for x |
| `fixef()$zi["(Intercept)"]` | `zi_Intercept` | Zero-inflation intercept |
| `VarCorr()$cond$subject_id["(Intercept)"]` | `sd_subject_id__Intercept` | Random intercept SD |
| `VarCorr()$cond$subject_id["x"]` | `sd_subject_id__x` | Random slope SD |
| `sigma()` | `phi` | Precision parameter |

**Prior specifications:**

```r
# Fixed effects (conditional model)
prior(normal(μ, σ_inflated), class = Intercept)
prior(normal(μ, σ_inflated), class = b, coef = x)

# Zero-inflation
prior(normal(μ, σ_inflated), class = Intercept, dpar = zi)

# Random effects (half-normal for SDs)
prior(normal(0, σ_inflated), class = sd, group = subject_id)

# Precision
prior(gamma(α, β), class = phi)
```

### 5. Bayesian Refitting (`brms`)

Uses Stan via brms to fit:
```r
bf(
  y ~ x + (1 + x | subject_id),
  zi ~ 1,
  family = zero_inflated_beta()
)
```

**MCMC settings:**
- Chains: 4
- Iterations: 2000 per chain
- Warmup: 1000
- Diagnostics: Rhat, ESS automatically checked

## Conformal Prediction

### Split Conformal Prediction

**Algorithm:**

1. Split data into train/calibration/test
2. Fit model on training data
3. Compute nonconformity scores on calibration set:
   ```
   sᵢ = |yᵢ - ŷᵢ|
   ```
4. Find quantile: q = quantile(s, (n+1)(1-α)/n)
5. Prediction interval: [ŷ - q, ŷ + q] ∩ [0, 1]

**Advantages:**
- Fast (single model fit)
- Exact coverage guarantee
- Simple to implement

**Implementation:**
```r
split_conformal_prediction(
  model_brms,
  train_data,
  calibration_data,
  test_data,
  alpha = 0.1  # 90% coverage
)
```

### Jackknife+ Conformal Prediction

**Algorithm:**

1. For each training observation i:
   - Fit model on data excluding observation i
   - Predict on held-out observation i
   - Predict on all test observations
2. Compute LOO residuals
3. For each test point j:
   - Augment residuals with LOO prediction variations
   - Compute quantile
   - Form prediction interval

**Advantages:**
- Can provide tighter intervals
- Better for small samples
- More adaptive to local structure

**Disadvantages:**
- Computationally expensive (n model fits)
- Requires careful implementation

**Implementation:**
```r
jackknife_plus_conformal(
  formula = bf_brms,
  train_data,
  test_data,
  alpha = 0.1,
  family = zero_inflated_beta(),
  prior = priors_brms
)
```

### Coverage Guarantees

Under exchangeability assumption:
```
P(Y_new ∈ [L, U]) ≥ 1 - α
```

For finite samples:
```
Coverage ≥ ⌈(n+1)(1-α)⌉/(n+1)
```

## API Reference

### Main Functions

#### `generate_zibeta_data(n_subjects, n_obs_per_subject, zi_prob)`

**Parameters:**
- `n_subjects`: Integer, number of groups (default: 50)
- `n_obs_per_subject`: Integer, observations per group (default: 10)
- `zi_prob`: Numeric [0,1], zero-inflation probability (default: 0.2)

**Returns:** Data frame with columns:
- `y`: Response variable
- `x`: Predictor variable
- `subject_id`: Grouping factor

#### `extract_and_inflate_priors(model_freq, inflation_factor)`

**Parameters:**
- `model_freq`: Fitted glmmTMB model
- `inflation_factor`: Numeric, variance inflation factor (default: 2)

**Returns:** List of prior specifications with brms-compatible names

#### `split_conformal_prediction(model_brms, train_data, calibration_data, test_data, alpha)`

**Parameters:**
- `model_brms`: Fitted brms model
- `train_data`: Training data frame
- `calibration_data`: Calibration data frame
- `test_data`: Test data frame
- `alpha`: Miscoverage level (default: 0.1)

**Returns:** List containing:
- `predictions`: Point predictions
- `lower`: Lower bounds
- `upper`: Upper bounds
- `coverage`: Empirical coverage
- `target_coverage`: Target coverage
- `mean_width`: Mean interval width
- `conformal_quantile`: Conformity score quantile

#### `run_zibeta_workflow(use_jackknife)`

**Parameters:**
- `use_jackknife`: Logical, use Jackknife+ instead of split (default: FALSE)

**Returns:** List containing:
- `model_freq`: glmmTMB model object
- `model_brms`: brms model object
- `conformal_results`: Conformal prediction results
- `data`: List of train/calib/test data

## Performance Considerations

### Computational Complexity

- **glmmTMB fitting**: O(n × p²) for n observations, p parameters
- **brms fitting**: O(iterations × n × p²)
- **Split conformal**: O(n_calib + n_test)
- **Jackknife+**: O(n_train × iterations × p²)

### Memory Requirements

- **glmmTMB**: Minimal (~100MB for moderate datasets)
- **brms**: ~500MB-2GB (Stan compilation + sampling)
- **Posterior predictive**: ~(n_draws × n_test × 8 bytes)

### Scaling Recommendations

For large datasets (n > 10,000):
- Use split conformal instead of Jackknife+
- Reduce MCMC iterations (e.g., 1000 total)
- Use `cmdstanr` backend for brms
- Consider variational inference (`algorithm = "meanfield"`)

## References

1. **Zero-inflated beta regression:**
   - Ospina, R., & Ferrari, S. L. (2012). A general class of zero-or-one inflated beta regression models. *Computational Statistics & Data Analysis*, 56(6), 1609-1623.

2. **glmmTMB:**
   - Brooks, M. E., et al. (2017). glmmTMB balances speed and flexibility among packages for zero-inflated generalized linear mixed modeling. *The R Journal*, 9(2), 378-400.

3. **brms:**
   - Bürkner, P. C. (2017). brms: An R package for Bayesian multilevel models using Stan. *Journal of Statistical Software*, 80(1), 1-28.

4. **Conformal prediction:**
   - Vovk, V., Gammerman, A., & Shafer, G. (2005). *Algorithmic learning in a random world*. Springer.
   - Barber, R. F., et al. (2021). Predictive inference with the jackknife+. *The Annals of Statistics*, 49(1), 486-507.

5. **Empirical Bayes:**
   - Efron, B. (2012). *Large-scale inference: empirical Bayes methods for estimation, testing, and prediction*. Cambridge University Press.

## License

See main repository LICENSE file.
