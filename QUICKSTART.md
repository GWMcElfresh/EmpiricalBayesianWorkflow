# Quick Start Guide

## Prerequisites

1. Install R (>= 4.0.0)
2. Install required packages:

```r
install.packages(c("glmmTMB", "brms", "dplyr", "ggplot2", "tidyr"))
```

Note: Installing `brms` will also install Stan, which may take some time.

## Running the Workflow

### Option 1: Run the Complete Example

```bash
cd /path/to/EmpiricalBayesianWorkflow
Rscript example.R
```

### Option 2: Use in Your R Session

```r
# Source the workflow functions
source("zero_inflated_beta_workflow.R")

# Run with default settings (split conformal)
results <- run_zibeta_workflow(use_jackknife = FALSE)

# Access results
print(results$conformal_results$coverage)

# View the plot
# (will be displayed automatically)
```

### Option 3: Test the Implementation

```bash
Rscript test_workflow.R
```

## Customization

### Using Your Own Data

```r
source("zero_inflated_beta_workflow.R")

# Your data should have:
# - y: response variable (0 to 1, with possible zeros)
# - x: predictor(s)
# - subject_id: grouping variable for random effects

# Fit frequentist model
model_freq <- glmmTMB(
  y ~ x + (1 + x | subject_id),
  data = your_data,
  family = beta_family(),
  ziformula = ~ 1
)

# Extract priors
priors <- extract_and_inflate_priors(model_freq)

# Continue with Bayesian fitting and conformal prediction...
```

### Adjusting Coverage Level

```r
# For 95% coverage instead of 90%
cp_results <- split_conformal_prediction(
  model_brms = your_model,
  train_data = train,
  calibration_data = calib,
  test_data = test,
  alpha = 0.05  # 1 - 0.95
)
```

### Using Jackknife+ Instead of Split Conformal

```r
# More powerful but slower
results <- run_zibeta_workflow(use_jackknife = TRUE)
```

## Understanding the Output

The workflow returns a list with:

- `model_freq`: Fitted glmmTMB model (frequentist)
- `model_brms`: Fitted brms model (Bayesian)
- `conformal_results`: List containing:
  - `predictions`: Point predictions for test data
  - `lower`: Lower bounds of prediction intervals
  - `upper`: Upper bounds of prediction intervals
  - `coverage`: Empirical coverage on test set
  - `target_coverage`: Desired coverage level
  - `mean_width`: Average width of intervals

## Troubleshooting

### Stan/brms Installation Issues

If you encounter issues with Stan:

1. Make sure you have a C++ compiler installed
2. See: https://github.com/stan-dev/rstan/wiki/RStan-Getting-Started

### Memory Issues

For large datasets:
- Reduce the number of MCMC iterations
- Use fewer chains
- Consider using `cmdstanr` backend for brms

### Slow Jackknife+ 

Jackknife+ fits n models (where n = training set size):
- Use split conformal for faster results
- Reduce training set size
- Use parallel processing if available

## Next Steps

1. Read the main README.md for detailed documentation
2. Explore the code in `zero_inflated_beta_workflow.R`
3. Modify the workflow for your specific use case
4. Consider extending to other distributions or model structures

## Support

For issues or questions, please open an issue on GitHub.
