# EmpiricalBayesianWorkflow AI Instructions

## Package Architecture

This is an R package implementing an **empirical Bayesian workflow** with conformal prediction. The core pipeline: frequentist model fitting → prior extraction → Bayesian refitting → conformal prediction intervals.

### Key Components

- **`R/workflow.R`**: Main entry point via `run_workflow()` orchestrating the complete pipeline
- **`R/family_mapping.R`**: Critical parameter mapping system between glmmTMB and brms (see below)
- **`R/frequentist_fitting.R`**: glmmTMB backend for initial model fitting
- **`R/bayesian_fitting.R`**: brms/cmdstanr backend for Bayesian refitting
- **`R/prior_construction.R`**: Extracts frequentist estimates and inflates variances for priors
- **`R/conformal_prediction.R`**: Split conformal and Jackknife+ methods for prediction intervals
- **`R/data_generation.R`**: Synthetic data generation for testing (zero-inflated beta, etc.)

## Critical Design Patterns

### 1. Family Mapping System (Most Complex Component)

The package's core innovation is **exact parameter correspondence** between glmmTMB and brms across 8+ model families. Each family has a configuration in `get_supported_families()`:

```r
zero_inflated_beta = list(
  name = "zero_inflated_beta",
  glmmTMB_family = "beta_family()",
  brms_family = "zero_inflated_beta()",
  has_zi = TRUE,
  params = c("phi"),
  param_classes = list(phi = "phi"),
  zi_link = "logit"
)
```

**When adding new families:** Must specify parameter mappings for fixed effects (β), zero-inflation (zi), dispersion (phi/shape), and random effects (sd). Use `map_parameter_name()` to translate between frameworks.

### 2. Prior Construction Pattern

`extract_and_inflate_priors()` follows this sequence:
1. Extract fixed effects from conditional model (β)
2. Extract zero-inflation parameters (if `has_zi = TRUE`)
3. Extract random effect SDs and correlations
4. Extract family-specific parameters (phi, shape, etc.)
5. **Inflate variances** by factor (default 2×) to avoid overconfident priors
6. Map names to brms conventions: `b_Intercept`, `zi_Intercept`, `sd_subject_id__Intercept`

Always preserve this naming: brms expects `b_` prefix for fixed effects, `zi_` for zero-inflation, `sd_GROUP__TERM` for random effects.

#### Prior Family Mappings (Critical for CI/CD)

The following table documents the **exact prior structure** that brms expects for each supported family. These mappings are essential for the `extract_and_inflate_priors()` and `create_brms_priors()` functions.

**General Structure:**
- All families: Fixed effects get `class = "b"`, Intercepts get `class = "Intercept"`
- Random effects: Always `class = "sd"` with `group` specified
- Family-specific parameters: `class = "phi"`, `"shape"`, or `"sigma"`
- Zero-inflation: `class = "Intercept"`, `dpar = "zi"`

##### Gaussian Family
```r
# brms prior classes:
class = "b"          # Fixed effect slopes (flat prior)
class = "Intercept"  # Fixed effect intercept (student_t(3, μ_y, 2.5))
class = "sd"         # Random effect SDs (student_t(3, 0, 2.5))
class = "sigma"      # Residual SD (student_t(3, 0, 2.5))

# glmmTMB extracts:
- Fixed Effects (cond): "(Intercept)", "x", ...
- Dispersion (disp): "(Intercept)" → maps to sigma
```

##### Beta Family
```r
# brms prior classes:
class = "b"          # Fixed effect slopes (flat prior)
class = "Intercept"  # Fixed effect intercept (student_t(3, 0, 2.5))
class = "sd"         # Random effect SDs (student_t(3, 0, 2.5))
class = "phi"        # Precision parameter (gamma(0.01, 0.01))

# glmmTMB extracts:
- Fixed Effects (cond): "(Intercept)", "x", ...
- Dispersion (disp): "(Intercept)" → maps to phi
```

##### Gamma Family
```r
# brms prior classes:
class = "b"          # Fixed effect slopes (flat prior)
class = "Intercept"  # Fixed effect intercept (student_t(3, μ_y, 2.5))
class = "sd"         # Random effect SDs (student_t(3, 0, 2.5))
class = "shape"      # Shape parameter (gamma(0.01, 0.01))

# glmmTMB extracts:
- Fixed Effects (cond): "(Intercept)", "x", ...
- Dispersion (disp): "(Intercept)" → maps to shape
```

##### Binomial Family
```r
# brms prior classes:
class = "b"          # Fixed effect slopes (flat prior)
class = "Intercept"  # Fixed effect intercept (student_t(3, 0, 2.5))
class = "sd"         # Random effect SDs (student_t(3, 0, 2.5))

# glmmTMB extracts:
- Fixed Effects (cond): "(Intercept)", "x", ...
- NO dispersion parameter for binomial
```

##### Poisson Family
```r
# brms prior classes:
class = "b"          # Fixed effect slopes (flat prior)
class = "Intercept"  # Fixed effect intercept (student_t(3, log(μ_y), 2.5))
class = "sd"         # Random effect SDs (student_t(3, 0, 2.5))

# glmmTMB extracts:
- Fixed Effects (cond): "(Intercept)", "x", ...
- NO dispersion parameter for Poisson
```

##### Negative Binomial (nbinom2) Family
```r
# brms prior classes:
class = "b"          # Fixed effect slopes (flat prior)
class = "Intercept"  # Fixed effect intercept (student_t(3, log(μ_y), 2.5))
class = "sd"         # Random effect SDs (student_t(3, 0, 2.5))
class = "shape"      # Dispersion parameter (inv_gamma(0.4, 0.3))

# glmmTMB extracts:
- Fixed Effects (cond): "(Intercept)", "x", ...
- Dispersion (disp): "(Intercept)" → maps to shape
```

##### Zero-Inflated Beta Family
```r
# brms prior classes:
class = "b"          # Conditional fixed effect slopes (flat prior)
class = "Intercept"  # Conditional intercept (student_t(3, 0, 2.5))
class = "sd"         # Random effect SDs (student_t(3, 0, 2.5))
class = "phi"        # Precision parameter (gamma(0.01, 0.01))
class = "Intercept", dpar = "zi"  # ZI intercept (logistic(0, 1))

# glmmTMB extracts:
- Fixed Effects (cond): "(Intercept)", "x", ...
- Zero-Inflation (zi): "(Intercept)" → maps to zi_Intercept
- Dispersion (disp): "(Intercept)" → maps to phi
```

##### Zero-Inflated Negative Binomial (nbinom2) Family
```r
# brms prior classes:
class = "b"          # Conditional fixed effect slopes (flat prior)
class = "Intercept"  # Conditional intercept (student_t(3, log(μ_y), 2.5))
class = "sd"         # Random effect SDs (student_t(3, 0, 2.5))
class = "shape"      # Dispersion parameter (inv_gamma(0.4, 0.3))
class = "Intercept", dpar = "zi"  # ZI intercept (logistic(0, 1))

# glmmTMB extracts:
- Fixed Effects (cond): "(Intercept)", "x", ...
- Zero-Inflation (zi): "(Intercept)" → maps to zi_Intercept
- Dispersion (disp): "(Intercept)" → maps to shape
```

**Critical Naming Conventions:**

| Component | glmmTMB Name | brms Prior Name | Notes |
|-----------|--------------|-----------------|-------|
| Fixed effect intercept | `(Intercept)` | `b_Intercept` | Use `class = "Intercept"` in prior |
| Fixed effect slope | `x` | `b_x` | Use `class = "b", coef = "x"` in prior |
| Random intercept SD | `Std.Dev.subject_id.(Intercept)` | `sd_subject_id__Intercept` | Note double underscore |
| Random slope SD | `Std.Dev.subject_id.x` | `sd_subject_id__x` | Note double underscore |
| ZI intercept | `(Intercept)` in zi component | `zi_Intercept` | Use `dpar = "zi"` in prior |
| Dispersion (beta) | `(Intercept)` in disp | `phi` | Use `class = "phi"` |
| Dispersion (gamma/nbinom) | `(Intercept)` in disp | `shape` | Use `class = "shape"` |
| Residual SD (gaussian) | `(Intercept)` in disp | `sigma` | Use `class = "sigma"` |

**Prior Construction Workflow:**

1. **Extract from glmmTMB:** Use `summary(model)$coefficients$cond` for fixed effects, `$zi` for zero-inflation, VarCorr for random effects
2. **Map parameter names:** Apply `map_parameter_name()` with correct component type
3. **Set prior values:** Use extracted estimates as location, inflate variance
4. **Create brms prior:** Use `set_prior()` with correct class, coef, group, dpar
5. **Validate structure:** Check against `brms::get_prior()` output for the family

**Common CI/CD Failures:**

- **Missing phi/shape priors:** Beta/gamma/nbinom families REQUIRE dispersion priors
- **Wrong ZI prior format:** Must specify `dpar = "zi"` for zero-inflation component
- **Incorrect random effect naming:** Must use double underscore `sd_GROUP__TERM`
- **Missing coef specification:** Slope priors need `coef = "x"` to target specific coefficients

### 3. Dual Backend Support

The package defaults to **cmdstanr** (faster) but falls back to rstan:

```r
backend <- match.arg(backend, c("cmdstanr", "rstan"))
```

When modifying Bayesian fitting functions, always check both backends are handled. cmdstanr is preferred; add it to GitHub Actions workflows.

### 4. Data Splitting Convention

Standard split in `run_workflow()`:
- **Train**: 50% (fit frequentist model)
- **Calibration**: 25% (compute conformal scores)
- **Test**: 25% (evaluate coverage)

The calibration set is crucial—never use training data for conformal calibration (invalidates coverage guarantees).

## Development Workflows

### Running Tests

```r
# Full test suite
devtools::test()

# Specific test file
testthat::test_file("tests/testthat/test-mtcars_workflow.R")

# Run end-to-end validation
Rscript validate_workflow.R
```

**Important:** Tests use `skip_on_cran()` for computationally intensive checks. End-to-end tests reduce chains/iterations (e.g., `chains = 1, iter = 1000`).

### Building Documentation

```r
# Regenerate .Rd files from roxygen2 comments
roxygen2::roxygenise()

# Check documentation consistency
devtools::check()
```

Always use roxygen2 comments (`#'`) in R/ files. Export public functions with `@export`. Private helpers should not have `@export`.

### Local Package Development

```r
# Load all functions without installing
devtools::load_all()

# Install locally
devtools::install()

# Check package (like R-CMD-check CI)
devtools::check()
```

### GitHub Actions CI

See `.github/workflows/R-CMD-check.yaml`. Key steps:
1. Install system dependencies (Linux: libcurl, libssl)
2. Install R packages via `r-lib/actions/setup-r-dependencies@v2`
3. **Install cmdstan** explicitly (not automatic with cmdstanr)
4. Run `roxygen2::roxygenise()` before checks
5. Run `rcmdcheck::rcmdcheck()`

**When adding dependencies:** Update both `DESCRIPTION` (Suggests/Imports) and the `extra-packages` section in the workflow.

## Package-Specific Conventions

### Function Naming

- `fit_*`: Model fitting functions (frequentist/Bayesian)
- `extract_*`: Parameter extraction from fitted models
- `create_*`: Object construction (families, priors)
- `posterior_predict_*`: Generate predictions from Bayesian models
- `conformal_prediction_*`: Conformal methods (split/jackknife)
- `get_*`: Accessors/queries (families, configs)
- `run_*`: High-level workflow orchestrators

### Return Objects

Main functions return S3 objects with custom classes:
- `run_workflow()` → `workflow_results` (list with `model_freq`, `model_bayes`, `conformal_results`, `data`, `priors`)
- Conformal methods → `conformal_prediction` (list with `predictions`, `lower`, `upper`, `coverage`, `mean_width`)

Define `print.workflow_results()` and `print.conformal_prediction()` methods for clean output.

### Formula Handling

Mixed model formulas follow lme4/glmmTMB syntax:
```r
y ~ x + (1 + x | subject_id)  # Random intercept + slope
```

Zero-inflation formulas are separate:
```r
zi_formula = ~1  # ZI probability depends only on intercept
```

Use `brms::bf()` for complex formulas with multiple components (conditional + ZI).

### Error Handling in Workflow

`run_workflow()` uses `tryCatch()` sparingly. Prefer early validation:
```r
if (!family %in% names(get_supported_families())) {
  stop("Unsupported family: ", family)
}
```

For Bayesian fitting failures, catch and provide informative messages about priors/data issues.

## Integration Points

### glmmTMB → brms Pipeline

1. Fit glmmTMB model: `fit_frequentist_model()`
2. Extract parameters: `extract_parameters(model_freq)` (generic S3 method)
3. Map to brms names: `map_parameter_name(param, family_config, component)`
4. Construct priors: `create_brms_priors(prior_list)`
5. Fit brms: `fit_bayesian_model(..., prior = priors, backend = "cmdstanr")`

**Critical:** The `extract_parameters()` function is S3 generic with methods for different model classes. Always return a list with `beta_fixef`, `zi_fixef`, `family_config`, `random_effects`.

### External Dependencies

- **glmmTMB**: Frequentist fitting (zero-inflation support)
- **brms**: Bayesian wrapper around Stan (handles formula translation)
- **cmdstanr**: Fast Stan interface (requires manual cmdstan installation)
- **dplyr/tidyr**: Data manipulation in examples only (not core functions)

The package is **agnostic by design**—can theoretically swap backends (e.g., replace glmmTMB with lme4 for non-ZI models).

## Common Pitfalls

1. **Parameter name mismatches**: Always use `map_parameter_name()` when moving between glmmTMB and brms. Manual name mapping will break.

2. **Zero-inflation detection**: Check `family_config$has_zi` before accessing `zi_fixef`. Not all families support ZI.

3. **Random effect naming**: brms uses double underscore: `sd_subject_id__Intercept`. glmmTMB uses `Std.Dev.subject_id.(Intercept)`.

4. **Conformal data leakage**: Never compute nonconformity scores on training data. Always use separate calibration set.

5. **cmdstanr not installed**: Test for cmdstanr availability and fall back to rstan gracefully:
   ```r
   if (!requireNamespace("cmdstanr", quietly = TRUE)) {
     warning("cmdstanr not available, using rstan")
     backend <- "rstan"
   }
   ```

6. **Family-specific parameters**: Beta has `phi`, gamma has `shape`, negative binomial has `shape`. These are extracted differently from glmmTMB—see `extract_family_parameters()`.

## Adding New Features

### Adding a New Family

1. Add entry to `get_supported_families()` with correct parameter mappings
2. Update `create_glmmTMB_family()` to handle the family
3. Update `create_brms_family()` to map to brms equivalent
4. Add test in `tests/testthat/test-family_mapping.R`
5. Update README.md with the new family

### Adding a New Conformal Method

1. Create function `conformal_prediction_<method>()` in `R/conformal_prediction.R`
2. Follow return structure: `list(predictions, lower, upper, coverage, target_coverage, mean_width)`
3. Add class `"conformal_prediction"` to return object
4. Update `run_workflow()` to accept new method in `conformal_method` argument
5. Add unit test

## File Organization

- **R/**: All package functions (use roxygen2 comments)
- **man/**: Auto-generated .Rd files (don't edit manually)
- **tests/testthat/**: Unit tests (prefix with `test-`)
- **tests/**: Integration tests (e.g., `test_workflow.R`)
- **Root scripts**: Examples (`example.R`, `zero_inflated_beta_workflow.R`) for demonstration only

## References

- [TECHNICAL.md](TECHNICAL.md): Mathematical details of zero-inflated beta models
- [IMPLEMENTATION_SUMMARY.md](IMPLEMENTATION_SUMMARY.md): Complete feature checklist
- [QUICKSTART.md](QUICKSTART.md): User-facing quick start guide
