# GasolineYield Example Workflow
# Full empirical Bayesian workflow with zero-inflated beta regression

library(EmpiricalBayesianWorkflow)
library(betareg)
library(ggplot2)

# ============================================================
# 1. Load and Prepare Data
# ============================================================

# Load GasolineYield dataset from betareg package
data("GasolineYield", package = "betareg")

# Introduce structural zeros for zero-inflated modeling
set.seed(12345)
gasoline <- GasolineYield
gasoline$yield_zi <- gasoline$yield

# Set bottom 10% to zero to simulate structural zeros
threshold <- quantile(gasoline$yield, 0.10)
gasoline$yield_zi[gasoline$yield < threshold] <- 0

# Inspect the data
cat("Dataset summary:\n")
summary(gasoline)
cat("\nNumber of zeros:", sum(gasoline$yield_zi == 0), "/", nrow(gasoline), "\n")

# ============================================================
# 2. Fit Frequentist Model
# ============================================================

library(glmmTMB)

# Define model components:
# - Mean (mu): conditional on temperature with batch random effects
# - Precision (phi): conditional on temperature (distributional regression)
# - Zero-inflation: conditional on pressure

model_freq <- glmmTMB(
  formula = yield_zi ~ temp + (1 | batch),
  data = gasoline,
  family = beta_family(),
  ziformula = ~ pressure,
  dispformula = ~ temp,  # Precision depends on temperature
  REML = FALSE
)

cat("\n============================================================\n")
cat("Frequentist Model Summary\n")
cat("============================================================\n")
summary(model_freq)

# ============================================================
# 3. Extract and Inflate Priors
# ============================================================

# Setup family configuration
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

# Extract priors with inflation factor of 2.5
priors <- extract_and_inflate_priors(model_freq, inflation_factor = 2.5)

cat("\n============================================================\n")
cat("Extracted Priors\n")
cat("============================================================\n")
cat("Prior parameters:", names(priors), "\n\n")

# Show example priors
cat("Temperature effect on mean (b_temp):\n")
print(priors$b_temp)
cat("\nPressure effect on zero-inflation (zi_pressure):\n")
print(priors$zi_pressure)

# ============================================================
# 4. Fit Bayesian Model
# ============================================================

library(brms)

# Convert priors to brms format
brms_priors <- create_brms_priors(priors)

# Define distributional formula
brms_formula <- bf(
  yield_zi ~ temp + (1 | batch),    # Mean model
  phi ~ temp,                        # Precision model
  zi ~ pressure                      # Zero-inflation model
)

# Fit Bayesian model
cat("\n============================================================\n")
cat("Fitting Bayesian Model (this may take a few minutes...)\n")
cat("============================================================\n")

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

cat("\n============================================================\n")
cat("Bayesian Model Summary\n")
cat("============================================================\n")
summary(model_bayes)

# Plot posterior distributions
plot(model_bayes)

# ============================================================
# 5. Compare Frequentist vs Bayesian Estimates
# ============================================================

cat("\n============================================================\n")
cat("Parameter Comparison: Frequentist vs Bayesian\n")
cat("============================================================\n")

# Fixed effects
freq_fixef <- fixef(model_freq)$cond
bayes_fixef <- fixef(model_bayes)

comparison_df <- data.frame(
  Parameter = c("Intercept", "temp"),
  Frequentist = freq_fixef[c("(Intercept)", "temp")],
  Bayesian = bayes_fixef[c("Intercept", "temp"), "Estimate"],
  Bayes_SE = bayes_fixef[c("Intercept", "temp"), "Est.Error"]
)

print(comparison_df)

# Zero-inflation comparison
freq_zi <- fixef(model_freq)$zi
bayes_zi <- fixef(model_bayes, component = "zi")

cat("\nZero-Inflation Parameters:\n")
comparison_zi <- data.frame(
  Parameter = c("ZI Intercept", "ZI pressure"),
  Frequentist = freq_zi[c("(Intercept)", "pressure")],
  Bayesian = bayes_zi[c("Intercept", "pressure"), "Estimate"],
  Bayes_SE = bayes_zi[c("Intercept", "pressure"), "Est.Error"]
)
print(comparison_zi)

# ============================================================
# 6. Split Conformal Prediction
# ============================================================

cat("\n============================================================\n")
cat("Performing Split Conformal Prediction\n")
cat("============================================================\n")

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

# Refit Bayesian model on training data only
cat("Refitting model on training data...\n")

model_freq_train <- glmmTMB(
  formula = yield_zi ~ temp + (1 | batch),
  data = train_data,
  family = beta_family(),
  ziformula = ~ pressure,
  dispformula = ~ temp,
  REML = FALSE
)

attr(model_freq_train, "family_config") <- family_config
class(model_freq_train) <- c("freq_glmmTMB_zero_inflated_beta", 
                             "freq_glmmTMB", 
                             class(model_freq_train))

priors_train <- extract_and_inflate_priors(model_freq_train, inflation_factor = 2.5)
brms_priors_train <- create_brms_priors(priors_train)

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
  refresh = 0,
  seed = 12345
)

# Perform conformal prediction
# Note: conformal_prediction_split expects 'y' as response variable
calib_data$y <- calib_data$yield_zi
test_data$y <- test_data$yield_zi

cp_results <- conformal_prediction_split(
  model = model_bayes_train,
  calibration_data = calib_data,
  test_data = test_data,
  alpha = 0.1  # 90% target coverage
)

cat("\n============================================================\n")
cat("Conformal Prediction Results\n")
cat("============================================================\n")
print(cp_results)

# ============================================================
# 7. Visualize Results
# ============================================================

cat("\n============================================================\n")
cat("Creating Visualizations\n")
cat("============================================================\n")

# Prepare results dataframe
results_df <- data.frame(
  observed = test_data$yield_zi,
  predicted = cp_results$predictions,
  lower = cp_results$lower,
  upper = cp_results$upper,
  covered = (test_data$yield_zi >= cp_results$lower) & 
            (test_data$yield_zi <= cp_results$upper),
  temp = test_data$temp,
  pressure = test_data$pressure
)

# Plot 1: Prediction intervals
p1 <- ggplot(results_df, aes(x = 1:nrow(results_df))) +
  geom_point(aes(y = observed, color = covered), size = 3) +
  geom_point(aes(y = predicted), shape = 4, size = 4, color = "blue") +
  geom_errorbar(aes(ymin = lower, ymax = upper), alpha = 0.3, width = 0.2) +
  scale_color_manual(
    values = c("red", "darkgreen"),
    labels = c("Miscovered", "Covered")
  ) +
  labs(
    x = "Test Observation",
    y = "Gasoline Yield",
    title = "90% Conformal Prediction Intervals",
    subtitle = sprintf("Coverage: %.1f%% (target: 90%%)", 100 * cp_results$coverage),
    color = "Status"
  ) +
  theme_minimal() +
  theme(legend.position = "bottom")

print(p1)
ggsave("gasoline_conformal_intervals.png", p1, width = 10, height = 6)

# Plot 2: Observed vs Predicted
p2 <- ggplot(results_df, aes(x = predicted, y = observed)) +
  geom_point(aes(color = covered), size = 3) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray") +
  scale_color_manual(
    values = c("red", "darkgreen"),
    labels = c("Miscovered", "Covered")
  ) +
  labs(
    x = "Predicted Yield",
    y = "Observed Yield",
    title = "Observed vs Predicted Values",
    color = "Status"
  ) +
  theme_minimal() +
  theme(legend.position = "bottom")

print(p2)
ggsave("gasoline_obs_vs_pred.png", p2, width = 8, height = 6)

# Plot 3: Coverage by temperature
p3 <- ggplot(results_df, aes(x = temp, y = observed)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2, fill = "blue") +
  geom_point(aes(color = covered), size = 3) +
  geom_line(aes(y = predicted), color = "blue", size = 1) +
  scale_color_manual(
    values = c("red", "darkgreen"),
    labels = c("Miscovered", "Covered")
  ) +
  labs(
    x = "Temperature",
    y = "Gasoline Yield",
    title = "Predictions vs Temperature",
    subtitle = "With 90% conformal intervals",
    color = "Status"
  ) +
  theme_minimal() +
  theme(legend.position = "bottom")

print(p3)
ggsave("gasoline_vs_temperature.png", p3, width = 10, height = 6)

# ============================================================
# 8. Summary Statistics
# ============================================================

cat("\n============================================================\n")
cat("Summary Statistics\n")
cat("============================================================\n")

cat(sprintf("Coverage: %.2f (target: 0.90)\n", cp_results$coverage))
cat(sprintf("Mean interval width: %.4f\n", cp_results$mean_width))
cat(sprintf("Conformal quantile: %.4f\n", cp_results$conformal_quantile))

# Analyze miscoverage
n_miscovered <- sum(!results_df$covered)
cat(sprintf("\nMiscovered observations: %d / %d (%.1f%%)\n", 
            n_miscovered, nrow(results_df), 100 * n_miscovered / nrow(results_df)))

# Summary by coverage status
cat("\nMiscovered observations characteristics:\n")
miscovered_df <- results_df[!results_df$covered, ]
if (nrow(miscovered_df) > 0) {
  cat(sprintf("  Mean observed: %.4f\n", mean(miscovered_df$observed)))
  cat(sprintf("  Mean predicted: %.4f\n", mean(miscovered_df$predicted)))
  cat(sprintf("  Mean abs error: %.4f\n", mean(abs(miscovered_df$observed - miscovered_df$predicted))))
}

cat("\nWorkflow complete! All plots saved to working directory.\n")
