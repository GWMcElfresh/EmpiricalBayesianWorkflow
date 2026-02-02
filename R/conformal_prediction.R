#' Split Conformal Prediction for Zero-Inflated Beta Models
#'
#' Performs split conformal prediction to generate prediction intervals with
#' guaranteed coverage under exchangeability assumptions.
#'
#' @param model Fitted model object (Bayesian or frequentist).
#' @param calibration_data Calibration data for computing nonconformity scores.
#' @param test_data Test data for predictions.
#' @param alpha Miscoverage level (1 - coverage). Default is 0.1 for 90% coverage.
#' @param method Method for computing predictions: "posterior_median" or "mean". Default is "posterior_median".
#' @param ndraws Number of posterior draws for Bayesian models. Default is 100.
#' @param score_function Function to compute nonconformity scores. Default is absolute residual.
#'
#' @return A list containing:
#'   \item{predictions}{Point predictions for test data}
#'   \item{lower}{Lower bounds of prediction intervals}
#'   \item{upper}{Upper bounds of prediction intervals}
#'   \item{coverage}{Empirical coverage on test set}
#'   \item{target_coverage}{Target coverage level}
#'   \item{mean_width}{Average width of prediction intervals}
#'   \item{conformal_quantile}{Computed conformal quantile}
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Split data
#' n <- nrow(data)
#' idx <- sample(1:n)
#' train_idx <- idx[1:floor(0.5*n)]
#' calib_idx <- idx[(floor(0.5*n)+1):floor(0.75*n)]
#' test_idx <- idx[(floor(0.75*n)+1):n]
#'
#' # Fit model and predict
#' model <- fit_bayesian_zibeta(formula, data = data[train_idx, ])
#' cp_results <- conformal_prediction_split(
#'   model,
#'   calibration_data = data[calib_idx, ],
#'   test_data = data[test_idx, ],
#'   alpha = 0.1
#' )
#' }
conformal_prediction_split <- function(model,
                                       calibration_data,
                                       test_data,
                                       alpha = 0.1,
                                       method = c("posterior_median", "mean"),
                                       ndraws = 100,
                                       score_function = NULL) {
  method <- match.arg(method)
  
  # Default score function: absolute residual
  if (is.null(score_function)) {
    score_function <- function(y, y_pred) abs(y - y_pred)
  }
  
  # Generate predictions for calibration data
  if (inherits(model, "bayes_brms") || inherits(model, "brmsfit") || 
      inherits(model, "zibeta_bayes_brms")) {
    pp_calibration <- posterior_predict_model(model, newdata = calibration_data, ndraws = ndraws, allow_new_levels = TRUE)
    
    if (method == "posterior_median") {
      pred_calibration <- apply(pp_calibration, 2, stats::median)
    } else {
      pred_calibration <- apply(pp_calibration, 2, mean)
    }
  } else {
    stop("Model class not supported for conformal prediction")
  }
  
  # Calculate nonconformity scores
  scores <- score_function(calibration_data$y, pred_calibration)
  
  # Calculate conformal quantile
  n_calib <- length(scores)
  q_level <- ceiling((n_calib + 1) * (1 - alpha)) / n_calib
  q_conformal <- stats::quantile(scores, probs = q_level, na.rm = TRUE)
  
  # Generate predictions for test data
  if (inherits(model, "bayes_brms") || inherits(model, "brmsfit") || 
      inherits(model, "zibeta_bayes_brms")) {
    pp_test <- posterior_predict_model(model, newdata = test_data, ndraws = ndraws, allow_new_levels = TRUE)
    
    if (method == "posterior_median") {
      pred_test <- apply(pp_test, 2, stats::median)
    } else {
      pred_test <- apply(pp_test, 2, mean)
    }
  }
  
  # Construct prediction intervals (bounded to [0, 1])
  lower <- pmax(pred_test - q_conformal, 0)
  upper <- pmin(pred_test + q_conformal, 1)
  
  # Calculate empirical coverage on test set
  coverage <- mean(test_data$y >= lower & test_data$y <= upper, na.rm = TRUE)
  
  # Calculate interval widths
  widths <- upper - lower
  
  results <- list(
    predictions = pred_test,
    lower = lower,
    upper = upper,
    coverage = coverage,
    target_coverage = 1 - alpha,
    mean_width = mean(widths),
    conformal_quantile = q_conformal,
    method = "split_conformal"
  )
  
  class(results) <- c("conformal_prediction", "list")
  
  return(results)
}

#' Jackknife+ Conformal Prediction
#'
#' Performs Jackknife+ conformal prediction, which can provide tighter intervals
#' but requires fitting n leave-one-out models.
#'
#' @param formula Model formula.
#' @param train_data Training data.
#' @param test_data Test data for predictions.
#' @param alpha Miscoverage level. Default is 0.1.
#' @param fit_function Function to fit model. Should accept formula and data.
#' @param predict_function Function to predict. Should accept model and newdata.
#' @param ... Additional arguments passed to fit_function.
#'
#' @return A list with the same structure as \code{conformal_prediction_split}.
#'
#' @export
conformal_prediction_jackknife <- function(formula,
                                           train_data,
                                           test_data,
                                           alpha = 0.1,
                                           fit_function,
                                           predict_function,
                                           ...) {
  n_train <- nrow(train_data)
  n_test <- nrow(test_data)
  
  # Store LOO predictions
  loo_predictions <- numeric(n_train)
  test_predictions <- matrix(NA, nrow = n_test, ncol = n_train)
  
  message("Running Jackknife+ conformal prediction...")
  
  # Leave-one-out fitting
  for (i in 1:n_train) {
    if (i %% 10 == 0) {
      message(sprintf("  Fitting LOO model %d/%d", i, n_train))
    }
    
    # Leave out observation i
    train_loo <- train_data[-i, ]
    
    # Fit model
    model_loo <- fit_function(formula, data = train_loo, ...)
    
    # Predict on left-out observation
    loo_predictions[i] <- predict_function(model_loo, train_data[i, , drop = FALSE])
    
    # Predict on test set
    test_predictions[, i] <- predict_function(model_loo, test_data)
  }
  
  # Calculate residuals (nonconformity scores)
  residuals <- abs(train_data$y - loo_predictions)
  
  # For each test point, calculate prediction interval
  lower <- numeric(n_test)
  upper <- numeric(n_test)
  
  for (j in 1:n_test) {
    # Jackknife+ augmented residuals
    pred_j_median <- stats::median(test_predictions[j, ])
    augmented_residuals <- c(residuals, abs(test_predictions[j, ] - pred_j_median))
    
    # Calculate quantile
    q_level <- ceiling((n_train + 1) * (1 - alpha)) / (n_train + 1)
    q_val <- stats::quantile(augmented_residuals, probs = q_level, na.rm = TRUE)
    
    # Prediction interval
    lower[j] <- max(pred_j_median - q_val, 0)
    upper[j] <- min(pred_j_median + q_val, 1)
  }
  
  # Calculate coverage
  pred_test <- apply(test_predictions, 1, stats::median)
  coverage <- mean(test_data$y >= lower & test_data$y <= upper, na.rm = TRUE)
  widths <- upper - lower
  
  results <- list(
    predictions = pred_test,
    lower = lower,
    upper = upper,
    coverage = coverage,
    target_coverage = 1 - alpha,
    mean_width = mean(widths),
    loo_residuals = residuals,
    method = "jackknife_plus"
  )
  
  class(results) <- c("conformal_prediction", "list")
  
  return(results)
}

#' Print Method for Conformal Prediction Results
#'
#' @param x Conformal prediction results object
#' @param ... Additional arguments (ignored)
#' @export
print.conformal_prediction <- function(x, ...) {
  cat("Conformal Prediction Results\n")
  cat("=============================\n")
  cat(sprintf("Method: %s\n", x$method))
  cat(sprintf("Target coverage: %.1f%%\n", x$target_coverage * 100))
  cat(sprintf("Empirical coverage: %.1f%%\n", x$coverage * 100))
  cat(sprintf("Mean interval width: %.4f\n", x$mean_width))
  cat(sprintf("Number of predictions: %d\n", length(x$predictions)))
  invisible(x)
}
