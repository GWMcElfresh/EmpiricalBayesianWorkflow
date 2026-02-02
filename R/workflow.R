#' Run Complete Mixed Model Workflow with Flexible Family
#'
#' Executes the complete empirical Bayesian workflow including:
#' 1. Data generation or splitting
#' 2. Frequentist model fitting
#' 3. Prior extraction and inflation
#' 4. Bayesian model refitting with cmdstanr
#' 5. Conformal prediction
#'
#' Supports multiple model families with exact parameter mapping.
#'
#' @param data Optional data frame. If NULL, synthetic data is generated.
#' @param formula Model formula. Required if data is provided.
#' @param family Family specification. Default is "zero_inflated_beta".
#' @param zi_formula Zero-inflation formula. Default is \code{~1}.
#' @param n_subjects Number of subjects for synthetic data. Default is 30.
#' @param n_obs_per_subject Observations per subject for synthetic data. Default is 20.
#' @param zi_prob Zero-inflation probability for synthetic data. Default is 0.15.
#' @param train_prop Proportion of data for training. Default is 0.5.
#' @param calib_prop Proportion of data for calibration. Default is 0.25.
#' @param inflation_factor Prior inflation factor. Default is 2.
#' @param conformal_method Conformal method: "split" or "jackknife". Default is "split".
#' @param alpha Miscoverage level for conformal prediction. Default is 0.1.
#' @param backend Stan backend: "cmdstanr" (default) or "rstan".
#' @param verbose Whether to print progress messages. Default is TRUE.
#' @param ... Additional arguments passed to fitting functions.
#'
#' @return A list containing:
#'   \item{model_freq}{Fitted frequentist model}
#'   \item{model_bayes}{Fitted Bayesian model}
#'   \item{conformal_results}{Conformal prediction results}
#'   \item{data}{List of train/calibration/test data}
#'   \item{priors}{Prior specifications}
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Run with synthetic data (zero-inflated beta)
#' results <- run_workflow(verbose = TRUE)
#'
#' # Run with different family
#' results <- run_workflow(family = "gamma", verbose = TRUE)
#'
#' # Run with your own data
#' results <- run_workflow(
#'   data = my_data,
#'   formula = y ~ x1 + x2 + (1 + x1 | subject_id),
#'   family = "beta"
#' )
#' }
run_workflow <- function(data = NULL,
                         formula = NULL,
                         family = "zero_inflated_beta",
                         zi_formula = ~1,
                         n_subjects = 30,
                         n_obs_per_subject = 20,
                         zi_prob = 0.15,
                         train_prop = 0.5,
                         calib_prop = 0.25,
                         inflation_factor = 2,
                         conformal_method = c("split", "jackknife"),
                         alpha = 0.1,
                         backend = c("cmdstanr", "rstan"),
                         chains = 4,
                         iter = 2000,
                         warmup = 1000,
                         cores = 4,
                         seed = 123,
                         verbose = TRUE,
                         ...) {
  conformal_method <- match.arg(conformal_method)
  backend <- match.arg(backend)

  if (verbose) {
    message(rep("=", 60))
    message("Mixed Model Workflow with Conformal Prediction")
    message(sprintf("Family: %s", family))
    message(rep("=", 60))
  }

  # Step 1: Data preparation
  if (is.null(data)) {
    if (verbose) message("\n[1] Generating synthetic data...")

    full_data <- generate_zibeta_data(
      n_subjects = n_subjects,
      n_obs_per_subject = n_obs_per_subject,
      zi_prob = zi_prob,
      seed = seed
    )

    # Default formula for synthetic data
    if (is.null(formula)) {
      formula <- y ~ x + (1 + x | subject_id)
    }

    if (verbose) {
      message(sprintf(
        "  Generated %d observations from %d subjects",
        nrow(full_data), length(unique(full_data$subject_id))
      ))
      message(sprintf("  Zero proportion: %.2f%%", mean(full_data$y == 0) * 100))
    }
  } else {
    if (is.null(formula)) {
      stop("Formula must be provided when using custom data")
    }
    full_data <- data
    if (verbose) message("\n[1] Using provided data...")
  }

  # Split data
  n <- nrow(full_data)
  if (!is.null(seed)) set.seed(seed)
  idx <- sample(1:n)

  n_train <- floor(train_prop * n)
  n_calib <- floor(calib_prop * n)

  train_idx <- idx[1:n_train]
  calib_idx <- idx[(n_train + 1):(n_train + n_calib)]
  test_idx <- idx[(n_train + n_calib + 1):n]

  train_data <- full_data[train_idx, ]
  calib_data <- full_data[calib_idx, ]
  test_data <- full_data[test_idx, ]

  if (verbose) {
    message(sprintf(
      "  Train: %d, Calibration: %d, Test: %d",
      nrow(train_data), nrow(calib_data), nrow(test_data)
    ))
  }

  # Step 2: Fit frequentist model
  if (verbose) message("\n[2] Fitting frequentist model...")

  model_freq <- fit_frequentist_model(
    formula = formula,
    data = train_data,
    family = family,
    zi_formula = zi_formula,
    ...
  )

  if (verbose) message("  Model fitted successfully")

  # Step 3: Extract and inflate priors
  if (verbose) message("\n[3] Extracting and inflating parameters for priors...")

  prior_list <- extract_and_inflate_priors(model_freq, inflation_factor = inflation_factor)
  priors_brms <- create_brms_priors(prior_list)

  if (verbose) {
    message(sprintf("  Extracted %d prior specifications", length(prior_list)))
  }

  # Step 4: Fit Bayesian model
  if (verbose) message(sprintf("\n[4] Fitting Bayesian model with %s backend...", backend))

  model_bayes <- fit_bayesian_model(
    formula = formula,
    data = train_data,
    family = family,
    zi_formula = zi_formula,
    prior = priors_brms,
    backend = backend,
    chains = chains,
    iter = iter,
    warmup = warmup,
    cores = cores,
    seed = seed,
    ...
  )

  if (verbose) message("  Bayesian model fitted successfully")

  # Step 5: Conformal prediction
  if (conformal_method == "split") {
    if (verbose) message("\n[5] Performing split conformal prediction...")

    cp_results <- conformal_prediction_split(
      model = model_bayes,
      calibration_data = calib_data,
      test_data = test_data,
      alpha = alpha
    )
  } else {
    if (verbose) {
      message("\n[5] Performing Jackknife+ conformal prediction...")
      message("  WARNING: This may take several minutes")
    }

    # Define fit and predict functions for jackknife+
    fit_fn <- function(formula, data, ...) {
      fit_bayesian_model(
        formula = formula,
        data = data,
        family = family,
        zi_formula = zi_formula,
        prior = priors_brms,
        backend = backend,
        chains = chains,
        iter = iter,
        warmup = warmup,
        cores = cores,
        seed = seed,
        ...
      )
    }

    predict_fn <- function(model, newdata) {
      pp <- posterior_predict_model(model, newdata, ndraws = 50)
      apply(pp, 2, stats::median)
    }

    cp_results <- conformal_prediction_jackknife(
      formula = formula,
      train_data = train_data,
      test_data = test_data,
      alpha = alpha,
      fit_function = fit_fn,
      predict_function = predict_fn,
      family = family
    )
  }

  if (verbose) {
    message("\n  Conformal Prediction Results:")
    message(sprintf("    Target coverage: %.1f%%", cp_results$target_coverage * 100))
    message(sprintf("    Empirical coverage: %.1f%%", cp_results$coverage * 100))
    message(sprintf("    Mean interval width: %.4f", cp_results$mean_width))
  }

  if (verbose) {
    message("\n", rep("=", 60))
    message("Workflow completed successfully!")
    message(rep("=", 60))
  }

  # Return results
  results <- list(
    model_freq = model_freq,
    model_bayes = model_bayes,
    conformal_results = cp_results,
    data = list(
      train = train_data,
      calibration = calib_data,
      test = test_data
    ),
    priors = prior_list
  )

  class(results) <- c("workflow_results", "list")

  return(results)
}

#' Run Zero-Inflated Beta Workflow (Deprecated)
#'
#' @param ... Arguments passed to run_workflow
#' @export
run_zibeta_workflow <- function(...) {
  .Deprecated("run_workflow")
  run_workflow(..., family = "zero_inflated_beta")
}

#' Print Method for Workflow Results
#'
#' @param x Workflow results object
#' @param ... Additional arguments (ignored)
#' @export
print.workflow_results <- function(x, ...) {
  cat("Mixed Model Workflow Results\n")
  cat("============================\n\n")
  cat("Frequentist model: ", class(x$model_freq)[1], "\n")
  cat("Bayesian model: ", class(x$model_bayes)[1], "\n")
  cat("\nConformal Prediction:\n")
  print(x$conformal_results)
  invisible(x)
}

#' @export
print.zibeta_workflow <- function(x, ...) {
  print.workflow_results(x, ...)
}
