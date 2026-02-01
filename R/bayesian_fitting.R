#' Fit Bayesian Mixed Model with Flexible Family
#'
#' Fits a mixed model using Bayesian methods (Stan/brms) with cmdstanr backend.
#' Supports multiple families including zero-inflation.
#'
#' @param formula A brms formula object or standard formula.
#' @param data Data frame containing the variables.
#' @param family Family specification. Can be:
#'   - Character: "beta", "zero_inflated_beta", "gamma", etc.
#'   - brms family object
#'   Default is "zero_inflated_beta".
#' @param zi_formula Formula for zero-inflation component (for ZI families).
#' @param prior Optional brms prior specifications.
#' @param backend Stan backend: "cmdstanr" (default) or "rstan".
#' @param chains Number of MCMC chains. Default is 4.
#' @param iter Total iterations per chain. Default is 2000.
#' @param warmup Warmup iterations per chain. Default is 1000.
#' @param cores Number of cores for parallel processing. Default is 4.
#' @param seed Random seed for reproducibility.
#' @param ... Additional arguments passed to brms::brm.
#'
#' @return A fitted Bayesian model object.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # With informative priors from frequentist model
#' model_freq <- fit_frequentist_model(
#'   y ~ x + (1 + x | subject_id),
#'   data = data,
#'   family = "zero_inflated_beta"
#' )
#' prior_list <- extract_and_inflate_priors(model_freq)
#' priors <- create_brms_priors(prior_list)
#'
#' # Fit Bayesian model with cmdstanr
#' model_bayes <- fit_bayesian_model(
#'   y ~ x + (1 + x | subject_id),
#'   data = data,
#'   family = "zero_inflated_beta",
#'   prior = priors,
#'   backend = "cmdstanr"
#' )
#'
#' # Different family
#' model_gamma <- fit_bayesian_model(
#'   y ~ x + (1 | subject_id),
#'   data = data,
#'   family = "gamma",
#'   prior = priors
#' )
#' }
fit_bayesian_model <- function(formula,
                               data,
                               family = "zero_inflated_beta",
                               zi_formula = ~1,
                               prior = NULL,
                               backend = c("cmdstanr", "rstan"),
                               chains = 4,
                               iter = 2000,
                               warmup = 1000,
                               cores = 4,
                               seed = 123,
                               ...) {
  backend <- match.arg(backend)
  
  if (!requireNamespace("brms", quietly = TRUE)) {
    stop("brms package is required but not installed.")
  }
  
  # Get family object
  if (is.character(family)) {
    family_config <- get_family_config(family)
    family_obj <- create_brms_family(family_config)
  } else if (inherits(family, "brmsfamily") || inherits(family, "family")) {
    family_obj <- family
    family_config <- NULL
  } else {
    stop("family must be a character string or brms/stats family object")
  }
  
  # Ensure formula has family specification
  if (!inherits(formula, "brmsformula")) {
    # Check if family has ZI component
    has_zi <- FALSE
    if (!is.null(family_config)) {
      has_zi <- family_config$has_zi
    } else if (is.character(family) && grepl("zero_inflated", family)) {
      has_zi <- TRUE
    }
    
    # Convert to brmsformula
    if (has_zi) {
      formula <- brms::bf(formula, zi = zi_formula, family = family_obj)
    } else {
      formula <- brms::bf(formula, family = family_obj)
    }
  }
  
  # Set up backend
  if (backend == "cmdstanr") {
    if (!requireNamespace("cmdstanr", quietly = TRUE)) {
      warning("cmdstanr not available, falling back to rstan")
      backend <- "rstan"
    }
  }
  
  # Fit model
  model <- brms::brm(
    formula = formula,
    data = data,
    prior = prior,
    chains = chains,
    iter = iter,
    warmup = warmup,
    cores = cores,
    seed = seed,
    backend = backend,
    refresh = 0,
    silent = 2,
    ...
  )
  
  # Add class attribute
  family_name <- if (!is.null(family_config)) family_config$name else "custom"
  class(model) <- c(
    paste0("bayes_brms_", family_name),
    "bayes_brms",
    class(model)
  )
  
  # Store family config as attribute
  if (!is.null(family_config)) {
    attr(model, "family_config") <- family_config
  }
  
  return(model)
}

#' Fit Bayesian Zero-Inflated Beta Mixed Model (Deprecated)
#'
#' This function is deprecated. Use \code{fit_bayesian_model} instead.
#'
#' @param ... Arguments passed to \code{fit_bayesian_model}
#' @export
fit_bayesian_zibeta <- function(...) {
  .Deprecated("fit_bayesian_model")
  fit_bayesian_model(..., family = "zero_inflated_beta")
}

#' Posterior Predictive Samples from Bayesian Model
#'
#' Generate posterior predictive samples for new data.
#'
#' @param model Fitted Bayesian model object.
#' @param newdata Data frame for predictions.
#' @param ndraws Number of posterior draws. Default is 100.
#' @param ... Additional arguments.
#'
#' @return Matrix of posterior predictive samples (ndraws × nrow(newdata)).
#'
#' @export
posterior_predict_model <- function(model, newdata, ndraws = 100, ...) {
  UseMethod("posterior_predict_model")
}

#' @export
posterior_predict_model.bayes_brms <- function(model, newdata, ndraws = 100, ...) {
  if (!requireNamespace("brms", quietly = TRUE)) {
    stop("brms package is required.")
  }
  
  brms::posterior_predict(model, newdata = newdata, ndraws = ndraws, ...)
}

#' @export
posterior_predict_model.zibeta_bayes_brms <- function(model, newdata, ndraws = 100, ...) {
  # Backward compatibility
  posterior_predict_model.bayes_brms(model, newdata, ndraws, ...)
}

#' @export
posterior_predict_model.default <- function(model, newdata, ndraws = 100, ...) {
  # Try standard posterior_predict if available
  if (requireNamespace("brms", quietly = TRUE) && inherits(model, "brmsfit")) {
    return(brms::posterior_predict(model, newdata = newdata, ndraws = ndraws, ...))
  }
  
  stop("No method available for posterior prediction with this model class.")
}

#' Posterior Predictive Samples (Deprecated)
#'
#' @param ... Arguments passed to posterior_predict_model
#' @export
posterior_predict_zibeta <- function(...) {
  .Deprecated("posterior_predict_model")
  posterior_predict_model(...)
}
