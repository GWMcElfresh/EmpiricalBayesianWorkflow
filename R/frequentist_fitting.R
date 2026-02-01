#' Fit Frequentist Mixed Model with Flexible Family
#'
#' Fits a mixed model using a frequentist approach with support for multiple
#' families including zero-inflation and hurdle models.
#' Supports multiple backends, with glmmTMB as the preferred method.
#'
#' @param formula Formula for the conditional model, including random effects.
#' @param data Data frame containing the variables.
#' @param family Family specification. Can be:
#'   - Character: "beta", "zero_inflated_beta", "gamma", "binomial", "poisson", etc.
#'   - Family object from glmmTMB or stats
#'   Default is "zero_inflated_beta".
#' @param zi_formula Formula for the zero-inflation component. Default is \code{~1}.
#'   Only used for zero-inflated families.
#' @param method Backend to use: "glmmTMB" (preferred), "auto". Default is "auto".
#' @param ... Additional arguments passed to the backend fitting function.
#'
#' @return A fitted model object with class attribute indicating the backend and family used.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Zero-inflated beta (default)
#' data <- generate_zibeta_data(n_subjects = 20, n_obs_per_subject = 10)
#' model <- fit_frequentist_model(
#'   y ~ x + (1 + x | subject_id),
#'   data = data,
#'   family = "zero_inflated_beta"
#' )
#'
#' # Regular beta
#' model <- fit_frequentist_model(
#'   y ~ x + (1 | subject_id),
#'   data = data,
#'   family = "beta"
#' )
#'
#' # Poisson
#' model <- fit_frequentist_model(
#'   count ~ x + (1 | subject_id),
#'   data = count_data,
#'   family = "poisson"
#' )
#' }
fit_frequentist_model <- function(formula, 
                                  data, 
                                  family = "zero_inflated_beta",
                                  zi_formula = ~1,
                                  method = c("auto", "glmmTMB"),
                                  ...) {
  method <- match.arg(method)
  
  # Auto-detect available method
  if (method == "auto") {
    if (requireNamespace("glmmTMB", quietly = TRUE)) {
      method <- "glmmTMB"
    } else {
      stop("No supported backend found. Please install 'glmmTMB'.")
    }
  }
  
  # Get family configuration if character
  if (is.character(family)) {
    family_config <- get_family_config(family)
    family_obj <- create_glmmTMB_family(family_config)
  } else {
    # Assume it's already a family object
    family_obj <- family
    family_config <- NULL
  }
  
  # Fit using selected method
  if (method == "glmmTMB") {
    if (!requireNamespace("glmmTMB", quietly = TRUE)) {
      stop("glmmTMB package is required but not installed.")
    }
    
    # Determine if we need zero-inflation
    use_zi <- FALSE
    if (!is.null(family_config)) {
      use_zi <- family_config$has_zi
    } else if (is.character(family) && grepl("zero_inflated", family)) {
      use_zi <- TRUE
    }
    
    # Build model call
    if (use_zi) {
      model <- glmmTMB::glmmTMB(
        formula = formula,
        data = data,
        family = family_obj,
        ziformula = zi_formula,
        ...
      )
    } else {
      model <- glmmTMB::glmmTMB(
        formula = formula,
        data = data,
        family = family_obj,
        ...
      )
    }
    
    # Add class attributes for method dispatch
    family_name <- if (!is.null(family_config)) family_config$name else "custom"
    class(model) <- c(
      paste0("freq_glmmTMB_", family_name),
      "freq_glmmTMB",
      class(model)
    )
    
    # Store family config as attribute
    if (!is.null(family_config)) {
      attr(model, "family_config") <- family_config
    }
  }
  
  return(model)
}

#' Fit Frequentist Zero-Inflated Beta Mixed Model (Deprecated)
#'
#' This function is deprecated. Use \code{fit_frequentist_model} instead.
#'
#' @param ... Arguments passed to \code{fit_frequentist_model}
#' @export
fit_frequentist_zibeta <- function(...) {
  .Deprecated("fit_frequentist_model")
  fit_frequentist_model(..., family = "zero_inflated_beta")
}

#' Extract Parameter Estimates from Frequentist Model
#'
#' Generic function to extract parameter estimates from different model backends.
#'
#' @param model Fitted frequentist model object.
#' @param ... Additional arguments.
#'
#' @return List containing fixed effects and variance components.
#' @export
extract_parameters <- function(model, ...) {
  UseMethod("extract_parameters")
}

#' @export
extract_parameters.freq_glmmTMB <- function(model, ...) {
  # Get family configuration
  family_config <- attr(model, "family_config")
  
  # Extract fixed effects
  fixed_effects <- glmmTMB::fixef(model)
  
  # Extract variance components
  vc <- glmmTMB::VarCorr(model)
  
  # Extract random effects standard deviations
  re_sd <- NULL
  if (!is.null(vc$cond) && length(vc$cond) > 0) {
    re_sd <- attr(vc$cond[[1]], "stddev")
  }
  
  # Extract family-specific parameters
  family_params <- if (!is.null(family_config)) {
    extract_family_parameters(model, family_config)
  } else {
    list()
  }
  
  # Build result
  result <- list(
    beta_fixef = fixed_effects$cond,
    re_sd = re_sd,
    method = "glmmTMB"
  )
  
  # Add zero-inflation effects if present
  if (!is.null(fixed_effects$zi)) {
    result$zi_fixef <- fixed_effects$zi
  }
  
  # Add dispersion effects if present
  if (!is.null(fixed_effects$disp)) {
    result$disp_fixef <- fixed_effects$disp
  }
  
  # Add family-specific parameters
  result <- c(result, family_params)
  
  # Add family config
  if (!is.null(family_config)) {
    result$family_config <- family_config
  }
  
  return(result)
}

#' @export
extract_parameters.zibeta_freq_glmmTMB <- function(model, ...) {
  # Backward compatibility
  # Add family config
  if (is.null(attr(model, "family_config"))) {
    attr(model, "family_config") <- get_family_config("zero_inflated_beta")
  }
  
  # Add the freq_glmmTMB class
  class(model) <- c("freq_glmmTMB", class(model))
  
  # Call the generic method
  extract_parameters.freq_glmmTMB(model, ...)
}

#' @export
extract_parameters.default <- function(model, ...) {
  stop("No method available for this model class. Supported: glmmTMB")
}
