#' Extract and Inflate Prior Specifications for Bayesian Fitting
#'
#' Converts frequentist parameter estimates into prior specifications for
#' Bayesian models. Inflates variances to account for estimation uncertainty
#' and avoid overconfident priors. Works with multiple model families.
#'
#' @param model Fitted frequentist model object.
#' @param inflation_factor Factor by which to inflate prior variances. Default is 2.
#'
#' @return A list of prior specifications with brms-compatible naming.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' model_freq <- fit_frequentist_model(
#'   y ~ x + (1 + x | subject_id),
#'   data = data,
#'   family = "zero_inflated_beta"
#' )
#' priors <- extract_and_inflate_priors(model_freq, inflation_factor = 2)
#' }
extract_and_inflate_priors <- function(model, inflation_factor = 2) {
  # Extract parameters using generic method
  params <- extract_parameters(model)
  
  # Get family config if available
  family_config <- params$family_config
  
  priors <- list()
  
  # Fixed effects priors for conditional model
  if (!is.null(params$beta_fixef)) {
    for (i in seq_along(params$beta_fixef)) {
      param_name <- names(params$beta_fixef)[i]
      mean_val <- params$beta_fixef[i]
      # Inflate the variance
      sd_val <- abs(mean_val) * inflation_factor + 0.5
      
      # Map parameter name
      if (!is.null(family_config)) {
        brms_name <- map_parameter_name(param_name, family_config, "cond")
      } else {
        brms_name <- if (param_name == "(Intercept)") "Intercept" else param_name
      }
      
      # Store with b_ prefix for brms
      if (brms_name == "Intercept") {
        priors[["b_Intercept"]] <- c(mean = mean_val, sd = sd_val)
      } else {
        priors[[paste0("b_", brms_name)]] <- c(mean = mean_val, sd = sd_val)
      }
    }
  }
  
  # Fixed effects priors for zero-inflation part
  if (!is.null(params$zi_fixef)) {
    for (i in seq_along(params$zi_fixef)) {
      param_name <- names(params$zi_fixef)[i]
      mean_val <- params$zi_fixef[i]
      sd_val <- abs(mean_val) * inflation_factor + 0.5
      
      # Map parameter name
      if (!is.null(family_config)) {
        brms_name <- map_parameter_name(param_name, family_config, "zi")
      } else {
        brms_name <- if (param_name == "(Intercept)") "zi_Intercept" else paste0("zi_", param_name)
      }
      
      priors[[brms_name]] <- c(mean = mean_val, sd = sd_val)
    }
  }
  
  # Fixed effects priors for dispersion part (if present)
  if (!is.null(params$disp_fixef)) {
    for (i in seq_along(params$disp_fixef)) {
      param_name <- names(params$disp_fixef)[i]
      mean_val <- params$disp_fixef[i]
      sd_val <- abs(mean_val) * inflation_factor + 0.5
      
      brms_name <- if (param_name == "(Intercept)") "disp_Intercept" else paste0("disp_", param_name)
      priors[[brms_name]] <- c(mean = mean_val, sd = sd_val)
    }
  }
  
  # Random effects priors
  if (!is.null(params$re_sd)) {
    for (i in seq_along(params$re_sd)) {
      param_name <- names(params$re_sd)[i]
      mean_val <- params$re_sd[i]
      sd_val <- mean_val * inflation_factor
      
      priors[[paste0("sd_", param_name)]] <- c(mean = 0, sd = sd_val)
    }
  }
  
  # Family-specific parameters
  # Phi parameter (beta, zero-inflated beta)
  if (!is.null(params$phi)) {
    priors[["phi"]] <- c(mean = params$phi, sd = params$phi * inflation_factor)
  }
  
  # Shape parameter (gamma, negative binomial)
  if (!is.null(params$shape)) {
    priors[["shape"]] <- c(mean = params$shape, sd = params$shape * inflation_factor)
  }
  
  attr(priors, "inflation_factor") <- inflation_factor
  attr(priors, "source_method") <- params$method
  if (!is.null(family_config)) {
    attr(priors, "family") <- family_config$name
  }
  
  return(priors)
}

#' Create brms Prior Objects from Prior List
#'
#' Converts a list of prior specifications into brms prior objects.
#' Handles different families and their specific parameters.
#'
#' @param prior_list List of prior specifications from \code{extract_and_inflate_priors}.
#' @param formula Optional formula to match priors against.
#'
#' @return A vector of brms prior specifications or a data frame of priors.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' prior_list <- extract_and_inflate_priors(model_freq)
#' brms_priors <- create_brms_priors(prior_list)
#' }
create_brms_priors <- function(prior_list, formula = NULL) {
  if (!requireNamespace("brms", quietly = TRUE)) {
    stop("brms package is required for creating brms priors.")
  }
  
  prior_specs <- list()
  
  # Get family from attributes if available
  family_name <- attr(prior_list, "family")
  
  for (param_name in names(prior_list)) {
    if (param_name %in% c("inflation_factor", "source_method", "family")) {
      next  # Skip attributes
    }
    
    prior_vals <- prior_list[[param_name]]
    
    if (grepl("^sd_", param_name)) {
      # Standard deviation priors (half-normal)
      param_clean <- gsub("^sd_", "", param_name)
      prior_specs[[length(prior_specs) + 1]] <- brms::prior_string(
        sprintf("normal(0, %.4f)", prior_vals["sd"]),
        class = "sd"
      )
    } else if (param_name == "b_Intercept") {
      # Intercept
      prior_specs[[length(prior_specs) + 1]] <- brms::prior_string(
        sprintf("normal(%.4f, %.4f)", prior_vals["mean"], prior_vals["sd"]),
        class = "Intercept"
      )
    } else if (grepl("^b_", param_name)) {
      # Fixed effect
      coef_name <- gsub("^b_", "", param_name)
      prior_specs[[length(prior_specs) + 1]] <- brms::prior_string(
        sprintf("normal(%.4f, %.4f)", prior_vals["mean"], prior_vals["sd"]),
        class = "b",
        coef = coef_name
      )
    } else if (param_name == "zi_Intercept") {
      # ZI intercept
      prior_specs[[length(prior_specs) + 1]] <- brms::prior_string(
        sprintf("normal(%.4f, %.4f)", prior_vals["mean"], prior_vals["sd"]),
        class = "Intercept",
        dpar = "zi"
      )
    } else if (grepl("^zi_", param_name) && param_name != "zi_Intercept") {
      # ZI coefficient
      coef_name <- gsub("^zi_", "", param_name)
      prior_specs[[length(prior_specs) + 1]] <- brms::prior_string(
        sprintf("normal(%.4f, %.4f)", prior_vals["mean"], prior_vals["sd"]),
        class = "b",
        coef = coef_name,
        dpar = "zi"
      )
    } else if (param_name == "disp_Intercept") {
      # Dispersion intercept
      prior_specs[[length(prior_specs) + 1]] <- brms::prior_string(
        sprintf("normal(%.4f, %.4f)", prior_vals["mean"], prior_vals["sd"]),
        class = "Intercept",
        dpar = "disp"
      )
    } else if (param_name == "phi") {
      # Precision parameter (beta families)
      shape <- prior_vals["mean"]^2 / prior_vals["sd"]^2
      rate <- prior_vals["mean"] / prior_vals["sd"]^2
      prior_specs[[length(prior_specs) + 1]] <- brms::prior_string(
        sprintf("gamma(%.4f, %.4f)", shape, rate),
        class = "phi"
      )
    } else if (param_name == "shape") {
      # Shape parameter (gamma, negative binomial)
      shape_alpha <- prior_vals["mean"]^2 / prior_vals["sd"]^2
      shape_beta <- prior_vals["mean"] / prior_vals["sd"]^2
      prior_specs[[length(prior_specs) + 1]] <- brms::prior_string(
        sprintf("gamma(%.4f, %.4f)", shape_alpha, shape_beta),
        class = "shape"
      )
    }
  }
  
  # Combine into brms prior object
  if (length(prior_specs) > 0) {
    do.call(c, prior_specs)
  } else {
    NULL
  }
}
