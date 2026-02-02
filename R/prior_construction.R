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
        priors[["b_Intercept"]] <- c(mean = unname(mean_val), sd = unname(sd_val))
      } else {
        priors[[paste0("b_", brms_name)]] <- c(mean = unname(mean_val), sd = unname(sd_val))
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

      priors[[brms_name]] <- c(mean = unname(mean_val), sd = unname(sd_val))
    }
  }

  # Fixed effects priors for dispersion part (if present)
  # parsing: parameter name matching relying on map_parameter_name
  if (!is.null(params$disp_fixef)) {
    for (i in seq_along(params$disp_fixef)) {
      param_name <- names(params$disp_fixef)[i]
      mean_val <- params$disp_fixef[i]
      sd_val <- abs(mean_val) * inflation_factor + 0.5

      # Map to brms name (e.g. sigma_Intercept, phi_x)
      if (!is.null(family_config)) {
         brms_name <- map_parameter_name(param_name, family_config, "disp")
         priors[[brms_name]] <- c(mean = unname(mean_val), sd = unname(sd_val))
      } else {
         # Fallback if no config
         brms_name <- if (param_name == "(Intercept)") "disp_Intercept" else paste0("disp_", param_name)
         priors[[brms_name]] <- c(mean = unname(mean_val), sd = unname(sd_val))
      }
    }
  }

  # Random effects priors
  if (!is.null(params$re_sd)) {
    for (i in seq_along(params$re_sd)) {
      param_name <- names(params$re_sd)[i]
      mean_val <- params$re_sd[i]
      sd_val <- mean_val * inflation_factor

      priors[[paste0("sd_", param_name)]] <- c(mean = 0, sd = unname(sd_val))
    }
  }

  # Family-specific parameters (SCALAR only)
  # Only add if NOT already covered by disp_fixef (which implies a formula).
  # glmmTMB puts scalar dispersion params in disp_fixef in recent versions?
  # Or does it use fit$sigma?
  # If we have disp_fixef, we prioritize that distribution regression over the scalar fallback.
  
  has_disp_model <- !is.null(params$disp_fixef) && length(params$disp_fixef) > 0

  if (!has_disp_model) {
      if (!is.null(params$phi)) {
        priors[["phi"]] <- c(mean = unname(params$phi), sd = unname(params$phi * inflation_factor))
      }
      if (!is.null(params$shape)) {
        priors[["shape"]] <- c(mean = unname(params$shape), sd = unname(params$shape * inflation_factor))
      }
      if (!is.null(params$sigma)) {
        priors[["sigma"]] <- c(mean = unname(params$sigma), sd = unname(params$sigma * inflation_factor))
      }
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

  # Helper to create normal prior string
  p_norm <- function(m, s) sprintf("normal(%.4f, %.4f)", m, s)

  for (param_name in names(prior_list)) {
    if (param_name %in% c("inflation_factor", "source_method", "family")) next

    prior_vals <- prior_list[[param_name]]
    
    # regex matches
    
    if (grepl("^sd_", param_name)) {
      # Random effects SD
      param_clean <- gsub("^sd_", "", param_name)
      prior_specs[[length(prior_specs) + 1]] <- brms::prior_string(
        sprintf("normal(0, %.4f)", prior_vals["sd"]),
        class = "sd"
      )
      
    } else if (param_name == "b_Intercept") {
       prior_specs[[length(prior_specs) + 1]] <- brms::prior_string(
        p_norm(prior_vals["mean"], prior_vals["sd"]),
        class = "Intercept"
      )
      
    } else if (grepl("^b_", param_name)) {
       coef_name <- gsub("^b_", "", param_name)
       prior_specs[[length(prior_specs) + 1]] <- brms::prior_string(
        p_norm(prior_vals["mean"], prior_vals["sd"]),
        class = "b", coef = coef_name
      )
      
    } else if (grepl("^zi_", param_name)) {
       # Zero-inflation
       is_intercept <- grepl("_Intercept$", param_name)
       if (is_intercept) {
           prior_specs[[length(prior_specs) + 1]] <- brms::prior_string(
            p_norm(prior_vals["mean"], prior_vals["sd"]),
            class = "Intercept", dpar = "zi"
          )
       } else {
           coef_name <- gsub("^zi_", "", param_name)
           prior_specs[[length(prior_specs) + 1]] <- brms::prior_string(
            p_norm(prior_vals["mean"], prior_vals["sd"]),
            class = "b", coef = coef_name, dpar = "zi"
          )
       }
       
    } else if (grepl("^(sigma|phi|shape|disp)_", param_name)) {
       # Distributional parameters (sigma, phi, shape, or generic disp)
       parts <- strsplit(param_name, "_")[[1]]
       dpar_name <- parts[1]
       suffix <- paste(parts[-1], collapse="_")
       
       if (suffix == "Intercept") {
           prior_specs[[length(prior_specs) + 1]] <- brms::prior_string(
            p_norm(prior_vals["mean"], prior_vals["sd"]),
            class = "Intercept", dpar = dpar_name
          )
       } else {
           prior_specs[[length(prior_specs) + 1]] <- brms::prior_string(
            p_norm(prior_vals["mean"], prior_vals["sd"]),
            class = "b", coef = suffix, dpar = dpar_name
          )
       }

    } else if (param_name %in% c("phi", "shape", "sigma")) {
       # Scalar parameters (no formula model)
       # Use Gamma prior for these positive parameters
       shape_alpha <- prior_vals["mean"]^2 / prior_vals["sd"]^2
       shape_beta <- prior_vals["mean"] / prior_vals["sd"]^2
       prior_specs[[length(prior_specs) + 1]] <- brms::prior_string(
         sprintf("gamma(%.4f, %.4f)", shape_alpha, shape_beta),
         class = param_name
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
