#' Family Mapping and Configuration
#'
#' This module handles the mapping between glmmTMB and brms families,
#' ensuring exact parameter correspondence.

#' Get Supported Families
#'
#' Returns a list of families supported by the workflow with their
#' parameter mappings between glmmTMB and brms.
#'
#' @return A list of family configurations
#' @export
#'
#' @examples
#' families <- get_supported_families()
#' names(families)
get_supported_families <- function() {
  list(
    beta = list(
      name = "beta",
      glmmTMB_family = "beta_family()",
      brms_family = "Beta()",
      has_zi = FALSE,
      has_hurdle = FALSE,
      params = c("phi"),
      param_classes = list(phi = "phi"),
      description = "Beta distribution for (0,1) responses"
    ),
    gaussian = list(
      name = "gaussian",
      glmmTMB_family = "gaussian()",
      brms_family = "gaussian()",
      has_zi = FALSE,
      has_hurdle = FALSE,
      params = c("sigma"),
      param_classes = list(sigma = "sigma"),
      description = "Gaussian distribution for continuous responses"
    ),
    zero_inflated_beta = list(
      name = "zero_inflated_beta",
      glmmTMB_family = "beta_family()",
      brms_family = "zero_inflated_beta()",
      has_zi = TRUE,
      has_hurdle = FALSE,
      params = c("phi"),
      param_classes = list(phi = "phi"),
      zi_link = "logit",
      description = "Zero-inflated beta for [0,1) with structural zeros"
    ),
    gamma = list(
      name = "gamma",
      glmmTMB_family = "Gamma(link='log')",
      brms_family = "Gamma(link='log')",
      has_zi = FALSE,
      has_hurdle = FALSE,
      params = c("shape"),
      param_classes = list(shape = "shape"),
      description = "Gamma distribution for positive continuous responses"
    ),
    binomial = list(
      name = "binomial",
      glmmTMB_family = "binomial()",
      brms_family = "binomial()",
      has_zi = FALSE,
      has_hurdle = FALSE,
      params = character(0),
      param_classes = list(),
      description = "Binomial for binary or proportion responses"
    ),
    poisson = list(
      name = "poisson",
      glmmTMB_family = "poisson()",
      brms_family = "poisson()",
      has_zi = FALSE,
      has_hurdle = FALSE,
      params = character(0),
      param_classes = list(),
      description = "Poisson for count data"
    ),
    zero_inflated_poisson = list(
      name = "zero_inflated_poisson",
      glmmTMB_family = "poisson()",
      brms_family = "zero_inflated_poisson()",
      has_zi = TRUE,
      has_hurdle = FALSE,
      params = character(0),
      param_classes = list(),
      zi_link = "logit",
      description = "Zero-inflated Poisson for count data with excess zeros"
    ),
    nbinom2 = list(
      name = "nbinom2",
      glmmTMB_family = "nbinom2()",
      brms_family = "negbinomial()",
      has_zi = FALSE,
      has_hurdle = FALSE,
      params = c("shape"),
      param_classes = list(shape = "shape"),
      description = "Negative binomial (NB2 parameterization) for overdispersed counts"
    ),
    zero_inflated_nbinom2 = list(
      name = "zero_inflated_nbinom2",
      glmmTMB_family = "nbinom2()",
      brms_family = "zero_inflated_negbinomial()",
      has_zi = TRUE,
      has_hurdle = FALSE,
      params = c("shape"),
      param_classes = list(shape = "shape"),
      zi_link = "logit",
      description = "Zero-inflated negative binomial for overdispersed counts with excess zeros"
    )
  )
}

#' Get Family Configuration
#'
#' Retrieves the configuration for a specific family, including parameter mappings.
#'
#' @param family_name Name of the family (e.g., "beta", "zero_inflated_beta")
#'
#' @return A list with family configuration
#' @export
get_family_config <- function(family_name) {
  families <- get_supported_families()

  if (!family_name %in% names(families)) {
    stop(sprintf(
      "Family '%s' not supported. Supported families: %s",
      family_name, paste(names(families), collapse = ", ")
    ))
  }

  families[[family_name]]
}

#' Map glmmTMB Family to brms Family
#'
#' Creates the appropriate brms family object based on family configuration.
#'
#' @param family_config Family configuration from get_family_config
#' @param link Optional link function override
#'
#' @return A brms family object
#' @export
create_brms_family <- function(family_config, link = NULL) {
  if (!requireNamespace("brms", quietly = TRUE)) {
    stop("brms package required")
  }

  family_name <- family_config$name

  # Create brms family object using brmsfamily() for correct API
  brms_family <- switch(family_name,
    "beta" = brms::brmsfamily("Beta", link = link %||% "logit"),
    "gaussian" = brms::brmsfamily("gaussian", link = link %||% "identity"),
    "zero_inflated_beta" = brms::brmsfamily("zero_inflated_beta", link = link %||% "logit"),
    "gamma" = brms::brmsfamily("Gamma", link = link %||% "log"),
    "binomial" = brms::brmsfamily("binomial", link = link %||% "logit"),
    "poisson" = brms::brmsfamily("poisson", link = link %||% "log"),
    "zero_inflated_poisson" = brms::brmsfamily("zero_inflated_poisson", link = link %||% "log"),
    "nbinom2" = brms::brmsfamily("negbinomial", link = link %||% "log"),
    "zero_inflated_nbinom2" = brms::brmsfamily("zero_inflated_negbinomial", link = link %||% "log"),
    stop(sprintf("Family '%s' not implemented", family_name))
  )

  return(brms_family)
}

#' Create glmmTMB Family
#'
#' Creates the appropriate glmmTMB family object.
#'
#' @param family_config Family configuration from get_family_config
#'
#' @return A glmmTMB family object
#' @export
create_glmmTMB_family <- function(family_config) {
  if (!requireNamespace("glmmTMB", quietly = TRUE)) {
    stop("glmmTMB package required")
  }

  family_name <- family_config$name

  # Create glmmTMB family
  glmmTMB_family <- switch(family_name,
    "beta" = glmmTMB::beta_family(),
    "gaussian" = stats::gaussian(),
    "zero_inflated_beta" = glmmTMB::beta_family(), # ZI handled via ziformula
    "gamma" = stats::Gamma(link = "log"),
    "binomial" = stats::binomial(),
    "poisson" = stats::poisson(),
    "zero_inflated_poisson" = stats::poisson(), # ZI handled via ziformula
    "nbinom2" = glmmTMB::nbinom2(),
    "zero_inflated_nbinom2" = glmmTMB::nbinom2(), # ZI handled via ziformula
    stop(sprintf("Family '%s' not implemented", family_name))
  )

  return(glmmTMB_family)
}

#' Map Parameter Names from glmmTMB to brms
#'
#' Ensures exact parameter name correspondence between packages.
#'
#' @param param_name Parameter name from glmmTMB
#' @param family_config Family configuration
#' @param component Component: "cond" (conditional), "zi" (zero-inflation), or "disp" (dispersion)
#'
#' @return Corresponding brms parameter name
#' @export
map_parameter_name <- function(param_name, family_config, component = "cond") {
  # Determine family-specific dispersion parameter name (sigma, phi, shape)
  # Default to "disp" if not found
  disp_name <- "disp"
  if (!is.null(family_config$param_classes)) {
     # Usually the first parameter in params list is the main dispersion/shape param
     # params = c("sigma") or c("phi") or c("shape")
     if (length(family_config$params) > 0) {
       disp_name <- family_config$params[1]
     }
  }

  prefix <- ""
  suffix <- param_name
  
  if (component == "zi") {
      prefix <- "zi"
  } else if (component == "disp") {
      prefix <- disp_name
  } else {
      # Conditional model
      if (param_name == "(Intercept)") return("Intercept")
      return(param_name)
  }

  # Clean up Intercept name
  if (param_name == "(Intercept)") {
      suffix <- "Intercept"
  }
  
  return(paste0(prefix, "_", suffix))
}

#' Extract Family-Specific Parameters
#'
#' Extracts dispersion/shape parameters specific to each family.
#'
#' @param model glmmTMB model object
#' @param family_config Family configuration
#'
#' @return Named list of family-specific parameters
#' @export
extract_family_parameters <- function(model, family_config) {
  params <- list()

  family_name <- family_config$name

  # Beta and zero-inflated beta: phi (precision)
  if (family_name %in% c("beta", "zero_inflated_beta")) {
    params$phi <- stats::sigma(model)
  }

  # Gaussian: sigma parameter
  else if (family_name == "gaussian") {
    params$sigma <- stats::sigma(model)
  }

  # Gamma: shape parameter
  else if (family_name == "gamma") {
    # glmmTMB returns dispersion (1/shape) as sigma
    s <- tryCatch(stats::sigma(model), error = function(e) NULL)

    if (!is.null(s)) {
      params$shape <- 1 / s
    } else {
      params$shape <- 1 # Default
    }
  }

  # Negative binomial: shape/size parameter
  else if (family_name %in% c("nbinom2", "zero_inflated_nbinom2")) {
    # glmmTMB returns theta (shape) as sigma
    s <- tryCatch(stats::sigma(model), error = function(e) NULL)

    if (!is.null(s)) {
      params$shape <- s
    } else {
      params$shape <- 1 # Default
    }
  }

  # No additional parameters for binomial, poisson, etc.

  return(params)
}

# Null-coalescing operator (internal use only)
# Not exported to avoid conflicts with rlang::`%||%`
`%||%` <- function(x, y) {
  if (is.null(x)) y else x
}
