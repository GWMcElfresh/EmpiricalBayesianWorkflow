#' Generate Example Data for Zero-Inflated Beta Distribution
#'
#' Creates synthetic data with zero-inflated beta distributed responses,
#' random effects for subjects (intercepts and slopes), and a continuous predictor.
#'
#' @param n_subjects Number of subjects (groups). Default is 50.
#' @param n_obs_per_subject Number of observations per subject. Default is 10.
#' @param zi_prob Zero-inflation probability. Default is 0.2.
#' @param seed Random seed for reproducibility. If NULL, uses current RNG state.
#'
#' @return A data frame with columns:
#'   \item{y}{Response variable, bounded in [0, 1] with possible zeros}
#'   \item{x}{Predictor variable (continuous)}
#'   \item{subject_id}{Grouping factor for random effects}
#'
#' @export
#'
#' @examples
#' # Generate data with default parameters
#' data <- generate_zibeta_data()
#' head(data)
#'
#' # Custom parameters
#' data <- generate_zibeta_data(n_subjects = 20, n_obs_per_subject = 15, zi_prob = 0.1)
#' summary(data)
generate_zibeta_data <- function(n_subjects = 50, 
                                 n_obs_per_subject = 10, 
                                 zi_prob = 0.2,
                                 seed = 123) {
  if (!is.null(seed)) {
    set.seed(seed)
  }
  
  # Input validation
  stopifnot(n_subjects > 0)
  stopifnot(n_obs_per_subject > 0)
  stopifnot(zi_prob >= 0 && zi_prob <= 1)
  
  # Create grouping structure
  subject_id <- rep(1:n_subjects, each = n_obs_per_subject)
  n_total <- n_subjects * n_obs_per_subject
  
  # Generate predictor
  x <- stats::rnorm(n_total)
  
  # Random effects for intercept and slope
  re_intercept <- stats::rnorm(n_subjects, 0, 0.5)
  re_slope <- stats::rnorm(n_subjects, 0, 0.3)
  
  # Linear predictor for beta mean (logit scale)
  # Fixed effects: intercept = 0.5, slope = 0.8
  mu_logit <- 0.5 + 0.8 * x + re_intercept[subject_id] + re_slope[subject_id] * x
  mu <- stats::plogis(mu_logit)
  
  # Precision parameter (phi)
  phi <- 10
  
  # Generate beta distributed values
  shape1 <- mu * phi
  shape2 <- (1 - mu) * phi
  y_beta <- stats::rbeta(n_total, shape1, shape2)
  
  # Add zero-inflation
  zi_indicator <- stats::rbinom(n_total, 1, zi_prob)
  y <- ifelse(zi_indicator == 1, 0, y_beta)
  
  # Create data frame
  data <- data.frame(
    y = y,
    x = x,
    subject_id = factor(subject_id)
  )
  
  return(data)
}
