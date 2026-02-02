# Script to inspect and record default priors and parameter mappings
# Saves output to tests/prior_reference.txt

library(EmpiricalBayesianWorkflow)
library(brms)
library(glmmTMB)
library(dplyr)

output_file <- "tests/prior_reference.txt"
sink(output_file)

cat("Prior Inspection and Mapping Reference\n")
cat("=====================================\n\n")

# Helper function to inspect a family
inspect_family <- function(family_config_name, formula, data, zi_formula = NULL, brms_fam_args = list(), glmmTMB_fam_args = list()) {
    cat(sprintf("### Family: %s\n", family_config_name))
    cat(rep("-", 40), "\n")

    # Get configuration
    config <- get_family_config(family_config_name)
    cat("Config Description: ", config$description, "\n\n")

    # 1. Create family objects
    brms_fam <- tryCatch(create_brms_family(config), error = function(e) e$message)
    glmmTMB_fam <- tryCatch(create_glmmTMB_family(config), error = function(e) e$message)

    # 2. brms::get_prior
    cat("--- brms::get_prior ---\n")

    # Handle brmsformula vs standard formula
    brms_formula <- formula

    priors <- tryCatch(
        {
            get_prior(brms_formula, data = data, family = brms_fam)
        },
        error = function(e) paste("Error:", e$message)
    )

    if (is.data.frame(priors)) {
        # Print relevant columns that exist
        cols_to_print <- c("prior", "class", "coef", "group", "resp", "dpar", "nlpar", "bound")
        valid_cols <- intersect(cols_to_print, names(priors))
        # Ensure we print valid columns only if we have rows or at least columns
        if (length(valid_cols) > 0) {
            print(priors[, valid_cols, drop = FALSE])
        } else {
            print(priors)
        }
    } else {
        print(priors)
    }
    cat("\n")

    # 3. glmmTMB parameters
    cat("--- glmmTMB Parameters ---\n")
    fit <- tryCatch(
        {
            # glmmTMB expects standard formula + ziformula
            # We separate main and ZI formulas.

            glmm_main_formula <- formula

            # Default ZI formula from argument
            glmm_zi_formula <- ~0
            if (!is.null(zi_formula)) {
                glmm_zi_formula <- zi_formula
            }

            if (inherits(formula, "brmsformula")) {
                glmm_main_formula <- formula$formula
                # If no explicit ZI formula passed, try to extract from bf
                if (is.null(zi_formula) && !is.null(formula$pforms$zi)) {
                    glmm_zi_formula <- formula$pforms$zi
                }
            }

            # Fit with minimal iterations to get parameter names structure
            glmmTMB(glmm_main_formula,
                data = data,
                family = glmmTMB_fam,
                ziformula = glmm_zi_formula,
                control = glmmTMBControl(optCtrl = list(maxit = 1))
            )
        },
        error = function(e) paste("Error setting up model:", e$message)
    )

    if (!is.character(fit)) {
        # Extract fixed effects using fixef()
        fe <- fixef(fit)

        cat("Fixed Effects (cond):\n")
        print(names(fe$cond))

        if (length(fe$zi) > 0) {
            cat("Zero-Inflation (zi):\n")
            print(names(fe$zi))
        }

        if (length(fe$disp) > 0) {
            cat("Dispersion (disp):\n")
            print(names(fe$disp))
        }

        # Check for auxiliary params (sigma/shape often in family specific or VarCorr?)
        # For glmmTMB, dispersion/shape matches often appear in sigma(fit) or family specific extraction
        # print(sigma(fit)) # dependent on family
    } else {
        print(fit)
    }
    cat("\n\n")
}

# --- Data Gen ---

set.seed(123)
N <- 50
data_gen <- data.frame(
    y_gauss = rnorm(N, 10, 2),
    y_beta = rbeta(N, 2, 5),
    y_gamma = rgamma(N, shape = 2, rate = 0.5),
    y_bin = rbinom(N, 1, 0.5),
    y_count = rpois(N, 5),
    y_nb = rnbinom(N, size = 2, mu = 5),
    x = rnorm(N),
    g = sample(letters[1:5], N, replace = TRUE)
)
# Zero inflated versions
data_gen$y_zibeta <- ifelse(runif(N) < 0.2, 0, rbeta(N, 2, 5))
data_gen$y_zinb <- ifelse(runif(N) < 0.2, 0, rnbinom(N, size = 2, mu = 5))

# --- Inspections ---

# 1. Gaussian
inspect_family("gaussian",
    formula = y_gauss ~ x + (1 | g),
    data = data_gen
)

# 2. Beta
inspect_family("beta",
    formula = y_beta ~ x + (1 | g),
    data = data_gen
)

# 3. Gamma
inspect_family("gamma",
    formula = y_gamma ~ x + (1 | g),
    data = data_gen
)

# 4. Binomial (Logistic)
inspect_family("binomial",
    formula = y_bin ~ x + (1 | g),
    data = data_gen
)

# 5. Poisson
inspect_family("poisson",
    formula = y_count ~ x + (1 | g),
    data = data_gen
)

# 6. Negative Binomial (nbinom2)
inspect_family("nbinom2",
    formula = y_nb ~ x + (1 | g),
    data = data_gen
)

# 7. Zero-Inflated Beta
# Note: ZI formula is usually ~1 by default in run_workflow
# Need to specify ZI formula for brms get_prior if we want to see zi parameters
# We pass explicit zi_formula for clarity and to simulate how workflow handles it internally if split.
inspect_family("zero_inflated_beta",
    formula = bf(y_zibeta ~ x + (1 | g), zi ~ 1),
    zi_formula = ~1,
    data = data_gen
)

# 8. Zero-Inflated Negative Binomial
inspect_family("zero_inflated_nbinom2",
    formula = bf(y_zinb ~ x + (1 | g), zi ~ 1),
    zi_formula = ~1,
    data = data_gen
)


sink()
cat("Inspection complete. Results in", output_file, "\n")
