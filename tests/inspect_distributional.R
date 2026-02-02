
library(glmmTMB)
library(brms)

# Helper to inspect distributional models
inspect_dist <- function(name, fam, gen_data_fn) {
    cat(sprintf("=== %s Distributional Model ===\n", name))
    
    set.seed(123)
    data <- gen_data_fn()
    
    # Fit with glmmTMB using dispformula
    # Note: glmmTMB uses dispformula for the dispersion parameter.
    # For Gaussian, this creates a model for log(sigma).
    # For Beta, it models phi.
    
    cat("Fitting glmmTMB...\n")
    # Using explicit dispformula ~ x
    fit <- tryCatch(
        glmmTMB(y ~ x, dispformula = ~ x, data = data, family = fam),
        error = function(e) e
    )
    
    if (inherits(fit, "glmmTMB")) {
        fe <- fixef(fit)
        cat("Fixed Effects (cond):\n")
        print(fe$cond)
        
        cat("Fixed Effects (disp):\n")
        print(fe$disp) # This should contain intercept and x
        
        # Check what brms expects for this family
        cat("\nBrms default priors for comparison:\n")
        # Define brms equivalent
        bf_obj <- switch(name,
            "gaussian" = bf(y ~ x, sigma ~ x),
            "beta" = bf(y ~ x, phi ~ x),
            "gamma" = bf(y ~ x, shape ~ x)
        )
        
        # Approximate brms family
        brms_fam <- switch(name,
            "gaussian" = brmsfamily("gaussian"),
            "beta" = brmsfamily("Beta"),
            "gamma" = brmsfamily("Gamma", link="log")
        )

        priors <- get_prior(bf_obj, data=data, family=brms_fam)
        # Filter for relevant dpar
        dpars <- c("sigma", "phi", "shape")
        print(priors[priors$dpar %in% dpars, c("prior", "class", "coef", "dpar")])
        
    } else {
        cat("glmmTMB fit failed:", fit$message, "\n")
    }
    cat("\n")
}

# Generators
gen_linear <- function() {
    N <- 100
    x <- rnorm(N)
    # Heteroscedastic gaussian
    sigma <- exp(0.5 + 0.5 * x)
    y <- rnorm(N, mean = 2 + 1*x, sd = sigma)
    data.frame(y=y, x=x)
}

gen_beta <- function() {
    N <- 200
    x <- rnorm(N)
    mu <- plogis(0.5 + 0.5 * x)
    phi <- exp(2 + 0.5 * x)
    # rbeta param
    y <- rbeta(N, mu*phi, (1-mu)*phi)
    # Avoid 0/1
    y <- (y * (N - 1) + 0.5) / N
    data.frame(y=y, x=x)
}

gen_gamma <- function() {
    N <- 100
    x <- rnorm(N)
    mu <- exp(1 + 0.5 * x)
    # shape modeled as function of x
    shape <- exp(1 + 0.5 * x) 
    # gamma var = mu^2 / shape. glmmTMB disp = 1/shape? 
    # actually glmmTMB dispformula models the dispersion parameter.
    # For Gamma, Var = mu^2 * dispersion_param.  dispersion_param = 1/shape.
    # So if we use dispformula, we are modeling log(1/shape) = -log(shape).
    # brms models log(shape).
    # So coefficients should be inverted (sign flip) if modeling shape ~ x vs disp ~ x.
    
    y <- rgamma(N, shape = shape, scale = mu/shape)
    data.frame(y=y, x=x)
}

inspect_dist("gaussian", gaussian(), gen_linear)
inspect_dist("beta", beta_family(), gen_beta)
inspect_dist("gamma", Gamma(link="log"), gen_gamma)
