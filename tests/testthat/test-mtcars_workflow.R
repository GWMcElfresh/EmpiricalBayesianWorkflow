test_that("End-to-end workflow runs with Gaussian family on mtcars", {
    skip_on_cran()

    # Prepare data
    data(mtcars)
    # Extract manufacturer from rownames (first word)
    mtcars$manufacturer <- factor(gsub(" .*", "", rownames(mtcars)))

    # Use subset of data for speed if needed, but mtcars is small (32 rows)
    # Ensure we have enough data for splitting
    n <- nrow(mtcars)

    # Run workflow
    # Predicting mpg from hp and wt, with random intercept for manufacturer
    # Using gaussian family

    # We need to handle the issue where some splits might have missing random effect levels
    # but run_workflow handles splitting internally.

    # Set seed for reproducibility
    set.seed(123)

    results <- tryCatch(
        {
            run_workflow(
                data = mtcars,
                formula = mpg ~ hp + wt + (1 | manufacturer),
                family = "gaussian",
                n_subjects = length(unique(mtcars$manufacturer)), # Ignored when data is provided
                n_obs_per_subject = 1, # Ignored
                verbose = FALSE,
                chains = 1, # Reduce chains/iter for speed in test
                iter = 1000,
                warmup = 500
            )
        },
        error = function(e) {
            fail(paste("Workflow failed with error:", e$message))
            NULL
        }
    )

    skip_if(is.null(results))

    # Verify structure
    expect_s3_class(results, "workflow_results")
    expect_true(!is.null(results$model_freq))
    expect_true(!is.null(results$model_bayes))
    expect_true(!is.null(results$conformal_results))

    # Verify Frequentist model
    # glmmTMB gaussian family
    expect_s3_class(results$model_freq, "glmmTMB")
    expect_equal(family(results$model_freq)$family, "gaussian")

    # Verify Bayesian model
    expect_s3_class(results$model_bayes, "brmsfit")
    expect_equal(family(results$model_bayes)$family, "gaussian")

    # Verify Conformal Prediction
    cp <- results$conformal_results
    expect_s3_class(cp, "conformal_prediction")
    expect_true(is.numeric(cp$coverage))
    expect_true(cp$coverage >= 0 && cp$coverage <= 1)

    # Check coverage is reasonable (should be around 0.9 for alpha=0.1)
    # With small sample size (mtcars n=32), split conformal might vary
    # But shouldn't be 0 or 1 completely (unless lucky/unlucky)
    expect_true(cp$coverage > 0.5)

    # Check priors were constructed
    priors <- results$priors
    expect_type(priors, "list")
    expect_true("sigma" %in% names(priors))
    expect_true("b_Intercept" %in% names(priors))
})
