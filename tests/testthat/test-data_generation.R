# Tests for data generation functions

test_that("generate_zibeta_data creates correct structure", {
  data <- generate_zibeta_data(n_subjects = 5, n_obs_per_subject = 4, zi_prob = 0.2)
  
  expect_s3_class(data, "data.frame")
  expect_equal(nrow(data), 20)
  expect_equal(ncol(data), 3)
  expect_true(all(c("y", "x", "subject_id") %in% names(data)))
})

test_that("generate_zibeta_data produces valid responses", {
  data <- generate_zibeta_data(n_subjects = 10, n_obs_per_subject = 10)
  
  expect_true(all(data$y >= 0))
  expect_true(all(data$y <= 1))
  expect_true(is.numeric(data$y))
  expect_true(is.numeric(data$x))
  expect_s3_class(data$subject_id, "factor")
})

test_that("generate_zibeta_data respects zero-inflation", {
  # With high ZI probability
  data_high_zi <- generate_zibeta_data(n_subjects = 20, n_obs_per_subject = 20, zi_prob = 0.5, seed = 123)
  prop_zeros <- mean(data_high_zi$y == 0)
  
  # Should have some zeros
  expect_true(prop_zeros > 0)
  
  # With zero ZI probability
  data_no_zi <- generate_zibeta_data(n_subjects = 20, n_obs_per_subject = 20, zi_prob = 0, seed = 456)
  prop_zeros_no_zi <- mean(data_no_zi$y == 0)
  
  # Should have fewer or no structural zeros
  expect_true(prop_zeros_no_zi < prop_zeros)
})

test_that("generate_zibeta_data seed works", {
  data1 <- generate_zibeta_data(n_subjects = 5, n_obs_per_subject = 4, seed = 789)
  data2 <- generate_zibeta_data(n_subjects = 5, n_obs_per_subject = 4, seed = 789)
  
  expect_equal(data1$y, data2$y)
  expect_equal(data1$x, data2$x)
})

test_that("generate_zibeta_data validates inputs", {
  expect_error(generate_zibeta_data(n_subjects = 0))
  expect_error(generate_zibeta_data(n_obs_per_subject = -1))
  expect_error(generate_zibeta_data(zi_prob = 1.5))
  expect_error(generate_zibeta_data(zi_prob = -0.1))
})
