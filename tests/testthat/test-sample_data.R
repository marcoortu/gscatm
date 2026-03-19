test_that("sample data generator creates valid output", {
  sample_data <- generate_sample_data()

  # Check dimensions
  expect_equal(dim(sample_data$dtm), c(10, 50))
  expect_equal(dim(sample_data$covariates), c(10, 3))
  expect_equal(dim(sample_data$true_topics), c(3, 50))
  expect_equal(dim(sample_data$true_theta), c(10, 3))

  # Check probability constraints
  expect_equal(rowSums(sample_data$true_topics), rep(1, 3))
  expect_equal(rowSums(sample_data$true_theta), rep(1, 10))

  # Check data types
  expect_s4_class(sample_data$dtm, "dgCMatrix")
  expect_true(is.data.frame(sample_data$covariates))


  # Generate sample data
  set.seed(42)
  sample_data <- generate_sample_data(
    n_docs = 10,
    n_terms = 50,
    n_topics = 3,
    n_covariates = 3
  )

  # Fit model with default reference topic (last topic)
  model1 <- fit_gsca_topic_model(
    dtm = sample_data$dtm,
    covariates = sample_data$covariates,
    K = 3,
    model = "zero-inflated"
  )

  # Check model components
  expect_s3_class(model1, "gsca_topic_model")

  # Fit model with custom reference topic
  model2 <- fit_gsca_topic_model(
    dtm = sample_data$dtm,
    covariates = sample_data$covariates,
    K = 3,
    model = "logistic-normal",
    reference_topic = 1  # Use first topic as reference
  )

  # Check model components
  expect_s3_class(model2, "gsca_topic_model")

})
