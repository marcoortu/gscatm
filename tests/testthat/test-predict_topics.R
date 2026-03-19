test_that("predict_topics works correctly", {
  library(gscatm)
  library(Matrix)
  library(quanteda)

  # Prepare sample data
  docs <- c("This is a sample document.",
            "This document is another example.",
            "Topic modeling is fun.",
            "We are integrating GSCA with topic modeling.")
  dfm <- dfm(tokens(docs))
  dtm <- as(dfm, "dgCMatrix")
  covariates <- data.frame(length = rowSums(dtm))

  # Fit the model
  K <- 2
  model <- fit_gsca_topic_model(dtm, covariates, K = K)

  # Create new documents for prediction
  new_docs <- c("This is a new document.",
                "Another new example.")

  # Instead of using dfm_lookup with dtm, create dfm with same features as original
  new_dfm <- dfm(tokens(new_docs))
  # Ensure vocabulary alignment with original dfm
  new_dfm <- dfm_match(new_dfm, featnames(dfm))
  new_dtm <- as(new_dfm, "dgCMatrix")
  new_covariates <- data.frame(length = rowSums(new_dtm))

  # Make predictions
  theta_pred <- predict_topics(model, new_dtm, new_covariates)

  # Check basic properties
  expect_true(is.matrix(theta_pred))
  expect_equal(dim(theta_pred), c(nrow(new_dtm), K))
  expect_true(all(abs(rowSums(theta_pred) - 1) < 1e-6))
  expect_true(all(theta_pred >= 0 & theta_pred <= 1))
})

test_that("predict_topics handles errors correctly", {
  library(gscatm)
  library(Matrix)
  library(quanteda)

  # Prepare sample data with more documents to avoid empty cases
  docs <- c("This is a sample document.",
            "This document is another example.",
            "Topic modeling is fun.",
            "We are integrating GSCA with topic modeling.")
  dfm <- dfm(tokens(docs))
  dtm <- as(dfm, "dgCMatrix")
  covariates <- data.frame(length = rowSums(dtm))

  # Fit model
  K <- 2
  model <- fit_gsca_topic_model(dtm, covariates, K = K)

  # Test with mismatched dimensions
  new_docs <- c("This is a new document.")
  # Create dfm with same features as original
  new_dfm <- dfm(tokens(new_docs))
  new_dfm <- dfm_match(new_dfm, featnames(dfm))
  new_dtm <- as(new_dfm, "dgCMatrix")
  bad_covariates <- data.frame(length = c(1, 2))  # Wrong number of rows

  expect_error(predict_topics(model, new_dtm, bad_covariates))

  # Test with wrong model class
  fake_model <- list(K = 2)
  class(fake_model) <- "fake_model"
  expect_error(predict_topics(fake_model, new_dtm, covariates))
})
