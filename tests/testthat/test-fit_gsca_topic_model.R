# tests/testthat/test-fit_gsca_topic_model.R

test_that("fit_gsca_topic_model works correctly", {
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

  # Create path matrix
  P <- ncol(covariates)
  K <- 2
  path_matrix <- matrix(0, P + K, P + K)
  path_matrix[1:P, (P+1):(P+K)] <- 1  # Paths from covariates to topics

  # Fit the model
  model <- fit_gsca_topic_model(dtm, covariates, K = K,
                                path_matrix = path_matrix,
                                max_iter = 10)

  # Basic model checks
  expect_s3_class(model, "gsca_topic_model")

  # Check dimensions of topic model components
  expect_equal(dim(model$phi), c(K, ncol(dtm)))
  expect_equal(dim(model$theta), c(nrow(dtm), K))

  # Check dimensions of GSCA components
  expect_equal(dim(model$gamma), c(nrow(dtm), P + K))
  expect_equal(dim(model$path_matrix), c(P + K, P + K))
  expect_equal(dim(model$weight_matrix), c(P + K, ncol(dtm)))

  # Check probability constraints
  expect_true(all(abs(rowSums(model$phi) - 1) < 1e-6))
  expect_true(all(abs(rowSums(model$theta) - 1) < 1e-6))

  # Check structural validity
  expect_true(all(model$gamma[, 1:P] - scale(as.matrix(covariates)) < 1e-6))
  expect_true(all(model$path_matrix[1:P, 1:P] == 0))  # No paths between covariates
})

test_that("fit_gsca_topic_model converges correctly", {
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

  # Parameters
  P <- ncol(covariates)
  K <- 2

  # Create path matrix
  path_matrix <- matrix(0, P + K, P + K)
  path_matrix[1:P, (P+1):(P+K)] <- 1

  # Test convergence
  model <- fit_gsca_topic_model(dtm, covariates, K = K,
                                path_matrix = path_matrix,
                                max_iter = 10)

  # Check model components
  expect_s3_class(model, "gsca_topic_model")
  expect_equal(dim(model$phi), c(K, ncol(dtm)))
  expect_equal(dim(model$theta), c(nrow(dtm), K))
  expect_equal(dim(model$gamma), c(nrow(dtm), P + K))

  # Check probability constraints
  expect_true(all(abs(rowSums(model$phi) - 1) < 1e-6))
  expect_true(all(abs(rowSums(model$theta) - 1) < 1e-6))

  # Check structural constraints
  expect_true(all(model$path_matrix[1:P, 1:P] == 0))
  expect_true(all(model$path_matrix[(P+1):(P+K), 1:P] == 0))

  # Test max_iter boundary
  model_quick <- fit_gsca_topic_model(dtm, covariates, K = K,
                                      path_matrix = path_matrix,
                                      max_iter = 1)
  expect_s3_class(model_quick, "gsca_topic_model")
})

test_that("fit_gsca_topic_model handles invalid inputs", {
  library(gscatm)
  library(Matrix)
  library(quanteda)

  # Prepare sample data
  docs <- c("This is a sample document.",
            "This document is another example.")

  dfm <- dfm(tokens(docs))
  dtm <- as(dfm, "dgCMatrix")
  covariates <- data.frame(length = rowSums(dtm))

  # Test invalid K
  expect_error(fit_gsca_topic_model(dtm, covariates, K = 0))

  # Test invalid path matrix
  P <- ncol(covariates)
  K <- 2
  bad_path_matrix <- matrix(0, P, P)  # Wrong dimensions
  expect_error(fit_gsca_topic_model(dtm, covariates, K = K,
                                    path_matrix = bad_path_matrix))

})

test_that("fit_gsca_topic_model handles handles zero-inflated topic distribution", {
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

  # Parameters
  P <- ncol(covariates)
  K <- 2

  # Create path matrix
  path_matrix <- matrix(0, P + K, P + K)
  path_matrix[1:P, (P+1):(P+K)] <- 1

  # Test convergence
  model <- fit_gsca_topic_model(dtm, covariates, K = K,
                                path_matrix = path_matrix,
                                model = "zero-inflated",
                                max_iter = 10)
  # Check model components
  expect_s3_class(model, "gsca_topic_model")
  expect_equal(dim(model$phi), c(K, ncol(dtm)))
  expect_equal(dim(model$theta), c(nrow(dtm), K))
  expect_equal(dim(model$gamma), c(nrow(dtm), P + K))

})

test_that("fit_gsca_topic_model handles handles dirichlet topic distribution", {
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

  # Parameters
  P <- ncol(covariates)
  K <- 2

  # Create path matrix
  path_matrix <- matrix(0, P + K, P + K)
  path_matrix[1:P, (P+1):(P+K)] <- 1

  # Test convergence
  model <- fit_gsca_topic_model(dtm, covariates, K = K,
                                path_matrix = path_matrix,
                                model = "dirichlet",
                                max_iter = 10)
  # Check model components
  expect_s3_class(model, "gsca_topic_model")
  expect_equal(dim(model$phi), c(K, ncol(dtm)))
  expect_equal(dim(model$theta), c(nrow(dtm), K))
  expect_equal(dim(model$gamma), c(nrow(dtm), P + K))
})

test_that("fit_gsca_topic_model effects table is correctly structured", {
  library(gscatm)
  set.seed(1)
  sample_data <- generate_sample_data(n_docs = 12, n_terms = 30,
                                      n_topics = 3, n_covariates = 2)
  model <- fit_gsca_topic_model(
    sample_data$dtm, sample_data$covariates, K = 3, max_iter = 10
  )
  eff <- model$effects
  P   <- model$P       # 2 covariates
  K   <- model$K       # 3 topics => K-1 = 2 non-ref topics

  # Dimensions: P * (K-1) rows
  expect_equal(nrow(eff), P * (K - 1))

  # Topic column must repeat each non-ref topic P times
  expect_equal(
    eff$Topic,
    rep(paste0("Topic_", setdiff(1:K, K)), each = P)
  )

  # Covariate column must cycle through all covariate names K-1 times
  expect_equal(
    eff$Covariate,
    rep(colnames(model$covariates), K - 1)
  )

  # All numeric fields must be finite
  expect_true(all(is.finite(eff$Log_Odds_Coefficient)))
  expect_true(all(is.finite(eff$Std_Error)))
  expect_true(all(eff$Odds_Ratio > 0))
  expect_true(all(eff$CI_Lower <= eff$Odds_Ratio))
  expect_true(all(eff$CI_Upper >= eff$Odds_Ratio))

  # No bootstrap columns before bootstrap_ci = TRUE
  expect_false(model$bootstrap_ci_computed)
  expect_false("Boot_SE" %in% names(eff))
})

test_that("parametric bootstrap adds correct columns to effects table", {
  library(gscatm)
  set.seed(42)
  sample_data <- generate_sample_data(n_docs = 15, n_terms = 30,
                                      n_topics = 3, n_covariates = 2)
  model_boot <- fit_gsca_topic_model(
    sample_data$dtm, sample_data$covariates, K = 3,
    max_iter = 10, bootstrap_ci = TRUE, n_bootstrap = 5
  )

  expect_true(model_boot$bootstrap_ci_computed)

  eff <- model_boot$effects
  # Bootstrap columns present
  expect_true(all(c("Boot_SE", "Boot_CI_Lower", "Boot_CI_Upper") %in% names(eff)))

  # Boot_SE must be non-negative
  expect_true(all(eff$Boot_SE >= 0, na.rm = TRUE))

  # Boot CIs on OR scale: lower <= OR <= upper
  expect_true(all(eff$Boot_CI_Lower <= eff$Odds_Ratio, na.rm = TRUE))
  expect_true(all(eff$Boot_CI_Upper >= eff$Odds_Ratio, na.rm = TRUE))

  # bootstrap list components
  boot <- model_boot$bootstrap
  expect_true(all(c("boot_coefs", "boot_se", "boot_mean", "n_valid") %in% names(boot)))
  expect_equal(nrow(boot$boot_coefs), 5L)
  expect_equal(ncol(boot$boot_coefs), model_boot$P * (model_boot$K - 1))
})

test_that("fit_gsca_topic_model handles wrong topic model specificaiton", {
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

  # Parameters
  P <- ncol(covariates)
  K <- 2

  # Create path matrix
  path_matrix <- matrix(0, P + K, P + K)
  path_matrix[1:P, (P+1):(P+K)] <- 1

  # Test convergence
  expect_error(fit_gsca_topic_model(dtm, covariates, K = K,
                                path_matrix = path_matrix,
                                model = "none",
                                max_iter = 10))
})
