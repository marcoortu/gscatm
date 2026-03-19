# R/predict_topics.R
#' Predict Topic Proportions
#'
#' @param model Fitted GSCA topic model
#' @param new_dtm New document-term matrix
#' @param new_covariates New covariates
#' @return Matrix of predicted topic proportions
#' @export
predict_topics <- function(model, new_dtm, new_covariates) {
  if (!inherits(model, "gsca_topic_model")) {
    stop("Model must be of class 'gsca_topic_model'.")
  }

  # Convert covariates consistently with fit_gsca_topic_model (model.matrix)
  if (is.data.frame(new_covariates)) {
    new_covariates <- model.matrix(~ . - 1, data = new_covariates)
  }

  # Get dimensions
  N_new <- nrow(new_dtm)
  P <- ncol(new_covariates)
  K <- model$K

  # Input validation
  if (nrow(new_covariates) != nrow(new_dtm)) {
    stop(sprintf("Number of rows in new_covariates (%d) must match number of documents in new_dtm (%d)",
                 nrow(new_covariates), nrow(new_dtm)))
  }

  # Check covariates dimension matches the model
  expected_P <- nrow(model$path_matrix) - model$K
  if (P != expected_P) {
    stop(sprintf("Number of columns in new_covariates (%d) must match original number of covariates (%d)",
                 P, expected_P))
  }

  # Initialize component scores for new data
  gamma_new <- matrix(0, N_new, P + K)

  # Set covariate components: standardize using training data statistics
  train_cov <- model$covariates
  cov_center <- colMeans(train_cov)
  cov_scale  <- apply(train_cov, 2, sd)
  cov_scale[cov_scale == 0] <- 1  # avoid division by zero for constant columns
  gamma_new[, 1:P] <- scale(new_covariates, center = cov_center, scale = cov_scale)

  # Compute topic components using path coefficients and weights
  for(k in 1:K) {
    topic_idx <- P + k
    predictors <- which(model$path_matrix[, topic_idx] != 0)

    if(length(predictors) > 0) {
      # First apply structural paths
      gamma_new[, topic_idx] <- gamma_new[, predictors, drop = FALSE] %*%
        model$path_matrix[predictors, topic_idx, drop = FALSE]
    }
  }

  # Convert to topic proportions based on model type
  topic_scores <- gamma_new[, (P+1):(P+K)]
  model_type   <- model$model_type

  if (is.null(model_type) || model_type == "logistic-normal") {
    theta <- exp(topic_scores)
    theta <- theta / rowSums(theta)

  } else if (model_type == "dirichlet") {
    if (!requireNamespace("MCMCpack", quietly = TRUE))
      stop("The 'MCMCpack' package is required for Dirichlet prediction.")
    dirichlet_params <- exp(topic_scores)
    theta <- t(apply(dirichlet_params, 1, function(a) MCMCpack::rdirichlet(1, a)))

  } else if (model_type == "zero-inflated") {
    sparsity_threshold <- if (!is.null(model$sparsity_threshold))
      model$sparsity_threshold else 1e-5
    theta <- exp(topic_scores)
    theta[theta < sparsity_threshold] <- 0
    theta <- theta / rowSums(theta)
  }

  # Final update using word probabilities
  if(ncol(new_dtm) == ncol(model$phi)) {
    log_probs <- as.matrix(new_dtm %*% t(log(model$phi + 1e-10)))
    theta <- theta * exp(log_probs)
    theta <- theta / rowSums(theta)
  } else {
    warning("Vocabulary mismatch between new documents and model. Using only structural predictions.")
  }

  return(theta)
}
