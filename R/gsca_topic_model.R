# R/fit_gsca_topic_model.R
#' Fit GSCA Topic Model
#'
#' Fits a Generalized Structured Component Analysis (GSCA) Topic Model to document-term matrix
#' with covariates. The model combines topic modeling with structural equation modeling to
#' analyze relationships between document covariates and topic proportions.
#'
#' @param dtm Document-term matrix (dgCMatrix)
#' @param covariates Data frame or matrix of covariates
#' @param K Number of topics to extract
#' @param model Model type for topic distributions. One of "logistic-normal", "dirichlet",
#'        or "zero-inflated"
#' @param path_matrix Optional matrix specifying the structural paths between components.
#'        If NULL, defaults to paths from all covariates to all topics
#' @param weight_matrix Optional initial weight matrix for GSCA estimation.
#'        If NULL, initialized randomly
#' @param max_iter Maximum number of iterations for the EM algorithm. Default is 500
#' @param tol Convergence tolerance for parameter updates. Default is 1e-3
#' @param alpha Smoothing parameter for topic-word distributions. Default is 0.1
#' @param reference_topic Optional numeric value specifying which topic to use as reference
#'        category for effect estimation. Must be between 1 and K.
#'        If NULL, defaults to last topic (K)
#' @param bootstrap_ci Logical. If TRUE, run a parametric bootstrap to compute
#'        uncertainty-propagated standard errors and confidence intervals for the
#'        covariate effects (paper Section 4.2). Default is FALSE.
#' @param n_bootstrap Number of bootstrap replications. Default is 200.
#'        Ignored when bootstrap_ci = FALSE.
#'
#' @return A list containing the fitted model parameters:
#' \describe{
#'   \item{phi}{Matrix of topic-word distributions}
#'   \item{theta}{Matrix of document-topic proportions}
#'   \item{gamma}{Matrix of component scores}
#'   \item{path_matrix}{Final structural path coefficients}
#'   \item{weight_matrix}{Final weight matrix}
#'   \item{effects}{Data frame of covariate effects on topic proportions.
#'     When bootstrap_ci = TRUE, three additional columns are appended:
#'     Boot_SE, Boot_CI_Lower, Boot_CI_Upper.}
#'   \item{fit}{List of model fit statistics}
#'   \item{convergence}{List of convergence diagnostics}
#'   \item{bootstrap_ci_computed}{Logical indicating whether bootstrap CIs were computed}
#' }
#'
#' @export
#'
#' @examples
#' # Generate sample data
#' set.seed(42)
#' sample_data <- generate_sample_data(
#'   n_docs = 10,
#'   n_terms = 50,
#'   n_topics = 3,
#'   n_covariates = 3
#' )
#'
#' # Fit model with default reference topic (last topic)
#' model1 <- fit_gsca_topic_model(
#'   dtm = sample_data$dtm,
#'   covariates = sample_data$covariates,
#'   K = 3,
#'   model = "logistic-normal"
#' )
#'
#' # Fit model with bootstrap CIs (B = 50 for illustration)
#' \dontrun{
#' model2 <- fit_gsca_topic_model(
#'   dtm = sample_data$dtm,
#'   covariates = sample_data$covariates,
#'   K = 3,
#'   model = "logistic-normal",
#'   bootstrap_ci = TRUE,
#'   n_bootstrap = 50
#' )
#' model2$effects[, c("Covariate","Topic","Log_Odds_Coefficient",
#'                    "Boot_SE","Boot_CI_Lower","Boot_CI_Upper")]
#' }
fit_gsca_topic_model <- function(dtm, covariates, K,
                                 model = c("logistic-normal", "dirichlet", "zero-inflated"),
                                 path_matrix = NULL,
                                 weight_matrix = NULL,
                                 max_iter = 500,
                                 tol = 1e-3,
                                 alpha = 0.1,
                                 reference_topic = NULL,
                                 bootstrap_ci = FALSE,
                                 n_bootstrap = 200) {

  # Match arguments
  model <- match.arg(model)

  # Preserve original covariates data frame for bootstrap refitting
  covariates_df <- covariates

  # Basic setup and validation
  if (nrow(dtm) < 2) {
    stop("Need at least 2 documents to fit the model")
  }

  # Validate reference_topic
  if (is.null(reference_topic)) {
    reference_topic <- K  # Default to last topic
  } else {
    if (!is.numeric(reference_topic) || reference_topic < 1 || reference_topic > K) {
      stop(sprintf("reference_topic must be a number between 1 and %d", K))
    }
  }

  N <- nrow(dtm)
  V <- ncol(dtm)


  # covariates <- as.matrix(covariates)
  covariates <- model.matrix(~ . - 1, data = covariates)

  P <- ncol(covariates)

  # Initialize path matrix if not provided
  if (is.null(path_matrix)) {
    path_matrix <- matrix(0, K + P, K + P)
    path_matrix[1:P, (P+1):(P+K)] <- 1
  }

  # Initialize weight matrix if not provided
  if (is.null(weight_matrix)) {
    W_cov <- matrix(0, P, V)
    W_top <- matrix(runif(K * V), K, V)
    W_top <- W_top / rowSums(W_top)
    weight_matrix <- rbind(W_cov, W_top)
  }

  # Initialize topics
  phi <- matrix(runif(K * V), nrow = K, ncol = V)
  phi <- phi / rowSums(phi)

  # Initialize component scores
  gamma <- matrix(0, N, K + P)
  gamma[, 1:P] <- scale(covariates, center = TRUE, scale = TRUE)

  # Initialize theta
  topic_scores <- gamma[, (P+1):(P+K)]
  theta_new <- exp(topic_scores)
  theta <- theta_new / rowSums(theta_new)

  # Function to calculate perplexity
  calculate_perplexity <- function(dtm, phi, theta) {
    log_prob <- log(theta %*% phi + 1e-10)
    log_prob[is.infinite(log_prob)] <- -700
    log_per_word_prob <- sum(dtm * log_prob, na.rm = TRUE)
    total_words <- sum(dtm)
    perplexity <- exp(-log_per_word_prob / total_words)
    return(perplexity)
  }

  # Store perplexity history
  perplexity_history <- numeric(max_iter)
  best_perplexity <- Inf
  best_phi <- phi
  best_theta <- theta

  # EM Algorithm with GSCA
  for(iter in 1:max_iter) {
    # Store old values
    phi_old <- phi
    theta_old <- theta
    gamma_old <- gamma

    # E-Step
    # Compute theta from gamma
    # topic_scores <- gamma[, (P+1):(P+K)]
    # theta_new <- exp(topic_scores)
    # theta <- theta_new / rowSums(theta_new)
    topic_scores <- gamma[, (P+1):(P+K)]  # Calculate topic scores

    # Replace theta update logic here
    if (model == "logistic-normal") {
      theta_new <- exp(topic_scores)
      theta <- theta_new / rowSums(theta_new)
    } else if (model == "dirichlet") {
      if (!requireNamespace("MCMCpack", quietly = TRUE)) {
        stop("The 'MCMCpack' package is required for Dirichlet sampling.")
      }
      dirichlet_params <- exp(topic_scores)  # Ensure non-negativity
      theta <- t(apply(dirichlet_params, 1, MCMCpack::rdirichlet, n = 1))
    } else if (model == "zero-inflated") {
      sparsity_threshold <- 1e-5
      theta_new <- exp(topic_scores)
      theta_new[theta_new < sparsity_threshold] <- 0
      theta <- theta_new / rowSums(theta_new)
    } else {
      stop("Invalid model specified. Choose 'logistic-normal', 'dirichlet', or 'zero-inflated'.")
    }

    # Compute log-likelihood
    log_theta <- log(theta + 1e-10)
    log_phi <- log(phi + 1e-10)
    log_likelihood <- log_theta + dtm %*% t(log_phi)

    # Numerical stability
    log_likelihood <- sweep(log_likelihood, 1, apply(log_likelihood, 1, max), "-")

    # Compute posterior probabilities
    posterior <- exp(log_likelihood)
    posterior <- posterior / rowSums(posterior)

    # M-Step
    # Update phi
    for(k in 1:K) {
      weighted_counts <- dtm * posterior[, k]
      phi[k, ] <- colSums(weighted_counts) + alpha
      phi[k, ] <- phi[k, ] / sum(phi[k, ])
    }

    # Update regression parameters (GSCA Step)
    # Compute expected topic counts per document
    topic_counts <- posterior

    lambda <- 0.1 * exp(-iter/100)  # Decreasing regularization

    # Ridge regression for GSCA
    for(k in 1:K) {
      topic_idx <- P + k
      predictors <- which(path_matrix[, topic_idx] != 0)

      if(length(predictors) > 0) {
        X <- gamma[, predictors, drop = FALSE]
        y <- topic_counts[, k]

        # Ridge regression
        XtX <- t(X) %*% X + diag(lambda, ncol(X))
        Xty <- t(X) %*% y
        beta <- try(solve(XtX, Xty), silent = TRUE)

        # if(!inherits(beta, "try-error")) {
        path_matrix[predictors, topic_idx] <- beta
        gamma[, topic_idx] <- X %*% beta
        # }
      }
    }

    # Calculate perplexity
    current_perplexity <- tryCatch({
      calculate_perplexity(dtm, phi, theta)
    }, error = function(e) {
      return(Inf)
    })

    perplexity_history[iter] <- current_perplexity

    # Store best model if perplexity improves
    if(current_perplexity < best_perplexity) {
      best_perplexity <- current_perplexity
      best_phi <- phi
      best_theta <- theta
    }

    # Check convergence
    delta_phi <- mean(abs(phi - phi_old), na.rm = TRUE)
    delta_theta <- mean(abs(theta - theta_old), na.rm = TRUE)
    delta_gamma <- mean(abs(gamma - gamma_old), na.rm = TRUE)

    # Check perplexity stability
    perplexity_window <- perplexity_history[max(1, iter-10):iter]
    perplexity_stable <- iter > 10 &&
      !any(is.na(perplexity_window)) &&
      !any(is.infinite(perplexity_window)) &&
      sd(perplexity_window, na.rm = TRUE) < tol

    # Convergence criteria
    deltas_ok <- !is.na(delta_phi) && !is.na(delta_theta) && !is.na(delta_gamma) &&
      max(delta_phi, delta_theta, delta_gamma) < tol

    if(deltas_ok || perplexity_stable) {
      message("Converged at iteration ", iter)
      break
    }
  }

  # Use best model if it didn't converge
  if(iter == max_iter) {
    message("Using best model found (lowest perplexity)")
    phi <- best_phi
    theta <- best_theta
  }


  # Effect Estimation:
  coefficients <- path_matrix[1:P, (P + 1):(P + K)]
  if(is.vector(coefficients)){
    coefficients <- matrix(coefficients, nrow=1)
  }
  rownames(coefficients) <- colnames(gamma)[1:P]
  colnames(coefficients) <- paste0("Topic_", 1:K)

  # Use last topic as reference category
  # Use specified topic as reference category
  # Reorder columns to put reference topic last for log ratios
  topic_cols <- 1:K
  topic_cols <- c(topic_cols[-reference_topic], reference_topic)
  theta_ordered <- theta[, topic_cols, drop=FALSE]
  # Calculate log ratios against reference topic
  log_ratios <- log(theta_ordered[, -K, drop=FALSE] / theta_ordered[, K])

  # Initialize a list to store models
  glm_models <- list()

  # Initialize a data frame to store results
  non_ref_topics <- setdiff(1:K, reference_topic)
  summary_table <- data.frame(
    Covariate = rep(colnames(covariates), K-1),
    Topic = paste0("Topic_", non_ref_topics),
    Reference_Topic = paste0("Topic_", reference_topic),
    Log_Odds_Coefficient = NA,
    Std_Error = NA,
    z_value = NA,
    p_value = NA,
    Odds_Ratio = NA,
    CI_Lower = NA,
    CI_Upper = NA
  )

  # Counter for the rows in summary_table
  counter <- 1

  # Loop over each topic (except reference topic)
  for (k in 1:(K-1)) {
    # Response variable: log ratio of topic k vs reference topic
    y <- log_ratios[, k]

    # Create a data frame for modeling
    glm_data <- data.frame(y = y, scale(covariates))

    # Calculate weights for heteroscedasticity
    weights <- theta[, k] * theta[, K]  # optimal weights for multinomial

    # Fit weighted logistic regression using gaussian family for log-ratios
    glm_fit <- glm(y ~ ., data = glm_data, weights = weights, family = gaussian())

    # Store the model
    glm_models[[k]] <- glm_fit

    # Extract summary statistics
    summary_glm <- summary(glm_fit)
    coef_table <- summary_glm$coefficients

    # Loop over covariates
    for (p in 1:P) {
      summary_table$Log_Odds_Coefficient[counter] <- coef_table[p, "Estimate"]
      summary_table$Std_Error[counter] <- coef_table[p, "Std. Error"]
      summary_table$z_value[counter] <- coef_table[p, "t value"]  # t-value is z-value in this context
      summary_table$p_value[counter] <- coef_table[p, "Pr(>|t|)"]

      # Calculate odds ratio and confidence intervals
      log_coef <- coef_table[p, "Estimate"]
      se <- coef_table[p, "Std. Error"]

      summary_table$Odds_Ratio[counter] <- exp(log_coef)
      summary_table$CI_Lower[counter] <- exp(log_coef - 1.96 * se)
      summary_table$CI_Upper[counter] <- exp(log_coef + 1.96 * se)

      counter <- counter + 1
    }
  }
  summary_table$adjusted_p_value <- p.adjust(summary_table$p_value, method = "BH")

  # Add overall model fit statistics
  model_fit <- list()
  for (k in 1:(K-1)) {
    model_fit[[k]] <- list(
      AIC = AIC(glm_models[[k]]),
      BIC = BIC(glm_models[[k]]),
      Deviance = deviance(glm_models[[k]]),
      Null_Deviance = glm_models[[k]]$null.deviance,
      McFadden_R2 = 1 - (deviance(glm_models[[k]]) / glm_models[[k]]$null.deviance)
    )
  }


  # Prepare return object
  model_obj <- list(
    phi = phi,
    theta = theta,
    gamma = gamma,
    path_matrix = path_matrix,
    weight_matrix = weight_matrix,
    K = K,
    P=P,
    vocab = colnames(dtm),
    model_type = model,
    covariates = covariates,
    effects = summary_table,
    fit = model_fit,
    convergence = list(
      iterations = iter,
      final_delta = max(delta_phi, delta_theta, delta_gamma),
      perplexity_history = perplexity_history[1:iter],
      final_perplexity = perplexity_history[iter],
      best_perplexity = best_perplexity
    ),
    bootstrap_ci_computed = FALSE
  )

  class(model_obj) <- "gsca_topic_model"

  # ----- Parametric bootstrap (paper Section 4.2) -----
  # Propagates uncertainty from the full GSCA-TM estimation through to the
  # ALR-WLS regression.  For each bootstrap replication b:
  #   1. Draw W^(b)_i ~ Multinomial(L_i, q_hat_i), q_hat_i = sum_k theta_ik phi_k
  #   2. Refit GSCA-TM on W^(b) to obtain theta^(b), phi^(b)
  #   3. Align bootstrap topics to original via Hungarian algorithm (L1 cost)
  #   4. Compute ALR-WLS effects on aligned theta^(b)
  # Bootstrap SE = sd of the B bootstrap log-odds coefficients.
  if (isTRUE(bootstrap_ci)) {
    if (n_bootstrap < 2)
      stop("n_bootstrap must be at least 2")

    message("Running parametric bootstrap (B = ", n_bootstrap,
            "). This may take a while...")

    n_effects  <- P * (K - 1L)
    boot_coefs <- matrix(NA_real_, nrow = n_bootstrap, ncol = n_effects)

    # Fitted word probabilities for bootstrap DTM generation
    q_hat <- as.matrix(theta %*% phi)     # N x V
    L     <- rowSums(dtm)                  # document lengths

    for (b in seq_len(n_bootstrap)) {
      tryCatch({
        # Step 1: generate bootstrap DTM
        dtm_boot <- Matrix::Matrix(0L, nrow = N, ncol = V, sparse = TRUE)
        colnames(dtm_boot) <- colnames(dtm)
        for (i in seq_len(N)) {
          if (L[i] > 0L) {
            q_i <- pmax(q_hat[i, ], 0)
            s   <- sum(q_i)
            if (s > 0) {
              q_i <- q_i / s
              dtm_boot[i, ] <- as.integer(rmultinom(1L, L[i], q_i))
            }
          }
        }
        dtm_boot <- as(dtm_boot, "dgCMatrix")

        # Step 2: refit GSCA-TM on bootstrap DTM (no nested bootstrap)
        model_b <- suppressMessages(
          fit_gsca_topic_model(
            dtm             = dtm_boot,
            covariates      = covariates_df,
            K               = K,
            model           = model,
            max_iter        = max_iter,
            tol             = tol,
            alpha           = alpha,
            reference_topic = reference_topic,
            bootstrap_ci    = FALSE
          )
        )

        # Step 3: align bootstrap topics to original (L1 Hungarian)
        cost_mat <- matrix(0, K, K)
        for (j in seq_len(K))
          for (l in seq_len(K))
            cost_mat[j, l] <- sum(abs(phi[j, ] - model_b$phi[l, ]))

        if (requireNamespace("clue", quietly = TRUE)) {
          perm <- as.integer(clue::solve_LSAP(cost_mat))
        } else {
          # Greedy fallback
          perm      <- integer(K)
          available <- seq_len(K)
          for (j in seq_len(K)) {
            best      <- which.min(cost_mat[j, available])
            perm[j]   <- available[best]
            available  <- available[-best]
          }
        }
        theta_b_aligned <- model_b$theta[, perm, drop = FALSE]

        # Step 4: compute ALR-WLS effects on aligned theta
        # Reorder for reference topic
        topic_cols_b <- c(setdiff(1:K, reference_topic), reference_topic)
        theta_b_ord  <- theta_b_aligned[, topic_cols_b, drop = FALSE]
        log_ratios_b <- log(theta_b_ord[, -K, drop = FALSE] / theta_b_ord[, K])

        coefs_b <- numeric(n_effects)
        idx_b   <- 1L
        for (k_b in 1:(K - 1)) {
          y_b       <- log_ratios_b[, k_b]
          glm_data_b <- data.frame(y = y_b, scale(covariates))
          weights_b <- theta_b_aligned[, k_b] * theta_b_aligned[, K]
          fit_b     <- glm(y ~ ., data = glm_data_b, weights = weights_b,
                           family = gaussian())
          ct_b      <- summary(fit_b)$coefficients
          for (p_b in 1:P) {
            coefs_b[idx_b] <- ct_b[p_b, "Estimate"]
            idx_b <- idx_b + 1L
          }
        }
        boot_coefs[b, ] <- coefs_b

      }, error = function(e) NULL)   # silently skip failed replications
    }

    # Compute bootstrap SEs and CIs
    n_valid  <- colSums(!is.na(boot_coefs))
    boot_se  <- apply(boot_coefs, 2L, sd,   na.rm = TRUE)
    boot_mean <- apply(boot_coefs, 2L, mean, na.rm = TRUE)

    # Append bootstrap SEs and CIs to effects table (on OR scale)
    model_obj$effects$Boot_SE       <- boot_se
    model_obj$effects$Boot_CI_Lower <- exp(
      model_obj$effects$Log_Odds_Coefficient - 1.96 * boot_se
    )
    model_obj$effects$Boot_CI_Upper <- exp(
      model_obj$effects$Log_Odds_Coefficient + 1.96 * boot_se
    )

    model_obj$bootstrap <- list(
      boot_coefs = boot_coefs,
      boot_se    = boot_se,
      boot_mean  = boot_mean,
      n_valid    = n_valid
    )
    model_obj$bootstrap_ci_computed <- TRUE
  }

  return(model_obj)
}

#' Search for Optimal Number of Topics
#'
#' This function searches for the optimal number of topics by fitting the GSCA Topic Model
#' across a range of topic numbers and computing model diagnostics.
#'
#' @param dtm Document-term matrix (dgCMatrix)
#' @param covariates Data frame or matrix of covariates
#' @param K_seq Sequence of topic numbers to evaluate
#' @param n Number of top K values to return
#' @param ... Additional arguments to pass to fit_gsca_topic_model
#'
#' @return A list containing:
#' \describe{
#'   \item{top_K}{Top n K values based on Semantic Coherence}
#'   \item{metrics}{A list of computed metrics for each K}
#'   \item{models}{List of fitted models for the top n K values}
#' }
#'
#' @export
search_optimal_topics <- function(dtm, covariates, K_seq, n=5, ...) {
  residuals_list <- numeric(length(K_seq))
  perplexity_list <- numeric(length(K_seq))
  semantic_coherence_list <- numeric(length(K_seq))
  lower_bound_list <- numeric(length(K_seq))
  models_list <- list()

  for (i in seq_along(K_seq)) {
    K <- K_seq[i]
    cat("Fitting model with K =", K, "\n")
    model <- fit_gsca_topic_model(dtm, covariates, K, ...)
    models_list[[i]] <- model

    # Get final perplexity as an approximation of Held-Out Likelihood
    perplexity <- model$convergence$final_perplexity
    perplexity_list[i] <- perplexity

    # Compute Residuals (sum of squared differences between observed and expected counts)
    residuals <- compute_residuals(model, dtm)
    residuals_list[i] <- residuals

    # Compute Semantic Coherence
    semantic_coherence <- compute_semantic_coherence(model, dtm)
    semantic_coherence_list[i] <- semantic_coherence

    # Approximate Lower Bound using negative perplexity
    lower_bound_list[i] <- -perplexity
  }

  # Plot the metrics
  par(mfrow=c(2,2))
  plot(K_seq, residuals_list, type='b', xlab='Number of Topics (K)', ylab='Residuals', main='Residuals')
  plot(K_seq, -perplexity_list, type='b', xlab='Number of Topics (K)', ylab='Held-Out Likelihood (approx.)', main='Held-Out Likelihood')
  plot(K_seq, semantic_coherence_list, type='b', xlab='Number of Topics (K)', ylab='Semantic Coherence', main='Semantic Coherence')
  plot(K_seq, lower_bound_list, type='b', xlab='Number of Topics (K)', ylab='Lower Bound (approx.)', main='Lower Bound')

  # Return list of top n K values based on Semantic Coherence
  top_n_indices <- order(semantic_coherence_list, decreasing=TRUE)[1:n]
  top_n_K <- K_seq[top_n_indices]

  return(list(top_K=top_n_K, metrics=list(residuals=residuals_list,
                                          perplexity=perplexity_list,
                                          semantic_coherence=semantic_coherence_list,
                                          lower_bound=lower_bound_list),
              models=models_list[top_n_indices]))
}

# Helper function to compute residuals
compute_residuals <- function(model, dtm) {
  theta <- model$theta  # N x K
  phi <- model$phi      # K x V
  N <- nrow(dtm)
  V <- ncol(dtm)
  expected_counts <- theta %*% phi  # N x V
  total_words <- rowSums(dtm)       # N x 1
  # Scale expected_counts by total words in each document
  expected_counts <- expected_counts * total_words
  # Compute residuals
  residuals <- sum((dtm - expected_counts)^2)
  return(residuals)
}

# Helper function to compute semantic coherence
compute_semantic_coherence <- function(model, dtm, top_n=10) {
  phi <- model$phi      # K x V matrix
  K <- nrow(phi)
  vocab <- model$vocab
  N_d <- nrow(dtm)
  term_counts <- colSums(dtm > 0)  # Number of documents each term appears in
  semantic_coherence <- numeric(K)
  for (k in 1:K) {
    # Get top_n words for topic k
    top_words_idx <- order(phi[k,], decreasing=TRUE)[1:top_n]
    # Initialize coherence for this topic
    coherence <- 0
    # Loop over word pairs
    for (m in 2:top_n) {
      for (l in 1:(m-1)) {
        word_m <- top_words_idx[m]
        word_l <- top_words_idx[l]
        D_wl <- term_counts[word_l]
        # Compute D(wm, wl): number of documents containing both words
        D_wm_wl <- sum(dtm[, word_m] > 0 & dtm[, word_l] > 0)
        # Add to coherence
        coherence <- coherence + log((D_wm_wl + 1) / D_wl)
      }
    }
    semantic_coherence[k] <- coherence
  }
  # Return average semantic coherence across topics
  mean_coherence <- mean(semantic_coherence)
  return(mean_coherence)
}
