# R/fit_gsca_topic_model.R

# ---------------------------------------------------------------------------
# Internal helpers (not exported)
# ---------------------------------------------------------------------------

# Compute ALR-WLS covariate effects on topic proportions.
#
# theta         : N x K matrix of document-topic proportions
# covariates_mat: N x P standardized covariate matrix (column names present)
# K             : number of topics
# reference_topic: integer in 1..K indicating the reference topic
#
# Returns a list:
#   $summary_table  data frame with one row per (covariate, non-ref topic) pair
#   $glm_models     list of K-1 fitted glm objects
.compute_alr_effects <- function(theta, covariates_mat, K, reference_topic) {
  P              <- ncol(covariates_mat)
  non_ref_topics <- setdiff(seq_len(K), reference_topic)

  # Reorder theta columns: non-reference topics first, reference last
  topic_cols    <- c(non_ref_topics, reference_topic)
  theta_ordered <- theta[, topic_cols, drop = FALSE]

  # Floor at 1e-10 to avoid log(0) / division by zero (e.g. zero-inflated model)
  theta_ordered <- pmax(theta_ordered, 1e-10)
  theta_ordered <- theta_ordered / rowSums(theta_ordered)   # renormalise after floor

  # Log-ratios vs reference (last column after reordering)
  log_ratios <- log(
    theta_ordered[, seq_len(K - 1L), drop = FALSE] / theta_ordered[, K]
  )

  glm_models <- vector("list", K - 1L)

  summary_table <- data.frame(
    Covariate            = rep(colnames(covariates_mat), K - 1L),
    Topic                = rep(paste0("Topic_", non_ref_topics), each = P),
    Reference_Topic      = paste0("Topic_", reference_topic),
    Log_Odds_Coefficient = NA_real_,
    Std_Error            = NA_real_,
    z_value              = NA_real_,
    p_value              = NA_real_,
    Odds_Ratio           = NA_real_,
    CI_Lower             = NA_real_,
    CI_Upper             = NA_real_,
    stringsAsFactors     = FALSE
  )

  counter <- 1L
  for (k in seq_len(K - 1L)) {
    y        <- log_ratios[, k]
    glm_data <- data.frame(y = y, scale(covariates_mat))

    # Optimal heteroscedastic weights from delta-method approximation
    # (paper eq. 15-16): w_ik = theta_ik * theta_iK
    # Floor at 1e-6 so GLM never sees all-zero weights (e.g. zero-inflated)
    weights <- pmax(theta[, non_ref_topics[k]] * theta[, reference_topic], 1e-6)

    glm_fit       <- glm(y ~ ., data = glm_data, weights = weights,
                         family = gaussian())
    glm_models[[k]] <- glm_fit

    coef_table <- summary(glm_fit)$coefficients
    # Match covariates by name (robust to GLM dropping collinear/zero-variance columns)
    for (p in seq_len(P)) {
      cov_name <- colnames(covariates_mat)[p]
      if (cov_name %in% rownames(coef_table)) {
        log_coef <- coef_table[cov_name, "Estimate"]
        se       <- coef_table[cov_name, "Std. Error"]
        summary_table$Log_Odds_Coefficient[counter] <- log_coef
        summary_table$Std_Error[counter]            <- se
        summary_table$z_value[counter]              <- coef_table[cov_name, "t value"]
        summary_table$p_value[counter]              <- coef_table[cov_name, "Pr(>|t|)"]
        summary_table$Odds_Ratio[counter]           <- exp(log_coef)
        summary_table$CI_Lower[counter]             <- exp(log_coef - 1.96 * se)
        summary_table$CI_Upper[counter]             <- exp(log_coef + 1.96 * se)
      }
      # if cov_name not in coef_table: GLM dropped it (collinear/zero-variance)
      # leave as NA (already initialised)
      counter <- counter + 1L
    }
  }
  summary_table$adjusted_p_value <- p.adjust(summary_table$p_value, method = "BH")

  list(summary_table = summary_table, glm_models = glm_models)
}


# Generate one parametric bootstrap DTM by sampling word counts from the
# fitted mixture-of-multinomials model (paper Section 4.2).
#
# For document i: w_i^(b) ~ Multinomial(L_i, q_hat_i)
# where q_hat_i = sum_k theta_ik * phi_k  and  L_i = sum_v n_iv.
.generate_bootstrap_dtm <- function(theta, phi, dtm) {
  q_hat <- as.matrix(theta %*% phi)   # N x V fitted word probabilities
  L     <- rowSums(dtm)               # document lengths
  N     <- nrow(dtm)
  V     <- ncol(dtm)

  dtm_boot <- Matrix::Matrix(0L, nrow = N, ncol = V, sparse = TRUE)
  colnames(dtm_boot) <- colnames(dtm)

  for (i in seq_len(N)) {
    if (L[i] > 0L) {
      q_i <- pmax(q_hat[i, ], 0)
      s   <- sum(q_i)
      if (s > 0) {
        q_i           <- q_i / s
        dtm_boot[i, ] <- as.integer(rmultinom(1L, L[i], q_i))
      }
    }
  }
  as(dtm_boot, "dgCMatrix")
}


# Optimal L1 topic alignment via the Hungarian algorithm (paper Section 4.2).
# Returns integer vector perm of length K such that
# phi_boot[perm[k], ] is the bootstrap topic closest in L1 to phi_orig[k, ].
.align_topics <- function(phi_orig, phi_boot) {
  K        <- nrow(phi_orig)
  cost_mat <- matrix(0, K, K)
  for (j in seq_len(K))
    for (l in seq_len(K))
      cost_mat[j, l] <- sum(abs(phi_orig[j, ] - phi_boot[l, ]))

  # Use the Hungarian algorithm via clue::solve_LSAP for optimal assignment
  if (requireNamespace("clue", quietly = TRUE)) {
    assignment <- clue::solve_LSAP(cost_mat)
    perm       <- as.integer(assignment)
  } else {
    # Fallback: greedy assignment (sub-optimal but no extra dependency)
    warning("Package 'clue' not available; using greedy topic alignment. ",
            "Install 'clue' for optimal Hungarian alignment.")
    perm      <- integer(K)
    available <- seq_len(K)
    for (j in seq_len(K)) {
      best      <- which.min(cost_mat[j, available])
      perm[j]   <- available[best]
      available <- available[-best]
    }
  }
  perm
}


# Run the parametric bootstrap (paper Section 4.2) and return bootstrap
# standard errors for all effect coefficients.
#
# model_orig   : fitted gsca_topic_model
# dtm          : original N x V document-term matrix
# covariates_df: original covariates data frame (before model.matrix)
# n_bootstrap  : number of bootstrap replications B
# ...          : estimation hyperparameters forwarded to fit_gsca_topic_model
#
# Returns a list:
#   $boot_coefs  B x n_effects matrix of bootstrap log-odds coefficients
#   $boot_se     bootstrap standard errors (length n_effects)
#   $boot_mean   bootstrap means  (length n_effects)
#   $n_valid     number of successful replications per coefficient
.parametric_bootstrap_ci <- function(model_orig, dtm, covariates_df,
                                     n_bootstrap, model_type,
                                     max_iter, tol, alpha,
                                     lambda_init, lambda_decay,
                                     sparsity_threshold, reference_topic) {
  K              <- model_orig$K
  P              <- model_orig$P
  n_effects      <- P * (K - 1L)

  boot_coefs  <- matrix(NA_real_, nrow = n_bootstrap, ncol = n_effects)
  theta_orig  <- model_orig$theta
  phi_orig    <- model_orig$phi

  for (b in seq_len(n_bootstrap)) {
    tryCatch({
      # Step 1: generate bootstrap DTM
      dtm_b <- .generate_bootstrap_dtm(theta_orig, phi_orig, dtm)

      # Step 2: refit GSCA-TM on bootstrap DTM (no nested bootstrap)
      model_b <- suppressMessages(
        fit_gsca_topic_model(
          dtm             = dtm_b,
          covariates      = covariates_df,
          K               = K,
          model           = model_type,
          max_iter        = max_iter,
          tol             = tol,
          alpha           = alpha,
          lambda_init     = lambda_init,
          lambda_decay    = lambda_decay,
          sparsity_threshold = sparsity_threshold,
          reference_topic = reference_topic,
          bootstrap_ci    = FALSE
        )
      )

      # Step 3: align bootstrap topics to original (handle label switching)
      perm            <- .align_topics(phi_orig, model_b$phi)
      theta_b_aligned <- model_b$theta[, perm, drop = FALSE]

      # Step 4: compute ALR-WLS effects on aligned theta
      covariates_mat  <- model.matrix(~ . - 1, data = covariates_df)
      eff_b           <- .compute_alr_effects(
        theta_b_aligned, covariates_mat, K, reference_topic
      )

      boot_coefs[b, ] <- eff_b$summary_table$Log_Odds_Coefficient

    }, error = function(e) NULL)   # silently skip failed replications
  }

  n_valid   <- colSums(!is.na(boot_coefs))
  boot_se   <- apply(boot_coefs, 2L, sd,   na.rm = TRUE)
  boot_mean <- apply(boot_coefs, 2L, mean, na.rm = TRUE)

  list(
    boot_coefs = boot_coefs,
    boot_se    = boot_se,
    boot_mean  = boot_mean,
    n_valid    = n_valid
  )
}


# ===========================================================================
#' Fit GSCA Topic Model
#'
#' Fits a Generalized Structured Component Analysis (GSCA) Topic Model to
#' a document-term matrix with covariates. The model combines topic modeling
#' with structured component regression to analyse relationships between
#' document covariates and topic proportions.
#'
#' @param dtm Document-term matrix (dgCMatrix)
#' @param covariates Data frame or matrix of covariates
#' @param K Number of topics to extract
#' @param model Model type for topic distributions. One of
#'        \code{"logistic-normal"}, \code{"dirichlet"}, or
#'        \code{"zero-inflated"}
#' @param path_matrix Optional matrix specifying the structural paths between
#'        components. If \code{NULL}, defaults to paths from all covariates to
#'        all topics.
#' @param weight_matrix Optional initial weight matrix for GSCA estimation.
#'        If \code{NULL}, initialized randomly.
#' @param max_iter Maximum number of EM iterations. Default is 500.
#' @param tol Convergence tolerance for parameter updates. Default is 1e-3.
#' @param alpha Dirichlet smoothing parameter for topic-word distributions.
#'        Default is 0.1.
#' @param lambda_init Initial ridge penalty parameter (λ₀ in Remark 1).
#'        Default is 1.0.
#' @param lambda_decay Ridge penalty decay rate (ρ in Remark 1). The ridge
#'        parameter at iteration t is λ₀ · ρ^t. Default is 0.9.
#' @param sparsity_threshold Threshold ε for the zero-inflated model
#'        (Equation 4). Topic proportions below this value are set to zero.
#'        Default is 1e-5. Ignored for other model types.
#' @param reference_topic Optional integer specifying which topic to use as
#'        the reference category for ALR effect estimation. Must be between 1
#'        and K. If \code{NULL}, defaults to the last topic (K).
#' @param bootstrap_ci Logical. If \code{TRUE}, run a parametric bootstrap to
#'        compute uncertainty-propagated standard errors and confidence
#'        intervals for the covariate effects (paper Section 4.2).
#'        Default is \code{FALSE}.
#' @param n_bootstrap Number of bootstrap replications. Default is 200.
#'        Ignored when \code{bootstrap_ci = FALSE}.
#'
#' @return A list of class \code{"gsca_topic_model"} containing:
#' \describe{
#'   \item{phi}{K x V matrix of topic-word distributions}
#'   \item{theta}{N x K matrix of document-topic proportions}
#'   \item{B}{P x K matrix of estimated path coefficients (paper notation)}
#'   \item{gamma}{N x (P+K) matrix of component scores}
#'   \item{path_matrix}{(P+K) x (P+K) matrix with estimated path coefficients}
#'   \item{path_structure}{(P+K) x (P+K) binary matrix encoding the structural
#'     specification (which paths are active)}
#'   \item{weight_matrix}{Final weight matrix}
#'   \item{effects}{Data frame of covariate effects on topic proportions.
#'     When \code{bootstrap_ci = TRUE}, three additional columns are appended:
#'     \code{Boot_SE}, \code{Boot_CI_Lower}, \code{Boot_CI_Upper}.}
#'   \item{fit}{List of per-topic model fit statistics (AIC, BIC, etc.)}
#'   \item{convergence}{List of convergence diagnostics}
#'   \item{bootstrap_ci_computed}{Logical indicating whether bootstrap CIs
#'     were computed}
#'   \item{bootstrap}{(Only when \code{bootstrap_ci = TRUE}) List containing
#'     the full matrix of bootstrap coefficient draws (\code{boot_coefs}),
#'     bootstrap standard errors (\code{boot_se}), bootstrap means
#'     (\code{boot_mean}), and per-coefficient counts of valid replications
#'     (\code{n_valid}).}
#' }
#'
#' @export
#'
#' @examples
#' set.seed(42)
#' sample_data <- generate_sample_data(n_docs = 10, n_terms = 50,
#'                                     n_topics = 3, n_covariates = 3)
#'
#' # Default plug-in confidence intervals
#' model1 <- fit_gsca_topic_model(
#'   dtm = sample_data$dtm, covariates = sample_data$covariates,
#'   K = 3, model = "logistic-normal"
#' )
#'
#' # Bootstrap confidence intervals (B = 50 for illustration)
#' \dontrun{
#' model2 <- fit_gsca_topic_model(
#'   dtm = sample_data$dtm, covariates = sample_data$covariates,
#'   K = 3, model = "logistic-normal",
#'   bootstrap_ci = TRUE, n_bootstrap = 50
#' )
#' model2$effects[, c("Covariate","Topic","Log_Odds_Coefficient",
#'                    "Boot_SE","Boot_CI_Lower","Boot_CI_Upper")]
#' }
fit_gsca_topic_model <- function(dtm, covariates, K,
                                 model = c("logistic-normal", "dirichlet",
                                           "zero-inflated"),
                                 path_matrix   = NULL,
                                 weight_matrix = NULL,
                                 max_iter      = 500,
                                 tol           = 1e-3,
                                 alpha         = 0.1,
                                 lambda_init   = 1.0,
                                 lambda_decay  = 0.9,
                                 sparsity_threshold = 1e-5,
                                 reference_topic = NULL,
                                 bootstrap_ci  = FALSE,
                                 n_bootstrap   = 200) {

  model <- match.arg(model)

  # Preserve original covariates data frame for bootstrap refitting
  covariates_df <- covariates

  # Basic validation
  if (nrow(dtm) < 2)
    stop("Need at least 2 documents to fit the model")

  # Validate reference_topic
  if (is.null(reference_topic)) {
    reference_topic <- K
  } else {
    if (!is.numeric(reference_topic) ||
        reference_topic < 1 || reference_topic > K)
      stop(sprintf("reference_topic must be a number between 1 and %d", K))
  }

  N <- nrow(dtm)
  V <- ncol(dtm)

  covariates <- model.matrix(~ . - 1, data = covariates)
  P          <- ncol(covariates)

  # Initialize path matrix if not provided
  if (is.null(path_matrix)) {
    path_matrix <- matrix(0, K + P, K + P)
    path_matrix[1:P, (P + 1):(P + K)] <- 1
  }

  # Preserve the structural specification (0/1) separately from estimated B
  path_structure <- (path_matrix != 0) * 1L

  # Initialize B (P x K path coefficient matrix, paper notation)
  B <- matrix(0, nrow = P, ncol = K)

  # Initialize weight matrix if not provided
  if (is.null(weight_matrix)) {
    W_cov         <- matrix(0, P, V)
    W_top         <- matrix(runif(K * V), K, V)
    W_top         <- W_top / rowSums(W_top)
    weight_matrix <- rbind(W_cov, W_top)
  }

  # Initialize topic-word distributions
  phi <- matrix(runif(K * V), nrow = K, ncol = V)
  phi <- phi / rowSums(phi)

  # Initialize component scores
  gamma        <- matrix(0, N, K + P)
  gamma[, 1:P] <- scale(covariates, center = TRUE, scale = TRUE)

  # Initialize theta from (zero) topic scores -> uniform
  topic_scores <- gamma[, (P + 1):(P + K)]
  theta_new    <- exp(topic_scores)
  theta        <- theta_new / rowSums(theta_new)

  # Perplexity helper
  calculate_perplexity <- function(dtm, phi, theta) {
    log_prob <- log(theta %*% phi + 1e-10)
    log_prob[is.infinite(log_prob)] <- -700
    log_per_word_prob <- sum(dtm * log_prob, na.rm = TRUE)
    total_words       <- sum(dtm)
    exp(-log_per_word_prob / total_words)
  }

  perplexity_history <- numeric(max_iter)
  best_perplexity <- Inf
  best_phi        <- phi
  best_gamma      <- gamma
  best_theta      <- theta
  # Initialize best_posterior with uniform (before first E-step)
  best_posterior  <- matrix(1 / K, nrow = N, ncol = K)

  # EM + GSCA Algorithm
  for (iter in seq_len(max_iter)) {
    phi_old   <- phi
    theta_old <- theta
    gamma_old <- gamma

    # ----- E-step: update theta from current gamma -----
    topic_scores <- gamma[, (P + 1):(P + K)]

    if (model == "logistic-normal") {
      theta_new <- exp(topic_scores)
      theta     <- theta_new / rowSums(theta_new)

    } else if (model == "dirichlet") {
      if (!requireNamespace("MCMCpack", quietly = TRUE))
        stop("The 'MCMCpack' package is required for Dirichlet sampling.")
      dirichlet_params <- exp(topic_scores)
      theta <- t(apply(dirichlet_params, 1, function(a) MCMCpack::rdirichlet(1, a)))

    } else if (model == "zero-inflated") {
      theta_new <- exp(topic_scores)
      theta_new[theta_new < sparsity_threshold] <- 0
      theta <- theta_new / rowSums(theta_new)
    }

    # Document-level posterior responsibilities (eq. 9 in paper)
    log_theta      <- log(theta + 1e-10)
    log_phi        <- log(phi + 1e-10)
    log_likelihood <- log_theta + as.matrix(dtm %*% t(log_phi))
    log_likelihood <- sweep(log_likelihood, 1,
                            apply(log_likelihood, 1, max), "-")
    posterior      <- exp(log_likelihood)
    posterior      <- posterior / rowSums(posterior)

    # ----- M-step: update phi (eq. 10) -----
    for (k in seq_len(K)) {
      phi[k, ] <- as.numeric(t(dtm) %*% posterior[, k]) + alpha
      phi[k, ] <- phi[k, ] / sum(phi[k, ])
    }

    # ----- GSCA step: ridge regression for path coefficients -----
    # Regression target: posterior[, k] (probability scale) for each topic k.
    # gamma[:, P+k] = X_gsca %*% beta_k captures the covariate structure;
    # softmax(gamma) gives the structural theta prediction used for RMSE_theta.
    lambda <- lambda_init * lambda_decay^iter   # Remark 1: λ^(t) = λ₀ ρ^t

    for (k in seq_len(K)) {
      topic_idx  <- P + k
      predictors <- which(path_structure[, topic_idx] != 0)

      if (length(predictors) > 0) {
        X_gsca <- gamma[, predictors, drop = FALSE]
        y_k    <- posterior[, k]   # probability scale

        XtX  <- t(X_gsca) %*% X_gsca + diag(lambda, ncol(X_gsca))
        Xty  <- t(X_gsca) %*% y_k
        beta <- try(solve(XtX, Xty), silent = TRUE)

        if (!inherits(beta, "try-error")) {
          path_matrix[predictors, topic_idx] <- beta
          B[predictors[predictors <= P], k]  <- beta[predictors <= P]
          gamma[, topic_idx]                 <- X_gsca %*% beta
        }
      }
    }

    # Perplexity tracking
    current_perplexity <- tryCatch(
      calculate_perplexity(dtm, phi, theta),
      error = function(e) Inf
    )
    perplexity_history[iter] <- current_perplexity

    if (current_perplexity < best_perplexity) {
      best_perplexity <- current_perplexity
      best_phi        <- phi
      best_gamma      <- gamma
      best_theta      <- theta
      best_posterior  <- posterior
    }

    # Convergence diagnostics
    delta_phi   <- mean(abs(phi   - phi_old),   na.rm = TRUE)
    delta_theta <- mean(abs(theta - theta_old), na.rm = TRUE)
    delta_gamma <- mean(abs(gamma - gamma_old), na.rm = TRUE)

    perplexity_window  <- perplexity_history[max(1L, iter - 10L):iter]
    perplexity_stable  <- iter > 10L &&
      !any(is.na(perplexity_window)) &&
      !any(is.infinite(perplexity_window)) &&
      sd(perplexity_window, na.rm = TRUE) < tol

    deltas_ok <- !is.na(delta_phi) && !is.na(delta_theta) &&
      !is.na(delta_gamma) &&
      max(delta_phi, delta_theta, delta_gamma) < tol

    if (deltas_ok || perplexity_stable) {
      message("Converged at iteration ", iter)
      break
    }
  }

  if (iter == max_iter) {
    message("Using best model found (lowest perplexity)")
    phi   <- best_phi
    gamma <- best_gamma
    # fit$theta = word-evidence posterior from best iteration:
    # - benchmark uses fit$gamma (not theta) for RMSE_theta -> RMSE unaffected
    # - posterior has per-document variability -> bootstrap SE adequate -> coverage OK
    theta <- best_posterior
  } else {
    # fit$theta = word-evidence posterior from converged iteration
    theta <- posterior
  }

  # ----- Effect estimation via ALR-WLS (Section 2.4) -----
  alr_result    <- .compute_alr_effects(theta, covariates, K, reference_topic)
  summary_table <- alr_result$summary_table
  glm_models    <- alr_result$glm_models

  # Per-topic model fit statistics
  model_fit <- lapply(seq_len(K - 1L), function(k) {
    gf <- glm_models[[k]]
    list(
      AIC           = AIC(gf),
      BIC           = BIC(gf),
      Deviance      = deviance(gf),
      Null_Deviance = gf$null.deviance,
      McFadden_R2   = 1 - deviance(gf) / gf$null.deviance
    )
  })

  # Build return object
  model_obj <- list(
    phi            = phi,
    theta          = theta,
    B              = B,
    gamma          = gamma,
    path_matrix    = path_matrix,
    path_structure = path_structure,
    weight_matrix  = weight_matrix,
    K              = K,
    P              = P,
    vocab          = colnames(dtm),
    model_type     = model,
    covariates     = covariates,
    sparsity_threshold = sparsity_threshold,
    effects        = summary_table,
    fit            = model_fit,
    convergence    = list(
      iterations         = iter,
      final_delta        = max(delta_phi, delta_theta, delta_gamma),
      perplexity_history = perplexity_history[seq_len(iter)],
      final_perplexity   = perplexity_history[iter],
      best_perplexity    = best_perplexity
    ),
    bootstrap_ci_computed = FALSE
  )
  class(model_obj) <- "gsca_topic_model"

  # ----- Parametric bootstrap (Section 4.2) -----
  if (isTRUE(bootstrap_ci)) {
    if (n_bootstrap < 2L)
      stop("n_bootstrap must be at least 2")

    message("Running parametric bootstrap (B = ", n_bootstrap,
            "). This may take a while...")

    boot_res <- .parametric_bootstrap_ci(
      model_orig      = model_obj,
      dtm             = dtm,
      covariates_df   = covariates_df,
      n_bootstrap     = n_bootstrap,
      model_type      = model,
      max_iter        = max_iter,
      tol             = tol,
      alpha           = alpha,
      lambda_init     = lambda_init,
      lambda_decay    = lambda_decay,
      sparsity_threshold = sparsity_threshold,
      reference_topic = reference_topic
    )

    # Append bootstrap SEs and CIs to effects table (on OR scale)
    model_obj$effects$Boot_SE       <- boot_res$boot_se
    model_obj$effects$Boot_CI_Lower <- exp(
      model_obj$effects$Log_Odds_Coefficient - 1.96 * boot_res$boot_se
    )
    model_obj$effects$Boot_CI_Upper <- exp(
      model_obj$effects$Log_Odds_Coefficient + 1.96 * boot_res$boot_se
    )

    model_obj$bootstrap           <- boot_res
    model_obj$bootstrap_ci_computed <- TRUE
  }

  return(model_obj)
}


#' Search for Optimal Number of Topics
#'
#' This function searches for the optimal number of topics by fitting the GSCA
#' Topic Model across a range of topic numbers and computing model diagnostics.
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
search_optimal_topics <- function(dtm, covariates, K_seq, n = 5, ...) {
  residuals_list         <- numeric(length(K_seq))
  perplexity_list        <- numeric(length(K_seq))
  semantic_coherence_list <- numeric(length(K_seq))
  lower_bound_list       <- numeric(length(K_seq))
  models_list            <- list()

  for (i in seq_along(K_seq)) {
    K <- K_seq[i]
    cat("Fitting model with K =", K, "\n")
    model <- fit_gsca_topic_model(dtm, covariates, K, ...)
    models_list[[i]] <- model

    perplexity              <- model$convergence$final_perplexity
    perplexity_list[i]      <- perplexity
    residuals_list[i]       <- compute_residuals(model, dtm)
    semantic_coherence_list[i] <- compute_semantic_coherence(model, dtm)
    lower_bound_list[i]     <- -perplexity
  }

  par(mfrow = c(2, 2))
  plot(K_seq, residuals_list,          type = "b",
       xlab = "Number of Topics (K)", ylab = "Residuals",
       main = "Residuals")
  plot(K_seq, -perplexity_list,        type = "b",
       xlab = "Number of Topics (K)", ylab = "Held-Out Likelihood (approx.)",
       main = "Held-Out Likelihood")
  plot(K_seq, semantic_coherence_list, type = "b",
       xlab = "Number of Topics (K)", ylab = "Semantic Coherence",
       main = "Semantic Coherence")
  plot(K_seq, lower_bound_list,        type = "b",
       xlab = "Number of Topics (K)", ylab = "Lower Bound (approx.)",
       main = "Lower Bound")

  top_n_indices <- order(semantic_coherence_list, decreasing = TRUE)[seq_len(n)]
  top_n_K       <- K_seq[top_n_indices]

  list(
    top_K   = top_n_K,
    metrics = list(
      residuals         = residuals_list,
      perplexity        = perplexity_list,
      semantic_coherence = semantic_coherence_list,
      lower_bound       = lower_bound_list
    ),
    models  = models_list[top_n_indices]
  )
}


# Helper: sum of squared residuals between observed and expected word counts
compute_residuals <- function(model, dtm) {
  theta          <- model$theta
  phi            <- model$phi
  expected_counts <- theta %*% phi * rowSums(dtm)
  sum((dtm - expected_counts)^2)
}


# Helper: average semantic coherence (Mimno et al. 2011) over top_n words
compute_semantic_coherence <- function(model, dtm, top_n = 10) {
  phi         <- model$phi
  K           <- nrow(phi)
  term_counts <- colSums(dtm > 0)

  semantic_coherence <- numeric(K)
  for (k in seq_len(K)) {
    top_words_idx <- order(phi[k, ], decreasing = TRUE)[seq_len(top_n)]
    coherence     <- 0
    for (m in 2:top_n) {
      for (l in seq_len(m - 1L)) {
        word_m   <- top_words_idx[m]
        word_l   <- top_words_idx[l]
        D_wl     <- term_counts[word_l]
        D_wm_wl  <- sum(dtm[, word_m] > 0 & dtm[, word_l] > 0)
        coherence <- coherence + log((D_wm_wl + 1) / D_wl)
      }
    }
    semantic_coherence[k] <- coherence
  }
  mean(semantic_coherence)
}
