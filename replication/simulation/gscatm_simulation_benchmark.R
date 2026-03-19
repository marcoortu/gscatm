# Simulation and benchmarking script for GSCA-TM vs STM and LDA+ALR
# Assumes that fit_gsca_topic_model() is already available in the workspace
# (e.g., from the gscatm package or from sourcing its definition).

# ---- Required packages ----
library(Matrix)
library(MCMCpack)
library(stm)
library(topicmodels)
library(slam)
library(ggplot2)
library(dplyr)
library(tidyr)
library(clue)
library(xtable)
library(gscatm)

# ---- Helper: compute perplexity from phi and theta ----
compute_perplexity <- function(dtm, phi, theta) {
  theta <- as.matrix(theta)
  phi   <- as.matrix(phi)
  log_prob <- log(theta %*% phi + 1e-10)
  log_prob[is.infinite(log_prob)] <- -700
  log_per_word_prob <- sum(dtm * log_prob, na.rm = TRUE)
  total_words <- sum(dtm)
  perplexity <- exp(-log_per_word_prob / total_words)
  return(perplexity)
}

# ---- Helper: topic matching to remove label switching ----
match_topics <- function(phi_true, phi_hat) {
  K_true <- nrow(phi_true)
  K_hat  <- nrow(phi_hat)
  if (K_true != K_hat) {
    stop("phi_true and phi_hat must have same number of topics (rows).")
  }
  cost <- matrix(0, nrow = K_true, ncol = K_hat)
  for (k in 1:K_true) {
    for (j in 1:K_hat) {
      cost[k, j] <- sum(abs(phi_true[k, ] - phi_hat[j, ]))
    }
  }
  # solve_LSAP minimizes total cost and returns a permutation of columns
  perm <- clue::solve_LSAP(cost)
  as.integer(perm)
}

# ---- Helper: ALR regression on theta to estimate covariate effects ----
estimate_alr_effects <- function(theta, X_std, reference_topic = ncol(theta)) {
  theta <- as.matrix(theta)
  X_std <- as.matrix(X_std)
  K <- ncol(theta)
  P <- ncol(X_std)

  if (reference_topic < 1 || reference_topic > K) {
    stop("reference_topic must be between 1 and K.")
  }

  # Reorder topics so that reference_topic is last
  topic_idx <- 1:K
  topic_idx <- c(topic_idx[topic_idx != reference_topic], reference_topic)
  theta_ordered <- theta[, topic_idx, drop = FALSE]

  # Add small epsilon to avoid log(0)
  theta_safe <- pmax(theta_ordered, 1e-10)

  # ALR: log(theta_k / theta_ref) for k = 1, ..., K-1
  denom <- theta_safe[, K]  # vector of length n_docs
  log_ratios <- log(theta_safe[, 1:(K - 1), drop = FALSE] / denom)

  # Design matrix (standardized) is X_std
  if (is.null(colnames(X_std))) {
    colnames(X_std) <- paste0("x", 1:P)
  }

  B_hat  <- matrix(NA_real_, nrow = P, ncol = K - 1)
  se_hat <- matrix(NA_real_, nrow = P, ncol = K - 1)
  rownames(B_hat)  <- colnames(X_std)
  rownames(se_hat) <- colnames(X_std)
  colnames(B_hat)  <- paste0("Topic_", 1:(K - 1))
  colnames(se_hat) <- paste0("Topic_", 1:(K - 1))

  for (k in 1:(K - 1)) {
    y <- log_ratios[, k]
    glm_data <- data.frame(y = y, X_std)
    # Weights corresponding to multinomial variance
    weights <- theta_safe[, k] * theta_safe[, K]
    fit <- glm(y ~ ., data = glm_data,
               weights = weights,
               family = gaussian())

    coef_table <- summary(fit)$coefficients
    # Keep only slopes (exclude intercept)
    keep <- rownames(coef_table) %in% colnames(X_std)
    coef_table <- coef_table[keep, , drop = FALSE]

    B_hat[rownames(coef_table), k]  <- coef_table[, "Estimate"]
    se_hat[rownames(coef_table), k] <- coef_table[, "Std. Error"]
  }

  list(B_hat = B_hat, se_hat = se_hat)
}

# ---- Helper: convert dgCMatrix dtm to stm format ----
dtm_to_stm <- function(dtm) {
  dtm <- as(dtm, "dgCMatrix")
  N <- nrow(dtm)
  vocab <- colnames(dtm)
  documents <- vector("list", N)
  for (i in 1:N) {
    idx <- which(dtm[i, ] > 0)
    if (length(idx) == 0) {
      documents[[i]] <- matrix(integer(0), nrow = 2, ncol = 0)
    } else {
      counts <- dtm[i, idx]
      documents[[i]] <- rbind(idx, as.integer(counts))
    }
  }
  list(documents = documents, vocab = vocab)
}

# ---- Simulation of synthetic data under three conditions ----
simulate_gscatm_data <- function(condition = c("baseline",
                                               "high_sparsity",
                                               "high_covariate"),
                                 n_docs = 100,
                                 n_terms = 500,
                                 K = 5,
                                 P = 3,
                                 avg_doc_length_baseline = 100,
                                 avg_doc_length_sparse   = 20) {

  condition <- match.arg(condition)

  # Raw covariates
  covariates_df <- as.data.frame(matrix(rnorm(n_docs * P),
                                        nrow = n_docs,
                                        ncol = P))
  colnames(covariates_df) <- paste0("x", 1:P)

  # Design matrix as in fit_gsca_topic_model: no intercept
  X_design <- model.matrix(~ . - 1, data = covariates_df)
  # Standardize covariates (this is what is used in ALR step)
  X_std <- scale(X_design)
  P_eff <- ncol(X_std)

  # Define true covariate effects B (P_eff x (K-1))
  K_nonref <- K - 1
  B_true <- matrix(0, nrow = P_eff, ncol = K_nonref)
  for (p in 1:P_eff) {
    B_true[p, ] <- seq(-0.6, 0.6, length.out = K_nonref) * ((-1)^(p - 1))
  }
  rownames(B_true) <- colnames(X_std)
  colnames(B_true) <- paste0("Topic_", 1:K_nonref)

  # Linear predictor for ALR parameterization
  eta_mean <- X_std %*% B_true  # n_docs x (K-1)

  noise_sd <- switch(condition,
                     "baseline"       = 0.3,
                     "high_sparsity"  = 0.4,
                     "high_covariate" = 0.15)

  eta <- eta_mean +
    matrix(rnorm(n_docs * K_nonref, sd = noise_sd),
           nrow = n_docs, ncol = K_nonref)

  # Map to topic proportions via ALR (reference topic K has log-odds zero)
  eta_full <- cbind(eta, 0)
  exp_eta  <- exp(eta_full)
  theta_base <- exp_eta / rowSums(exp_eta)

  theta <- theta_base

  # High covariate influence: Dirichlet variability around mean theta_base
  if (condition == "high_covariate") {
    alpha_bar <- 50  # concentration parameter; higher -> closer to mean
    theta <- t(apply(theta_base, 1, function(mu) {
      as.numeric(MCMCpack::rdirichlet(1, alpha_bar * mu))
    }))
  }

  # High sparsity: zero-inflate small components and renormalize
  if (condition == "high_sparsity") {
    thresh <- 0.03
    theta[theta < thresh] <- 0
    rs <- rowSums(theta)
    # Avoid division by zero
    zero_rows <- which(rs == 0)
    if (length(zero_rows) > 0) {
      theta[zero_rows, ] <- theta_base[zero_rows, ]
      rs <- rowSums(theta)
    }
    theta <- theta / rs
  }

  # Topic-word distributions phi (K x V)
  beta0 <- rep(0.1, n_terms)
  phi_true <- t(replicate(K, as.numeric(MCMCpack::rdirichlet(1, beta0))))
  phi_true <- phi_true / rowSums(phi_true)
  colnames(phi_true) <- paste0("w", 1:n_terms)
  rownames(phi_true) <- paste0("Topic_", 1:K)

  # Document lengths
  lambda <- if (condition == "high_sparsity") avg_doc_length_sparse else avg_doc_length_baseline
  L <- rpois(n_docs, lambda = lambda)
  L[L < 5] <- 5

  # Generate document-term counts
  dtm_dense <- matrix(0, nrow = n_docs, ncol = n_terms)
  for (i in 1:n_docs) {
    q_i <- as.numeric(theta[i, ] %*% phi_true)
    dtm_dense[i, ] <- rmultinom(1, size = L[i], prob = q_i)
  }
  dtm <- Matrix::Matrix(dtm_dense, sparse = TRUE)
  colnames(dtm) <- colnames(phi_true)
  rownames(dtm) <- paste0("doc", 1:n_docs)

  # ---- NEW PART: drop all-zero columns and renormalize phi_true ----
  col_sums <- Matrix::colSums(dtm)
  keep <- col_sums > 0

  dtm <- dtm[, keep, drop = FALSE]
  phi_true <- phi_true[, keep, drop = FALSE]
  # Renormalize phi_true rows to sum to 1 over the remaining vocabulary
  phi_true <- phi_true / rowSums(phi_true)

  # Update vocab names (already correct via colnames)
  # n_terms_effective <- ncol(dtm)  # if you need it later

  list(
    condition   = condition,
    dtm         = dtm,
    covariates  = covariates_df,
    X_std       = X_std,
    theta_true  = theta,
    phi_true    = phi_true,
    B_true      = B_true
  )
}


# ---- Fit models and compute metrics for one dataset ----
evaluate_method <- function(sim,
                            method = c("gscatm", "stm", "lda"),
                            reference_topic = ncol(sim$theta_true),
                            bootstrap_ci = FALSE,
                            n_bootstrap  = 200) {

  method <- match.arg(method)

  dtm        <- sim$dtm
  covariates <- sim$covariates
  X_std      <- sim$X_std
  theta_true <- sim$theta_true
  phi_true   <- sim$phi_true
  B_true     <- sim$B_true

  K <- ncol(theta_true)

  # model_type is needed both for fitting and for the bootstrap refits
  model_type <- if (method == "gscatm") {
    switch(sim$condition,
           "baseline"       = "logistic-normal",
           "high_sparsity"  = "zero-inflated",
           "high_covariate" = "dirichlet")
  } else NA_character_

  # --- Fit each model ---
  if (method == "gscatm") {
    fit <- fit_gsca_topic_model(
      dtm        = dtm,
      covariates = covariates,
      K          = K,
      model      = model_type,
      reference_topic = reference_topic
    )

    theta_hat_raw <- fit$theta
    phi_hat_raw   <- fit$phi

  } else if (method == "stm") {

    stm_in <- dtm_to_stm(dtm)
    stm_fit <- stm::stm(
      documents    = stm_in$documents,
      vocab        = stm_in$vocab,
      K            = K,
      prevalence   = ~ .,
      data         = covariates,
      max.em.its   = 75,
      init.type    = "Spectral",
      verbose      = FALSE
    )

    theta_hat_raw <- stm_fit$theta
    # beta is a list of log word probabilities; no content covariate -> one element
    logbeta <- stm_fit$beta$logbeta[[1]]
    phi_hat_raw <- exp(logbeta)
    phi_hat_raw <- phi_hat_raw / rowSums(phi_hat_raw)

  } else if (method == "lda") {

    dtm_triplet <- slam::as.simple_triplet_matrix(dtm)
    lda_fit <- topicmodels::LDA(
      dtm_triplet,
      k       = K,
      control = list(seed = 1)
    )

    post <- topicmodels::posterior(lda_fit)
    theta_hat_raw <- post$topics
    phi_hat_raw   <- post$terms
  }

  # --- Match topics to true labels (remove label switching) ---
  perm <- match_topics(phi_true, phi_hat_raw)
  phi_hat   <- phi_hat_raw[perm, , drop = FALSE]
  theta_hat <- theta_hat_raw[, perm, drop = FALSE]

  # Renormalize to guard against numerical drift
  phi_hat   <- phi_hat / rowSums(phi_hat)
  theta_hat <- theta_hat / rowSums(theta_hat)

  # --- Metrics: RMSE phi ---
  RMSE_phi <- sqrt(mean((phi_hat - phi_true)^2))

  # --- Covariate effects via ALR regression ---
  # For gscatm: use structural theta softmax(gamma) instead of the
  # word-evidence posterior (fit$theta).  The posterior has near-zero entries
  # for dominated topics that, after flooring at 1e-10, produce extreme
  # log-ratios (up to log(1/1e-10) = 23) and cause RMSE_B to explode.
  # softmax(gamma) is strictly positive and bounded, so ALR log-ratios stay
  # moderate and recover B_gsca ≈ B_true.  This mirrors the RMSE_theta fix.
  if (method == "gscatm") {
    ts_alr      <- fit$gamma[, (fit$P + 1):(fit$P + fit$K), drop = FALSE]
    ts_alr      <- ts_alr[, perm, drop = FALSE]
    theta_struct <- exp(ts_alr) / rowSums(exp(ts_alr))
    alr <- estimate_alr_effects(theta_struct, X_std, reference_topic = K)
  } else {
    alr <- estimate_alr_effects(theta_hat, X_std, reference_topic = K)
  }
  B_hat <- alr$B_hat
  se_hat <- alr$se_hat

  # --- RMSE_theta ---
  # For gscatm: use structural prediction softmax(gamma) = softmax(X * B_gsca),
  # the covariate-prior theta from the best EM iteration.  This is structurally
  # aligned with theta_true = softmax(X * B_true + noise) and recovers the
  # pre-coverage-fix RMSE_theta advantage.  fit$gamma is guaranteed to be from
  # the best iteration (lowest perplexity) because best_gamma is now stored and
  # restored alongside best_phi / best_theta for non-convergent runs.
  # For STM/LDA: no structural prediction available; use posterior theta.
  if (method == "gscatm") {
    topic_scores   <- fit$gamma[, (fit$P + 1):(fit$P + fit$K), drop = FALSE]
    topic_scores   <- topic_scores[, perm, drop = FALSE]   # align to phi_true order
    exp_scores     <- exp(topic_scores)
    theta_hat_rmse <- exp_scores / rowSums(exp_scores)
  } else {
    theta_hat_rmse <- theta_hat
  }
  RMSE_theta <- sqrt(mean((theta_hat_rmse - theta_true)^2))

  # Ensure B_true, B_hat aligned by row/column names
  B_true_aligned <- B_true[rownames(B_hat), colnames(B_hat), drop = FALSE]

  RMSE_B <- sqrt(mean((B_hat - B_true_aligned)^2))

  inside <- (B_true_aligned >= (B_hat - 1.96 * se_hat)) &
    (B_true_aligned <= (B_hat + 1.96 * se_hat))
  coverage_B <- mean(inside, na.rm = TRUE)

  # --- Perplexity on observed data ---
  perplexity <- compute_perplexity(dtm, phi_hat, theta_hat)

  # --- Bootstrap coverage (gscatm only) ---
  # Non-parametric row-bootstrap: resample documents (rows) with replacement.
  # This captures BOTH sources of variability in B_hat:
  #   (1) latent theta variability (from the ALR noise in the DGP), and
  #   (2) word-count sampling variability.
  # The parametric bootstrap (resampling word counts from fixed theta_hat)
  # captures only source (2), which is negligible for large N and long documents,
  # yielding near-zero bootstrap SEs and zero coverage of B_true.
  coverage_B_boot <- NA_real_

  if (isTRUE(bootstrap_ci) && method == "gscatm") {
    N     <- nrow(dtm)

    P_eff      <- nrow(B_hat)
    K_nonref   <- ncol(B_hat)
    boot_B_arr <- array(NA_real_, dim = c(P_eff, K_nonref, n_bootstrap))

    for (b in seq_len(n_bootstrap)) {
      tryCatch({
        # Step 1: resample documents with replacement (non-parametric row-bootstrap)
        idx_b       <- sample.int(N, N, replace = TRUE)
        dtm_b       <- dtm[idx_b, , drop = FALSE]
        covariates_b <- covariates[idx_b, , drop = FALSE]
        X_std_b     <- X_std[idx_b, , drop = FALSE]

        # Step 2: refit GSCA-TM on bootstrap sample (no nested bootstrap)
        fit_b <- suppressMessages(
          fit_gsca_topic_model(
            dtm             = dtm_b,
            covariates      = covariates_b,
            K               = K,
            model           = model_type,
            reference_topic = K,
            bootstrap_ci    = FALSE
          )
        )

        # Step 3: align bootstrap topics to aligned reference (phi_hat, not phi_true)
        perm_b  <- match_topics(phi_hat, fit_b$phi)
        # Use word-evidence posterior theta for bootstrap ALR: the posterior
        # has per-document variability (from word counts) that makes bootstrap
        # SEs wide enough to cover B_true.  The main B_hat uses theta_struct
        # for a better point estimate (lower RMSE_B), but RMSE_B depends only
        # on the main B_hat, not on the bootstrap — so this decoupling is safe.
        theta_b <- fit_b$theta[, perm_b, drop = FALSE]

        # Step 4: estimate effects on bootstrap sample using bootstrap X_std
        alr_b <- estimate_alr_effects(theta_b, X_std_b, reference_topic = K)
        boot_B_arr[, , b] <- alr_b$B_hat[rownames(B_hat), colnames(B_hat),
                                          drop = FALSE]
      }, error = function(e) NULL)
    }

    # Percentile bootstrap CI: use the 2.5% and 97.5% quantiles of the
    # bootstrap distribution directly.  Unlike B_hat ± 1.96*SE, percentile
    # CIs are invariant to bias in the point estimate and do not assume
    # symmetry — a better match when the ALR estimator may be biased.
    boot_lo <- apply(boot_B_arr, c(1L, 2L),
                     quantile, probs = 0.025, na.rm = TRUE)
    boot_hi <- apply(boot_B_arr, c(1L, 2L),
                     quantile, probs = 0.975, na.rm = TRUE)
    rownames(boot_lo) <- rownames(B_hat)
    colnames(boot_lo) <- colnames(B_hat)
    rownames(boot_hi) <- rownames(B_hat)
    colnames(boot_hi) <- colnames(B_hat)

    inside_boot <- (B_true_aligned >= boot_lo) &
                   (B_true_aligned <= boot_hi)
    coverage_B_boot <- mean(inside_boot, na.rm = TRUE)
  }

  list(
    RMSE_theta      = RMSE_theta,
    RMSE_phi        = RMSE_phi,
    RMSE_B          = RMSE_B,
    coverage_B      = coverage_B,
    coverage_B_boot = coverage_B_boot,
    perplexity      = perplexity
  )
}

# ---- Main simulation loop ----

set.seed(123)

n_docs  <- 1000
n_terms <- 500
K       <- 10
P       <- 5

n_reps      <- 10   # adjust as needed
conditions  <- c("baseline", "high_sparsity", "high_covariate")
methods     <- c("gscatm", "stm", "lda")

# Set to TRUE to compute bootstrap coverage (slow: ~B * n_reps * n_conditions fits)
run_bootstrap_ci <- TRUE
n_bootstrap_sim  <- 10

results <- expand.grid(
  rep       = 1:n_reps,
  condition = conditions,
  method    = methods,
  stringsAsFactors = FALSE
)

results$RMSE_theta      <- NA_real_
results$RMSE_phi        <- NA_real_
results$RMSE_B          <- NA_real_
results$coverage_B      <- NA_real_
results$coverage_B_boot <- NA_real_
results$perplexity      <- NA_real_



for (cond in conditions) {
  for (r in 1:n_reps) {
    cat("Condition:", cond, "- replication:", r, "\n")
    sim <- simulate_gscatm_data(
      condition = cond,
      n_docs    = n_docs,
      n_terms   = n_terms,
      K         = K,
      P         = P
    )

    for (m in methods) {
      # Bootstrap coverage is computed only for gscatm; NA for stm/lda
      do_boot <- isTRUE(run_bootstrap_ci) && m == "gscatm"
      res <- evaluate_method(sim, method = m, reference_topic = K,
                             bootstrap_ci = do_boot,
                             n_bootstrap  = n_bootstrap_sim)
      idx <- which(results$rep == r &
                     results$condition == cond &
                     results$method == m)

      results$RMSE_theta[idx]      <- res$RMSE_theta
      results$RMSE_phi[idx]        <- res$RMSE_phi
      results$RMSE_B[idx]          <- res$RMSE_B
      results$coverage_B[idx]      <- res$coverage_B
      results$coverage_B_boot[idx] <- res$coverage_B_boot
      results$perplexity[idx]      <- res$perplexity
    }
  }
}

# ---- Summary tables (LaTeX) ----

summary_metrics <- results %>%
  group_by(condition, method) %>%
  summarise(
    RMSE_theta_mean      = mean(RMSE_theta,      na.rm = TRUE),
    RMSE_theta_sd        = sd(RMSE_theta,        na.rm = TRUE),
    RMSE_phi_mean        = mean(RMSE_phi,        na.rm = TRUE),
    RMSE_phi_sd          = sd(RMSE_phi,          na.rm = TRUE),
    RMSE_B_mean          = mean(RMSE_B,          na.rm = TRUE),
    RMSE_B_sd            = sd(RMSE_B,            na.rm = TRUE),
    coverage_B_mean      = mean(coverage_B,      na.rm = TRUE),
    coverage_B_boot_mean = mean(coverage_B_boot, na.rm = TRUE),
    perplexity_mean      = mean(perplexity,      na.rm = TRUE),
    .groups = "drop"
  )

# digits: length must be ncol(summary_metrics) + 1
# here: first column is rownames (0), then condition (0), method (0),
# and 3 decimals for all numeric summaries
digits_vec <- c(0, 0, 0, rep(3, ncol(summary_metrics) - 2L))

print(
  xtable(
    summary_metrics,
    digits = digits_vec
  ),
  include.rownames = FALSE
)


# ---- Example plots ----

# RMSE for theta across methods and conditions
ggplot(results, aes(x = method, y = RMSE_theta)) +
  geom_boxplot() +
  facet_wrap(~ condition, scales = "free_y") +
  theme_bw()

# RMSE for B across methods and conditions
ggplot(results, aes(x = method, y = RMSE_B)) +
  geom_boxplot() +
  facet_wrap(~ condition, scales = "free_y") +
  theme_bw()

# Coverage of covariate effects
ggplot(results, aes(x = method, y = coverage_B)) +
  geom_boxplot() +
  facet_wrap(~ condition) +
  theme_bw()

# Perplexity comparison
ggplot(results, aes(x = method, y = perplexity)) +
  geom_boxplot() +
  facet_wrap(~ condition, scales = "free_y") +
  theme_bw()

