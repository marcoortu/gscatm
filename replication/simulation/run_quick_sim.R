pkg <- "c:/Users/Utente/OneDrive - Università/OneDrive - Università di Cagliari/UNIVERSITA/workspace/R/gscatm"
suppressWarnings(devtools::load_all(pkg, quiet = TRUE))
suppressPackageStartupMessages({
  library(Matrix); library(MCMCpack); library(clue)
  library(topicmodels); library(slam); library(dplyr); library(xtable)
})

# ========== Helpers (self-contained, no source()) ==========

compute_perplexity <- function(dtm, phi, theta) {
  theta <- as.matrix(theta); phi <- as.matrix(phi)
  lp <- log(theta %*% phi + 1e-10); lp[is.infinite(lp)] <- -700
  exp(-sum(dtm * lp, na.rm = TRUE) / sum(dtm))
}

match_topics <- function(phi_true, phi_hat) {
  K <- nrow(phi_true)
  cost <- matrix(0, K, K)
  for (k in seq_len(K)) for (j in seq_len(K))
    cost[k, j] <- sum(abs(phi_true[k, ] - phi_hat[j, ]))
  as.integer(clue::solve_LSAP(cost))
}

estimate_alr_effects <- function(theta, X_std, reference_topic = ncol(theta)) {
  theta <- as.matrix(theta); X_std <- as.matrix(X_std)
  K <- ncol(theta); P <- ncol(X_std)
  topic_idx <- c(setdiff(seq_len(K), reference_topic), reference_topic)
  th <- pmax(theta[, topic_idx, drop = FALSE], 1e-10)
  lr <- log(th[, seq_len(K - 1), drop = FALSE] / th[, K])
  nrt <- topic_idx[seq_len(K - 1)]
  if (is.null(colnames(X_std))) colnames(X_std) <- paste0("x", seq_len(P))
  B_hat  <- matrix(NA_real_, P, K - 1,
                   dimnames = list(colnames(X_std), paste0("Topic_", seq_len(K - 1))))
  se_hat <- B_hat
  for (k in seq_len(K - 1)) {
    y <- lr[, k]; w <- th[, k] * th[, K]
    fit <- glm(y ~ ., data = data.frame(y = y, X_std), weights = w, family = gaussian())
    ct  <- summary(fit)$coefficients
    ct  <- ct[setdiff(rownames(ct), "(Intercept)"), , drop = FALSE]
    for (j in seq_len(nrow(ct))) {
      nm <- rownames(ct)[j]
      if (nm %in% rownames(B_hat)) {
        B_hat[nm, k]  <- ct[j, "Estimate"]
        se_hat[nm, k] <- ct[j, "Std. Error"]
      }
    }
  }
  list(B_hat = B_hat, se_hat = se_hat)
}

simulate_gscatm_data <- function(condition = c("baseline", "high_sparsity", "high_covariate"),
                                 n_docs = 200, n_terms = 100, K = 3, P = 2,
                                 avg_len_base = 100, avg_len_sparse = 20) {
  condition <- match.arg(condition)
  cov_df <- as.data.frame(matrix(rnorm(n_docs * P), n_docs, P))
  colnames(cov_df) <- paste0("x", seq_len(P))
  X_std <- scale(model.matrix(~ . - 1, data = cov_df))
  K_nonref <- K - 1
  B_true <- matrix(0, P, K_nonref,
                   dimnames = list(colnames(X_std), paste0("Topic_", seq_len(K_nonref))))
  for (p in seq_len(P))
    B_true[p, ] <- seq(-0.6, 0.6, length.out = K_nonref) * ((-1)^(p - 1))
  eta_mean <- X_std %*% B_true
  noise_sd <- switch(condition, baseline = 0.3, high_sparsity = 0.4, high_covariate = 0.15)
  eta  <- eta_mean + matrix(rnorm(n_docs * K_nonref, sd = noise_sd), n_docs, K_nonref)
  ee   <- exp(cbind(eta, 0)); theta <- ee / rowSums(ee)
  if (condition == "high_covariate")
    theta <- t(apply(theta, 1, function(mu)
      as.numeric(MCMCpack::rdirichlet(1, 50 * mu))))
  if (condition == "high_sparsity") {
    theta[theta < 0.03] <- 0; rs <- rowSums(theta)
    zr <- which(rs == 0)
    if (length(zr) > 0) theta[zr, ] <- ee[zr, ] / rowSums(ee[zr, , drop = FALSE])
    theta <- theta / rowSums(theta)
  }
  phi_true <- t(replicate(K, as.numeric(MCMCpack::rdirichlet(1, rep(0.1, n_terms)))))
  phi_true <- phi_true / rowSums(phi_true)
  colnames(phi_true) <- paste0("w", seq_len(n_terms))
  rownames(phi_true) <- paste0("Topic_", seq_len(K))
  lambda <- if (condition == "high_sparsity") avg_len_sparse else avg_len_base
  L <- pmax(rpois(n_docs, lambda), 5L)
  dm <- matrix(0L, n_docs, n_terms)
  for (i in seq_len(n_docs))
    dm[i, ] <- rmultinom(1L, L[i], theta[i, ] %*% phi_true)
  dtm <- Matrix::Matrix(dm, sparse = TRUE)
  colnames(dtm) <- colnames(phi_true); rownames(dtm) <- paste0("doc", seq_len(n_docs))
  keep <- Matrix::colSums(dtm) > 0
  dtm <- dtm[, keep]; phi_true <- phi_true[, keep]
  phi_true <- phi_true / rowSums(phi_true)
  list(condition = condition, dtm = dtm, covariates = cov_df, X_std = X_std,
       theta_true = theta, phi_true = phi_true, B_true = B_true)
}

bootstrap_coverage <- function(theta_g, phi_g, dtm, covariates, X_std,
                               B_true, K, model_type, n_bootstrap, max_iter) {
  N <- nrow(dtm); L <- rowSums(dtm)
  q_hat   <- as.matrix(theta_g %*% phi_g)
  P_eff   <- nrow(B_true); K_nonref <- ncol(B_true)
  boot_B  <- array(NA_real_, c(P_eff, K_nonref, n_bootstrap))
  for (b in seq_len(n_bootstrap)) {
    tryCatch({
      dtm_b <- Matrix::Matrix(0L, N, ncol(dtm), sparse = TRUE)
      colnames(dtm_b) <- colnames(dtm)
      for (i in seq_len(N)) {
        if (L[i] > 0L) {
          q_i <- pmax(q_hat[i, ], 0); q_i <- q_i / sum(q_i)
          dtm_b[i, ] <- as.integer(rmultinom(1L, L[i], q_i))
        }
      }
      dtm_b <- as(dtm_b, "dgCMatrix")
      fit_b <- suppressMessages(
        fit_gsca_topic_model(dtm_b, covariates, K = K, model = model_type,
                             reference_topic = K, max_iter = max_iter, bootstrap_ci = FALSE)
      )
      perm_b  <- match_topics(phi_g, fit_b$phi)
      theta_b <- fit_b$theta[, perm_b, drop = FALSE]
      theta_b <- theta_b / rowSums(theta_b)
      alr_b   <- estimate_alr_effects(theta_b, X_std, K)
      boot_B[, , b] <- alr_b$B_hat[rownames(B_true), colnames(B_true), drop = FALSE]
    }, error = function(e) NULL)
  }
  boot_se <- apply(boot_B, c(1L, 2L), sd, na.rm = TRUE)
  mean((B_true >= (rowSums(boot_B, dims = 1L) * 0 + # reuse B_hat from outer scope
         boot_B[, , 1L] * NA) | TRUE), na.rm = TRUE) # placeholder — computed below
  list(boot_se = boot_se, boot_B = boot_B)
}

evaluate_one <- function(sim, n_bootstrap = 10, max_iter = 50) {
  dtm        <- sim$dtm; K <- ncol(sim$theta_true)
  covariates <- sim$covariates; X_std <- sim$X_std
  phi_true   <- sim$phi_true; B_true <- sim$B_true; theta_true <- sim$theta_true
  model_type <- switch(sim$condition,
    baseline = "logistic-normal", high_sparsity = "zero-inflated",
    high_covariate = "dirichlet")

  # --- GSCA-TM ---
  fit_g   <- fit_gsca_topic_model(dtm, covariates, K = K, model = model_type,
                                  reference_topic = K, max_iter = max_iter,
                                  bootstrap_ci = FALSE)
  perm_g  <- match_topics(phi_true, fit_g$phi)
  # theta_g   = softmax(X*B): covariate prediction — used for RMSE_theta
  # theta_g_p = posterior (word-evidence): used for B_hat and coverage
  # Separation avoids reference-topic alignment artefacts in estimate_alr_effects
  # when softmax(X*B) is used (rigid structure → wrong ALR after permutation).
  theta_g   <- fit_g$theta[,           perm_g, drop = FALSE]
  theta_g   <- theta_g / rowSums(theta_g)
  theta_g_p <- fit_g$theta_posterior[, perm_g, drop = FALSE]
  theta_g_p <- theta_g_p / rowSums(theta_g_p)
  phi_g     <- fit_g$phi[perm_g, , drop = FALSE]; phi_g <- phi_g / rowSums(phi_g)
  alr_g   <- estimate_alr_effects(theta_g_p, X_std, K)
  B_hat_g <- alr_g$B_hat[rownames(B_true), colnames(B_true), drop = FALSE]
  se_g    <- alr_g$se_hat[rownames(B_true), colnames(B_true), drop = FALSE]

  cov_g  <- mean((B_true >= (B_hat_g - 1.96 * se_g)) &
                 (B_true <= (B_hat_g + 1.96 * se_g)), na.rm = TRUE)

  # Bootstrap coverage for GSCA-TM — non-parametric row-bootstrap.
  # Resamples documents (rows) with replacement to capture BOTH sources of
  # variability: latent-theta noise AND word-count sampling.
  # Both B_hat_g and boot_B are computed from theta_posterior (consistent scale).
  N <- nrow(dtm)
  boot_B <- array(NA_real_, c(nrow(B_hat_g), ncol(B_hat_g), n_bootstrap))
  for (b in seq_len(n_bootstrap)) {
    tryCatch({
      idx_b        <- sample.int(N, N, replace = TRUE)
      dtm_b        <- dtm[idx_b, , drop = FALSE]
      covariates_b <- covariates[idx_b, , drop = FALSE]
      X_std_b      <- X_std[idx_b, , drop = FALSE]
      fit_b <- suppressMessages(
        fit_gsca_topic_model(dtm_b, covariates_b, K = K, model = model_type,
                             reference_topic = K, max_iter = max_iter, bootstrap_ci = FALSE)
      )
      perm_b    <- match_topics(phi_g, fit_b$phi)
      theta_b_p <- fit_b$theta_posterior[, perm_b, drop = FALSE]
      theta_b_p <- theta_b_p / rowSums(theta_b_p)
      alr_b   <- estimate_alr_effects(theta_b_p, X_std_b, K)
      boot_B[, , b] <- alr_b$B_hat[rownames(B_hat_g), colnames(B_hat_g), drop = FALSE]
    }, error = function(e) NULL)
  }
  boot_se <- apply(boot_B, c(1L, 2L), sd, na.rm = TRUE)
  cov_boot_g <- mean((B_true >= (B_hat_g - 1.96 * boot_se)) &
                     (B_true <= (B_hat_g + 1.96 * boot_se)), na.rm = TRUE)

  RMSE_theta_g <- sqrt(mean((theta_g   - theta_true)^2))   # theta_gsca vs theta_true
  RMSE_phi_g   <- sqrt(mean((phi_g     - phi_true)^2))
  RMSE_B_g     <- sqrt(mean((B_hat_g   - B_true)^2, na.rm = TRUE))
  perp_g       <- compute_perplexity(dtm, phi_g, theta_g)

  # --- LDA ---
  dtm_trip <- slam::as.simple_triplet_matrix(dtm)
  lda_fit  <- topicmodels::LDA(dtm_trip, k = K, control = list(seed = 1))
  post     <- topicmodels::posterior(lda_fit)
  theta_l_raw <- post$topics; phi_l_raw <- post$terms
  perm_l   <- match_topics(phi_true, phi_l_raw)
  theta_l  <- theta_l_raw[, perm_l, drop = FALSE]; theta_l <- theta_l / rowSums(theta_l)
  phi_l    <- phi_l_raw[perm_l, , drop = FALSE];   phi_l   <- phi_l   / rowSums(phi_l)
  alr_l    <- estimate_alr_effects(theta_l, X_std, K)
  B_hat_l  <- alr_l$B_hat[rownames(B_true), colnames(B_true), drop = FALSE]
  se_l     <- alr_l$se_hat[rownames(B_true), colnames(B_true), drop = FALSE]
  cov_l    <- mean((B_true >= (B_hat_l - 1.96 * se_l)) &
                   (B_true <= (B_hat_l + 1.96 * se_l)), na.rm = TRUE)
  RMSE_theta_l <- sqrt(mean((theta_l - theta_true)^2))
  RMSE_phi_l   <- sqrt(mean((phi_l   - phi_true)^2))
  RMSE_B_l     <- sqrt(mean((B_hat_l - B_true)^2, na.rm = TRUE))
  perp_l       <- compute_perplexity(dtm, phi_l, theta_l)

  list(
    gscatm = list(RMSE_theta = RMSE_theta_g, RMSE_phi = RMSE_phi_g, RMSE_B = RMSE_B_g,
                  coverage_B = cov_g, coverage_B_boot = cov_boot_g, perplexity = perp_g),
    lda    = list(RMSE_theta = RMSE_theta_l, RMSE_phi = RMSE_phi_l, RMSE_B = RMSE_B_l,
                  coverage_B = cov_l, coverage_B_boot = NA_real_,   perplexity = perp_l)
  )
}

# ========== Quick simulation ==========
set.seed(123)
N_REPS    <- 10     # 3 conditions x 10 reps = 30 datasets
N_DOCS    <- 300
N_TERMS   <- 150
K_SIM     <- 3
P_SIM     <- 2
B_SIM     <- 50     # bootstrap replications per fit
MAXITER   <- 50
conditions <- c("baseline", "high_sparsity", "high_covariate")
methods    <- c("gscatm", "lda")

records <- list()
for (cond in conditions) {
  cat("Condition:", cond, "\n")
  for (r in seq_len(N_REPS)) {
    cat("  rep", r, "/", N_REPS, "\n")
    sim <- simulate_gscatm_data(cond, N_DOCS, N_TERMS, K_SIM, P_SIM)
    res <- evaluate_one(sim, n_bootstrap = B_SIM, max_iter = MAXITER)
    for (m in methods)
      records[[length(records) + 1]] <- data.frame(
        condition = cond, method = m, rep = r,
        RMSE_theta      = res[[m]]$RMSE_theta,
        RMSE_phi        = res[[m]]$RMSE_phi,
        RMSE_B          = res[[m]]$RMSE_B,
        coverage_B      = res[[m]]$coverage_B,
        coverage_B_boot = res[[m]]$coverage_B_boot,
        perplexity      = res[[m]]$perplexity,
        stringsAsFactors = FALSE)
  }
}
results <- do.call(rbind, records)
saveRDS(results, file.path(pkg, "scripts/quick_sim_results.rds"))
cat("\nResults saved.\n")

# ========== LaTeX Table 1 ==========
summary_tbl <- results %>%
  group_by(condition, method) %>%
  summarise(
    RMSE_theta_mean = mean(RMSE_theta, na.rm = TRUE),
    RMSE_theta_sd   = sd(RMSE_theta,   na.rm = TRUE),
    RMSE_B_mean     = mean(RMSE_B,     na.rm = TRUE),
    RMSE_B_sd       = sd(RMSE_B,       na.rm = TRUE),
    coverage_B      = mean(coverage_B,      na.rm = TRUE),
    coverage_B_boot = mean(coverage_B_boot, na.rm = TRUE),
    perplexity      = mean(perplexity,      na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    Method    = toupper(method),
    Condition = gsub("_", " ", condition),
    `RMSE(Theta)` = sprintf("%.3f (%.3f)", RMSE_theta_mean, RMSE_theta_sd),
    `RMSE(B)`     = sprintf("%.3f (%.3f)", RMSE_B_mean,     RMSE_B_sd),
    `Cov(plug-in)`   = sprintf("%.3f", coverage_B),
    `Cov(bootstrap)` = ifelse(is.na(coverage_B_boot), "---",
                              sprintf("%.3f", coverage_B_boot)),
    Perplexity = sprintf("%.1f", perplexity)
  ) %>%
  select(Condition, Method, `RMSE(Theta)`, `RMSE(B)`,
         `Cov(plug-in)`, `Cov(bootstrap)`, Perplexity)

cat("\n\n========== LaTeX TABLE 1 ==========\n\n")
print(
  xtable(as.data.frame(summary_tbl),
         caption = paste0("Simulation results (mean and SD over ", N_REPS,
                          " replications; $N=", N_DOCS, "$, $V=", N_TERMS,
                          "$, $K=", K_SIM, "$, $P=", P_SIM,
                          "$, $B=", B_SIM, "$). ",
                          "Coverage target: 0.950. ",
                          "``---'' = bootstrap not implemented for LDA."),
         label   = "tab:simulation"),
  include.rownames = FALSE,
  booktabs         = TRUE,
  sanitize.text.function = identity,
  floating         = TRUE,
  table.placement  = "t"
)

cat("\n\n========== RAW NUMBERS ==========\n")
print(as.data.frame(summary_tbl), row.names = FALSE)
