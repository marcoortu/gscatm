# Simulation and benchmarking script for GSCA-TM vs STM and LDA+ALR
# Parallelized version using foreach + doParallel
#
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
library(parallel)
library(pbapply)    # progress bars for parallel apply
# ---- Load gscatm ----
# Option 1 (recommended): install once, then library(gscatm) works everywhere
#   devtools::install(".", upgrade = "never", quiet = TRUE)
# Option 2: load from source for the main process only
pkg_root <- normalizePath(file.path(getwd(), "../.."), mustWork = FALSE)
if (!file.exists(file.path(pkg_root, "DESCRIPTION"))) {
  pkg_root <- normalizePath(".", mustWork = FALSE)
}
cat("Loading gscatm from:", pkg_root, "\n")
devtools::load_all(pkg_root)

# Install the package so workers can use library(gscatm) without devtools
# This avoids 22+ workers all doing devtools::load_all() on the same directory
cat("Installing gscatm for parallel workers...\n")
devtools::install(pkg_root, upgrade = FALSE, quiet = TRUE)

# ---- Parallelization setup ----
n_cores <- parallel::detectCores() - 1L   # leave 1 core free
n_cores <- max(n_cores, 1L)
cat("Using", n_cores, "cores\n")

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

  topic_idx <- 1:K
  topic_idx <- c(topic_idx[topic_idx != reference_topic], reference_topic)
  theta_ordered <- theta[, topic_idx, drop = FALSE]

  theta_safe <- pmax(theta_ordered, 1e-10)

  denom <- theta_safe[, K]
  log_ratios <- log(theta_safe[, 1:(K - 1), drop = FALSE] / denom)

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
    weights <- theta_safe[, k] * theta_safe[, K]
    fit <- glm(y ~ ., data = glm_data, weights = weights, family = gaussian())

    coef_table <- summary(fit)$coefficients
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

# ---- Helper: generate one parametric bootstrap DTM ----
generate_bootstrap_dtm <- function(theta_hat, phi_hat, dtm) {
  q_hat <- as.matrix(theta_hat %*% phi_hat)
  L     <- rowSums(dtm)
  N     <- nrow(dtm)
  V     <- ncol(dtm)

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
  as(dtm_boot, "dgCMatrix")
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

  covariates_df <- as.data.frame(matrix(rnorm(n_docs * P),
                                        nrow = n_docs, ncol = P))
  colnames(covariates_df) <- paste0("x", 1:P)

  X_design <- model.matrix(~ . - 1, data = covariates_df)
  X_std <- scale(X_design)
  P_eff <- ncol(X_std)

  K_nonref <- K - 1
  B_true <- matrix(0, nrow = P_eff, ncol = K_nonref)
  for (p in 1:P_eff) {
    B_true[p, ] <- seq(-0.6, 0.6, length.out = K_nonref) * ((-1)^(p - 1))
  }
  rownames(B_true) <- colnames(X_std)
  colnames(B_true) <- paste0("Topic_", 1:K_nonref)

  eta_mean <- X_std %*% B_true
  noise_sd <- switch(condition,
                     "baseline"       = 0.3,
                     "high_sparsity"  = 0.4,
                     "high_covariate" = 0.15)
  eta <- eta_mean + matrix(rnorm(n_docs * K_nonref, sd = noise_sd),
                           nrow = n_docs, ncol = K_nonref)

  eta_full <- cbind(eta, 0)
  exp_eta  <- exp(eta_full)
  theta_base <- exp_eta / rowSums(exp_eta)
  theta <- theta_base

  if (condition == "high_covariate") {
    alpha_bar <- 50
    theta <- t(apply(theta_base, 1, function(mu)
      as.numeric(MCMCpack::rdirichlet(1, alpha_bar * mu))))
  }

  if (condition == "high_sparsity") {
    thresh <- 0.03
    theta[theta < thresh] <- 0
    rs <- rowSums(theta)
    zero_rows <- which(rs == 0)
    if (length(zero_rows) > 0) {
      theta[zero_rows, ] <- theta_base[zero_rows, ]
      rs <- rowSums(theta)
    }
    theta <- theta / rs
  }

  beta0 <- rep(0.1, n_terms)
  phi_true <- t(replicate(K, as.numeric(MCMCpack::rdirichlet(1, beta0))))
  phi_true <- phi_true / rowSums(phi_true)
  colnames(phi_true) <- paste0("w", 1:n_terms)
  rownames(phi_true) <- paste0("Topic_", 1:K)

  lambda <- if (condition == "high_sparsity") avg_doc_length_sparse else avg_doc_length_baseline
  L <- rpois(n_docs, lambda = lambda)
  L[L < 5] <- 5

  dtm_dense <- matrix(0, nrow = n_docs, ncol = n_terms)
  for (i in 1:n_docs) {
    q_i <- as.numeric(theta[i, ] %*% phi_true)
    dtm_dense[i, ] <- rmultinom(1, size = L[i], prob = q_i)
  }
  dtm <- Matrix::Matrix(dtm_dense, sparse = TRUE)
  colnames(dtm) <- colnames(phi_true)
  rownames(dtm) <- paste0("doc", 1:n_docs)

  col_sums <- Matrix::colSums(dtm)
  keep <- col_sums > 0
  dtm <- dtm[, keep, drop = FALSE]
  phi_true <- phi_true[, keep, drop = FALSE]
  phi_true <- phi_true / rowSums(phi_true)

  list(condition = condition, dtm = dtm, covariates = covariates_df,
       X_std = X_std, theta_true = theta, phi_true = phi_true, B_true = B_true)
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

  model_type <- if (method == "gscatm") {
    switch(sim$condition,
           "baseline"       = "logistic-normal",
           "high_sparsity"  = "zero-inflated",
           "high_covariate" = "dirichlet")
  } else NA_character_

  # --- Fit each model ---
  if (method == "gscatm") {
    fit <- fit_gsca_topic_model(dtm = dtm, covariates = covariates,
                                K = K, model = model_type,
                                reference_topic = reference_topic)
    theta_hat_raw <- fit$theta
    phi_hat_raw   <- fit$phi

  } else if (method == "stm") {
    stm_in <- dtm_to_stm(dtm)
    stm_fit <- stm::stm(documents = stm_in$documents, vocab = stm_in$vocab,
                         K = K, prevalence = ~ ., data = covariates,
                         max.em.its = 75, init.type = "Spectral", verbose = FALSE)
    theta_hat_raw <- stm_fit$theta
    logbeta <- stm_fit$beta$logbeta[[1]]
    phi_hat_raw <- exp(logbeta)
    phi_hat_raw <- phi_hat_raw / rowSums(phi_hat_raw)

  } else if (method == "lda") {
    dtm_triplet <- slam::as.simple_triplet_matrix(dtm)
    lda_fit <- topicmodels::LDA(dtm_triplet, k = K, control = list(seed = 1))
    post <- topicmodels::posterior(lda_fit)
    theta_hat_raw <- post$topics
    phi_hat_raw   <- post$terms
  }

  # --- Match topics to true labels ---
  perm <- match_topics(phi_true, phi_hat_raw)
  phi_hat   <- phi_hat_raw[perm, , drop = FALSE]
  theta_hat <- theta_hat_raw[, perm, drop = FALSE]
  phi_hat   <- phi_hat / rowSums(phi_hat)
  theta_hat <- theta_hat / rowSums(theta_hat)

  # --- Metrics (UNCHANGED from original) ---
  RMSE_theta <- sqrt(mean((theta_hat - theta_true)^2))
  RMSE_phi   <- sqrt(mean((phi_hat   - phi_true)^2))

  alr    <- estimate_alr_effects(theta_hat, X_std, reference_topic = K)
  B_hat  <- alr$B_hat
  se_hat <- alr$se_hat
  B_true_aligned <- B_true[rownames(B_hat), colnames(B_hat), drop = FALSE]
  RMSE_B <- sqrt(mean((B_hat - B_true_aligned)^2))

  inside <- (B_true_aligned >= (B_hat - 1.96 * se_hat)) &
    (B_true_aligned <= (B_hat + 1.96 * se_hat))
  coverage_B <- mean(inside, na.rm = TRUE)

  perplexity <- compute_perplexity(dtm, phi_hat, theta_hat)

  # --- Bootstrap coverage for GSCATM only (NEW) ---
  coverage_B_boot <- NA_real_

  if (isTRUE(bootstrap_ci) && method == "gscatm") {
    P_eff    <- nrow(B_hat)
    K_nonref <- ncol(B_hat)
    boot_B_arr <- array(NA_real_, dim = c(P_eff, K_nonref, n_bootstrap))

    for (b in seq_len(n_bootstrap)) {
      tryCatch({
        dtm_boot <- generate_bootstrap_dtm(theta_hat, phi_hat, dtm)
        fit_b <- suppressMessages(
          fit_gsca_topic_model(dtm_boot, covariates, K = K, model = model_type,
                               reference_topic = K, bootstrap_ci = FALSE))
        perm_b  <- match_topics(phi_hat, fit_b$phi)
        theta_b <- fit_b$theta[, perm_b, drop = FALSE]
        theta_b <- theta_b / rowSums(theta_b)
        alr_b   <- estimate_alr_effects(theta_b, X_std, reference_topic = K)
        boot_B_arr[, , b] <- alr_b$B_hat[rownames(B_hat), colnames(B_hat), drop = FALSE]
      }, error = function(e) NULL)
    }

    boot_se_mat <- apply(boot_B_arr, c(1L, 2L), sd, na.rm = TRUE)
    inside_boot <- (B_true_aligned >= (B_hat - 1.96 * boot_se_mat)) &
                   (B_true_aligned <= (B_hat + 1.96 * boot_se_mat))
    coverage_B_boot <- mean(inside_boot, na.rm = TRUE)
  }

  list(RMSE_theta = RMSE_theta, RMSE_phi = RMSE_phi, RMSE_B = RMSE_B,
       coverage_B = coverage_B, coverage_B_boot = coverage_B_boot,
       perplexity = perplexity)
}


# ---- Worker: one (condition, replication) pair ----
run_one_replication <- function(task, n_docs, n_terms, K, P,
                                methods, run_bootstrap_ci, n_bootstrap_sim) {
  cond <- task$cond
  r    <- task$r
  sim  <- simulate_gscatm_data(condition = cond, n_docs = n_docs, n_terms = n_terms,
                               K = K, P = P, avg_doc_length_baseline = 100)
  rows <- list()
  for (m in methods) {
    do_boot <- isTRUE(run_bootstrap_ci) && m == "gscatm"
    res <- evaluate_method(sim, method = m, reference_topic = K,
                           bootstrap_ci = do_boot, n_bootstrap = n_bootstrap_sim)
    rows[[m]] <- data.frame(
      rep = r, condition = cond, method = m,
      RMSE_theta = res$RMSE_theta, RMSE_phi = res$RMSE_phi,
      RMSE_B = res$RMSE_B, coverage_B = res$coverage_B,
      coverage_B_boot = res$coverage_B_boot, perplexity = res$perplexity,
      stringsAsFactors = FALSE)
  }
  do.call(rbind, rows)
}


# ---- Output directory ----

output_dir <- normalizePath(
  file.path(pkg_root, "replication", "simulation", "output", "parametric_sims"),
  mustWork = FALSE
)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
cat("Output directory:", output_dir, "\n")

# ---- Main simulation loop (PARALLEL) ----

set.seed(123)

n_docs  <- 100
n_terms <- 500
K       <- 10
P       <- 8

n_reps      <- 100
conditions  <- c("baseline", "high_sparsity", "high_covariate")
methods     <- c("gscatm", "stm", "lda")

run_bootstrap_ci <- TRUE
n_bootstrap_sim  <- 200

# Build task list: one entry per (condition, replication) pair
tasks <- expand.grid(cond = conditions, r = 1:n_reps, stringsAsFactors = FALSE)
task_list <- split(tasks, seq_len(nrow(tasks)))

cat("Running", length(task_list), "tasks across", n_cores, "cores...\n")
t_start <- Sys.time()

# Create cluster and export everything workers need
cl <- makeCluster(n_cores)
on.exit(stopCluster(cl), add = TRUE)

# Load packages and gscatm on each worker
clusterExport(cl, varlist = "pkg_root", envir = environment())
clusterEvalQ(cl, {
  library(Matrix)
  library(MCMCpack)
  library(stm)
  library(topicmodels)
  library(slam)
  library(clue)
  devtools::load_all(pkg_root)
})

# Export helper functions defined in this script + simulation parameters
clusterExport(cl, varlist = c(
  "simulate_gscatm_data", "evaluate_method", "run_one_replication",
  "generate_bootstrap_dtm", "compute_perplexity", "match_topics",
  "estimate_alr_effects", "dtm_to_stm",
  "n_docs", "n_terms", "K", "P", "methods",
  "run_bootstrap_ci", "n_bootstrap_sim"
), envir = environment())

# Reproducible parallel RNG
clusterSetRNGStream(cl, iseed = 123)

results_list <- pblapply(task_list, function(task) {
  run_one_replication(task, n_docs = n_docs, n_terms = n_terms,
                      K = K, P = P, methods = methods,
                      run_bootstrap_ci = run_bootstrap_ci,
                      n_bootstrap_sim = n_bootstrap_sim)
}, cl = cl)

stopCluster(cl)
on.exit(NULL)   # cancel the on.exit now that we stopped manually

results <- do.call(rbind, results_list)
rownames(results) <- NULL

t_end <- Sys.time()
cat("Total time:", format(t_end - t_start), "\n")


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

digits_vec <- c(0, 0, 0, rep(3, ncol(summary_metrics) - 2))

print(
  xtable(summary_metrics, digits = digits_vec),
  include.rownames = FALSE
)

# ---- Save results ----

saveRDS(results,         file.path(output_dir, "results.rds"))
saveRDS(summary_metrics, file.path(output_dir, "summary_metrics.rds"))
write.csv(results,         file.path(output_dir, "results.csv"),         row.names = FALSE)
write.csv(summary_metrics, file.path(output_dir, "summary_metrics.csv"), row.names = FALSE)

latex_file <- file.path(output_dir, "summary_table.tex")
sink(latex_file)
print(xtable(summary_metrics, digits = digits_vec), include.rownames = FALSE)
sink()

cat("Results saved to:", output_dir, "\n")

# ---- Plots ----

p_rmse_theta <- ggplot(results, aes(x = method, y = RMSE_theta)) +
  geom_boxplot() +
  facet_wrap(~ condition, scales = "free_y") +
  theme_bw() +
  labs(title = "RMSE theta by method and condition")

p_rmse_B <- ggplot(results, aes(x = method, y = RMSE_B)) +
  geom_boxplot() +
  facet_wrap(~ condition, scales = "free_y") +
  theme_bw() +
  labs(title = "RMSE B (covariate effects) by method and condition")

p_coverage <- ggplot(results, aes(x = method, y = coverage_B)) +
  geom_boxplot() +
  facet_wrap(~ condition) +
  geom_hline(yintercept = 0.95, linetype = "dashed", color = "red") +
  theme_bw() +
  labs(title = "Analytical (Wald) coverage of B_true")

p_coverage_boot <- ggplot(results %>% filter(method == "gscatm"),
       aes(x = condition, y = coverage_B_boot)) +
  geom_boxplot() +
  geom_hline(yintercept = 0.95, linetype = "dashed", color = "red") +
  theme_bw() +
  labs(title = "Parametric bootstrap coverage for GSCATM (target = 0.95)")

p_perplexity <- ggplot(results, aes(x = method, y = perplexity)) +
  geom_boxplot() +
  facet_wrap(~ condition, scales = "free_y") +
  theme_bw() +
  labs(title = "Perplexity by method and condition")

print(p_rmse_theta)
print(p_rmse_B)
print(p_coverage)
print(p_coverage_boot)
print(p_perplexity)

ggsave(file.path(output_dir, "rmse_theta.pdf"),    p_rmse_theta,    width = 9, height = 4)
ggsave(file.path(output_dir, "rmse_B.pdf"),        p_rmse_B,        width = 9, height = 4)
ggsave(file.path(output_dir, "coverage_B.pdf"),    p_coverage,      width = 9, height = 4)
ggsave(file.path(output_dir, "coverage_boot.pdf"), p_coverage_boot, width = 7, height = 4)
ggsave(file.path(output_dir, "perplexity.pdf"),    p_perplexity,    width = 9, height = 4)

cat("Plots saved to:", output_dir, "\n")






