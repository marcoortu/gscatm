library(gscatm)
library(quanteda)
library(quanteda.textstats)
library(dplyr)
library(ggplot2)
library(dplyr)
library(tidyr)
library(reshape2)
library(lubridate)  # For working with dates



# Load necessary libraries
library(dplyr)
library(readr)

# File paths
files <- c("../data/rawtext_DonaldTrump.tsv",
           "../data/rawtext_JoeBiden.tsv",
           "../data/rawtext_KamalaHarris.tsv",
           "../data/rawtext_MikePence.tsv")

# Function to read TSV files
read_tsv_files <- function(file) {
  read_tsv(file, col_types = cols())  # Automatically detects column types
}

# Read and merge all files
merged_data <- files %>%
  lapply(read_tsv_files) %>%
  bind_rows()

merged_data <- merged_data %>%
  mutate(
    party = case_when(
      POTUS %in% c("Donald Trump", "Mike Pence") ~ "republican",
      POTUS %in% c("Joe Biden", "Kamala Harris") ~ "democratic",
      TRUE ~ "unknown"  # Catch any unexpected cases
    ),
    candidate_type = case_when(
      (POTUS == "Donald Trump" & year(Date) < 2024) ~ "president",
      (POTUS == "Joe Biden" & year(Date) < 2024) ~ "president",
      (POTUS == "Kamala Harris" & year(Date) >= 2024) ~ "president",
      (POTUS == "Kamala Harris" & year(Date) < 2024) ~ "vice president",
      POTUS == "Mike Pence" ~ "vice president",
    ),
    year = year(ymd(Date))  # Extract year from the Date column
  )

# Display the merged dataframe structure
glimpse(merged_data)

# Write the dataset to a CSV file
write.csv(merged_data, "usa_election_speeches_2024.csv", row.names = FALSE)



df <- read.csv("usa_election_speeches_2024.csv")



# Create a sample dataframe with texts and covariates
df_sotu <- df %>%
  select(c("RawText", "year", "party", "candidate_type"))

# Sample if needed (optional)
# df_sample <- sample_n(df_sotu, size = 100)

# Select covariates
covariates <- df_sotu %>%
  select(c("year", "party", "candidate_type"))

# Process the documents
docs <- df_sotu$RawText
tokens <- tokens(docs, remove_punct = TRUE)
tokens <- tokens_remove(tokens, stopwords("en"))
# tokens <- tokens_wordstem(tokens)

# Create the DFM
dfm <- dfm(tokens)

# Identify the most frequent tokens (e.g., top 10)
freq <- textstat_frequency(dfm)
most_freq_terms <- freq$feature[1:10]  # Replace 10 with the number of tokens you want to remove

# Option 1: Remove tokens from the tokens object
tokens <- tokens_remove(tokens, pattern = most_freq_terms)

# Recreate the DFM after removing the tokens
# Recreate dfm and proceed
dfm <- dfm(tokens)
dtm <- as(dfm, "dgCMatrix")


K_seq <- 2:30
result <- search_optimal_topics(dtm, covariates, K_seq, n=3)
# View the top K values
print(result$top_K)

# Set number of topics
n_topics = 9

# Fit the model
model <- fit_gsca_topic_model(
  dtm          = dtm,
  covariates   = covariates,
  K            = n_topics,
  model        = "zero-inflated", #"logistic-normal", "dirichlet", "zero-inflated"
  max_iter     = 100,
  tol          = 1e-4,
  bootstrap_ci = TRUE,
  n_bootstrap  = 200
)

topic_terms_distribution(model, top_n_terms = 10)




# Print the summary table
print(model$effects)

print(model$fit)

# Add confidence intervals to the summary table


# Plotting
ggplot(model$effects, aes(x = Covariate, y = Log_Odds_Coefficient, color = Topic)) +
  geom_point(position = position_dodge(width = 0.5)) +
  geom_errorbar(aes(ymin = CI_Lower, ymax = CI_Upper), width = 0.2,
                position = position_dodge(width = 0.5)) +
  facet_wrap(~ Topic, scales = "free_y") +
  theme_bw() +
  labs(title = "Effects of Covariates on Topics",
       y = "Log-Odds Coefficient",
       x = "Covariate") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Prepare data for heatmap
heatmap_data <- dcast(model$effects, Covariate ~ Topic, value.var = "Log_Odds_Coefficient")

# Convert to matrix
rownames(heatmap_data) <- heatmap_data$Covariate
heatmap_matrix <- as.matrix(heatmap_data[, -1])

# Plot heatmap
heatmap(heatmap_matrix, Rowv = NA, Colv = NA, scale = "none",
        col = colorRampPalette(c("blue", "white", "red"))(100),
        margins = c(5, 10))



# Extract the effects data frame
effects_df <- model$effects

# Convert the effects data frame to a long format suitable for plotting
effects_long <- effects_df %>%
  mutate(Significant = adjusted_p_value < 0.05) %>%
  select(Covariate, Topic, Log_Odds_Coefficient, Std_Error, Significant) %>%
  mutate(Topic = factor(Topic, levels = unique(Topic)))

# Basic coefficient plot
ggplot(effects_long, aes(x = Covariate, y = Log_Odds_Coefficient, color = Topic)) +
  geom_point(position = position_dodge(width = 0.5), size = 3) +
  geom_errorbar(aes(ymin = Log_Odds_Coefficient - 1.96 * Std_Error,
                    ymax = Log_Odds_Coefficient + 1.96 * Std_Error),
                position = position_dodge(width = 0.5), width = 0.2) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  theme_minimal() +
  labs(title = "Effects of Covariates on Topic Proportions",
       x = "Covariate",
       y = "Log-Odds Coefficient") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


# Faceted plot by Topic
ggplot(effects_long, aes(x = Covariate, y = Log_Odds_Coefficient, color = Significant)) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = Log_Odds_Coefficient - 1.96 * Std_Error,
                    ymax = Log_Odds_Coefficient + 1.96 * Std_Error),
                width = 0.2) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  facet_wrap(~ Topic, scales = "free_y") +
  theme_minimal() +
  labs(title = "Effects of Covariates on Topic Proportions by Topic",
       x = "Covariate",
       y = "Log-Odds Coefficient") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))



# Select a continuous covariate (e.g., "Covariate1")
covariate_name <- "year"

# Create a sequence of values for the covariate
covariate_seq <- seq(min(covariates[, covariate_name]),
                     max(covariates[, covariate_name]), length.out = 4)

# Function to predict topic proportions
predict_topic_proportions <- function(value) {
  new_covariates <- covariates
  new_covariates[, covariate_name] <- value
  # Recalculate gamma and theta
  gamma_new <- model.matrix(~ . - 1, data = new_covariates) %*% model$path_matrix[1:model$P, ]
  gamma_new <- cbind(gamma_new, model$gamma[, (model$P + 1):(model$P + model$K)])
  topic_scores <- gamma_new[, (model$P + 1):(model$P + model$K)]
  theta_new <- exp(topic_scores)
  theta_new <- theta_new / rowSums(theta_new)
  # Return mean topic proportions
  colMeans(theta_new)
}

# Predict topic proportions over the covariate range
predicted_proportions <- sapply(covariate_seq, predict_topic_proportions)

# Convert to data frame for plotting
predicted_df <- data.frame(
  Covariate_Value = covariate_seq,
  t(predicted_proportions)
)

# Reshape the data for ggplot2
predicted_long <- predicted_df %>%
  pivot_longer(cols = -Covariate_Value, names_to = "Topic", values_to = "Proportion")

# Plot the predicted topic proportions
ggplot(predicted_long, aes(x = Covariate_Value, y = Proportion, color = Topic)) +
  geom_line(size = 1) +
  theme_minimal() +
  labs(title = paste("Predicted Topic Proportions Across", covariate_name),
       x = covariate_name,
       y = "Topic Proportion")



















## Real–data application: GSCA–TM with robustness checks vs LDA+ALR and STM

library(gscatm)
library(quanteda)
library(quanteda.textstats)
library(dplyr)
library(ggplot2)
library(tidyr)
library(reshape2)
library(lubridate)
library(readr)
library(Matrix)
library(topicmodels)
library(stm)
library(slam)
library(xtable)

## ------------------------------------------------------------------
## 1. Data import and covariate construction
## ------------------------------------------------------------------

files <- c("data/rawtext_DonaldTrump.tsv",
           "data/rawtext_JoeBiden.tsv",
           "data/rawtext_KamalaHarris.tsv",
           "data/rawtext_MikePence.tsv")

read_tsv_files <- function(file) {
  read_tsv(file, col_types = cols())
}

merged_data <- files %>%
  lapply(read_tsv_files) %>%
  bind_rows()

merged_data <- merged_data %>%
  mutate(
    party = case_when(
      POTUS %in% c("Donald Trump", "Mike Pence") ~ "republican",
      POTUS %in% c("Joe Biden", "Kamala Harris") ~ "democratic",
      TRUE ~ NA_character_
    ),
    candidate_type = case_when(
      (POTUS == "Donald Trump" & year(Date) < 2024) ~ "president",
      (POTUS == "Joe Biden"   & year(Date) < 2024)  ~ "president",
      (POTUS == "Kamala Harris" & year(Date) >= 2024) ~ "president",
      (POTUS == "Kamala Harris" & year(Date) < 2024)  ~ "vice_president",
      POTUS == "Mike Pence" ~ "vice_president",
      TRUE ~ NA_character_
    ),
    year = year(ymd(Date))
  ) %>%
  filter(!is.na(RawText), !is.na(party), !is.na(candidate_type))

write.csv(merged_data, "data/usa_election_speeches_2024.csv", row.names = FALSE)

df <- read.csv("data/usa_election_speeches_2024.csv")

df_sotu <- df %>%
  select(RawText, year, party, candidate_type) %>%
  mutate(
    party          = factor(party),
    candidate_type = factor(candidate_type)
  )

covariates <- df_sotu %>%
  select(year, party, candidate_type)

## ------------------------------------------------------------------
## 2. Text preprocessing and DTM
## ------------------------------------------------------------------

docs <- df_sotu$RawText

tokens <- tokens(
  x             = docs,
  remove_punct  = TRUE,
  remove_symbols = TRUE
)
tokens <- tokens_remove(tokens, stopwords("en"))

dfm <- dfm(tokens)

## remove very frequent tokens to reduce boilerplate
freq <- textstat_frequency(dfm)
most_freq_terms <- freq$feature[1:10]
tokens <- tokens_remove(tokens, pattern = most_freq_terms)

dfm <- dfm(tokens)
dfm <- dfm_trim(dfm, min_termfreq = 5)  # remove very rare terms
dtm <- as(dfm, "dgCMatrix")

## ------------------------------------------------------------------
## 3. GSCA–TM fit (main model)
## ------------------------------------------------------------------

K_seq <- 2:30
search_result <- search_optimal_topics(dtm, covariates, K_seq, n = 3)
print(search_result$top_K)

n_topics <- 9

model_gscatm <- fit_gsca_topic_model(
  dtm          = dtm,
  covariates   = covariates,
  K            = n_topics,
  model        = "zero-inflated",
  max_iter     = 100,
  tol          = 1e-4,
  bootstrap_ci = TRUE,
  n_bootstrap  = 200
)

topic_terms_distribution(model_gscatm, top_n_terms = 10)

## ------------------------------------------------------------------
## 4. Helper functions for robustness checks
## ------------------------------------------------------------------

## Convert dgCMatrix to STM input
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

## ALR regression on estimated theta for any method
estimate_alr_effects <- function(theta, covariates, reference_topic = ncol(theta)) {
  theta <- as.matrix(theta)
  K <- ncol(theta)
  if (reference_topic < 1 || reference_topic > K) {
    stop("reference_topic must be between 1 and K.")
  }

  ## design matrix as in GSCA–TM
  X_design <- model.matrix(~ . - 1, data = covariates)
  X_std <- scale(X_design)
  cov_names <- colnames(X_std)

  ## reorder topics so that reference is last
  topic_idx <- 1:K
  topic_idx <- c(topic_idx[topic_idx != reference_topic], reference_topic)
  theta_ord  <- theta[, topic_idx, drop = FALSE]
  theta_safe <- pmax(theta_ord, 1e-10)

  denom <- theta_safe[, K]
  log_ratios <- log(theta_safe[, 1:(K - 1), drop = FALSE] / denom)

  non_ref_topics <- topic_idx[1:(K - 1)]

  results <- list()

  for (k in 1:(K - 1)) {
    y <- log_ratios[, k]
    glm_data <- data.frame(y = y, X_std)
    weights <- theta_safe[, k] * theta_safe[, K]

    fit <- glm(y ~ ., data = glm_data,
               weights = weights,
               family = gaussian())

    coef_tab <- summary(fit)$coefficients

    ## drop intercept if present
    coef_tab <- coef_tab[setdiff(rownames(coef_tab), "(Intercept)"),
                         , drop = FALSE]

    if (nrow(coef_tab) == 0) next

    for (j in seq_len(nrow(coef_tab))) {
      nm <- rownames(coef_tab)[j]
      beta_hat <- coef_tab[j, "Estimate"]
      se_hat   <- coef_tab[j, "Std. Error"]

      results[[length(results) + 1L]] <- data.frame(
        Covariate        = nm,
        Topic            = paste0("Topic_", non_ref_topics[k]),
        Reference_Topic  = paste0("Topic_", reference_topic),
        Log_Odds_Coefficient = beta_hat,
        Std_Error            = se_hat,
        z_value              = coef_tab[j, "t value"],
        p_value              = coef_tab[j, "Pr(>|t|)"],
        Odds_Ratio           = exp(beta_hat),
        CI_Lower             = exp(beta_hat - 1.96 * se_hat),
        CI_Upper             = exp(beta_hat + 1.96 * se_hat),
        stringsAsFactors     = FALSE
      )
    }
  }

  if (length(results) == 0L) {
    stop("No covariate effects could be estimated.")
  }

  out <- do.call(rbind, results)
  out$adjusted_p_value <- p.adjust(out$p_value, method = "BH")
  rownames(out) <- NULL
  out
}

## ------------------------------------------------------------------
## 5. Covariate effects from GSCA–TM (main analysis)
## ------------------------------------------------------------------

# Use effects already computed inside fit_gsca_topic_model (avoids redundant ALR call)
effects_gscatm <- model_gscatm$effects

## LaTeX table of GSCA–TM effects (simplified)
effects_gscatm_table <- effects_gscatm %>%
  arrange(Covariate, Topic)

print(
  xtable(effects_gscatm_table,
         caption = "GSCA--TM covariate effects (ALR coefficients).",
         label   = "tab:gscatm-effects"),
  include.rownames = FALSE
)

## Plot GSCA–TM effects with plug-in CIs
ggplot(effects_gscatm, aes(x = Covariate, y = Log_Odds_Coefficient)) +
  geom_point() +
  geom_errorbar(aes(ymin = CI_Lower, ymax = CI_Upper),
                width = 0.2) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  facet_wrap(~ Topic, scales = "free_y") +
  theme_minimal() +
  labs(title = "GSCA--TM: covariate effects on topic prevalence (plug-in 95% CI)",
       x = "Covariate",
       y = "Log-odds coefficient")

## Bootstrap vs plug-in CI comparison (only if bootstrap was computed)
if (isTRUE(model_gscatm$bootstrap_ci_computed)) {

  ci_compare_gscatm <- effects_gscatm %>%
    select(Covariate, Topic, Log_Odds_Coefficient,
           CI_Lower, CI_Upper,
           Boot_CI_Lower, Boot_CI_Upper) %>%
    tidyr::pivot_longer(
      cols      = c(CI_Lower, CI_Upper, Boot_CI_Lower, Boot_CI_Upper),
      names_to  = "bound",
      values_to = "value"
    ) %>%
    mutate(
      CI_type = if_else(grepl("^Boot", bound), "Bootstrap", "Plug-in"),
      bound   = if_else(grepl("Lower", bound), "lower", "upper")
    ) %>%
    tidyr::pivot_wider(names_from = bound, values_from = value)

  ggplot(ci_compare_gscatm,
         aes(x = Covariate, y = Log_Odds_Coefficient,
             ymin = lower, ymax = upper,
             colour = CI_type, linetype = CI_type)) +
    geom_point(size = 2, colour = "black") +
    geom_errorbar(width = 0.25, position = position_dodge(width = 0.4)) +
    geom_hline(yintercept = 0, linetype = "dashed", colour = "grey50") +
    scale_colour_manual(values = c("Plug-in"   = "steelblue",
                                   "Bootstrap" = "firebrick")) +
    scale_linetype_manual(values = c("Plug-in"   = "solid",
                                     "Bootstrap" = "dashed")) +
    facet_wrap(~ Topic, scales = "free_y") +
    theme_minimal() +
    labs(title = "GSCA-TM: plug-in vs bootstrap 95% CIs",
         subtitle = "Bootstrap propagates first-stage EM estimation uncertainty",
         x = "Covariate", y = "Log-odds coefficient",
         colour = "CI type", linetype = "CI type") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
}

## ------------------------------------------------------------------
## 6. Robustness check: LDA+ALR and STM
## ------------------------------------------------------------------

## LDA with same number of topics
dtm_triplet <- as.simple_triplet_matrix(dtm)
set.seed(123)
lda_fit <- LDA(dtm_triplet, k = n_topics, control = list(seed = 123))
lda_post <- posterior(lda_fit)
theta_lda <- lda_post$topics

effects_lda <- estimate_alr_effects(
  theta          = theta_lda,
  covariates     = covariates,
  reference_topic = n_topics
)

## STM with prevalence on the same covariates
stm_input <- dtm_to_stm(dtm)
set.seed(123)
stm_fit <- stm(
  documents  = stm_input$documents,
  vocab      = stm_input$vocab,
  K          = n_topics,
  prevalence = ~ party + candidate_type + year,
  data       = covariates,
  max.em.its = 75,
  init.type  = "Spectral",
  verbose    = FALSE
)
theta_stm <- stm_fit$theta

effects_stm <- estimate_alr_effects(
  theta          = theta_stm,
  covariates     = covariates,
  reference_topic = n_topics
)

## ------------------------------------------------------------------
## 7. Comparative summary across methods (robustness)
## ------------------------------------------------------------------

effects_all <- bind_rows(
  mutate(effects_gscatm, method = "GSCA-TM"),
  mutate(effects_lda,    method = "LDA+ALR"),
  mutate(effects_stm,    method = "STM")
)

## Label-invariant summary: distribution of significant effects by covariate
effects_summary <- effects_all %>%
  mutate(
    Direction = case_when(
      adjusted_p_value < 0.05 & Log_Odds_Coefficient > 0 ~ "Positive",
      adjusted_p_value < 0.05 & Log_Odds_Coefficient < 0 ~ "Negative",
      TRUE ~ "Non-significant"
    )
  ) %>%
  group_by(method, Covariate, Direction) %>%
  summarise(n = n(), .groups = "drop")

## LaTeX table for robustness summary
effects_summary_wide <- effects_summary %>%
  pivot_wider(names_from = Direction, values_from = n, values_fill = 0) %>%
  arrange(Covariate, method)

print(
  xtable(effects_summary_wide,
         caption = "Robustness check: number of positive, negative, and non-significant ALR effects by method and covariate.",
         label   = "tab:robustness-summary"),
  include.rownames = FALSE
)

## Comparative plot: distribution of significant effects by method and covariate
ggplot(effects_summary, aes(x = Covariate, y = n, fill = Direction)) +
  geom_bar(stat = "identity", position = "stack") +
  facet_wrap(~ method) +
  theme_minimal() +
  labs(title = "Robustness check: distribution of covariate effects across methods",
       x = "Covariate",
       y = "Number of topic-level effects")

## Optional: direct coefficient comparison by method (subject to topic relabeling)
# For GSCA-TM, show bootstrap CIs when available; plug-in CIs otherwise
effects_all_plot <- effects_all %>%
  mutate(
    ymin = if_else(method == "GSCA-TM" & isTRUE(model_gscatm$bootstrap_ci_computed),
                   Boot_CI_Lower, CI_Lower),
    ymax = if_else(method == "GSCA-TM" & isTRUE(model_gscatm$bootstrap_ci_computed),
                   Boot_CI_Upper, CI_Upper)
  )

ggplot(effects_all_plot, aes(x = Topic, y = Log_Odds_Coefficient,
                              colour = method, group = method)) +
  geom_point(position = position_dodge(width = 0.4)) +
  geom_errorbar(aes(ymin = ymin, ymax = ymax),
                position = position_dodge(width = 0.4), width = 0.2) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  facet_wrap(~ Covariate, scales = "free_y") +
  theme_minimal() +
  labs(title = "Topic-level covariate effects by method (ALR coefficients)",
       subtitle = "GSCA-TM shows bootstrap CIs; LDA+ALR and STM show plug-in CIs",
       x = "Topic",
       y = "Log-odds coefficient")

## ------------------------------------------------------------------
## 8. Example: predicted topic proportions over year (GSCA–TM)
## ------------------------------------------------------------------

covariate_name <- "year"

## Original design matrix and its scaling (used as reference for predictions)
X_orig <- model.matrix(~ . - 1, data = covariates)
X_orig_scaled <- scale(X_orig)
center_vec <- attr(X_orig_scaled, "scaled:center")
scale_vec  <- attr(X_orig_scaled, "scaled:scale")
X_cols     <- colnames(X_orig)

## Sequence of values for the covariate of interest
covariate_seq <- seq(
  from = min(covariates[[covariate_name]], na.rm = TRUE),
  to   = max(covariates[[covariate_name]], na.rm = TRUE),
  length.out = 4
)

## Prediction function
predict_topic_proportions <- function(value) {
  new_covariates <- covariates
  new_covariates[[covariate_name]] <- value

  ## Design matrix for new covariates, aligned with original columns
  X_new <- model.matrix(~ . - 1, data = new_covariates)
  X_new <- X_new[, X_cols, drop = FALSE]

  ## Standardize using the same center/scale as X_orig
  X_new_centered <- sweep(X_new, 2, center_vec, "-")
  X_new_std      <- sweep(X_new_centered, 2, scale_vec, "/")

  ## Reconstruct gamma and theta following the GSCA–TM structure
  gamma_cov <- X_new_std %*% model_gscatm$path_matrix[1:model_gscatm$P, ]
  gamma_new <- cbind(
    gamma_cov,
    model_gscatm$gamma[, (model_gscatm$P + 1):(model_gscatm$P + model_gscatm$K)]
  )
  topic_scores <- gamma_new[, (model_gscatm$P + 1):(model_gscatm$P + model_gscatm$K)]
  theta_new <- exp(topic_scores)
  theta_new <- theta_new / rowSums(theta_new)

  ## Return mean topic proportions across documents
  colMeans(theta_new)
}

## Predict topic proportions over the covariate range
predicted_proportions <- sapply(covariate_seq, predict_topic_proportions)

## Convert to data frame for plotting
predicted_df <- data.frame(
  Covariate_Value = covariate_seq,
  t(predicted_proportions)
)

predicted_long <- predicted_df %>%
  pivot_longer(cols = -Covariate_Value, names_to = "Topic", values_to = "Proportion")

ggplot(predicted_long, aes(x = Covariate_Value, y = Proportion, color = Topic)) +
  geom_line(size = 1) +
  theme_minimal() +
  labs(title = paste("GSCA--TM: predicted topic proportions across", covariate_name),
       x = covariate_name,
       y = "Mean topic proportion")




