library(gscatm)
library(quanteda)
library(quanteda.textstats)
library(dplyr)
library(ggplot2)
library(dplyr)
library(tidyr)
library(reshape2)

# Load the State of the Union corpus
data("data_corpus_inaugural")


# Create a sample dataframe with texts and covariates
df_sotu <- data.frame(
  text = as.character(data_corpus_inaugural),
  year = docvars(data_corpus_inaugural, "Year"),
  party = docvars(data_corpus_inaugural, "Party"),
  stringsAsFactors = FALSE
)

# Sample if needed (optional)
# df_sample <- sample_n(df_sotu, size = 100)

# Select covariates
covariates <- df_sotu %>%
  select(c("year", "party"))

# Process the documents
docs <- df_sotu$text
tokens <- tokens(docs, remove_punct = TRUE)
tokens <- tokens_remove(tokens, stopwords("en"))
# tokens <- tokens_wordstem(tokens)

# Create the DFM
dfm <- dfm(tokens)

# Identify the most frequent tokens (e.g., top 10)
freq <- textstat_frequency(dfm)
most_freq_terms <- freq$feature[1:25]  # Replace 10 with the number of tokens you want to remove

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
n_topics = 5

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

# Bootstrap CI vs plug-in CI comparison
# Only shown if bootstrap was computed (bootstrap_ci = TRUE above)
if (isTRUE(model$bootstrap_ci_computed)) {

  # Stack plug-in and bootstrap intervals into long format for overlay
  ci_compare <- model$effects %>%
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

  ggplot(ci_compare,
         aes(x = Covariate, y = Log_Odds_Coefficient,
             ymin = lower, ymax = upper,
             colour = CI_type, linetype = CI_type)) +
    geom_point(size = 2, colour = "black") +
    geom_errorbar(width = 0.25, position = position_dodge(width = 0.4)) +
    geom_hline(yintercept = 0, linetype = "dashed", colour = "grey50") +
    scale_colour_manual(values = c("Plug-in" = "steelblue",
                                   "Bootstrap" = "firebrick")) +
    scale_linetype_manual(values = c("Plug-in" = "solid",
                                     "Bootstrap" = "dashed")) +
    facet_wrap(~ Topic, scales = "free_y") +
    theme_minimal() +
    labs(title = "Plug-in vs Bootstrap 95% CIs for covariate effects",
         subtitle = "Bootstrap propagates first-stage EM estimation uncertainty",
         x = "Covariate", y = "Log-odds coefficient",
         colour = "CI type", linetype = "CI type") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
}



# Select a continuous covariate (e.g., "Covariate1")
covariate_name <- "year"

# Create a sequence of values for the covariate
covariate_seq <- seq(min(covariates[, covariate_name]),
                     max(covariates[, covariate_name]), length.out = 10)

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
