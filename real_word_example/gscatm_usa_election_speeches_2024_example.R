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
files <- c("rawtext_DonaldTrump.tsv",
           "rawtext_JoeBiden.tsv",
           "rawtext_KamalaHarris.tsv",
           "rawtext_MikePence.tsv")

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
  dtm = dtm,
  covariates = covariates,
  K = n_topics,
  model = "zero-inflated", #"logistic-normal", "dirichlet", "zero-inflated"
  max_iter = 100,
  tol = 1e-4
)

topic_terms_distribution(model, top_n_terms = 10)




# Print the summary table
print(model$effects)

print(model$fit)

# Add confidence intervals to the summary table


# Plotting
ggplot(model$effects, aes(x = Covariate, y = Log_Odds_Coefficient, color = Topic)) +
  geom_point(position = position_dodge(width = 0.5)) +
  geom_errorbar(aes(ymin = min(model$effects$CI_Lower), ymax = max(model$effects$CI_Upper)), width = 0.2,
                position = position_dodge(width = 0.5)) +
  facet_wrap(~ Topic, scales = "free_y") +
  theme_bw() +
  labs(title = "Effects of Covariates on Topics",
       y = "Odds_Ratio Estimate",
       x = "Covariate") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Prepare data for heatmap
heatmap_data <- dcast(model$effects, Log_Odds_Coefficient ~ Topic, value.var = "Log_Odds_Coefficient")

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


