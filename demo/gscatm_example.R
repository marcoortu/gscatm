# Example usage of GSCA Topic Model
library(gscatm)
# Load required packages
library(ggplot2)
library(dplyr)
library(plotly)

# Generate sample data
set.seed(42)
sample_data <- generate_sample_data(
  n_docs = 30,
  n_terms = 150,
  n_topics = 3,
  n_covariates = 4,
  seed=23
)

# Fit the model
model <- fit_gsca_topic_model(
  dtm = sample_data$dtm,
  covariates = sample_data$covariates,
  K = 3,
  model = "logistic-normal"
)

# Examine topic-word distributions
top_terms <- apply(model$phi, 1, function(x) {
  names(sort(x, decreasing = TRUE))[1:10]
})
print("Top terms per topic:")
print(top_terms)

# Compare estimated vs true topic proportions
plot(sample_data$true_theta[,1], model$theta[,1],
     xlab = "True proportions", ylab = "Estimated proportions",
     main = "Topic 1 Proportions")

# Basic topic distribution plot
plot_topic_distribution(model)

# Interactive plot with custom settings
plot_topic_distribution_interactive(model, top_n = 3, min_prob = 0.02)

# Print summary of top 5 topics with 10 words each
print_topic_summary(model, n_words = 5, n_topics = 3)
