﻿# Generalized Structured Component Analysis for Topic Modeling

## Overview

The `gscatm` package provides a novel statistical framework for topic modeling that integrates Generalized Structured Component Analysis (GSCA) with probabilistic topic models. This methodology enables the estimation of relationships between document-level covariates and latent topic distributions in a unified manner. The package offers flexibility in topic proportion distributions, supporting logistic-normal, Dirichlet, and zero-inflated specifications.

The framework is designed to bridge exploratory and confirmatory text analysis, allowing hypothesis-driven research on thematic variations in textual data. The package includes functions for model estimation, topic prediction, visualization, and evaluation.

## Features

- **GSCA-Enhanced Topic Modeling**: Integrates GSCA with topic modeling to incorporate document-level covariates.
- **Flexible Topic Proportion Models**: Supports logistic-normal, Dirichlet, and zero-inflated distributions.
- **Efficient Estimation Procedure**: Combines Expectation-Maximization (EM) with GSCA optimization.
- **Interpretability**: Outputs a structured path coefficient matrix linking covariates to topic proportions.
- **Visualization Tools**: Functions for topic distribution plots and interactive visualizations.
- **Simulation and Real-World Applications**: Supports synthetic data generation and empirical applications.

## Installation

To install the package, use:

```r
# Install from GitHub (if applicable)
devtools::install_github("https://github.com/marcoortu/gscatm")

# Or load locally
install.packages("gscatm")


library(gscatm)

# Generate synthetic data
sample_data <- generate_sample_data(n_docs = 100, n_terms = 500, n_topics = 5, n_covariates = 3)


# Fit GSCA Topic Model
model <- fit_gsca_topic_model(
  dtm = sample_data$dtm,
  covariates = sample_data$covariates,
  K = 5,
  model = "logistic-normal"
)

new_topic_probs <- predict_topics(model, new_dtm, new_covariates)


# Plot topic distribution
plot_topic_distribution(model, type = "prevalence")

# Interactive topic visualization
plot_topic_distribution_interactive(model, type = "prevalence")


optimal_topics <- search_optimal_topics(sample_data$dtm, sample_data$covariates, K_seq = 2:10)

```
