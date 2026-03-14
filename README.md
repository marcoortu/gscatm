# Generalized Structured Component Analysis for Topic Modeling

## Overview

The `gscatm` package provides a novel statistical framework for topic modeling that integrates Generalized Structured Component Analysis (GSCA) with probabilistic topic models. This methodology enables the estimation of relationships between document-level covariates and latent topic distributions in a unified manner. The package offers flexibility in topic proportion distributions, supporting logistic-normal, Dirichlet, and zero-inflated specifications.

The framework is designed to bridge exploratory and confirmatory text analysis, allowing hypothesis-driven research on thematic variations in textual data. Covariate effects on topic proportions are estimated via **Additive Log-Ratio Weighted Least Squares (ALR-WLS)** regression, with optional **parametric bootstrap** confidence intervals for uncertainty propagation.

## Features

- **GSCA-Enhanced Topic Modeling**: Integrates GSCA with topic modeling to incorporate document-level covariates.
- **Flexible Topic Proportion Models**: Supports logistic-normal, Dirichlet, and zero-inflated distributions.
- **Efficient Estimation Procedure**: Combines Expectation-Maximization (EM) with GSCA optimization and decaying ridge regularization.
- **Covariate Effect Estimation**: ALR-WLS regression with plug-in and parametric bootstrap confidence intervals.
- **Parametric Bootstrap**: Uncertainty-propagated SEs and CIs via multinomial resampling and L1 Hungarian topic alignment (Section 4.2).
- **Interpretability**: Outputs a structured path coefficient matrix linking covariates to topic proportions.
- **Visualization Tools**: Static (ggplot2) and interactive (plotly) topic distribution plots; top-term bar charts.
- **Simulation and Real-World Applications**: Supports synthetic data generation and empirical applications.

## Installation

```r
# From local source
devtools::install("path/to/gscatm")

# From GitHub (when available)
# devtools::install_github("marcoortu/gscatm")
```

**Required packages** (installed automatically): `Matrix`, `MCMCpack`, `ggplot2`, `dplyr`, `plotly`

**Optional**: `clue` for optimal Hungarian topic alignment in bootstrap (a greedy fallback is used if not available)

## Functions

| Function | Description |
|---|---|
| `fit_gsca_topic_model()` | Fit the GSCA-TM and estimate covariate effects |
| `predict_topics()` | Predict topic proportions for new documents |
| `generate_sample_data()` | Generate synthetic data for testing |
| `plot_topic_distribution()` | Static ggplot2 plot of topic prevalence or content |
| `plot_topic_distribution_interactive()` | Interactive plotly version of topic distribution plot |
| `print_topic_summary()` | Print top terms and summary for each topic |
| `topic_terms_distribution()` | Bar chart of top terms per topic |
| `search_optimal_topics()` | Search for optimal number of topics K |

## Usage

```r
library(gscatm)

# Generate synthetic data
sample_data <- generate_sample_data(n_docs = 100, n_terms = 500,
                                    n_topics = 5, n_covariates = 3)

# Fit GSCA Topic Model (plug-in CIs)
model <- fit_gsca_topic_model(
  dtm        = sample_data$dtm,
  covariates = sample_data$covariates,
  K          = 5,
  model      = "logistic-normal"
)

# Covariate effects on topic proportions
model$effects[, c("Covariate", "Topic", "Odds_Ratio", "CI_Lower", "CI_Upper",
                  "adjusted_p_value")]

# Top terms per topic
apply(model$phi, 1, function(x) names(sort(x, decreasing = TRUE))[1:5])

# Parametric bootstrap CIs (Section 4.2)
model_boot <- fit_gsca_topic_model(
  dtm          = sample_data$dtm,
  covariates   = sample_data$covariates,
  K            = 5,
  bootstrap_ci = TRUE,
  n_bootstrap  = 200
)

# Compare plug-in vs bootstrap CIs
model_boot$effects[, c("Covariate", "Topic", "Odds_Ratio",
                       "CI_Lower", "CI_Upper",
                       "Boot_CI_Lower", "Boot_CI_Upper")]

# Predict for new documents
new_data <- generate_sample_data(n_docs = 10, n_terms = 500,
                                 n_topics = 5, n_covariates = 3)
theta_new <- predict_topics(model, new_data$dtm, new_data$covariates)

# Plot topic distribution
plot_topic_distribution(model, type = "prevalence")

# Interactive topic visualization
plot_topic_distribution_interactive(model, type = "prevalence")

# Print topic summaries
print_topic_summary(model, n_words = 10)

# Search for optimal K
optimal_topics <- search_optimal_topics(sample_data$dtm, sample_data$covariates,
                                        K_seq = 2:10)
```

## License

GPL-3 © 2024 Marco Ortu
