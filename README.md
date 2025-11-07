# Generalized Structured Component Analysis for Topic Modeling

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

## Open Source License

This project is distributed under the MIT License, a permissive open source license that allows reuse with minimal restrictions
while requiring attribution and preservation of the license notice.

### Common Open Source Licenses

Below are several of the most widely adopted open source licenses:

- **MIT License** – Simple and permissive, commonly used for libraries and applications.
- **Apache License 2.0** – Permissive license with explicit patent grants and contribution terms.
- **GNU General Public License v3.0 (GPLv3)** – Strong copyleft license ensuring derivative works remain open source.
- **BSD 3-Clause License** – Permissive license similar to MIT with an explicit non-endorsement clause.
- **Mozilla Public License 2.0 (MPL 2.0)** – File-level copyleft license that balances openness with flexibility.

### MIT License

MIT License

Copyright (c) 2024 Marco Ortu

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
