# GSCA-TM: Generalized Structured Component Analysis for Topic Modeling

## Overview

The `gscatm` package provides a novel statistical framework for topic modeling that integrates Generalized Structured Component Analysis (GSCA) with probabilistic topic models. This methodology enables the estimation of relationships between document-level covariates and latent topic distributions in a unified manner. The package offers flexibility in topic proportion distributions, supporting logistic-normal, Dirichlet, and zero-inflated specifications.

The framework is designed to bridge exploratory and confirmatory text analysis, allowing hypothesis-driven research on thematic variations in textual data. Covariate effects on topic proportions are estimated via **Additive Log-Ratio Weighted Least Squares (ALR-WLS)** regression, with optional **parametric bootstrap** confidence intervals for uncertainty propagation.

## Repository Structure

```
gscatm/
├── R/                        # Package source code
├── man/                      # Documentation (Rd files)
├── data/                     # Package datasets
├── demo/                     # Package demo scripts
├── tests/                    # Unit tests (testthat)
├── vignettes/                # Package vignettes
├── replication/
│   ├── simulation/           # Simulation study (Section 5)
│   └── real_application/     # Real-world applications (Section 6)
├── DESCRIPTION
├── NAMESPACE
├── .Rbuildignore
├── README.md
└── gscatm_0.1.0.tar.gz      # Pre-built package tarball
```

## Installation

```r
# Install from the pre-built tarball
install.packages("gscatm_0.1.0.tar.gz", repos = NULL, type = "source")

# Or install from source
devtools::install(".")

# Or install from GitHub
devtools::install_github("marcoortu/gscatm")
```

**Required packages** (installed automatically): `Matrix`, `MCMCpack`, `ggplot2`, `dplyr`, `plotly`

**Optional**: `clue` for optimal Hungarian topic alignment in bootstrap (a greedy fallback is used if not available)

## Features

- **GSCA-Enhanced Topic Modeling**: Integrates GSCA with topic modeling to incorporate document-level covariates.
- **Flexible Topic Proportion Models**: Supports logistic-normal, Dirichlet, and zero-inflated distributions.
- **Efficient Estimation Procedure**: Combines Expectation-Maximization (EM) with GSCA optimization and decaying ridge regularization.
- **Covariate Effect Estimation**: ALR-WLS regression with plug-in and parametric bootstrap confidence intervals.
- **Parametric Bootstrap**: Uncertainty-propagated SEs and CIs via multinomial resampling and L1 Hungarian topic alignment.
- **Interpretability**: Outputs a structured path coefficient matrix linking covariates to topic proportions.
- **Visualization Tools**: Static (ggplot2) and interactive (plotly) topic distribution plots; top-term bar charts.

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

## Quick Start

```r
library(gscatm)

# Generate synthetic data
sample_data <- generate_sample_data(n_docs = 100, n_terms = 500,
                                    n_topics = 5, n_covariates = 3)

# Fit GSCA Topic Model
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

# Parametric bootstrap CIs
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
```

## Replication Instructions

The `replication/` folder contains all scripts and data needed to reproduce the results presented in the paper. Before running the replication scripts, install `gscatm` and the additional dependencies listed below.

### Prerequisites

```r
# Install gscatm
install.packages("gscatm_0.1.0.tar.gz", repos = NULL, type = "source")

# Additional packages required by the replication scripts
install.packages(c(
  "stm",           # Structural Topic Model (benchmark comparison)
  "topicmodels",   # LDA via topicmodels (benchmark comparison)
  "slam",          # Sparse matrix conversion for topicmodels
  "clue",          # Hungarian algorithm for topic alignment
  "xtable",        # LaTeX table generation
  "tidyr",         # Data reshaping
  "quanteda",      # Text preprocessing (real-world applications)
  "quanteda.textstats",
  "reshape2",
  "lubridate",
  "readr"
))
```

### Simulation Study

The simulation study compares GSCA-TM against STM and LDA+ALR across three conditions (baseline, high sparsity, high covariate influence), evaluating RMSE for topic proportions, word distributions, and covariate effects, as well as bootstrap coverage.

```r
# Set working directory to the repository root
setwd("path/to/gscatm")

# Run the full simulation benchmark
source("replication/simulation/gscatm_simulation_benchmark.R")
```

The script generates:
- LaTeX summary tables with RMSE, coverage, and perplexity metrics
- Boxplots comparing methods across conditions

Key parameters (adjustable at the top of the script):
- `n_reps`: number of Monte Carlo replications (default: 10)
- `n_bootstrap_sim`: number of bootstrap samples per replication (default: 10)
- `n_docs`, `n_terms`, `K`, `P`: simulation dimensions

A lighter quick-simulation script is also provided:

```r
source("replication/simulation/run_quick_sim.R")
```

### Real-World Applications

#### US Presidential Inaugural Speeches

Applies GSCA-TM to the inaugural speeches corpus from the `quanteda` package, with year and party as covariates.

```r
source("replication/real_application/gscatm_usa_president_speech_example.R")
```

#### US 2024 Election Speeches

Applies GSCA-TM to campaign speeches from the 2024 US presidential election (Trump, Biden, Harris, Pence). The dataset is included in `replication/real_application/usa_election_speeches_2024.csv`.

```r
source("replication/real_application/gscatm_usa_election_speeches_2024_example.R")
```

## License

GPL-3 © 2024 Marco Ortu
