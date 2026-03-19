# R/generate_sample_data.R
#' Generate Sample Data for GSCA Topic Model
#'
#' Creates a synthetic dataset with known topic structure and covariate relationships
#' for demonstration and testing purposes.
#'
#' @param n_docs Number of documents (default: 10)
#' @param n_terms Number of unique terms (default: 50)
#' @param n_topics Number of topics (default: 3)
#' @param n_covariates Number of covariates (default: 3)
#' @param doc_length Average document length (default: 100)
#' @param seed Random seed for reproducibility (default: 42)
#' @return A list containing:
#'   \item{dtm}{Document-term matrix (sparse dgCMatrix)}
#'   \item{covariates}{Matrix of document covariates}
#'   \item{true_topics}{True topic-word distributions}
#'   \item{true_theta}{True document-topic proportions}
#' @export
#' @examples
#' sample_data <- generate_sample_data()
#' model <- fit_gsca_topic_model(sample_data$dtm, sample_data$covariates, K = 3)
generate_sample_data <- function(n_docs = 10, n_terms = 50, n_topics = 3,
                                 n_covariates = 3, doc_length = 100, seed = 42) {
  set.seed(seed)

  # Generate true topic-word distributions
  true_topics <- matrix(rgamma(n_topics * n_terms, shape = 1, rate = 1),
                        nrow = n_topics)
  true_topics <- sweep(true_topics, 1, rowSums(true_topics), "/")

  # Generate covariates with some correlation structure
  covariates <- matrix(rnorm(n_docs * n_covariates), nrow = n_docs)
  colnames(covariates) <- paste0("cov_", 1:n_covariates)

  # Generate topic proportions influenced by covariates
  beta <- matrix(rnorm(n_covariates * n_topics, sd = 0.75),
                 nrow = n_covariates)
  logits <- covariates %*% beta
  true_theta <- exp(logits)
  true_theta <- sweep(true_theta, 1, rowSums(true_theta), "/")

  # Generate documents
  dtm_ <- Matrix::Matrix(0, nrow = n_docs, ncol = n_terms, sparse = TRUE)
  colnames(dtm_) <- paste0("term_", 1:n_terms)

  for(d in 1:n_docs) {
    doc_length_d <- rpois(1, doc_length)
    topics <- sample(1:n_topics, doc_length_d, prob = true_theta[d,],
                     replace = TRUE)
    for(w in 1:doc_length_d) {
      term <- sample(1:n_terms, 1, prob = true_topics[topics[w],])
      dtm_[d, term] <- dtm_[d, term] + 1
    }
  }
  dtm_ <- as(dtm_, "dgCMatrix")
  return(list(
    dtm = dtm_,
    covariates = as.data.frame(covariates),
    true_topics = true_topics,
    true_theta = true_theta
  ))
}
