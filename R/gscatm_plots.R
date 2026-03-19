# R/plot_topic_distribution.R

#' Plot Topic Distribution
#'
#' Creates a horizontal bar plot of topic distributions with topic labels
#'
#' @param model Fitted GSCA topic model object
#' @param type Character string specifying the type of plot ("prevalence" or "content")
#' @param top_n Integer specifying number of top words to show in labels (default: 3)
#' @param min_prob Minimum probability threshold for filtering topics (default: 0.01)
#' @param show_numbers Logical; whether to show topic numbers in labels (default: TRUE)
#' @return A ggplot2 object
#' @importFrom ggplot2 ggplot aes geom_col theme_minimal theme element_text element_blank labs
#' @importFrom dplyr %>%
#' @export
#'
#' @examples
#' \dontrun{
#' model <- fit_gsca_topic_model(dtm, covariates, K = 10)
#' plot_topic_distribution(model)
#' }
plot_topic_distribution <- function(model,
                                    type = c("prevalence", "content"),
                                    top_n = 3,
                                    min_prob = 0.01,
                                    show_numbers = TRUE) {

  type <- match.arg(type)

  # Extract topic-word distributions and create labels
  phi <- model$phi
  vocab <- model$vocab

  # Get top words for each topic
  top_words <- apply(phi, 1, function(x) {
    words <- vocab[order(x, decreasing = TRUE)][1:top_n]
    paste(words, collapse = ", ")
  })

  if (type == "prevalence") {
    # Calculate average topic proportions
    topic_props <- colMeans(model$theta)

    # Create data frame for plotting
    plot_data <- data.frame(
      topic = 1:length(topic_props),
      proportion = topic_props,
      label = paste("Topic", 1:length(topic_props), ":", top_words)
    )

  } else {
    # For content type, use word distribution entropy
    topic_props <- apply(phi, 1, function(x) -sum(x * log(x)))
    plot_data <- data.frame(
      topic = 1:length(topic_props),
      proportion = topic_props,
      label = paste("Topic", 1:length(topic_props), ":", top_words)
    )
  }

  # Filter by minimum probability if needed
  plot_data <- plot_data[plot_data$proportion >= min_prob, ]

  # Sort by proportion
  plot_data <- plot_data[order(plot_data$proportion), ]

  # Create labels
  if (!show_numbers) {
    plot_data$label <- top_words
  }

  # Convert labels to factors to preserve order
  plot_data$label <- factor(plot_data$label, levels = plot_data$label)

  # Create the plot
  p <- ggplot(plot_data, aes(x = proportion, y = label)) +
    geom_col(fill = "steelblue", alpha = 0.7) +
    theme_minimal() +
    theme(
      axis.text.y = element_text(hjust = 1),
      panel.grid.major.y = element_blank(),
      panel.grid.minor.y = element_blank()
    ) +
    labs(
      x = "Expected Topic Proportions",
      y = NULL,
      title = if(type == "prevalence") "Topic Prevalence" else "Topic Content Distribution"
    )

  return(p)
}

#' Create Interactive Topic Distribution Plot
#'
#' Creates an interactive visualization of topic distributions using plotly
#'
#' @param model Fitted GSCA topic model object
#' @param type Character string specifying the type of plot ("prevalence" or "content")
#' @param top_n Integer specifying number of top words to show in labels
#' @param min_prob Minimum probability threshold for filtering topics
#' @return A plotly object
#' @import plotly
#' @export
plot_topic_distribution_interactive <- function(model,
                                                type = c("prevalence", "content"),
                                                top_n = 3,
                                                min_prob = 0.01) {

  # Get the ggplot object
  p <- plot_topic_distribution(model, type, top_n, min_prob)

  # Convert to plotly
  p_interactive <- ggplotly(p, tooltip = c("x", "y")) %>%
    layout(margin = list(l = 150)) # Adjust left margin for labels

  return(p_interactive)
}

#' Print Topic Summary
#'
#' Prints a formatted summary of topics with their top words and proportions
#'
#' @param model Fitted GSCA topic model object
#' @param n_words Number of top words to show per topic
#' @param n_topics Number of topics to show (NULL for all)
#' @export
print_topic_summary <- function(model, n_words = 10, n_topics = NULL) {
  # Extract model components
  phi <- model$phi
  vocab <- model$vocab
  topic_props <- colMeans(model$theta)

  # Determine number of topics to show
  K <- nrow(phi)
  if (is.null(n_topics)) n_topics <- K

  # Get top words and their probabilities for each topic
  topic_summaries <- lapply(1:K, function(k) {
    word_probs <- sort(phi[k,], decreasing = TRUE)
    top_indices <- head(order(phi[k,], decreasing = TRUE), n_words)
    data.frame(
      topic = k,
      proportion = topic_props[k],
      words = paste(vocab[top_indices], collapse = ", "),
      stringsAsFactors = FALSE
    )
  })

  # Combine into single data frame
  summary_df <- do.call(rbind, topic_summaries)

  # Sort by proportion and select top topics
  summary_df <- summary_df[order(summary_df$proportion, decreasing = TRUE), ]
  summary_df <- head(summary_df, n_topics)

  # Print formatted summary
  cat("\nTopic Model Summary:\n")
  cat("Number of topics:", K, "\n\n")

  for(i in 1:nrow(summary_df)) {
    cat(sprintf("Topic %d (%.1f%%):\n",
                summary_df$topic[i],
                summary_df$proportion[i] * 100))
    cat(summary_df$words[i], "\n\n")
  }
}


#' Display Top Terms for Each Topic
#'
#' Prints the most probable terms for each topic in a fitted GSCA topic model,
#' along with their corresponding probabilities. This function helps interpret
#' the semantic meaning of each discovered topic.
#'
#' @param gsca_model A fitted GSCA topic model object (class 'gsca_topic_model')
#'        containing topic-term distributions and vocabulary
#' @param top_n_terms Integer specifying how many top terms to display for each topic.
#'        Default is 5
#'
#' @details
#' For each topic in the model, this function:
#' 1. Extracts the probability distribution over terms
#' 2. Identifies the terms with highest probability
#' 3. Prints these terms along with their probabilities
#'
#' The output is formatted as:
#' Topic k:
#'   term1 (probability1)
#'   term2 (probability2)
#'   ...
#'
#' @return No return value, called for side effects (printing)
#'
#' @examples
#' # Assuming you have a fitted GSCA topic model:
#' topic_terms_distribution(model)
#'
#' # Display top 10 terms instead of default 5:
#' topic_terms_distribution(model, top_n_terms = 10)
#'
#' @seealso \code{\link{fit_gsca_topic_model}} for fitting the GSCA topic model
#'
#' @export
topic_terms_distribution <- function(gsca_model, top_n_terms=5){

  # Assuming 'gsca_model' is your fitted model object
  phi <- gsca_model$phi      # Topic-term distributions (K x V matrix)
  vocab <- gsca_model$vocab  # Vocabulary (vector of terms)
  # For each topic, get top 10 terms
  for (k in 1:gsca_model$K) {
    # Get the probabilities for the k-th topic
    topic_probs <- phi[k, ]

    # Get indices of top 10 terms
    top_term_indices <- order(topic_probs, decreasing = TRUE)[1:top_n_terms]

    # Get the corresponding terms
    top_terms <- vocab[top_term_indices]

    # Get the probabilities of the top terms
    top_probs <- topic_probs[top_term_indices]

    # Print the topic number and top terms with their probabilities
    cat(sprintf("Topic %d:\n", k))
    for (i in 1:length(top_terms)) {
      cat(sprintf("  %s (%.4f)\n", top_terms[i], top_probs[i]))
    }
    cat("\n")
  }

}
