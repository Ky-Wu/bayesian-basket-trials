library(basket)

#' MEM Posterior Analysis
#'
#' @param n_b Vector of basket sample sizes.
#' @param y Vector of basket observed responses.
#' @param p0 Historical/control response rate.
#' @param prior_prob Prior exchangeability pair inclusion probability, default = 0.1.
#'
#' @return Vector of posterior efficacy probabilities under the MEM framework.
#' @export
#'
#' @examples
MEMEfficacy <- function(n_b, y, p0, prior_prob = 0.1) {
  res <- basket::basket(y, n_b, name = letters[seq_along(n_b)], p0 = p0,
                        method = "exact",
                        prior = diag(length(n_b)) * (1 - prior_prob) +
                          matrix(prior_prob, nrow = length(n_b), ncol = length(n_b)))
  res$basket$post_prob
}
