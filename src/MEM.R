library(basket)

MEMEfficacy <- function(n_b, y, p0, prior_prob = 0.1) {
  res <- basket::basket(y, n_b, name = letters[seq_along(n_b)], p0 = p0,
                        method = "exact",
                        prior = diag(length(n_b)) * (1 - prior_prob) +
                          matrix(prior_prob, nrow = length(n_b), ncol = length(n_b)))
  res$basket$post_prob
}
