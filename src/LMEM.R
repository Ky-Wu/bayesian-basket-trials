library(partitions)

#' LMEM Partition Posterior Probability Calculations
#'
#' @param n_b Vector of basket sample sizes.
#' @param y Vector of basket responses.
#' @param d1 Design prior hyperparameter controlling propensity towards borrowing.
#' @param d2 Analysis prior hyperparameter controlling degree of borrowing.
#' @param a Hyperparameter of beta prior on basket response rate.
#' @param b Hyperparameter of beta prior on basket response rate.
#'
#' @return
#' @export
#'
#' @examples
gridSearchLMEMPartition <- function(n_b, y, d1 = 0, d2 = 2, a = 1, b = 1) {
  stopifnot(length(y) <= 10, length(y) == length(n_b))
  parts <- listParts(length(y))
  log_weights <- vapply(parts, function(part) {
    ll <- sum(vapply(part, function(x) {
      ys <- y[x]
      ns <- n_b[x]
      basket_weight(sum(ns), sum(ys), a, b)
    }, numeric(1)))
    ll
  }, numeric(1))
  K <- vapply(parts, length, numeric(1))
  search_weights <- logsumexp_weights(log_weights + d1 * log(K))
  prob_weights <- logsumexp_weights(log_weights + d2 * log(K))
  chosen_partition <- sample(which(search_weights == max(search_weights)), 1)
  list(part = parts[[chosen_partition]],
       post_prob = prob_weights[chosen_partition] / sum(prob_weights))
}

#' LMEM Posterior Analysis
#'
#' @param n_b Vector of basket sample sizes.
#' @param y Vector of basket observed responses.
#' @param p0 Historical/control response rate.
#' @param a Hyperparameter of beta priors on basket response rates.
#' @param b Hyperparameter of beta priors on basket response rates.
#' @param d1 Design prior hyperparameter controlling propensity towards borrowing.
#' @param d2 Analysis prior hyperparameter controlling degree of borrowing.
#'
#' @return Vector of posterior efficacy probabilities under the LMEM framework.
#' @export
#'
#' @examples
LMEMBasketEfficacy <- function(n_b, y, p0, a = 1, b = 1, d1 = 0, d2 = 0) {
  res <- gridSearchLMEMPartition(n_b, y, d1 = d1, d2 = d2)
  part <- res$part
  pp <- res$post_prob
  out <- numeric(length(y))
  for (x in part) {
    ys <- y[x]
    ns <- n_b[x]
    Pb <- vapply(seq_along(ys), function(i) {
      a_new <- a + ys[i] + pp * sum(ys[-i])
      b_new <- b + ns[i] - ys[i] + pp * sum(ns[-i] - ys[-i])
      pbeta(p0, a_new, b_new, lower.tail = FALSE)
    }, numeric(1))
    out[x] <- Pb
  }
  out
}

partitionFormat <- function(part) {
  parts <- vapply(part, function(x) paste0("(", paste0(x, collapse = ","), ")"), character(1))
  paste0(parts, collapse = "")
}

evaluateLMEMPriorPart <- function(n_b, y, delta = 0) {
  parts <- listParts(length(y))
  parts_format <- vapply(parts, partitionFormat, character(1))
  K <- vapply(parts, length, numeric(1))
  log_weights <- vapply(parts, function(part) {
    temp <- vapply(part, function(x) {
      -lchoose(sum(n_b[x]), sum(y[x]))
    }, numeric(1))
    sum(temp)
  }, numeric(1))
  log_weights <- log_weights + delta * log(K)
  pd <- logsumexp_weights(log_weights)
  data.frame(part = parts_format,
             prior_density = pd,
             prior_logweights = log_weights)
}

evaluateLMEMPostPart <- function(n_b, y, d1 = 0, d2 = 0, a = 1, b = 1) {
  parts <- listParts(length(y))
  search_prior <- evaluateLMEMPriorPart(n_b, y, delta = d1)
  pool_prior <- evaluateLMEMPriorPart(n_b, y, delta = d2)
  lls <- vapply(parts, function(part) {
    ll <- sum(vapply(part, function(x) {
      ys <- y[x]
      ns <- n_b[x]
      basket_weight(sum(ns), sum(ys), a, b)
    }, numeric(1)))
    ll
  }, numeric(1))
  search_post <- logsumexp_weights(search_prior$prior_logweights + lls)
  pool_post <- logsumexp_weights(pool_prior$prior_logweights + lls)
  data.frame(part = search_prior$part,
             search_prior = search_prior$prior_density,
             search_post = search_post,
             pool_prior = pool_prior$prior_density,
             pool_post = pool_post)
}

#' Log-sum-exp trick for normalizing small weights
#'
#' @param x Vector of unnormalized weights.
#'
#' @return Vector of normalized weights.
#' @export
#'
#' @examples
logsumexp_weights <- function(x) {
  xmax <- max(x)
  d <- xmax + log(sum(exp(x - xmax)))
  exp(x - d)
}

#' Unnormalized Marginal Density
#'
#' @param n Vector of sample sizes.
#' @param yi Vector of observed responses.
#' @param a Hyperparameter of beta priors on basket response rates.
#' @param b Hyperparameter of beta priors on basket response rates.
#'
#' @return Vector of unnormalized marginal densities.
#' @export
#'
#' @examples
basket_weight <- function(n, yi, a, b) {
  a_new <- a + yi
  b_new <- b + n - yi
  lbeta(a_new, b_new) - lbeta(a, b)
}


