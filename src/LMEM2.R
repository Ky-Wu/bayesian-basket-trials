library(partitions)

#' Basket Posterior KL Divergences
#'
#' @param ni Target basket sample size.
#' @param yi Target basket observed responses.
#' @param ns Vector of external basket sample sizes.
#' @param ys Vector of external basket observed responses.
#' @param a Hyperparameter of beta prior on basket response rate
#' @param b Hyperparameter of beta prior on basket response rate
#'
#' @return Computed KL divergences of response rate posterior
#' distribution in target basket from that of external baskets.
#' @export
#'
computePosteriorKLDivergences <- function(ni, yi, ns, ys, a = 1, b = 1) {
  a1 <- a + yi
  b1 <- b + ni - yi
  vapply(seq_along(ns),
         function(i) {
           a2 <- a + ys[i]
           b2 <- b + ns[i] - ys[i]
           t1 <- lbeta(a1, b1) - lbeta(a2, b2)
           t2 <- (a2 - a1) * (digamma(a2) - digamma(a2 + b2))
           t3 <- (b2 - b1) * (digamma(b2) - digamma(a2 + b2))
           t1 + t2 + t3
         }, numeric(1))
}

#' LMEM2 Partition Posterior Probability Calculations
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
gridSearchLMEM2Partition <- function(n_b, y, d1 = 0, d2 = 2, a = 1, b = 1) {

  stopifnot(length(y) <= 10, length(y) == length(n_b))
  parts <- listParts(length(y))
  log_weights <- vapply(parts, function(part) {
    ll <- sum(vapply(part, function(x) {
      ys <- y[x]
      ns <- n_b[x]
      #log_dbinombeta(sum(ns), sum(ys), a, b)
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


#' LMEM2 Posterior Analysis
#'
#' @param n_b Vector of basket sample sizes.
#' @param y Vector of basket observed responses.
#' @param p0 Historical/control response rate.
#' @param a Hyperparameter of beta priors on basket response rates.
#' @param b Hyperparameter of beta priors on basket response rates.
#' @param d1 Design prior hyperparameter controlling propensity towards borrowing.
#' @param d2 Analysis prior hyperparameter controlling degree of borrowing.
#'
#' @return Vector of posterior efficacy probabilities under the LMEM2 framework.
#' @export
#'
#' @examples
LMEM2BasketEfficacy <- function(n_b, y, p0, a = 1, b = 1, d1 = 0, d2 = 2) {
  res <- gridSearchLMEM2Partition(n_b, y, d1 = d1, d2 = d2)
  part <- res$part
  pp <- res$post_prob
  out <- numeric(length(y))
  for (x in part) {
    ys <- y[x]
    ns <- n_b[x]
    Pb <- vapply(seq_along(ys), function(i) {
      if (length(ys) > 1) {
        p_klds <- computePosteriorKLDivergences(ns[i], ys[i], ns[-i], ys[-i])
        aij <- (exp(-p_klds) * pp)^(exp(d2))
        #aij <- pp^(exp(d2))
      } else {
        aij <- 0
      }
      a_new <- a + ys[i] + sum(aij * ys[-i])
      b_new <- b + ns[i] - ys[i] + sum(aij * (ns[-i] - ys[-i]))
      pbeta(p0, a_new, b_new, lower.tail = FALSE)
    }, numeric(1))
    out[x] <- Pb
  }
  out
}

partitionFormat <- function(part) {
  parts <- vapply(part, function(x) paste0("(", paste0(x, collapse = ","), ")"),
                  character(1))
  paste0(parts, collapse = "")
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
