library(partitions)

logsumexp_weights <- function(x) {
  xmax <- max(x)
  d <- xmax + log(sum(exp(x - xmax)))
  exp(x - d)
}

basket_weight <- function(n, yi, a, b) {
  a_new <- a + yi
  b_new <- b + n - yi
  lbeta(a_new, b_new) - lbeta(a, b)
}

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

LMEM2BasketEfficacy <- function(n_b, y, p0, a = 1, b = 1, d1 = 0, d2 = 2) {
  # n_b: numeric vector, patients in each basket
  # y: numeric vector, observed responses
  # p0: numeric scalar: historical/control response rate
  # a, b: parameters of beta-prior on basket response rates
  # d1: prior hyperparameter controlling propensity towards borrowing
  # d2: prior hyperparameter controlling degree of borrowing
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
  parts <- vapply(part, function(x) paste0("(", paste0(x, collapse = ","), ")"), character(1))
  paste0(parts, collapse = "")
}
