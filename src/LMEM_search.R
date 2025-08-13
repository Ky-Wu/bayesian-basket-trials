log_dbinombeta <- function(n, yi, a, b) {
  a_new <- a + yi
  b_new <- b + n - yi
  lbeta(a_new, b_new) + lchoose(n, yi) - lbeta(a, b)
}

stepPartition <- function(n_b, y, d1 = 0, d2 = 2, a = 1, b = 1, WW_method = FALSE) {
  stopifnot(length(y) <= 10, length(y) == length(n_b))
  parts <- listParts(length(y))
  log_weights <- vapply(parts, function(part) {
    ll <- sum(vapply(part, function(x) {
      ys <- y[x]
      ns <- n_b[x]
      if (WW_method) {
        log_dbinombeta(sum(ns), sum(ys), a, b) - lchoose(sum(ns), sum(ys))
      } else {
        log_dbinombeta(sum(ns), sum(ys), a, b)
      }
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
