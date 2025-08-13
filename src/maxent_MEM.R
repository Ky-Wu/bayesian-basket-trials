library(partitions)

computeBetaEntropy <- function(a, b) {
  x <- digamma(a + b)
  -(a - 1) * (digamma(a) - x) - (b - 1) * (digamma(b) - x) + lbeta(a, b)
}

computeEntropy <- function(a, b, n_x) {
  stopifnot(n_x >= 1, a > 0, b > 0)
  x <- seq(0, n_x, by = 1)
  lp <- lbeta(a + x, b + n_x - x) + lchoose(n_x, x) - lbeta(a, b)
  #mean(-lp) # old version, not correct way to compute entropy...
  -sum(exp(lp) * lp)
}

log_dbinombeta <- function(n, yi, a, b) {
  a_new <- a + yi
  b_new <- b + n - yi
  lbeta(a_new, b_new) + lchoose(n, yi) - lbeta(a, b)
}

computeCrossEntropy <- function(n, yi, a, b) {
  stopifnot(n >= yi, length(n) == 1, length(yi) == 1, a > 0, b > 0)
  -log_dbinombeta(n, yi, a, b)
}

computePoolingBF <- function(n_b, ys, exchangeability_prior = 0.5, equal_factor = TRUE) {
  #factor <- (length(n_b) - 1)
  factor <- 1
  if (equal_factor) {
    sat_ll <- exp(sum(log_dbinombeta(n_b, ys, 1 + sum(ys), 1 + sum(n_b) - sum(ys))))
    null_ll <- exp(sum(log_dbinombeta(n_b, ys, 1 + ys, 1 + n_b - ys)))
    out <- rep((sat_ll * exchangeability_prior) /
                 (sat_ll * exchangeability_prior + null_ll * (1 - exchangeability_prior)) ,
               length(n_b))
  } else {
    sat_ll <- sum(log_dbinombeta(n_b, ys, 1 + sum(ys), 1 + sum(n_b) - sum(ys)))
    null_ll <- vapply(seq_along(n_b), function(i) {
      others <- sum(log_dbinombeta(n_b[-i], ys[-i], 1 + sum(ys[-i]),
                                   1 + sum(n_b[-i]) - sum(ys[-i])))
      others + log_dbinombeta(n_b[i], ys[i], 1 + ys[i], 1 + n_b[i] - ys[i])
    }, numeric(1))
    out <- (exp(sat_ll) * exchangeability_prior) /
      (exp(sat_ll) * exchangeability_prior + exp(null_ll) * (1 - exchangeability_prior))
  }
  out / factor
}

computeBorrowerEntropy <- function(n_b, y, pd = 0.5, part, a0 = 1, b0 = 1,
                                   equal_factor = TRUE, BF_factor = TRUE) {
  entropy <- sum(vapply(part, function(x) {
    ys <- y[x]
    units <- n_b[x]
    if (BF_factor) {
      d <- computePoolingBF(units, ys, pd, equal_factor = equal_factor)
    } else {
      d <- rep(1, length(ys))
    }
    total_units <- sum(units)
    sum(vapply(seq_along(ys), function(i) {
      if (length(ys) > 1) {
        a <- a0 + d[i] * sum(ys[-i]) + ys[i]
        b <- b0 + d[i] * (total_units - units[i] - sum(ys[-i])) + units[i] - ys[i]
      } else {
        a <- a0 + ys[i]
        b <- a0 + units[i] - ys[i]
      }
      #computeEntropy(a, b, units[i])
      #computeBetaEntropy(a, b)
      #print(paste0("units:", units[i]))
      #print(paste0("response:", ys[i]))
      computeCrossEntropy(units[i], ys[i], a, b)
    }, numeric(1)))
  }, numeric(1)))
  entropy
}

gridSearchMEPartition <- function(n_b, y, exchangeability_prior = 0.5, d1 = 0,
                                  a0 = 1, b0 = 1, equal_factor = FALSE, BF_factor = FALSE) {
  stopifnot(length(y) <= 10, length(y) == length(n_b))
  parts <- listParts(length(y))
  ents <- vapply(parts, function(part) {
    #print(paste0("part:", part))
    computeBorrowerEntropy(n_b, y, pd = exchangeability_prior, part,
                           a0 = a0, b0 = b0, equal_factor = equal_factor,
                           BF_factor = BF_factor)
  }, numeric(1))
  K <- vapply(parts, length, numeric(1))
  weights <- ents + d1 * log(K)
  parts[[which.min(weights)]]
}

basketEfficacy <- function(n_b, y, p0, part, exchangeability_prior = 0.5,
                           a0 = 1, b0 = 1, equal_factor = FALSE) {
  out <- numeric(length(n_b))
  for (x in part) {
    ns <- n_b[x]
    ys <- y[x]
    if (length(ns) == 1) {
      out[x] <- pbeta(p0, a0 + ys, b0 + ns - ys, lower.tail = FALSE)
    } else {
      d <- computePoolingBF(ns, ys, exchangeability_prior, equal_factor)
      ef <- vapply(seq_along(ns), function(i) {
        a_new <- a0 + ys[i] + d[i] * sum(ys[-i])
        b_new <- b0 + ns[i] - ys[i] + d[i] * (sum(ns[-i]) - sum(ys[-i]))
        pbeta(p0, a_new, b_new, lower.tail = FALSE)
      }, numeric(1))
      out[x] <- ef
    }
  }
  out
}


ME_efficacy <- function(n_b, y, p0, pd = 0.5, d1 = 0, equal_factor = FALSE,
                        search_factor = FALSE) {
  part_ME <- gridSearchMEPartition(n_b, y, exchangeability_prior = pd, d1 = d1,
                                   a0 = 1, b0 = 1, equal_factor = equal_factor,
                                   BF_factor = search_factor)
  basketEfficacy(n_b, y, p0, part_ME, exchangeability_prior = pd,
                 a0 = 1, b0 = 1, equal_factor = equal_factor)
}

