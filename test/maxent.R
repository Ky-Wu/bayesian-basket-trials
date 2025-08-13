library(partitions)

p <- c(0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.2, 0.2)
n <- 19
y <- rbinom(length(p), n, p)
r <- signif(y / n, 3)

parts <- listParts(length(p))

computeEntropy <- function(a, b, x_n) {
  x <- seq(0, x_n, by = 1)
  lp <- lbeta(a + x, b + x_n - x) + log(choose(x_n, x)) - lbeta(a, b)
  mean(-lp)
}

ents <- vapply(parts, function(part) {
  sum(vapply(part, function(x) {
    n_b <- n * length(x)
    y_b <- sum(y[x])
    a <- 1 + y_b
    b <- 1 + n_b - y_b
    computeEntropy(a, b, n) * length(x)
  }, numeric(1)))
}, numeric(1))

#parts[[which.min(ents)]]
# Maximum entropy partition
part <- parts[[which.max(ents)]]
part
rapply(part, function(ii) r[ii], how="replace")

x <- part[[1]]

log_dbinombeta <- function(n, yi, a, b) {
  a_new <- a + yi
  b_new <- b + n - yi
  lbeta(a_new, b_new) + log(choose(n, yi)) - lbeta(a, b)
}

computePoolingBF <- function(n, ys, d = 0.5) {
  sat_ll <- exp(sum(log_dbinombeta(n, ys, 1 + sum(ys), 1 + n * length(ys) - sum(ys))))
  null_ll <- exp(sum(log_dbinombeta(n, ys, 1 + ys, 1 + n - ys)))
  (sat_ll * d) / (sat_ll * d + null_ll * (1 - d))
}


bf <- computePoolingBF(n, y[x])

computeBorrowerEntropy <- function(pd, n, y, part) {
  entropy <- sum(vapply(part, function(x) {
    ys <- y[x]
    d <- computePoolingBF(n, ys, pd)
    nb <- length(ys) * n
    sum(vapply(seq_along(ys), function(i) {
      a <- 1 + d * sum(ys[-i]) + ys[i]
      b <- 1 + d * (nb - n - sum(ys[-i])) + n - ys[i]
      computeEntropy(a, b, n)
    }, numeric(1)))
  }, numeric(1)))
  entropy
}
loss <- function(pd) {
  -computeBorrowerEntropy(pd, n, y, part = part)
}

optim(par = 0.1, loss, method = "Brent",
      lower = 0, upper = 1)

### maximum entropy always favors complete pooling
### what if we apply maximum entropy to local mem?

Kj <- vapply(parts, length, numeric(1))

computePartitionLL <- function(n, y, delta, part) {
  ll <- sum(vapply(part, function(x) {
    n_b <- n * length(x)
    y_b <- sum(y[x])
    a <- 1 + y_b
    b <- 1 + n_b - y_b
    sum(log_dbinombeta(n, y[x], a, b))
  }, numeric(1)))
  ll
}
delta <- 0
prior_weights <- Kj^delta / sum(Kj^delta)
LL <- vapply(parts, function(part) computePartitionLL(n, y, delta, part), numeric(1))
weights <- LL + prior_weights
weights <- exp(weights) / sum(exp(weights))

## maximum entropy but accounting for BF

ents2 <- vapply(parts, function(part) {
  computeBorrowerEntropy(n, y, pd = 0.5, part)
}, numeric(1))
part2 <- parts[[which.max(ents2)]]
part2
rapply(part2, function(ii) r[ii], how="replace")

vapply(part2, function(part_i) {
  computePoolingBF(n, y[part_i], d = 0.5)
}, numeric(1))


