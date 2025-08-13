PooledBasketEfficacy <- function(n_b, y, p0, a = 1, b = 1) {
  stopifnot(length(n_b) == length(y))
  o <- pbeta(p0, a + sum(y), b + sum(n_b - y), lower.tail = FALSE)
  rep(o, length(n_b))
}
